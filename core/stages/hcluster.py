import sys
import os
import signal
from pathlib import Path
from collections import OrderedDict
from typing import Union, List, Dict, Set, Any, Tuple, Optional
from tempfile import mkstemp
from concurrent.futures import ProcessPoolExecutor, as_completed
from threading import Lock
from core.stages.base_stage import BaseStage
from core.utils.constants import (
    FASTA_EXTENSIONS,
    CENTROIDS_SUFFIX,
    CLUSTERED_SUFFIX,
    UNCLUSTERED_SUFFIX,
    CLUSTER_ID_FORMAT
)
from core.utils.hcutil import (
    list_of_files_at_path,          # finds fasta files in a directory
    read_fasta,                     # loads sequences from a fasta file
    seq_records_to_dict,            # maps sequence IDs to their records
    write_fasta,                    # saves sequence dictionary to fasta file
    parse_newick_file,              # loads tree structure from newick format
    parse_nodes_dict_to_steps,      # converts tree nodes into clustering steps
    combine_results,                # merges clustering outputs from multiple runs
    make_cluster_sets,              # groups sequences by cluster assignments
    run_sub_hc_vsearch,             # clusters sequences using vsearch
    run_sub_hc_mmseqs2,             # clusters sequences using mmseqs2
    combine_text_files,             # merges multiple text files together
    node_to_newick,                 # serializes tree node to newick format
    create_node_mapping_text,       # generates human-readable node relationships
    calculate_node_levels,          # computes depth of each node in tree
    Node
)


def _worker_init():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def _process_group_worker(
    grp_name        : str,                      # name of the clustering group
    step_nodes      : List[str],                # nodes involved in this clustering step
    fastas_user     : List[str],                # fasta files to cluster
    nodes_fasta_dir : str,                      # directory with centroid files
    hcluster_id     : float,                    # identity threshold for clustering
    hcluster_iddef  : int,                      # identity definition for clustering
    hcluster_tool   : str,                      # clustering tool (vsearch or mmseqs2)
    threads         : int,                      # number of threads for clustering
    species_file_map: Dict[str, List[str]],     # maps species names to their fasta files
    id_to_stem      : Dict[str, str],           # maps node IDs to filename stems
    precombined_file: Optional[str],            # pre-combined fasta if available
) -> Tuple[str, Optional[Dict[str, Set[str]]], List[str], bool]:
    """
    Process a clustering group.
    """
    try:
        fastas_filt_user: List[str] = []
        fastas_filt_nodes: List[str] = []
        
        for node_name in step_nodes:
            if node_name in species_file_map:
                matching_files = species_file_map[node_name]
            else:
                if id_to_stem and node_name in id_to_stem:
                    lookup_name = id_to_stem[node_name]
                else:
                    lookup_name = node_name
                
                # Match files by stem, with fallback for BUSCO files
                matching_files = [
                    x for x in fastas_user
                    if (Path(x).stem == lookup_name or 
                        Path(x).name.split('.')[0] == lookup_name or
                        Path(x).stem.replace('_busco_renamed', '') == lookup_name or
                        Path(x).stem.replace('_busco', '') == lookup_name)
                ]
            
            fastas_filt_user.extend(matching_files)
        
        fastas_nodes, _ = list_of_files_at_path(str(nodes_fasta_dir), ".fasta")
        if fastas_nodes is None:
            fastas_nodes = []
            
        fastas_filt_nodes = [
            x for node_name in step_nodes
            for x in fastas_nodes 
            if x.endswith(f"{node_name}{CENTROIDS_SUFFIX}")
        ]
        
        fastas_filt: List[str] = fastas_filt_user + fastas_filt_nodes
        if not fastas_filt:
            return (grp_name, None, [], False)

        use_precombined = (precombined_file and 
                          set(fastas_filt_user) == set(fastas_user) and 
                          not fastas_filt_nodes)
        
        if hcluster_tool == 'mmseqs2':
            current_results: Dict[str, Set[str]] = run_sub_hc_mmseqs2(
                fastas_filt,
                out_dir          = str(nodes_fasta_dir),
                ident            = hcluster_id,
                iddef            = hcluster_iddef,
                grp_name         = grp_name,
                threads          = threads,
                precombined_file = precombined_file if use_precombined else None,
            )
        else:
            current_results: Dict[str, Set[str]] = run_sub_hc_vsearch(
                fastas_filt,
                out_dir          = str(nodes_fasta_dir),
                ident            = hcluster_id,
                iddef            = hcluster_iddef,
                grp_name         = grp_name,
                threads          = threads,
                precombined_file = precombined_file if use_precombined else None,
            )
        
        prev_nodes_used: List[str] = [
            Path(x).stem.replace(".centroids", "") 
            for x in fastas_filt_nodes
        ]
        
        return (grp_name, current_results, prev_nodes_used, True)
        
    except Exception as e:
        print(f"Error processing group {grp_name}: {str(e)}", file=sys.stderr)
        return (grp_name, None, [], False)


class HCluster(BaseStage):    
    def __init__(self,
                 dir_base      : str,
                 dir_treeforge : str,
                 dir_hcluster  : str,
                 files_fasta   : List[str],
                 file_tree     : str,
                 threads       : int,
                 log           : str,
                 hc            : str,
                 bc            : str,
                 hcluster_id   : float,
                 hcluster_iddef: int,
                 hcluster_tool : str,
                 species_file_map: Dict[str, List[str]] = None,
                 id_to_stem: Dict[str, str] = None,
                 shared_printClass=None):

        super().__init__(log, hc, bc, threads, shared_printClass=shared_printClass)
        self.dir_base            = Path(dir_base)
        self.dir_treeforge       = Path(dir_treeforge)
        self.dir_hcluster        = Path(dir_hcluster)
        self.files_fasta         = files_fasta
        self.file_tree           = file_tree
        self.hcluster_id         = hcluster_id
        self.hcluster_iddef      = hcluster_iddef
        self.hcluster_tool       = hcluster_tool
        self.species_file_map    = species_file_map if species_file_map is not None else {}
        self.id_to_stem          = id_to_stem if id_to_stem is not None else {}
        
        self.nodes_fasta_dir     = self.dir_hcluster / "vsearch"
        self.clustered_dir       = self.dir_hcluster / "clusters"
        self.unclustered_dir     = self.dir_hcluster / "clusters"
        
        self.cluster_count       = 0
        self.total_sequences     = 0
        self.clustered_sequences = 0
        self.guide_tree_depth    = 0
        self.guide_tree_clades   = 0
        
        self.fasta_cache: Dict[str, Dict[str, str]] = {}
        self.progress_lock = Lock()
        
        self.return_dict         = {
            'cluster_count'      : self.cluster_count,
            'total_sequences'    : self.total_sequences,
            'clustered_sequences': self.clustered_sequences,
            'nodes_fasta_dir'    : str(self.nodes_fasta_dir),
            'clustered_dir'      : str(self.clustered_dir),
            'unclustered_dir'    : str(self.unclustered_dir)
        }
    
    def get_fasta_records(self, fp: str) -> Dict[str, str]:
        """Get FASTA records with caching"""
        if fp not in self.fasta_cache:
            records = read_fasta(fp)
            self.fasta_cache[fp] = seq_records_to_dict(records)
        return self.fasta_cache[fp]
    
    def _create_precombined_file(self, fastas_user: List[str]) -> Optional[str]:
        try:
            try:
                _, tmpf = mkstemp(text=True, suffix='.fasta', prefix='hcluster_combined_')
            except TypeError:
                _, tmpf = mkstemp(suffix='.fasta', prefix='hcluster_combined_')
            combine_text_files(fastas_user, tmpf)
            # self.printout('metric', f'Pre-combined {len(fastas_user)} FASTA files for efficiency')
            return tmpf
        except Exception as e:
            self.printout('warning', f'Failed to pre-combine FASTA files: {str(e)}')
            return None
        
    def run(self) -> Dict[str, Union[int, str]]:
        """
        Run hierarchical clustering.
        """
        self.nodes_fasta_dir.mkdir(parents=True, exist_ok=True)
        self.clustered_dir.mkdir(parents=True, exist_ok=True)
        self.unclustered_dir.mkdir(parents=True, exist_ok=True)

        self.printout('metric', 'Parsing Newick tree')
        tree: Node = parse_newick_file(self.file_tree)
        self.printout('metric', 'Parse Nodes to Steps')
        steps: OrderedDict[str, List[str]] = parse_nodes_dict_to_steps(tree)[0]
        n_groups = len(steps)
        
        node_levels: Dict[str, int] = calculate_node_levels(tree)
        self.printout('metric', [('clustering_groups', n_groups)])

        if n_groups == 0:
            self.printout('error', 'No clustering groups found in species tree')
            return self.return_dict
        
        labeled_tree_path = self.clustered_dir / "guide_tree_labeled.nwk"
        node_mapping_path = self.clustered_dir / "guide_tree_node_mapping.txt"
        
        try:
            newick_str = node_to_newick(tree)
            with open(labeled_tree_path, 'w') as f:
                f.write(newick_str)
        except Exception as e:
            self.printout('error', f'Failed to save labeled tree: {str(e)}')
        
        try:
            mapping_text = create_node_mapping_text(tree)
            with open(node_mapping_path, 'w') as f:
                f.write(mapping_text)
        except Exception as e:
            self.printout('error', f'Failed to save node mapping: {str(e)}')

        
        fastas_user: List[str] = [str(f) for f in self.files_fasta]
        if not fastas_user:
            self.printout('error', "No FASTA files found")
            return self.return_dict

        self.printout('metric', [('files_processed', len(fastas_user))])

        self.total_sequences = sum(len(self.get_fasta_records(fp)) for fp in fastas_user)
        self.printout('metric', [('total_sequences', self.total_sequences)])
        
        precombined_file = self._create_precombined_file(fastas_user)
        
        results: Dict[str, Dict[str, Set[str]]] = {}
        total_groups = n_groups
        completed_groups = 0
        
        completed_percentage = f'{completed_groups/total_groups*100:.1f}%'
        self.printout('progress', f"HCluster progress: {completed_groups:0{len(str(total_groups))}d}/{total_groups} ({completed_percentage}) groups completed")

        levels_to_groups: Dict[int, List[str]] = {}
        for grp_name in steps.keys():
            level = node_levels.get(grp_name, 0)
            if level not in levels_to_groups:
                levels_to_groups[level] = []
            levels_to_groups[level].append(grp_name)
        
        max_level = max(levels_to_groups.keys()) if levels_to_groups else 0
        self.guide_tree_depth  = max_level + 1
        self.guide_tree_clades = n_groups
        self.printout('metric', [('guide_tree_clades', n_groups), ('guide_tree_depth', self.guide_tree_depth)])
        
        min_threads_per_worker = 4
        
        if self.threads >= min_threads_per_worker * 2:
            max_workers = max(1, min(self.threads // min_threads_per_worker, n_groups))
            threads_per_worker = max(min_threads_per_worker, self.threads // max_workers)
        else:
            max_workers = 1
            threads_per_worker = self.threads
        
        # self.printout('metric', f'Using {max_workers} parallel workers ({threads_per_worker} threads each, total: {max_workers * threads_per_worker})')
        
        executor = None
        try:
            for level in sorted(levels_to_groups.keys()):
                level_groups = levels_to_groups[level]
                self.printout('metric', f'Processing depth {level}: {len(level_groups)} clades')
                
                try:
                    executor = ProcessPoolExecutor(max_workers=max_workers, initializer=_worker_init)
                except TypeError:
                    executor = ProcessPoolExecutor(max_workers=max_workers)
                try:
                    future_to_group = {}
                    for grp_name in level_groups:
                        step_nodes = steps[grp_name]
                        future = executor.submit(
                            _process_group_worker,
                            grp_name,                       # name of the clustering group
                            step_nodes,                     # nodes involved in this clustering step
                            fastas_user,                    # fasta files to cluster
                            str(self.nodes_fasta_dir),      # directory with centroid files
                            self.hcluster_id,               # identity threshold for clustering
                            self.hcluster_iddef,            # identity definition for clustering
                            self.hcluster_tool,             # clustering tool (vsearch or mmseqs2)
                            threads_per_worker,             # number of threads for clustering
                            self.species_file_map,          # maps species names to their fasta files
                            self.id_to_stem,                # maps node IDs to filename stems
                            precombined_file,               # pre-combined fasta if available
                        )
                        future_to_group[future] = grp_name
                    
                    for future in as_completed(future_to_group):
                        grp_name, current_results, prev_nodes_used, success = future.result()
                        
                        if not success or current_results is None:
                            if current_results is None:
                                self.printout('info', f'No files matched for group {grp_name}')
                        else:
                            if not prev_nodes_used:
                                results[grp_name] = current_results
                            else:
                                results_to_combine: List[Dict[str, Set[str]]] = [current_results]
                                for pnu in prev_nodes_used:
                                    if pnu in results:
                                        results_to_combine.append(results[pnu])
                                results_combined: Dict[str, Set[str]] = combine_results(results_to_combine)
                                results[grp_name] = results_combined
                            
                            clusters: List[Set[str]] = make_cluster_sets(results[grp_name])
                            clusters.sort(key=len, reverse=True)
                            self.write_clusters(fastas_user, clusters, grp_name)
                        
                        with self.progress_lock:
                            completed_groups += 1
                            completed_percentage = f'{completed_groups/total_groups*100:.1f}%'
                            self.printout('progress', f"HCluster progress: {completed_groups:0{len(str(total_groups))}d}/{total_groups} ({completed_percentage}) groups completed")
                finally:
                    if executor is not None:
                        executor.shutdown(wait=True)
                        executor = None
        
        except KeyboardInterrupt:
            if executor is not None:
                try:
                    executor.shutdown(wait=False, cancel_futures=True)
                except:
                    pass
            self.printout('error', 'Process interrupted by user')
            raise
        
        finally:
            if precombined_file and Path(precombined_file).exists():
                try:
                    Path(precombined_file).unlink()
                except Exception as e:
                    self.printout('warning', f'Failed to clean up temporary file: {str(e)}')
                        
        self.return_dict.update({
            'cluster_count'      : self.cluster_count,
            'total_sequences'    : self.total_sequences,
            'clustered_sequences': self.clustered_sequences,
            'guide_tree_depth'   : self.guide_tree_depth,
            'guide_tree_clades'  : self.guide_tree_clades,
        })
        
        self.printout('metric', [
            ('cluster_count', self.cluster_count),
            ('clustered_sequences', self.clustered_sequences)
        ])
        return self.return_dict

    def write_clusters(self, fastas: List[str], clusters: List[Set[str]], grp_name: str) -> None:
        """
        Write clustered and unclustered sequences to FASTA files.
        """
        all_seq_recs: Dict[str, str] = {}
        for fp in fastas:
            try:
                recs = self.get_fasta_records(fp)
                all_seq_recs.update(recs)
            except Exception as e:
                self.printout('warning', f'Failed to read FASTA file {fp}: {str(e)}')
                continue
        
        centroid_files, _ = list_of_files_at_path(str(self.nodes_fasta_dir), ".fasta")
        if centroid_files:
            for fp in centroid_files:
                if fp.endswith(CENTROIDS_SUFFIX):
                    try:
                        all_seq_recs.update(self.get_fasta_records(fp))
                    except Exception as e:
                        self.printout('warning', f'Failed to read centroid file {fp}: {str(e)}')
                        continue
            
        clustered_seq_recs  : Dict[str, str] = {}
        unclustered_seq_recs: Dict[str, str] = all_seq_recs.copy()
        
        for i, seq_name_set in enumerate(clusters):
            self.cluster_count            += 1
            cluster_id: str                = CLUSTER_ID_FORMAT.format(i + 1)
            seq_names : List[str]          = sorted(list(seq_name_set))
            
            seq_recs: Dict[str, str] = {}
            missing_seqs = []
            for sn in seq_names:
                if sn in all_seq_recs:
                    seq_recs[f"{cluster_id}||{sn}"] = all_seq_recs[sn]
                else:
                    missing_seqs.append(sn)
            
            if missing_seqs:
                self.printout('warning', f'Group {grp_name}: {len(missing_seqs)} sequences not found in input files')
                if len(missing_seqs) <= 5:
                    for sn in missing_seqs:
                        self.printout('warning', f'  Missing: {sn}')
            
            if not seq_recs:
                self.printout('warning', f'No sequences found for cluster {cluster_id} in group {grp_name}')
                continue
            
            clustered_seq_recs.update(seq_recs)
            self.clustered_sequences += len(seq_recs)
            
            for sn in seq_names:
                unclustered_seq_recs.pop(sn, None)
                
        out_fp_clustered: str = str(self.clustered_dir / f"{grp_name}{CLUSTERED_SUFFIX}")
        write_fasta(clustered_seq_recs, out_fp_clustered)
        
        out_fp_unclustered: str = str(self.unclustered_dir / f"{grp_name}{UNCLUSTERED_SUFFIX}")
        write_fasta(unclustered_seq_recs, out_fp_unclustered)