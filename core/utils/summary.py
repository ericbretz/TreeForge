import csv
import json
from pathlib import Path
from typing import Dict, Any, List, Optional

import ete3

from core.utils.decor        import convert_paths
from core.utils.iqtree_parser import parse_iqtree_directory
from core.utils.astral_parser import parse_verbose_csv


class MetricsSummary:
    """generate summary reports from master_dict"""

    def __init__(self, master_dict: Dict[str, Any], dir_metrics: Path,
                 hcluster_enabled: bool = False, output_super_matrix: bool = False,
                 run_timestamp: str = '',
                 dir_prank: Optional[str] = None,
                 dir_results: Optional[str] = None):
        self.master_dict         = master_dict
        self.dir_metrics         = Path(dir_metrics)
        self.hcluster_enabled    = hcluster_enabled
        self.output_super_matrix = output_super_matrix
        self.run_timestamp       = run_timestamp
        self.dir_prank           = Path(dir_prank)   if dir_prank   else None
        self.dir_results         = Path(dir_results) if dir_results else None

    def generate_all_reports(self) -> Dict[str, Path]:
        reports = {}

        reports['json'] = self.generate_json()

        reports['iteration_flow']     = self.generate_iteration_flow_csv()
        reports['pipeline_summary']   = self.generate_pipeline_summary_csv()
        reports['clustering_summary'] = self.generate_clustering_summary_csv()
        reports['mafft_summary']      = self.generate_mafft_summary_csv()
        reports['tree_processing']    = self.generate_tree_processing_csv()
        reports['prune_metrics']      = self.generate_prune_metrics_csv()
        reports['final_alignment']    = self.generate_final_alignment_csv()

        reports['gene_tree_stats']    = self.generate_gene_tree_stats_csv()
        reports['alignment_quality']  = self.generate_alignment_quality_csv()
        reports['species_coverage']   = self.generate_species_coverage_csv()
        reports['astral_support']     = self.generate_astral_support_csv()
        reports['astral_summary']     = self.dir_metrics / 'astral_summary.csv'

        if self.hcluster_enabled:
            reports['hcluster_metrics'] = self.generate_hcluster_csv()

        return reports

    def generate_json(self) -> Path:
        json_file = self.dir_metrics / 'run_metrics.json'
        with open(json_file, 'w') as f:
            json.dump(convert_paths(self.master_dict), f, indent=2)
        return json_file

    def generate_iteration_flow_csv(self) -> Path:
        csv_file = self.dir_metrics / 'iteration_flow.csv'
        rows     = []

        mcl_count = self.master_dict.get('mcl', {}).get('mcl_count', 0)

        iter_keys = sorted([k for k in self.master_dict.keys() if k.startswith('iter_')],
                           key=lambda x: int(x.split('_')[1]))

        for i, iter_key in enumerate(iter_keys):
            iter_data = self.master_dict[iter_key]

            mafft_data   = iter_data.get('mafft', {}).get('mafft', {})
            after_mafft  = len(mafft_data.get('aln_files', []))
            after_pxclsq = len(mafft_data.get('cln_files', []))
            after_iqtree = len(mafft_data.get('iqtree_files', []))

            if i == 0:
                input_count = mcl_count
            else:
                input_count = after_mafft

            tree_data       = iter_data.get('tree', {})
            trim_data       = tree_data.get('trim_tips', {})
            mask_data       = tree_data.get('mask_tips', {})
            write_data      = tree_data.get('write_tree', {})
            after_trimtips  = trim_data.get('processed_count', after_iqtree)
            after_masktips  = mask_data.get('processing_results', {}).get('trees_processed', after_trimtips)
            to_next         = len(write_data.get('written_files', []))
            to_prune        = write_data.get('files_moved_to_prune', 0)
            after_cutbranch = to_prune + to_next

            row = {
                'Iteration'      : iter_key.replace('_', '_'),
                'Input'          : input_count,
                'After_MAFFT'    : after_mafft,
                'After_pxclsq'   : after_pxclsq,
                'After_IQTree'   : after_iqtree,
                'After_TrimTips' : after_trimtips,
                'After_MaskTips' : after_masktips,
                'After_CutBranch': after_cutbranch,
                'To_Prune'       : to_prune,
                'To_Next_Iter'   : to_next if to_next > 0 else ''
            }
            rows.append(row)

        prune_data  = self.master_dict.get('prune', {}).get('prune', {})
        prune_input = len(prune_data.get('ortho1to1', []))

        rows.append({
            'Iteration'      : 'Prune',
            'Input'          : prune_input,
            'After_MAFFT'    : '',
            'After_pxclsq'   : '',
            'After_IQTree'   : '',
            'After_TrimTips' : '',
            'After_MaskTips' : '',
            'After_CutBranch': '',
            'To_Prune'       : '',
            'To_Next_Iter'   : ''
        })

        fieldnames = ['Iteration', 'Input', 'After_MAFFT', 'After_pxclsq', 'After_IQTree',
                      'After_TrimTips', 'After_MaskTips', 'After_CutBranch',
                      'To_Prune', 'To_Next_Iter']

        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        return csv_file

    def generate_pipeline_summary_csv(self) -> Path:
        csv_file = self.dir_metrics / 'pipeline_summary.csv'

        iter_keys        = [k for k in self.master_dict.keys() if k.startswith('iter_')]
        total_iterations = len(iter_keys)

        total_iter_time = 0.0
        for iter_key in iter_keys:
            iter_data        = self.master_dict.get(iter_key, {})
            mafft_time       = self._parse_time(iter_data.get('mafft', {}).get('metrics', {}).get('time', '0 s'))
            tree_time        = self._parse_time(iter_data.get('tree', {}).get('metrics', {}).get('time', '0 s'))
            total_iter_time += mafft_time + tree_time

        blast_data    = self.master_dict.get('blast', {})
        mcl_data      = self.master_dict.get('mcl', {})
        prune_data    = self.master_dict.get('prune', {}).get('prune', {})
        prune_metrics = prune_data.get('metrics', {})
        prank_data    = self.master_dict.get('prank', {})
        astral_data   = self.master_dict.get('astral', {})

        blast_time    = self._parse_time(blast_data.get('metrics', {}).get('time', '0 s'))
        mcl_time      = self._parse_time(mcl_data.get('metrics', {}).get('time', '0 s'))
        prune_time    = self._parse_time(self.master_dict.get('prune', {}).get('metrics', {}).get('time', '0 s'))
        prank_time    = self._parse_time(prank_data.get('metrics', {}).get('time', '0 s'))
        astral_time   = self._parse_time(astral_data.get('metrics', {}).get('time', '0 s'))

        total_pipeline_time = blast_time + mcl_time + total_iter_time + prune_time + prank_time + astral_time

        mean_ll = mean_len = mean_pi = ''
        if self.dir_prank and self.dir_prank.is_dir():
            gene_stats = parse_iqtree_directory(str(self.dir_prank))
            ll_vals    = [s['log_likelihood']            for s in gene_stats if s['log_likelihood']            is not None]
            len_vals   = [s['total_tree_length']         for s in gene_stats if s['total_tree_length']         is not None]
            pi_vals    = [s['pct_parsimony_informative'] for s in gene_stats if s['pct_parsimony_informative']  is not None]
            if ll_vals:  mean_ll  = f'{sum(ll_vals)  / len(ll_vals):.4f}'
            if len_vals: mean_len = f'{sum(len_vals) / len(len_vals):.6f}'
            if pi_vals:  mean_pi  = f'{sum(pi_vals)  / len(pi_vals):.2f}'

        astral_mean_support = astral_pct_gt_095 = ''
        verbose_csv = self._find_verbose_csv()
        if verbose_csv:
            ap = parse_verbose_csv(verbose_csv)
            if ap['num_branches']:
                astral_mean_support = f"{ap['mean_support']:.6f}"
                astral_pct_gt_095   = f"{ap['pct_branches_gt_0_95']:.2f}"

        row = {
            'run_timestamp'                 : self.run_timestamp,
            'contig_count'                  : blast_data.get('contig_count', 0),
            'mcl_count'                     : mcl_data.get('mcl_count', 0),
            'mcl_inflation'                 : mcl_data.get('mcl_inflation', ''),
            'mcl_hit_frac'                  : mcl_data.get('mcl_hit_frac', ''),
            'mcl_eval_power'                : mcl_data.get('mcl_eval_power', ''),
            'total_iterations'              : total_iterations,
            'ortho1to1_count'               : prune_metrics.get('ortho1to1_count', 0),
            'ortho_mi_count'                : prune_metrics.get('ortho_mi_count', 0),
            'num_trees'                     : astral_data.get('num_trees', 0),
            'finaltree'                     : astral_data.get('finaltree'),
            'super'                         : self.output_super_matrix,
            'supermatrix_taxa'              : astral_data.get('supermatrix_taxa', ''),
            'supermatrix_bp'                : astral_data.get('supermatrix_bp', ''),
            'bes_tree_count'                : self._resolve_bes_tree_count(astral_data),
            'mean_gene_tree_log_likelihood' : mean_ll,
            'mean_gene_tree_length'         : mean_len,
            'mean_pct_parsimony_informative': mean_pi,
            'astral_mean_branch_support'    : astral_mean_support,
            'astral_pct_branches_gt_0_95'   : astral_pct_gt_095,
            'blast_time'                    : f'{blast_time:.2f}',
            'mcl_time'                      : f'{mcl_time:.2f}',
            'iteration_time'                : f'{total_iter_time:.2f}',
            'prune_time'                    : f'{prune_time:.2f}',
            'prank_time'                    : f'{prank_time:.2f}',
            'astral_time'                   : f'{astral_time:.2f}',
            'total_time'                    : f'{total_pipeline_time:.2f}'
        }

        fieldnames = [
            'run_timestamp', 'contig_count', 'mcl_count',
            'mcl_inflation', 'mcl_hit_frac', 'mcl_eval_power',
            'total_iterations', 'ortho1to1_count', 'ortho_mi_count',
            'num_trees', 'finaltree', 'super',
            'supermatrix_taxa', 'supermatrix_bp', 'bes_tree_count',
            'mean_gene_tree_log_likelihood', 'mean_gene_tree_length',
            'mean_pct_parsimony_informative',
            'astral_mean_branch_support', 'astral_pct_branches_gt_0_95',
            'blast_time', 'mcl_time', 'iteration_time',
            'prune_time', 'prank_time', 'astral_time', 'total_time'
        ]

        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(row)

        return csv_file

    def generate_clustering_summary_csv(self) -> Path:
        csv_file   = self.dir_metrics / 'clustering_summary.csv'
        rows       = []

        blast_data = self.master_dict.get('blast', {})
        mcl_data   = self.master_dict.get('mcl', {})

        blast_row = {
            'stage'         : 'BLAST',
            'contig_count'  : blast_data.get('contig_count', 0),
            'raw_blast'     : blast_data.get('files', {}).get('raw_blast'),
            'concatenated'  : blast_data.get('files', {}).get('concatenated'),
            'mcl_count'     : '',
            'mcl_out'       : '',
            'ident_out'     : '',
            'blast_mcl_out' : '',
            'mcl_inflation' : '',
            'mcl_hit_frac'  : '',
            'mcl_eval_power': '',
            'time'          : self._parse_time(blast_data.get('metrics', {}).get('time', '0 s')),
        }
        rows.append(blast_row)

        mcl_row = {
            'stage'         : 'MCL',
            'contig_count'  : '',
            'raw_blast'     : '',
            'concatenated'  : '',
            'mcl_count'     : mcl_data.get('mcl_count', 0),
            'mcl_out'       : mcl_data.get('mcl_out'),
            'ident_out'     : mcl_data.get('ident_out'),
            'blast_mcl_out' : mcl_data.get('blast_mcl_out'),
            'mcl_inflation' : mcl_data.get('mcl_inflation', ''),
            'mcl_hit_frac'  : mcl_data.get('mcl_hit_frac', ''),
            'mcl_eval_power': mcl_data.get('mcl_eval_power', ''),
            'time'          : self._parse_time(mcl_data.get('metrics', {}).get('time', '0 s')),
        }
        rows.append(mcl_row)

        fieldnames = ['stage', 'contig_count', 'raw_blast', 'concatenated',
                      'mcl_count', 'mcl_out', 'ident_out', 'blast_mcl_out',
                      'mcl_inflation', 'mcl_hit_frac', 'mcl_eval_power', 'time']

        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        return csv_file

    def generate_mafft_summary_csv(self) -> Path:
        csv_file = self.dir_metrics / 'mafft_summary.csv'
        rows     = []

        iter_keys = sorted([k for k in self.master_dict.keys() if k.startswith('iter_')],
                           key=lambda x: int(x.split('_')[1]))

        mcl_count = self.master_dict.get('mcl', {}).get('mcl_count', 0)

        for i, iter_key in enumerate(iter_keys):
            iter_data    = self.master_dict[iter_key]
            mafft_data   = iter_data.get('mafft', {}).get('mafft', {})

            aln_files    = mafft_data.get('aln_files', [])
            cln_files    = mafft_data.get('cln_files', [])
            iqtree_files = mafft_data.get('iqtree_files', [])

            if i == 0:
                input_count = mcl_count
            else:
                prev_iter       = iter_keys[i - 1]
                prev_write_data = self.master_dict[prev_iter].get('tree', {}).get('write_tree', {})
                input_count     = len(prev_write_data.get('written_files', []))

            aln_count  = len(aln_files)
            cln_count  = len(cln_files)
            tree_count = len(iqtree_files)

            success_rate = (tree_count / input_count * 100) if input_count > 0 else 0

            rows.append({
                'iteration'      : iter_key,
                'input_clusters' : input_count,
                'aln_files'      : aln_count,
                'cln_files'      : cln_count,
                'iqtree_files'   : tree_count,
                'time'           : self._parse_time(iter_data.get('mafft', {}).get('metrics', {}).get('time', '0 s')),
                'success_rate'   : f'{success_rate:.1f}'
            })

        fieldnames = ['iteration', 'input_clusters', 'aln_files', 'cln_files',
                      'iqtree_files', 'time', 'success_rate']

        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        return csv_file

    def generate_tree_processing_csv(self) -> Path:
        csv_file  = self.dir_metrics / 'tree_processing.csv'
        rows      = []

        iter_keys = sorted([k for k in self.master_dict.keys() if k.startswith('iter_')],
                           key=lambda x: int(x.split('_')[1]))

        for iter_key in iter_keys:
            iter_data = self.master_dict[iter_key]
            tree_data = iter_data.get('tree', {})

            # TrimTips
            trim_data = tree_data.get('trim_tips', {})
            trim_thresh = trim_data.get('trimming_thresholds', {})
            rows.append({
                'iteration'                  : iter_key,
                'sub_stage'                  : 'trim_tips',
                'total_trees'                : trim_data.get('total_trees', 0),
                'trees_processed'            : trim_data.get('processed_count', 0),
                'trees_skipped'              : trim_data.get('skipped_count', 0),
                'error_count'                : trim_data.get('error_count', 0),
                'cutting_threshold_relative' : trim_thresh.get('relative_cutoff', ''),
                'cutting_threshold_absolute' : trim_thresh.get('absolute_cutoff', ''),
                'written_files'              : '',
                'files_moved_to_prune'       : '',
                'total_sequences'            : '',
                'modified_files_count'       : '',
                'monophyletic_tips_masked'   : '',
                'paraphyletic_tips_masked'   : '',
                'time'                       : ''
            })

            # MaskTips
            mask_data    = tree_data.get('mask_tips', {})
            mask_results = mask_data.get('processing_results', {})
            mask_masking = mask_data.get('masking_results', {})
            rows.append({
                'iteration'                  : iter_key,
                'sub_stage'                  : 'mask_tips',
                'total_trees'                : mask_data.get('file_counts', {}).get('total_tree_files', 0),
                'trees_processed'            : mask_results.get('trees_processed', 0),
                'trees_skipped'              : mask_results.get('trees_skipped', 0),
                'error_count'                : (mask_results.get('file_read_errors', 0)
                                                + mask_results.get('tree_parse_errors', 0)),
                'cutting_threshold_relative' : '',
                'cutting_threshold_absolute' : '',
                'written_files'              : '',
                'files_moved_to_prune'       : '',
                'total_sequences'            : '',
                'modified_files_count'       : '',
                'monophyletic_tips_masked'   : mask_masking.get('monophyletic_tips_masked', 0),
                'paraphyletic_tips_masked'   : mask_masking.get('paraphyletic_tips_masked', 0),
                'time'                       : ''
            })

            # WriteTree / CutBranches
            write_data = tree_data.get('write_tree', {})
            to_next    = len(write_data.get('written_files', []))
            to_prune   = write_data.get('files_moved_to_prune', 0)
            rows.append({
                'iteration'                  : iter_key,
                'sub_stage'                  : 'cut_branches',
                'total_trees'                : to_next + to_prune,
                'trees_processed'            : to_next + to_prune,
                'trees_skipped'              : write_data.get('skipped_count', 0),
                'error_count'                : write_data.get('error_count', 0),
                'cutting_threshold_relative' : '',
                'cutting_threshold_absolute' : '',
                'written_files'              : to_next,
                'files_moved_to_prune'       : to_prune,
                'total_sequences'            : write_data.get('total_sequences', ''),
                'modified_files_count'       : write_data.get('modified_files_count', ''),
                'monophyletic_tips_masked'   : '',
                'paraphyletic_tips_masked'   : '',
                'time'                       : self._parse_time(tree_data.get('metrics', {}).get('time', '0 s'))
            })

        fieldnames = [
            'iteration', 'sub_stage', 'total_trees', 'trees_processed', 'trees_skipped',
            'error_count', 'cutting_threshold_relative', 'cutting_threshold_absolute',
            'written_files', 'files_moved_to_prune', 'total_sequences', 'modified_files_count',
            'monophyletic_tips_masked', 'paraphyletic_tips_masked', 'time'
        ]

        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        return csv_file

    def generate_prune_metrics_csv(self) -> Path:
        csv_file     = self.dir_metrics / 'prune_metrics.csv'

        prune_data   = self.master_dict.get('prune', {}).get('prune', {})
        metrics      = prune_data.get('metrics', {})

        total_trees  = metrics.get('total_trees', 0)
        ortho1to1    = metrics.get('ortho1to1_count', 0)
        pruned_trees = metrics.get('pruned_trees', 0)
        ortho_mi     = metrics.get('ortho_mi_count', 0)
        skipped      = metrics.get('skipped_trees', 0)

        success_rate = ((ortho1to1 + pruned_trees) / total_trees * 100) if total_trees > 0 else 0

        row = {
            'total_trees'     : total_trees,
            'ortho1to1_count' : ortho1to1,
            'pruned_trees'    : pruned_trees,
            'ortho_mi_count'  : ortho_mi,
            'skipped_trees'   : skipped,
            'success_rate'    : f'{success_rate:.1f}',
            'time'            : self._parse_time(
                self.master_dict.get('prune', {}).get('metrics', {}).get('time', '0 s'))
        }

        fieldnames = ['total_trees', 'ortho1to1_count', 'pruned_trees',
                      'ortho_mi_count', 'skipped_trees', 'success_rate', 'time']

        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(row)

        return csv_file

    def generate_final_alignment_csv(self) -> Path:
        csv_file           = self.dir_metrics / 'final_alignment.csv'
        rows               = []

        prank_data         = self.master_dict.get('prank', {})
        astral_data        = self.master_dict.get('astral', {})

        prune_data         = self.master_dict.get('prune', {}).get('prune', {})
        prune_metrics      = prune_data.get('metrics', {})
        input_trees        = prune_metrics.get('ortho1to1_count', 0)

        alignments         = prank_data.get('alignments', [])
        prank_trees        = prank_data.get('trees', [])
        cleaned_alignments = prank_data.get('cleaned_alignments', [])

        prank_row = {
            'stage'                : 'prank',
            'fasta_from_tree'      : input_trees,
            'alignments'           : len(alignments)         if isinstance(alignments, list)         else 0,
            'cleaned_alignments'   : len(cleaned_alignments) if isinstance(cleaned_alignments, list) else 0,
            'prank_trees_generated': len(prank_trees)        if isinstance(prank_trees, list)        else 0,
            'gene_trees'           : '',
            'supermatrix_taxa'     : '',
            'supermatrix_bp'       : '',
            'bes_tree_count'       : '',
            'time'                 : self._parse_time(prank_data.get('metrics', {}).get('time', '0 s'))
        }
        rows.append(prank_row)

        astral_row = {
            'stage'                : 'astral',
            'fasta_from_tree'      : '',
            'alignments'           : '',
            'cleaned_alignments'   : '',
            'prank_trees_generated': '',
            'gene_trees'           : astral_data.get('num_trees', 0),
            'supermatrix_taxa'     : astral_data.get('supermatrix_taxa', ''),
            'supermatrix_bp'       : astral_data.get('supermatrix_bp', ''),
            'bes_tree_count'       : self._resolve_bes_tree_count(astral_data),
            'time'                 : self._parse_time(astral_data.get('metrics', {}).get('time', '0 s'))
        }
        rows.append(astral_row)

        fieldnames = ['stage', 'fasta_from_tree', 'alignments', 'cleaned_alignments',
                      'prank_trees_generated', 'gene_trees',
                      'supermatrix_taxa', 'supermatrix_bp', 'bes_tree_count', 'time']

        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        return csv_file

    def generate_hcluster_csv(self, hcluster_metrics: Dict[str, Any] = None) -> Path:
        csv_file             = self.dir_metrics / 'hcluster_metrics.csv'

        hc_metrics           = hcluster_metrics if hcluster_metrics else {}

        hcluster_data        = hc_metrics.get('hierarchical_clustering', {})
        config               = hc_metrics.get('configuration', {})
        gene_cluster         = hc_metrics.get('gene_clustering', {})
        guide_tree           = hc_metrics.get('guide_tree', {})
        busco                = hc_metrics.get('busco_extraction', {})
        input_files          = hc_metrics.get('input_files', {})

        total_sequences      = hcluster_data.get('total_sequences', 0)
        cluster_count        = hcluster_data.get('cluster_count', 0)
        clustered_sequences  = hcluster_data.get('clustered_sequences', 0)

        clustering_tool      = config.get('hcluster_tool', '')
        identity_threshold   = config.get('hcluster_id', '')
        identity_def         = config.get('hcluster_iddef', '')
        threads              = config.get('threads', '')
        guide_tree_depth     = hcluster_data.get('guide_tree_depth', '')
        guide_tree_clades    = hcluster_data.get('guide_tree_clades', '')

        input_file_count     = input_files.get('count', 0)
        busco_count          = busco.get('busco_count', 0)
        gene_count           = gene_cluster.get('gene_count', 0)
        gene_trees_count     = guide_tree.get('gene_trees_count', 0)

        unclustered          = (total_sequences - clustered_sequences
                                    if clustered_sequences <= total_sequences else 0)
        efficiency           = (clustered_sequences / total_sequences * 100) if total_sequences > 0 else 0

        seq_per_cluster_mean = (round(clustered_sequences / cluster_count, 2)
                                if cluster_count > 0 else '')

        row = {
            'Run_Timestamp'                : self.run_timestamp,
            'Input_Files'                  : input_file_count,
            'BUSCO_Extraction_Count'       : busco_count,
            'Gene_Clusters_Created'        : gene_count,
            'Guide_Tree_Gene_Count'        : gene_trees_count,
            'Guide_Tree_Depth'             : guide_tree_depth,
            'Guide_Tree_Clades'            : guide_tree_clades,
            'Total_Input_Sequences'        : total_sequences,
            'Sequence_Clusters_Produced'   : cluster_count,
            'Total_Sequences_Clustered'    : clustered_sequences,
            'Sequences_Unclustered'        : unclustered,
            'Clustering_Efficiency_Percent': f'{efficiency:.1f}',
            'Sequences_Per_Cluster_Mean'   : seq_per_cluster_mean,
            'Clustering_Tool'              : clustering_tool,
            'Identity_Threshold'           : identity_threshold,
            'Identity_Definition'          : identity_def,
            'Threads_Used'                 : threads,
        }

        fieldnames = [
            'Run_Timestamp', 'Input_Files', 'BUSCO_Extraction_Count',
            'Gene_Clusters_Created', 'Guide_Tree_Gene_Count',
            'Guide_Tree_Depth', 'Guide_Tree_Clades',
            'Total_Input_Sequences', 'Sequence_Clusters_Produced',
            'Total_Sequences_Clustered', 'Sequences_Unclustered',
            'Clustering_Efficiency_Percent', 'Sequences_Per_Cluster_Mean',
            'Clustering_Tool', 'Identity_Threshold', 'Identity_Definition', 'Threads_Used'
        ]

        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(row)

        return csv_file

    def generate_gene_tree_stats_csv(self) -> Path:
        csv_file = self.dir_metrics / 'gene_tree_stats.csv'

        fieldnames = [
            'gene_id', 'num_sequences', 'num_sites',
            'num_constant_sites', 'pct_constant',
            'num_parsimony_informative', 'pct_parsimony_informative',
            'num_distinct_patterns', 'substitution_model',
            'log_likelihood', 'AIC', 'AICc', 'BIC',
            'gamma_alpha', 'total_tree_length',
            'internal_branch_length', 'pct_internal_length',
            'cpu_time_s', 'wall_time_s',
        ]

        rows = []
        if self.dir_prank and self.dir_prank.is_dir():
            rows = parse_iqtree_directory(str(self.dir_prank))

        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)

        return csv_file

    def generate_alignment_quality_csv(self) -> Path:
        csv_file = self.dir_metrics / 'alignment_quality.csv'

        fieldnames = [
            'gene_id', 'num_sequences', 'alignment_length_bp',
            'mean_gap_pct', 'max_gap_pct', 'pct_missing_data',
            'sites_before_pxclsq', 'sites_after_pxclsq', 'pct_sites_retained',
        ]

        rows = []
        if self.dir_prank and self.dir_prank.is_dir():
            rows = self._compute_alignment_quality(self.dir_prank)

        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)

        return csv_file

    def _compute_alignment_quality(self, prank_dir: Path) -> List[Dict[str, Any]]:
        rows = []
        cln_files = sorted(prank_dir.glob('*-cln'))

        for cln_path in cln_files:
            gene_id  = cln_path.name.replace('-cln', '')
            seqs     = self._read_fasta_seqs(cln_path)
            if not seqs:
                continue

            lengths = [len(s) for s in seqs.values()]
            aln_len = lengths[0] if lengths else 0
            n_seqs  = len(seqs)

            gap_pcts = []
            total_gap_chars = 0
            for seq in seqs.values():
                gaps = seq.count('-') + seq.count('?')
                total_gap_chars += gaps
                gap_pcts.append(gaps / aln_len * 100 if aln_len else 0)

            mean_gap = sum(gap_pcts) / len(gap_pcts) if gap_pcts else 0
            max_gap  = max(gap_pcts) if gap_pcts else 0
            pct_missing = total_gap_chars / (n_seqs * aln_len) * 100 if (n_seqs * aln_len) else 0

            aln_stem     = cln_path.name.replace('-cln', '')
            aln_path     = prank_dir / aln_stem
            sites_before = None
            pct_retained = None
            if aln_path.exists():
                aln_seqs    = self._read_fasta_seqs(aln_path)
                if aln_seqs:
                    raw_lengths  = [len(s) for s in aln_seqs.values()]
                    sites_before = raw_lengths[0] if raw_lengths else None
                    if sites_before and aln_len:
                        pct_retained = round(aln_len / sites_before * 100, 2)

            rows.append({
                'gene_id'             : gene_id,
                'num_sequences'       : n_seqs,
                'alignment_length_bp' : aln_len,
                'mean_gap_pct'        : round(mean_gap, 2),
                'max_gap_pct'         : round(max_gap, 2),
                'pct_missing_data'    : round(pct_missing, 2),
                'sites_before_pxclsq' : sites_before,
                'sites_after_pxclsq'  : aln_len,
                'pct_sites_retained'  : pct_retained,
            })

        return rows

    @staticmethod
    def _read_fasta_seqs(filepath: Path) -> Dict[str, str]:
        seqs: Dict[str, str] = {}
        current_id = None
        current_seq: List[str] = []
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.rstrip()
                    if line.startswith('>'):
                        if current_id is not None:
                            seqs[current_id] = ''.join(current_seq)
                        current_id  = line[1:].split()[0]
                        current_seq = []
                    else:
                        current_seq.append(line)
            if current_id is not None:
                seqs[current_id] = ''.join(current_seq)
        except (OSError, IOError):
            pass
        return seqs

    def generate_species_coverage_csv(self) -> Path:
        csv_file = self.dir_metrics / 'species_coverage.csv'

        fieldnames = [
            'species', 'gene_trees_present', 'total_gene_trees',
            'coverage_pct',
        ]

        rows = []
        gene_trees_dir = self._find_gene_trees_dir()
        if gene_trees_dir and gene_trees_dir.is_dir():
            rows = self._compute_species_coverage(gene_trees_dir)

        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)

        return csv_file

    def _compute_species_coverage(self, gene_trees_dir: Path) -> List[Dict[str, Any]]:
        tree_files  = sorted(gene_trees_dir.glob('*.tre'))
        total_trees = len(tree_files)
        if total_trees == 0:
            return []

        species_tree_sets: Dict[str, set] = {}

        for tree_path in tree_files:
            try:
                t = ete3.Tree(str(tree_path))
                seen_in_tree: set = set()
                for leaf in t.get_leaves():
                    species = leaf.name.split('@')[0] if '@' in leaf.name else leaf.name
                    seen_in_tree.add(species)
                for species in seen_in_tree:
                    if species not in species_tree_sets:
                        species_tree_sets[species] = set()
                    species_tree_sets[species].add(tree_path.name)
            except Exception:
                continue

        rows = []
        for species, tree_set in sorted(species_tree_sets.items()):
            count = len(tree_set)
            rows.append({
                'species'           : species,
                'gene_trees_present': count,
                'total_gene_trees'  : total_trees,
                'coverage_pct'      : round(count / total_trees * 100, 2),
            })

        return rows

    def generate_astral_support_csv(self) -> Path:
        branch_file  = self.dir_metrics / 'astral_branch_support.csv'
        summary_file = self.dir_metrics / 'astral_summary.csv'

        branch_fields = [
            'branch_topology', 'num_gene_trees',
            'mean_local_posterior', 'median_local_posterior',
            'min_local_posterior', 'max_local_posterior',
        ]
        summary_fields = [
            'num_branches',
            'mean_support', 'median_support', 'min_support', 'max_support',
            'branches_gt_0.95', 'branches_lt_0.50', 'pct_branches_gt_0.95',
        ]

        verbose_csv = self._find_verbose_csv()
        parsed      = parse_verbose_csv(verbose_csv) if verbose_csv else {}
        branches    = parsed.get('branches', [])

        with open(branch_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=branch_fields, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(branches)

        summary_row = {
            'num_branches'        : parsed.get('num_branches', ''),
            'mean_support'        : parsed.get('mean_support', ''),
            'median_support'      : parsed.get('median_support', ''),
            'min_support'         : parsed.get('min_support', ''),
            'max_support'         : parsed.get('max_support', ''),
            'branches_gt_0.95'    : parsed.get('branches_gt_0_95', ''),
            'branches_lt_0.50'    : parsed.get('branches_lt_0_50', ''),
            'pct_branches_gt_0.95': parsed.get('pct_branches_gt_0_95', ''),
        }
        with open(summary_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=summary_fields)
            writer.writeheader()
            writer.writerow(summary_row)

        return branch_file

    def _parse_time(self, time_str: str) -> float:
        if not time_str:
            return 0.0
        try:
            return float(str(time_str).replace(' s', '').strip())
        except (ValueError, AttributeError):
            return 0.0

    def _resolve_bes_tree_count(self, astral_data: dict) -> Any:
        val = astral_data.get('bes_tree_count')
        if val is not None:
            return val
        if self.dir_results:
            bes_dir = self.dir_results / 'species_trees' / 'BES_Trees'
            if bes_dir.is_dir():
                return len(list(bes_dir.iterdir()))
        return ''

    def _find_verbose_csv(self) -> Optional[str]:
        if self.dir_results:
            candidate = self.dir_results / 'species_trees' / 'SpeciesTree.verbose.csv'
            if candidate.exists():
                return str(candidate)
        return None

    def _find_gene_trees_dir(self) -> Optional[Path]:
        if self.dir_results:
            candidate = self.dir_results / 'gene_trees'
            if candidate.is_dir():
                return candidate
        return None
