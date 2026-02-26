import csv
import json
import shutil
import subprocess
import sys

import ete3

from collections import defaultdict
from datetime    import datetime
from glob        import glob
from pathlib     import Path

from core.utils.stage_config      import PipelineConfig
from core.stages.blast            import Blast
from core.stages.busco            import Busco
from core.stages.genecluster      import GeneCluster
from core.stages.hcluster         import HCluster
from core.stages.mafft            import Mafft
from core.stages.mcl              import MCL
from core.stages.prank            import Prank
from core.stages.prune            import Prune
from core.stages.tree             import Tree
from core.stages.final            import Astral
from core.utils.clutter           import cleanup_files, final_cleanup
from core.utils.contig            import Contig
from core.utils.constants         import FASTA_EXTENSIONS
from core.utils.decor             import bilge_crew, tree_crew, convert_paths, _truncate_for_terminal
from core.utils.printout          import PrintOut
from core.utils.sublogger         import get_subprocess_dir

class DataHub:
    def __init__(self, args):
        # Load Configuration
        self.config = PipelineConfig.from_args(args)

        # General Parameters
        self.iterations                 = self.config.iterations
        self.current_iter               = self.config.current_iter
        self.threads                    = self.config.threads
        self.hcluster_tree_override     = self.config.hcluster_tree_override

        # Config Sections
        self.blast_config               = self.config.blast
        self.mcl_config                 = self.config.mcl
        self.mafft_config               = self.config.mafft
        self.tree_config                = self.config.tree
        self.prune_config               = self.config.prune
        self.prank_config               = self.config.prank
        self.super_config               = self.config.super
        self.hcluster_config            = self.config.hcluster
        self.busco_config               = self.config.busco
        
        # BLAST
        self.blast_evalue               = self.blast_config.evalue
        self.blast_max_targets          = self.blast_config.max_targets
        
        # MCL
        self.hit_frac_cutoff            = self.mcl_config.hit_frac_cutoff
        self.minimum_taxa               = self.mcl_config.minimum_taxa
        self.mcl_inflation              = self.mcl_config.inflation
        self.perfect_identity           = self.mcl_config.perfect_identity
        self.coverage_threshold         = self.mcl_config.coverage_threshold
        self.min_seq_length             = self.mcl_config.min_seq_length
        
        # MAFFT
        self.mafft_maxiter              = self.mafft_config.maxiter
        self.pxclsq_threshold           = self.mafft_config.pxclsq_threshold
        self.thread_divisor             = self.mafft_config.thread_divisor
        
        # Tree
        self.relative_cutoff            = self.tree_config.relative_cutoff
        self.absolute_cutoff            = self.tree_config.absolute_cutoff
        self.branch_cutoff              = self.tree_config.branch_cutoff
        self.mask_paralogs              = self.tree_config.mask_paralogs
        self.outlier_ratio              = self.tree_config.outlier_ratio
        self.max_trim_iterations        = self.tree_config.max_trim_iterations
        self.min_subtree_taxa           = self.tree_config.min_subtree_taxa
        self.min_tree_leaves            = self.tree_config.min_leaves
        self.tree_start_from_prev       = self.tree_config.start_from_prev

        # Prune
        self.orthocutoff                = self.prune_config.orthocutoff
        self.prune_relative_cutoff      = self.prune_config.relative_cutoff
        self.prune_absolute_cutoff      = self.prune_config.absolute_cutoff
        self.prune_outlier_ratio        = self.prune_config.outlier_ratio
        self.prune_max_trim_iterations  = self.prune_config.max_trim_iterations
        self.prune_min_tree_leaves      = self.prune_config.min_tree_leaves
        
        # PRANK
        self.seqtype                    = self.prank_config.seqtype
        self.prank_pxclsq_threshold     = self.prank_config.pxclsq_threshold
        self.bootstrap_replicates       = self.prank_config.bootstrap
        
        # Super 
        self.super_bootstrap            = self.super_config.bootstrap
        self.bes_support                = self.super_config.bes_support
        self.output_super_matrix        = self.super_config.output_matrix
        
        # HCluster
        self.hcluster_enabled           = self.hcluster_config.enabled
        self.hcluster_id                = self.hcluster_config.id
        self.hcluster_iddef             = self.hcluster_config.iddef
        self.hcluster_tool              = self.hcluster_config.tool
        self.hcluster_use_busco         = self.hcluster_config.use_busco
        
        # BUSCO
        self.busco_evalue               = self.busco_config.evalue
        self.busco_max_targets          = self.busco_config.max_targets
        self.busco_coverage_threshold   = self.busco_config.coverage_threshold

        # Terminal
        self.hc                         = self.config.hc
        self.bc                         = self.config.bc
        self.log                        = self.config.log
        self.subprocess_logs            = self.config.subprocess_logs

        # Developer Flags
        self.clutter                    = self.config.clutter

        # Directories and Files
        self.dir_base                   = self.config.dir_base
        self.output_dir                 = self.config.output_dir
        self.dir_treeforge              = self.config.dir_treeforge
        self.files_fasta                = self.config.files_fasta

        self.run_timestamp              = datetime.now().strftime('%Y%m%d_%H%M%S')
        dir_logs                        = self.dir_treeforge / '04_metrics' / 'logs'
        log_file_path                   = dir_logs / self.run_timestamp / f"treeforge_{self.run_timestamp}.log"
        log_file_path.parent.mkdir(parents=True, exist_ok=True)
        self.log_file_path              = str(log_file_path)
        self.printClass                 = PrintOut(self.log, args.highlight_color, args.background_color, log_file=self.log_file_path)
        if getattr(args, 'nocolor', False):
            self.printClass.set_nocolor(True)
        self.printout                   = self.printClass.printout
        if self.subprocess_logs:
            self.subprocess_dir         = get_subprocess_dir(dir_logs, self.run_timestamp)
        else:
            self.subprocess_dir         = None

        self._prior_check()

        # master dictionary
        self.master_dict                = {}
        self._dict_setup()
        
        # Clutter management
        self.files_to_preserve          = set()
        
        self._make_dirs()

    def _detect_previous_modes(self):
        """Detect which modes have previously run in directory."""
        modes_run = set()
        
        if not self.dir_treeforge.exists():
            return modes_run
        
        if (self.dir_treeforge / '02_analysis' / 'mcl').exists():
            modes_run.add('normal')
        elif (self.dir_treeforge / '02_analysis' / 'prune').exists():
            modes_run.add('normal')
        
        hcluster_dir = self.dir_treeforge / '03_results' / 'hcluster'
        if hcluster_dir.exists():
            if any(hcluster_dir.glob('clusters/*')) or any(hcluster_dir.glob('vsearch/*')):
                modes_run.add('hcluster')
        
        return modes_run
    
    def _prior_check(self):
        if not Path(self.dir_treeforge).exists():
            return
        
        if not any(Path(self.dir_treeforge).iterdir()):
            return
        
        previous_modes = self._detect_previous_modes()
        current_mode   = 'hcluster' if self.hcluster_enabled else 'normal'
        
        if current_mode in previous_modes:
            mode_name = 'Hierarchical clustering' if current_mode == 'hcluster' else 'Normal'
            self.printout('error', f'{mode_name} mode already run in directory')
            self.printout('error', _truncate_for_terminal(self.dir_treeforge.as_posix()))
            self.printout('info', 'Use a different output directory')
            sys.exit(1)
        
        return

    def _dict_setup(self):
        """Setup master dictionary with directory structure."""
        
        self.dir_logs    = self.dir_treeforge / '04_metrics' / 'logs'
        self.dir_metrics = self.dir_treeforge / '04_metrics'

        self.master_dict = {
            # Input processing
            'base'      : {'fai': self.dir_treeforge / '01_input' / 'fai',
                           'temp_fasta': self.dir_treeforge / '01_input' / 'renamed_fasta'
                           },
            
            # Analysis stages
            'blast'     : {
                'dir'  : self.dir_treeforge / '02_analysis' / 'blast',
                'files': {
                    'concatenated': self.dir_treeforge / '02_analysis' / 'blast' / 'concatenated.fasta',
                    'raw_blast'   : self.dir_treeforge / '02_analysis' / 'blast' / 'raw.blast'
                }
            },
            'mcl'       : {
                'dir': self.dir_treeforge / '02_analysis' / 'mcl',
            },
            'prune'     : {'dir': self.dir_treeforge / '02_analysis' / 'prune',
                           'ortho1to1': self.dir_treeforge / '02_analysis' / 'prune' / 'orthologs'},
            'prank': {'dir': self.dir_treeforge / '02_analysis' / 'prank'},
            'super': {'dir': self.dir_treeforge / '02_analysis' / 'super'},
            
            # Results
            'results': {
                'species_trees': {
                    'dir'  : self.dir_treeforge / '03_results' / 'species_trees',
                    'files': {
                        'coalescent' : self.dir_treeforge / '03_results' / 'species_trees' / 'SpeciesTree.coalescent.tre',
                        'molecular'  : self.dir_treeforge / '03_results' / 'species_trees' / 'SpeciesTree.molecular.tre',
                        'supermatrix': self.dir_treeforge / '03_results' / 'species_trees' / 'SuperMatrix.tre'
                    }
                },
                'gene_trees': {
                    'dir'  : self.dir_treeforge / '03_results' / 'gene_trees',
                    'files': {
                        'concatenated'  : self.dir_treeforge / '03_results' / 'gene_trees' / 'concatenated.tre',
                        'individual_dir': self.dir_treeforge / '03_results' / 'gene_trees' / 'individual'
                    }
                },
                'hcluster': {
                    'dir': self.dir_treeforge / '03_results' / 'hcluster'
                }
            },
            
            # Metrics
            'metrics': {
                'dir'  : self.dir_treeforge / '04_metrics'
            },
            
            # Logs
            'logs': {'dir': self.dir_treeforge / '04_metrics' / 'logs'},
            
            'gene_trees'      : {'dir': self.dir_treeforge / '03_results' / 'gene_trees'},
            'hcluster_enabled': {'dir': self.dir_treeforge / '03_results' / 'hcluster'},
            
            # HCluster-specific blast
            'hcluster_blast'  : {
                'dir'  : self.dir_treeforge / '03_results' / 'hcluster' / 'blast',
                'files': {
                    'concatenated': self.dir_treeforge / '03_results' / 'hcluster' / 'blast' / 'concatenated.fasta',
                    'raw_blast'   : self.dir_treeforge / '03_results' / 'hcluster' / 'blast' / 'raw.blast'
                }
            }
        }

    def _iteration_dirs(self, iteration):
        """Make sure that the directory for an iteration exists."""
        iter_key = f'iter_{iteration}'
        if iter_key not in self.master_dict:
            iter_base = self.dir_treeforge / '02_analysis' / 'iterations' / f'iter_{iteration}'
            self.master_dict[iter_key] = {
                'mafft': {'dir': iter_base / 'mafft'},
                'tree': {
                    'dir'      : iter_base / 'tree',
                    'trimmed'  : iter_base / 'tree' / 'trimmed'
                }
            }
            self.master_dict[iter_key]['mafft']['dir'].mkdir(parents=True, exist_ok=True)
            self.master_dict[iter_key]['tree']['dir'].mkdir(parents=True, exist_ok=True)
            self.master_dict[iter_key]['tree']['trimmed'].mkdir(parents=True, exist_ok=True)


    def _make_dirs(self):
        """Make only the directories we need"""
        def create_dirs(d: dict):
            for key, value in d.items():
                if isinstance(value, dict):
                    create_dirs(value)
                elif isinstance(value, Path) and key == 'dir':
                    value.mkdir(parents=True, exist_ok=True)
        
        if 'base' in self.master_dict and 'temp_fasta' in self.master_dict['base']:
            self.master_dict['base']['temp_fasta'].mkdir(parents=True, exist_ok=True)
        
        if 'metrics' in self.master_dict:
            create_dirs(self.master_dict['metrics'])
        if 'logs' in self.master_dict:
            create_dirs(self.master_dict['logs'])
        
        if self.hcluster_enabled and self.hcluster_tree_override:
            if 'hcluster_enabled' in self.master_dict:
                create_dirs(self.master_dict['hcluster_enabled'])
            if 'hcluster_blast' in self.master_dict:
                create_dirs(self.master_dict['hcluster_blast'])
            return
        
        if self.hcluster_enabled:
            if 'hcluster_blast' in self.master_dict:
                create_dirs(self.master_dict['hcluster_blast'])
            if 'hcluster_enabled' in self.master_dict:
                create_dirs(self.master_dict['hcluster_enabled'])
            if 'results' in self.master_dict and 'hcluster' in self.master_dict['results']:
                create_dirs(self.master_dict['results']['hcluster'])
            return
        
        if 'results' in self.master_dict:
            if 'species_trees' in self.master_dict['results']:
                create_dirs(self.master_dict['results']['species_trees'])
            if 'gene_trees' in self.master_dict['results']:
                create_dirs(self.master_dict['results']['gene_trees'])
        
        if 'blast' in self.master_dict:
            create_dirs(self.master_dict['blast'])
        
        if 'mcl' in self.master_dict:
            create_dirs(self.master_dict['mcl'])
        
        if 'prune' in self.master_dict:
            create_dirs(self.master_dict['prune'])
        
        if 'prank' in self.master_dict:
            create_dirs(self.master_dict['prank'])
        
        if 'super' in self.master_dict:
            create_dirs(self.master_dict['super'])

    def get_result_path(self, category, file_type=None):
        if file_type:
            return self.master_dict['results'][category]['files'][file_type]
        return self.master_dict['results'][category]['dir']

    def get_analysis_path(self, stage, file_type=None):
        if file_type:
            return self.master_dict[stage]['files'][file_type]
        return self.master_dict[stage]['dir']

    def _cleanup_intermediate_files(self, stage, iter_override=None):
        cleanup_files(
            self.master_dict,
            iter_override if iter_override is not None else self.current_iter,
            stage,
            self.clutter,
            self.printout,
            hit_frac_cutoff = self.hit_frac_cutoff,
            mcl_inflation   = self.mcl_inflation,
            cleanup_type    = 'intermediate'
        )

    def _cleanup_preserved_files(self, stage):
        cleanup_files(
            self.master_dict,
            self.current_iter,
            stage,
            self.clutter,
            self.printout,
            hit_frac_cutoff = self.hit_frac_cutoff,
            mcl_inflation   = self.mcl_inflation,
            cleanup_type    = 'preserved'
        )

    def _clutter_check(self):
        if self.clutter:
            self.printout('metric', 'Performing final cleanup')
            
            if self.dir_treeforge.exists():
                for item in self.dir_treeforge.iterdir():
                    if item.is_dir() and item.name not in {'03_results', '04_metrics'}:
                        try:
                            if not any(item.rglob('*')):
                                item.rmdir()
                                self.printout('metric', f'Removed empty directory: {item.name}')
                        except OSError:
                            pass
                
                self.printout('metric', f'Final cleanup complete')
            else:
                self.printout('error', 'TreeForge directory not found')

    def _move_fai(self):
        """Move the fai files to the TreeForge/01_input/fai directory.
        Makes things look cleaner."""
        files_fai = list(self.dir_base.glob('*.fai'))
        if files_fai:
            dir_fai = self.master_dict['base']['fai']
            dir_fai.mkdir(parents=True, exist_ok=True)
            for fai in files_fai:
                fai_dest = dir_fai / fai.name
                shutil.move(str(fai), str(fai_dest))
        
        # please make these logs optional, phyx family
        phyx_log_parent = self.dir_treeforge.parent / 'phyx.logfile'
        if phyx_log_parent.exists():
            phyx_log_parent.unlink()
        
        phyx_log_base = self.dir_base / 'phyx.logfile'
        if phyx_log_base.exists():
            phyx_log_base.unlink()

    def _finish_pipeline(self, gene_tree_count=None):
        """Finish the pipeline and cleanup files."""
        self.printout('title', 'TreeForge Complete')
        self._move_fai()
        protected_files = {
            'SpeciesTree.coalescent.tre',
            'SpeciesTree.molecular.tre',
            'SuperMatrix.tre',
            'run_metrics.json',
            'iteration_flow.csv',
            'pipeline_summary.csv',
            'clustering_summary.csv',
            'mafft_summary.csv',
            'tree_processing.csv',
            'prune_metrics.csv',
            'final_alignment.csv',
            'hcluster_metrics.csv',
            'hcluster_metrics.json',
            'hcluster_guide.tre'
        }
        protected_dirs = {
            '03_results',
            '04_metrics'
        }
        
        final_cleanup(
            root_dir        = self.dir_treeforge,
            logs_dir        = self.dir_logs,
            printout        = self.printout,
            protected_files = protected_files,
            protected_dirs  = protected_dirs,
            protected_exts  = {'.yaml'}
        ) if self.clutter else None
        
        from core.utils.summary import MetricsSummary

        summary = MetricsSummary(
            master_dict         = self.master_dict,
            dir_metrics         = self.dir_metrics,
            hcluster_enabled    = self.hcluster_enabled,
            output_super_matrix = self.output_super_matrix,
            run_timestamp       = self.run_timestamp,
            dir_prank           = self.dir_treeforge / '02_analysis' / 'prank',
            dir_results         = self.dir_treeforge / '03_results',
        )

        reports = summary.generate_all_reports()
        
        if gene_tree_count is None:
            gene_tree_count = len(glob(str(self.get_result_path('gene_trees')) + '/*.tre'))
        
        out_dict = {'Species Tree'   : '',
                    '   Molecular'   : str(self.get_result_path('species_trees', 'molecular')),
                    '   Coalescent'  : str(self.get_result_path('species_trees', 'coalescent')),
                    'Gene Trees'     : '',
                    '   Directory'   : str(self.get_result_path('gene_trees')),
                    '   Count'       : gene_tree_count,
                    'Metrics'        : '',
                    '   Reports'     : str(self.dir_metrics),
                    '   JSON'        : str(reports['json']),
                    '   Logs'        : str(self.dir_metrics / 'logs' / self.run_timestamp),
                    }
        if self.output_super_matrix:
            out_dict['SuperMatrix'] = str(self.get_result_path('species_trees', 'supermatrix'))
        if self.hcluster_enabled:
            out_dict['HCluster'] = str(self.get_result_path('hcluster'))
            if 'hcluster_metrics' in reports:
                out_dict['   HCluster Metrics'] = str(reports['hcluster_metrics'])
        self.printout('final', out_dict)
        
        self.printClass.close()

    def _finish_hcluster_pipeline(self, busco_result=None, genecluster_result=None, 
                                   guide_tree_result=None, hcluster_result=None):
        self.printout('title', 'TreeForge HCluster Complete')
        self._move_fai()
        
        metrics = {
            'mode'         : 'hierarchical_clustering',
            'run_timestamp': self.run_timestamp,
            'input_files'  : {
                'count': len(getattr(self, 'original_files_fasta', self.files_fasta)),
                'files': [str(f) for f in getattr(self, 'original_files_fasta', self.files_fasta)]
            },
            'directories': {
                'base'    : str(self.dir_base),
                'output'  : str(self.dir_treeforge),
                'hcluster': str(self.get_result_path('hcluster'))
            }
        }
        
        if busco_result:
            metrics['busco_extraction']        = convert_paths(busco_result) if isinstance(busco_result, dict) else busco_result
        if genecluster_result:
            metrics['gene_clustering']         = convert_paths(genecluster_result) if isinstance(genecluster_result, dict) else genecluster_result
        if guide_tree_result:
            metrics['guide_tree']              = convert_paths(guide_tree_result) if isinstance(guide_tree_result, dict) else guide_tree_result
        if hcluster_result:
            metrics['hierarchical_clustering'] = convert_paths(hcluster_result) if isinstance(hcluster_result, dict) else hcluster_result
        
        metrics['configuration'] = {
            'hcluster_id'             : self.hcluster_id,
            'hcluster_iddef'          : self.hcluster_iddef,
            'hcluster_tool'           : self.hcluster_tool,
            'hcluster_use_busco'      : self.hcluster_use_busco,
            'threads'                 : self.threads,
            'busco_evalue'            : self.busco_evalue,
            'busco_max_targets'       : self.busco_max_targets,
            'busco_coverage_threshold': self.busco_coverage_threshold,
            'blast_evalue'            : self.blast_evalue,
            'blast_max_targets'       : self.blast_max_targets,
            'minimum_taxa'            : self.minimum_taxa
        }
        
        metrics['paths'] = convert_paths({
            'hcluster_dir'  : self.master_dict.get('hcluster_enabled', {}).get('dir'),
            'hcluster_blast': self.master_dict.get('hcluster_blast', {}),
            'base'          : self.master_dict.get('base', {}),
            'logs'          : self.master_dict.get('logs', {})
        })
        
        json_file = self.dir_metrics / 'hcluster_metrics.json'
        with open(json_file, 'w') as f:
            json.dump(metrics, f, indent=2)
        
        csv_file = None
        try:
            from core.utils.summary import MetricsSummary
            summary = MetricsSummary(
                master_dict         = self.master_dict,
                dir_metrics         = self.dir_metrics,
                hcluster_enabled    = True,
                output_super_matrix = False,
                run_timestamp       = self.run_timestamp,
                dir_prank           = self.dir_treeforge / '02_analysis' / 'prank',
                dir_results         = self.dir_treeforge / '03_results',
            )
            csv_file = summary.generate_hcluster_csv(hcluster_metrics=metrics)
        except Exception as e:
            self.printout('warning', f'Could not generate HCluster CSV: {str(e)}')
        
        cluster_count  = hcluster_result.get('cluster_count', 0) if hcluster_result else 0
        clustered_seqs = hcluster_result.get('clustered_sequences', 0) if hcluster_result else 0
        
        out_dict = {
            'HCluster Results'  : '',
            '   Clusters Dir'   : str(self.get_result_path('hcluster')),
            '   Clstr Count'    : cluster_count,
            '   Clstr Sequences': clustered_seqs,
            'Metrics'           : '',
            '   Reports'        : str(self.dir_metrics),
            '   JSON'           : str(json_file),
            '   Logs'           : str(self.dir_metrics / 'logs' / self.run_timestamp)
        }
        
        if csv_file:
            out_dict['   CSV'] = str(csv_file)
        
        if guide_tree_result and 'guide_tree_file' in guide_tree_result:
            out_dict['   Guide Tree'] = guide_tree_result['guide_tree_file']
        
        self.printout('final', out_dict)
        
        if self.clutter:
            protected_dirs  = {'03_results', '04_metrics'}
            protected_files = {
                'hcluster_metrics.json',
                'hcluster_metrics.csv',
                'hcluster_guide.tre',
                'SpeciesTree.coalescent.tre',
                'SpeciesTree.molecular.tre',
                'SuperMatrix.tre',
                'run_metrics.json'
            }
            final_cleanup(
                root_dir        = self.dir_treeforge,
                logs_dir        = self.dir_logs,
                printout        = self.printout,
                protected_files = protected_files,
                protected_dirs  = protected_dirs,
                protected_exts  = {'.yaml'}
            )
        
        phyx_log_parent = self.dir_treeforge.parent / 'phyx.logfile'
        if phyx_log_parent.exists():
            phyx_log_parent.unlink()
        
        phyx_log_base = self.dir_base / 'phyx.logfile'
        if phyx_log_base.exists():
            phyx_log_base.unlink()
        
        self.printClass.close()
        
        return hcluster_result

    def _generate_species_map(self, target_files):
        self.species_map = {} # ID to FilePath
        self.stem_to_id  = {} # Original stem to ID
        self.id_to_stem  = {} # ID to Original stem

        prefix_to_files  = defaultdict(list)
        for f in target_files:
            prefix = f.name.split('.')[0]
            prefix_to_files[prefix].append(f)
            
        for prefix, files in prefix_to_files.items():
            if len(files) == 1:
                f                            = files[0]
                species_id                   = prefix
                self.species_map[species_id] = str(f)
                self.stem_to_id[f.stem]      = species_id
                self.id_to_stem[species_id]  = f.stem
            else:
                self.printout('info', f"Multiple prefix '{prefix}', assigning unique species IDs.")
                
                for f in files:
                    unique_id          = f.stem.replace('.', '_')
                    original_unique_id = unique_id
                    counter            = 1
                    while unique_id in self.species_map:
                        unique_id = f"{original_unique_id}_{counter}"
                        counter  += 1
                        
                    self.species_map[unique_id] = str(f)
                    self.stem_to_id[f.stem]     = unique_id
                    self.id_to_stem[unique_id]  = f.stem
                    self.printout('info', f"  Mapped '{f.name}' -> ID '{unique_id}'")
        
    def _contig_renamer(self):
        """Rename the contigs/transcripts/scaffolds.
        This should take care of any naming schemes or weird characters.
        Prevents treeforge and the stuff it wraps from crashing"""
        
        temp_fasta_dir = self.master_dict['base']['temp_fasta']
        temp_fasta_dir.mkdir(parents=True, exist_ok=True)
        
        contig_renamer   = Contig(self.dir_base, temp_fasta_dir)
        result           = contig_renamer.rename_contigs(self.files_fasta)
        self.renamed_map = contig_renamer.name_mapping
        
        if not self.renamed_map:
            self.printout('error', 'Failed to create gene-to-species mapping')
        
        return result
    
    def run(self):
        original_files_fasta   = self.files_fasta.copy() if self.hcluster_enabled else None
        self.busco_files_fasta = None

        self._generate_species_map(self.files_fasta)

        if self.hcluster_enabled and self.hcluster_tree_override:
            species_tree_file = Path(self.hcluster_tree_override)
            if not species_tree_file.exists():
                self.printout('error', f'Custom HCluster guide tree not found: {species_tree_file}')
                self.printClass.close()
                return {'status': 'no_guide_tree'}
            
            self.original_files_fasta     = original_files_fasta or self.files_fasta.copy()
            self.hcluster_guide_tree_file = str(species_tree_file)
            
            hcluster_result = self.hcluster()
            self._cleanup_preserved_files('hcluster')
            return self._finish_hcluster_pipeline(
                busco_result       = None,
                genecluster_result = None,
                guide_tree_result  = {'guide_tree_file': str(species_tree_file), 'source': 'user_provided'},
                hcluster_result    = hcluster_result
            )
        
        busco_result       = None
        genecluster_result = None
        
        if self.hcluster_enabled:
            busco_result = self.busco_extract()
            if busco_result and 'busco_files' in busco_result and busco_result['busco_files']:
                self.busco_files_fasta = [Path(f) for f in busco_result['busco_files']]
                self.files_fasta = self.busco_files_fasta
                self.printout('metric', [('busco_files_for_blast', len(self.files_fasta))])
            else:
                self.printout('error', 'No BUSCO files generated, using original files')
                self.busco_files_fasta = None
        
        self._contig_renamer()

        self.blast()
        
        if self.hcluster_enabled:
            genecluster_result = self.genecluster_for_hcluster()
            if genecluster_result and 'gene_files' in genecluster_result and genecluster_result['gene_files']:
                self.files_fasta = [Path(f) for f in genecluster_result['gene_files']]
                self.printout('metric', [('gene_cluster_files_for_pipeline', len(self.files_fasta))])
            else:
                self.printout('error', 'No gene cluster files generated, using original files')
            
            self._cleanup_intermediate_files('blast')
        else:
            mcl_result = self.mcl()
            self._cleanup_intermediate_files('blast')
            
            if mcl_result and mcl_result.get('mcl_count', 0) == 0:
                self.printout('error', 'No clusters produced by MCL')
                self.printout('error', 'Try adjusting the MCL parameters')
                # self._finish_pipeline(gene_tree_count=0)
                return
        
        if self.hcluster_enabled:
            guide_tree_result = self.hcluster_guide_tree()
            if not guide_tree_result:
                self.printout('error', 'Failed to generate hierarchical clustering guide tree')
                return
            self.original_files_fasta = original_files_fasta
            hcluster_result           = self.hcluster()
            self._cleanup_preserved_files('hcluster')
            return self._finish_hcluster_pipeline(
                busco_result       = busco_result,
                genecluster_result = genecluster_result,
                guide_tree_result  = guide_tree_result,
                hcluster_result    = hcluster_result
            )

        continue_iterations = True
        
        for i in range(self.iterations):
            if not continue_iterations:
                self.printout('info', f'Skipping remaining iterations.')
                break
            
            self.mafft()
            
            if (f'iter_{self.current_iter}' in self.master_dict and 
                'mafft' in self.master_dict[f'iter_{self.current_iter}'] and
                'mafft' in self.master_dict[f'iter_{self.current_iter}']['mafft'] and
                self.master_dict[f'iter_{self.current_iter}']['mafft']['mafft'].get('status') == 'no_files_to_process'):
                self.printout('info', f'No files in iteration {int(self.current_iter) + 1}')
                continue_iterations = False
                continue
            
            tree_result = self.tree()
            self._cleanup_intermediate_files('mafft')
            self._cleanup_intermediate_files('mcl')
            
            if tree_result and 'write_tree' in tree_result:
                modified_count = tree_result['write_tree'].get('modified_files_count', 0)
                # self.printout('info', f'Iteration {self.current_iter}: {modified_count} files modified')
                
                if modified_count == 0:
                    self.printout('info', f'No files modified in iteration {int(self.current_iter) + 1}')
                    continue_iterations = False
            
            if self.current_iter > 0:
                prev_iter = self.current_iter - 1
                self._cleanup_intermediate_files('mafft', iter_override=prev_iter)
                self._cleanup_intermediate_files('tree', iter_override=prev_iter)
            
            self.current_iter += 1 if self.current_iter < self.iterations - 1 else 0
            
        self.prune()
        for i in range(self.iterations):
            self.current_iter = i
            if f'iter_{i}' in self.master_dict:
                self._cleanup_intermediate_files('tree')
                self._cleanup_preserved_files('prune')
            
        self.prank()
        self._cleanup_intermediate_files('prank')
        self._cleanup_preserved_files('prank')
        self._cleanup_preserved_files('mcl')
        self.astral()
        self._cleanup_preserved_files('astral')
        self._finish_pipeline()

    @bilge_crew()
    def blast(self):
        self.printout('subtitle', 'BLAST')
        
        temp_fasta_dir = self.master_dict['base']['temp_fasta']
        all_files = [Path(f) for f in glob(str(temp_fasta_dir) + '/*') if f.endswith(tuple(['.fasta', '.fa', '.fas', '.fna']))]
        
        if self.hcluster_enabled:
            self.files_fasta_renamed = [f for f in all_files if '_busco_renamed.fa' in f.name]
            blast_dict = self.master_dict['hcluster_blast']
        else:
            self.files_fasta_renamed = [f for f in all_files if '_renamed.' in f.name and '_busco_' not in f.name]
            blast_dict = self.master_dict['blast']
        
        blast = Blast(self.dir_base,                                                # Input Base Directory
                      self.dir_treeforge,                                           # TreeForge Directory
                      self.files_fasta_renamed,                                     # List of renamed FASTA files
                      blast_dict['files']['concatenated'],                          # Concatenated FASTA file
                      blast_dict['files']['raw_blast'],                             # Raw BLAST file
                      self.threads,                                                 # Number of threads
                      self.log,                                                     # Log level
                      self.hc,                                                      # Highlight color
                      self.bc,                                                      # Background color
                      self.blast_evalue,                                            # BLAST E-value threshold
                      self.blast_max_targets,                                       # BLAST max target sequences
                      self.subprocess_dir,                                          # Subprocess output directory
                      shared_printClass=self.printClass)                            # Shared PrintOut instance
        blast.run()
        return blast.return_dict
    
    @bilge_crew()
    def mcl(self):
        self.printout('subtitle', 'MCL')
        self._iteration_dirs(self.current_iter)
        
        mcl = MCL(self.dir_base,                                                        # Input Base Directory
                  self.dir_treeforge,                                                   # TreeForge Directory
                  self.master_dict['mcl']['dir'],                                       # MCL Directory
                  self.master_dict[f'iter_{self.current_iter}']['mafft']['dir'],        # MAFFT Directory
                  self.master_dict['blast']['files']['raw_blast'],                      # Raw BLAST file
                  self.master_dict['blast']['files']['concatenated'],                   # Concatenated FASTA file
                  self.hit_frac_cutoff,                                                 # Hit fraction cutoff
                  self.minimum_taxa,                                                    # Minimum taxa for MCL clustering
                  self.threads,                                                         # Number of threads
                  self.log,                                                             # Log level
                  self.hc,                                                              # Highlight color
                  self.bc,                                                              # Background color
                  self.mcl_inflation,                                                   # MCL inflation parameter
                  self.perfect_identity,                                                # Perfect identity threshold
                  self.coverage_threshold,                                              # Coverage threshold for identical sequences
                  self.min_seq_length,                                                  # Minimum sequence length
                  self.subprocess_dir,                                                  # Subprocess output directory
                  shared_printClass=self.printClass)                                    # Shared PrintOut instance
        mcl.run()
        return mcl.return_dict

    @bilge_crew()
    def mafft(self):
        self.printout('subtitle', f'MAFFT - {self.current_iter + 1}')
        self._iteration_dirs(self.current_iter)
        prev_trimmed_tree_dir   = None
        use_prev_trees_as_start = False
        if self.current_iter > 0 and self.tree_start_from_prev:
            prev_iter_key = f'iter_{self.current_iter - 1}'
            if prev_iter_key in self.master_dict and 'tree' in self.master_dict[prev_iter_key]:
                prev_trimmed_tree_dir   = self.master_dict[prev_iter_key]['tree'].get('trimmed')
                use_prev_trees_as_start = prev_trimmed_tree_dir is not None
        
        mafft = Mafft(self.dir_base,                                                    # Input Base Directory
                      self.master_dict[f'iter_{self.current_iter}']['mafft']['dir'],    # MAFFT Directory
                      self.threads,                                                     # Number of threads
                      self.log,                                                         # Log level
                      self.hc,                                                          # Highlight color
                      self.bc,                                                          # Background color
                      self.mafft_maxiter,                                               # MAFFT max iterations
                      self.pxclsq_threshold,                                            # pxclsq probability threshold
                      self.thread_divisor,                                              # Thread division factor
                      prev_trimmed_tree_dir=str(prev_trimmed_tree_dir) \
                        if prev_trimmed_tree_dir else None,                             # Previous trimmed tree directory
                      use_prev_trees_as_start=use_prev_trees_as_start,                  # Use previous trees as starting trees
                      subprocess_dir=self.subprocess_dir,                               # Subprocess output directory
                      shared_printClass=self.printClass)                                # Shared PrintOut instance
        mafft.run()
        return mafft.return_dict

    @tree_crew()
    def tree(self):
        self.printout('subtitle', f'Tree - {self.current_iter + 1}')
        self._iteration_dirs(self.current_iter)
        tree = Tree(self.dir_base,                                                      # Input Base Directory
                    self.dir_treeforge,                                                 # TreeForge Directory
                    self.master_dict[f'iter_{self.current_iter}']['tree']['dir'],       # Tree Directory
                    self.master_dict[f'iter_{self.current_iter}']['mafft']['dir'],      # MAFFT Directory
                    self.master_dict[f'iter_{self.current_iter}']['tree']['trimmed'],   # Trimmed Tree Directory
                    self.master_dict['prune']['dir'],                                   # Prune Directory
                    self.iterations,                                                    # Number of iterations
                    self.current_iter,                                                  # Current iteration
                    self.minimum_taxa,                                                  # Minimum taxa for MCL clustering
                    self.master_dict['blast']['files']['concatenated'],                 # Concatenated FASTA file
                    self.threads,                                                       # Number of threads
                    self.log,                                                           # Log level
                    self.hc,                                                            # Highlight color
                    self.bc,                                                            # Background color
                    self.relative_cutoff,                                               # Relative cutoff for trimming tips
                    self.absolute_cutoff,                                               # Absolute cutoff for trimming tips
                    self.branch_cutoff,                                                 # Branch cutoff for cutting branches
                    self.mask_paralogs,                                                 # Mask paraphyletic tips
                    self.outlier_ratio,                                                 # Ratio threshold for outlier detection
                    self.max_trim_iterations,                                           # Maximum trimming iterations
                    self.min_subtree_taxa,                                              # Minimum taxa for valid subtrees
                    self.min_tree_leaves,                                               # Minimum leaves for valid tree
                    clutter=self.clutter,                                               # Pass clutter flag
                    shared_printClass=self.printClass)                                  # Shared PrintOut instance
        tree.run()
        return tree.return_dict

    @bilge_crew()
    def prune(self):
        self.printout('subtitle', 'Prune')
        self._iteration_dirs(self.current_iter)
        prune = Prune(self.dir_base,                                                    # Input Base Directory
                      self.master_dict['prune']['dir'],                                 # Prune Directory
                      self.master_dict[f'iter_{self.current_iter}']['tree']['trimmed'], # Trimmed Tree Directory
                      self.master_dict['prune']['ortho1to1'],                           # 1to1 Ortholog Directory
                      self.orthocutoff,                                                 # Ortholog Minimum Taxa Cutoff
                      self.threads,                                                     # Number of threads
                      self.log,                                                         # Log level
                      self.hc,                                                          # Highlight color
                      self.bc,                                                          # Background color
                      self.prune_relative_cutoff,                                       # Relative tip cutoff for pruning
                      self.prune_absolute_cutoff,                                       # Absolute tip cutoff for pruning
                      self.prune_outlier_ratio,                                         # Outlier ratio for pruning
                      self.prune_max_trim_iterations,                                   # Max trim iterations for pruning
                      self.prune_min_tree_leaves,                                       # Min tree leaves for pruning
                      shared_printClass=self.printClass)                                # Shared PrintOut instance
        prune.run()
        return prune.return_dict

    @bilge_crew()
    def prank(self):
        self.printout('subtitle', 'PRANK')
        
        prank = Prank(self.master_dict['prune']['ortho1to1'],                           # 1to1 Ortholog Directory
                      self.master_dict['prank']['dir'],                                 # PRANK Directory
                      '.tre',                                                           # Tree File Extension
                      self.seqtype,                                                     # Sequence type
                      self.master_dict['blast']['files']['concatenated'],               # Concatenated FASTA file
                      self.threads,                                                     # Number of threads
                      self.log,                                                         # Log level
                      self.hc,                                                          # Highlight color
                      self.bc,                                                          # Background color
                      self.prank_pxclsq_threshold,                                      # pxclsq probability threshold
                      self.bootstrap_replicates,                                        # IQ-TREE bootstrap replicates
                      self.subprocess_dir,                                              # Subprocess output directory
                      shared_printClass=self.printClass)                                # Shared PrintOut instance
        prank.run()
        return prank.return_dict

    @bilge_crew()
    def astral(self):
        self.printout('subtitle', 'Astral')
        
        id_to_stem = getattr(self, 'id_to_stem', {})
        astral = Astral(self.master_dict['prank']['dir'],                               # PRANK Directory
                        self.master_dict['super']['dir'],                               # Supertree Directory
                        self.dir_treeforge,                                             # TreeForge Directory
                        self.renamed_map,                                               # Renamed Map
                        self.output_super_matrix,                                       # Output supermatrix
                        self.super_bootstrap,                                           # Supermatrix bootstrap replicates
                        self.threads,                                                   # Number of threads
                        self.log,                                                       # Log level
                        self.hc,                                                        # Highlight color
                        self.bc,                                                        # Background color
                        self.get_result_path('species_trees', 'coalescent'),            # Coalescent tree path
                        self.get_result_path('species_trees', 'molecular'),             # Molecular tree path
                        self.get_result_path('species_trees', 'supermatrix'),           # Supermatrix tree path
                        self.get_result_path('gene_trees'),                             # Gene trees directory
                        id_to_stem,                                                     # Species ID -> original name
                        self.subprocess_dir,                                            # Subprocess output directory
                        shared_printClass=self.printClass)                              # Shared PrintOut instance
        astral.run()
        return astral.return_dict

    @bilge_crew()
    def busco_extract(self):
        self.printout('subtitle', 'BUSCO Extraction')
        
        id_to_stem = getattr(self, 'id_to_stem', {})
        
        busco = Busco(
                        self.dir_base,
                        self.dir_treeforge,
                        self.master_dict['hcluster_enabled']['dir'],
                        self.master_dict['hcluster_enabled']['dir'],
                        self.files_fasta,
                        self.threads,
                        self.log,
                        self.hc,
                        self.bc,
                        self.busco_evalue,
                        self.busco_max_targets,
                        self.busco_coverage_threshold,
                        self.hcluster_id,
                        self.hcluster_iddef,
                        None,
                        id_to_stem,
                        self.subprocess_dir,
                        shared_printClass=self.printClass
        )
        result = busco.run()
        return result
    
    @bilge_crew()
    def genecluster_for_hcluster(self):
        self.printout('subtitle', 'Gene Clustering for HCluster')
        
        genecluster_dir = self.master_dict['hcluster_enabled']['dir'] / 'genecluster'
        
        genecluster = GeneCluster(
                        self.dir_base,
                        self.dir_treeforge,
                        genecluster_dir,
                        self.master_dict['hcluster_blast']['files']['raw_blast'],
                        self.master_dict['hcluster_blast']['files']['concatenated'],
                        self.minimum_taxa,
                        self.threads,
                        self.log,
                        self.hc,
                        self.bc,
                        shared_printClass=self.printClass
        )
        result = genecluster.run()
        return result
    
    @bilge_crew()
    def hcluster_guide_tree(self):
        self.printout('subtitle', 'HCluster Guide Tree Generation')
        
        self.current_iter = 0
        self._iteration_dirs(self.current_iter)
        
        genecluster_dir = self.master_dict['hcluster_enabled']['dir'] / 'genecluster'
        self.printout('metric', 'Running MAFFT')
        
        mafft = Mafft(
            self.dir_base,
            genecluster_dir,
            self.threads,
            self.log,
            self.hc,
            self.bc,
            self.mafft_maxiter,
            self.pxclsq_threshold,
            self.thread_divisor,
            fast_mode=True,
            subprocess_dir=self.subprocess_dir,
            shared_printClass=self.printClass
        )
        mafft_result = mafft.run()
        
        if not mafft_result or 'mafft' not in mafft_result:
            self.printout('error', 'MAFFT failed for guide tree')
            return None
        
        if 'iqtree_files' not in mafft_result['mafft'] or not mafft_result['mafft']['iqtree_files']:
            self.printout('error', 'No IQ-TREE files generated from MAFFT')
            return None
        
        treefile_files = list(Path(genecluster_dir).glob('*.treefile'))
        
        if not treefile_files:
            self.printout('error', 'No treefile files found from MAFFT')
            return None
        
        self.printout('metric', [('gene_trees_found', len(treefile_files))])
        
        hcluster_tree_dir = self.dir_treeforge / '03_results' / 'hcluster_guide_tree'
        hcluster_tree_dir.mkdir(parents=True, exist_ok=True)
        
        concat_tree_file = hcluster_tree_dir / 'concat.tre'
        gene_trees_dir   = hcluster_tree_dir / 'gene_trees'
        gene_trees_dir.mkdir(parents=True, exist_ok=True)
        
        self.printout('metric', 'Concatenating gene trees')
        processed_count = 0
        skipped_count   = 0
        
        with open(concat_tree_file, 'w') as outfile:
            for treefile in sorted(treefile_files):
                try:
                    tree = ete3.Tree(str(treefile))

                    for node in tree.traverse():
                        if node.is_leaf():
                            lookup_key = node.name.replace('_', '@')
                            if hasattr(self, 'renamed_map') and lookup_key in self.renamed_map:
                                _, source_file = self.renamed_map[lookup_key]
                                
                                source_stem = Path(source_file).stem
                                if source_stem.endswith('_busco'):
                                    original_stem = source_stem[:-6]
                                else:
                                    original_stem = source_stem
                                
                                if hasattr(self, 'stem_to_id') and original_stem in self.stem_to_id:
                                    node.name = self.stem_to_id[original_stem]
                                else:
                                    node.name = source_file.split('.')[0]
                            else:
                                seq_id = node.name
                                if '@' in seq_id:
                                    species_name = seq_id.split('@')[0]
                                elif '_' in seq_id:
                                    species_name = seq_id.split('_')[0]
                                else:
                                    species_name = seq_id.split('.')[0]
                                node.name = species_name

                    tree_str = tree.write(format=1)
                    outfile.write(tree_str + '\n')
                    
                    gene_tree_file = gene_trees_dir / treefile.name
                    with open(gene_tree_file, 'w') as gene_out:
                        gene_out.write(tree_str)
                    
                    processed_count += 1
                except Exception as e:
                    self.printout('warning', f'Failed to process tree file {treefile}: {str(e)}')
                    skipped_count += 1
                    continue
        
        if processed_count == 0:
            self.printout('error', 'No valid tree files processed')
            return None
        
        # self.printout('metric', f'Successfully processed {processed_count} tree files')
        if skipped_count > 0:
            self.printout('warning', f'Skipped {skipped_count} tree files due to errors')
        
        astral_tree_file = hcluster_tree_dir / 'hcluster_guide.tre'
        self.printout('metric', 'Running ASTRAL for guide tree')
        
        cmd         = f'astral -i {concat_tree_file} -o {astral_tree_file}'
        return_code = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        if return_code.returncode != 0:
            self.printout('error', 'ASTRAL failed for hierarchical clustering guide tree')
            if return_code.stderr:
                self.printout('error', return_code.stderr.decode('utf-8'))
            return None
        
        if not astral_tree_file.exists():
            self.printout('error', 'ASTRAL tree file was not created')
            return None
        try:
            tree_ete = ete3.Tree(str(astral_tree_file))
            
            max_n = -1
            for node in tree_ete.traverse():
                if not node.is_leaf() and node.name and node.name.startswith('N'):
                    try:
                        n_num = int(node.name[1:])
                        max_n = max(max_n, n_num)
                    except ValueError:
                        pass
            
            counter       = max_n + 1
            labeled_count = 0
            for node in tree_ete.traverse("postorder"):
                if not node.is_leaf():
                    if not node.name or node.name.strip() == "":
                        node.name      = f"N{counter}"
                        counter       += 1
                        labeled_count += 1
            
            tree_ete.write(outfile=str(astral_tree_file), format=1)
            self.printout('metric', [('annotated_internal_nodes', labeled_count)])
        except Exception as e:
            self.printout('warning', f'Failed to annotate internal nodes in guide tree: {str(e)}')

        self.printout('metric', [('hcluster_guide_tree_created', str(astral_tree_file))])
        
        self.hcluster_guide_tree_file = str(astral_tree_file)
        
        return {'guide_tree_file': str(astral_tree_file), 'gene_trees_count': processed_count}
    
    @bilge_crew()
    def hcluster(self):
        """Run hierarchical"""
        self.printout('subtitle', 'Hierarchical Clustering')
        
        if hasattr(self, 'hcluster_guide_tree_file'):
            species_tree_file = self.hcluster_guide_tree_file
            self.printout('metric', 'Using hierarchical clustering guide tree')
        else:
            species_tree_file = str(self.get_result_path('species_trees', 'coalescent'))
            self.printout('metric', 'Using main coalescent tree')
        
        if self.hcluster_use_busco and hasattr(self, 'busco_files_fasta') and self.busco_files_fasta:
            temp_fasta_dir = self.master_dict['base']['temp_fasta']
            all_files      = [Path(f) for f in glob(str(temp_fasta_dir) + '/*') if f.endswith(tuple(['.fasta', '.fa', '.fas', '.fna']))]
            hcluster_files = [f for f in all_files if '_busco_renamed.fa' in f.name]
            
            if not hcluster_files:
                self.printout('warning', 'No renamed BUSCO files found')
                hcluster_files = self.busco_files_fasta
            
            # self.printout('metric', [('hcluster_using_busco_files', len(hcluster_files))])
        elif hasattr(self, 'original_files_fasta') and self.original_files_fasta:
            hcluster_files = [Path(f) for f in self.original_files_fasta]
            # self.printout('metric', [('hcluster_using_original_files', len(hcluster_files))])
        else:
            hcluster_files = [f for f in self.dir_base.iterdir() 
                            if f.suffix in FASTA_EXTENSIONS]
            # self.printout('metric', [('hcluster_using_input_dir_files', len(hcluster_files))])
        
        if not hcluster_files:
            self.printout('error', 'No input files found for hierarchical clustering')
            return {'status': 'no_files'}
        
        hcluster_map = {}
        if self.hcluster_use_busco and hasattr(self, 'busco_files_fasta') and self.busco_files_fasta:
            # self.printout('metric', 'Building species map for BUSCO files')
            for f in hcluster_files:
                f_path = Path(f) if not isinstance(f, Path) else f
                stem = f_path.stem
                
                if '_busco_renamed' in stem:
                    species_name = stem.replace('_busco_renamed', '')
                elif '_busco' in stem:
                    species_name = stem.replace('_busco', '')
                else:
                    if '_' in stem:
                        species_name = stem.split('_')[0]
                    else:
                        species_name = stem.split('.')[0]
                
                if species_name not in hcluster_map:
                    hcluster_map[species_name] = []
                hcluster_map[species_name].append(str(f_path))
            
            self.printout('metric', [('busco_species_mapped', len(hcluster_map))])
        elif hasattr(self, 'species_map'):
            for id, filepath in self.species_map.items():
                hcluster_map[id] = [filepath]
        else:
            for f in hcluster_files:
                f_path = Path(f) if not isinstance(f, Path) else f
                species_name = f_path.name.split('.')[0]
                if species_name not in hcluster_map:
                    hcluster_map[species_name] = []
                hcluster_map[species_name].append(str(f_path))

        id_to_stem = getattr(self, 'id_to_stem', {})
        
        hcluster = HCluster(
            self.dir_base,
            self.dir_treeforge,
            self.master_dict['hcluster_enabled']['dir'],
            hcluster_files,
            species_tree_file,
            self.threads,
            self.log,
            self.hc,
            self.bc,
            self.hcluster_id,
            self.hcluster_iddef,
            self.hcluster_tool,
            species_file_map=hcluster_map,
            id_to_stem=id_to_stem,
            shared_printClass=self.printClass
        )
        
        hcluster_result = hcluster.run()
        self._cleanup_intermediate_files('hcluster')
        return hcluster_result