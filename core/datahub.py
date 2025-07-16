import os
import json
import csv
import shutil
from pathlib             import Path
from core.utils.decor    import bilge_crew, tree_crew
from core.utils.contig   import Contig
from glob                import glob
from core.stages.blast   import Blast
from core.stages.mcl     import MCL
from core.stages.mafft   import Mafft
from core.stages.tree    import Tree
from core.stages.prune   import Prune
from core.stages.prank   import Prank
from core.stages.final   import Astral
from datetime            import datetime
from core.utils.printout import PrintOut
from core.utils.clutter  import cleanup_files, remove_empty_dirs, final_cleanup

class DataHub:
    def __init__(self, args):
        # TreeForge
        self.iterations                 = args.iter                         # Number of MAFFT/TREE iterations
        self.current_iter               = 0                                 # Current iteration
        self.threads                    = args.threads                      # Number of threads

        # Tool
        self.hit_frac_cutoff            = args.mcl_hit_frac_cutoff          # Hit fraction cutoff
        self.minimum_taxa               = args.mcl_minimum_taxa             # Minimum taxa for MCL clustering
        self.orthocutoff                = args.prune_orthocutoff            # Ortholog Minimum Taxa Cutoff
        self.seqtype                    = args.prank_seqtype                # Sequence type
        
        # BLAST
        self.blast_evalue               = args.blast_evalue                 # BLAST E-value threshold
        self.blast_max_targets          = args.blast_max_targets            # BLAST max target sequences
        
        # MCL
        self.mcl_inflation              = args.mcl_inflation                # MCL inflation parameter
        self.perfect_identity           = args.mcl_perfect_identity         # Perfect identity threshold
        self.coverage_threshold         = args.mcl_coverage_threshold       # Coverage threshold for identical sequences
        self.min_seq_length             = args.mcl_min_seq_length           # Minimum sequence length
        
        # MAFFT
        self.mafft_maxiter              = args.mafft_maxiter                # MAFFT max iterations
        self.pxclsq_threshold           = args.mafft_pxclsq_threshold       # pxclsq probability threshold (MAFFT)
        self.thread_divisor             = args.mafft_thread_divisor         # Thread division factor
        
        # Tree
        self.relative_cutoff            = args.tree_relative_cutoff         # Relative cutoff for trimming tips
        self.absolute_cutoff            = args.tree_absolute_cutoff         # Absolute cutoff for trimming tips
        self.branch_cutoff              = args.tree_branch_cutoff           # Branch cutoff for cutting branches
        self.mask_paralogs              = args.tree_mask_paralogs           # Mask paraphyletic tips
        self.outlier_ratio              = args.tree_outlier_ratio           # Ratio threshold for outlier detection
        self.max_trim_iterations        = args.tree_max_trim_iterations     # Maximum trimming iterations
        self.min_subtree_taxa           = args.tree_min_subtree_taxa        # Minimum taxa for valid subtrees
        self.min_tree_leaves            = args.tree_min_leaves              # Minimum leaves for valid tree
        
        # Prune
        self.prune_relative_cutoff      = args.prune_relative_cutoff        # Relative tip cutoff for pruning
        self.prune_absolute_cutoff      = args.prune_absolute_cutoff        # Absolute tip cutoff for pruning
        self.prune_outlier_ratio        = args.prune_outlier_ratio          # Outlier ratio for pruning
        self.prune_max_trim_iterations  = args.prune_max_trim_iterations    # Max trim iterations for pruning
        self.prune_min_tree_leaves      = args.prune_min_tree_leaves         # Min tree leaves for pruning
        
        # PRANK
        self.prank_pxclsq_threshold     = args.prank_pxclsq_threshold       # pxclsq probability threshold (PRANK)
        self.bootstrap_replicates       = args.prank_bootstrap              # IQ-TREE bootstrap replicates
        
        # Super 
        self.super_bootstrap            = args.super_bootstrap              # Supermatrix bootstrap replicates
        
        # Terminal
        self.hc                         = args.highlight_color              # Highlight color
        self.bc                         = args.background_color             # Background color
        self.log                        = args.log                          # Log level

        # Developer Flags
        self.skip                       = args.SKIP if args.SKIP else '-'   # Skip flag
        self.save                       = args.SAVE                         # Save flag
        self.clutter                    = args.clutter                      # Clutter flag

        # Output
        self.output_super_matrix        = args.output_super_matrix          # Output supermatrix

        # Directories and Files
        self.dir_base                   = Path(args.input_dir)                                            # Input Base Directory
        self.output_dir                 = Path(args.output_dir) if args.output_dir else self.dir_base # Output Directory
        self.dir_temp_fasta             = Path(os.path.join(self.output_dir, 'temp_fasta'))           # Temp Fasta Directory
        self.dir_treeforge              = Path(os.path.join(self.output_dir, 'TreeForge'))            # TreeForge Directory
        self.files_fasta                = [Path(f) for f in glob(os.path.join(self.dir_base, '*')) if f.endswith(tuple(['.fasta', '.fa', '.fas', '.fna']))]                     # List of FASTA files
        self.files_fasta_renamed        = [Path(f) for f in glob(os.path.join(self.dir_treeforge, 'temp_fasta', '*')) if f.endswith(tuple(['.fasta', '.fa', '.fas', '.fna']))]  # List of renamed FASTA files
                
        # master dictionary and skip flags
        self.master_dict = {}   # Master dictionary
        self._skip_check()      # Skip flags
        self._dict_setup()      # Master dictionary
        
        # Clutter management
        self.files_to_preserve = set()  # Files that should not be deleted
        
        self._load_previous()   # Load previous run
        self._make_dirs()       # Make directories

        self.run_timestamp  = datetime.now().strftime('%Y%m%d_%H%M%S')                       # Timestamp

        self.printClass     = PrintOut(self.log,args.highlight_color, args.background_color) # Print class
        self.printout       = self.printClass.printout                                       # Print function

    def _skip_check(self):
        """Set skip flags for each step."""
        skip_mapping = {
            'b': 'blast',
            'm': 'mcl',
            'a': 'mafft',
            't': 'tree',
            'p': 'prune',
            'r': 'prank',
            's': 'astral'
        }
        
        for flag, step in skip_mapping.items():
            skip_value = flag in self.skip
            setattr(self, f'{step}_skip', skip_value)

    def _dict_setup(self):
        """Setup master dictionary with directory structure."""
        """Tried to find a middle ground between organization and ll ./x/y/z/ hell"""
        
        self.dir_logs = Path(os.path.join(self.dir_treeforge, 'logs'))

        self.master_dict = {
            'base'      : {'fai': Path(os.path.join(self.dir_treeforge, 'fai')),
                           'temp_fasta': Path(os.path.join(self.dir_treeforge, 'temp_fasta'))
                           },
            'blast'     : {
                'dir': Path(os.path.join(self.dir_treeforge, 'blast')),
                'files': {
                    'concatenated': Path(os.path.join(self.dir_treeforge, 'blast', 'concatenated.fasta')),
                    'raw_blast'   : Path(os.path.join(self.dir_treeforge, 'blast', 'raw.blast'))
                }
            },
            'mcl'       : {
                'dir': Path(os.path.join(self.dir_treeforge, 'mcl')),
            },
            'prune'     : {'dir': Path(os.path.join(self.dir_treeforge, 'prune')),
                           'ortho1to1': Path(os.path.join(self.dir_treeforge, 'prune', '1to1ortho'))},
            'prank'     : {'dir': Path(os.path.join(self.dir_treeforge, 'prank'))},
            'super'     : {'dir': Path(os.path.join(self.dir_treeforge, 'super'))},
            'logs'      : {'dir': Path(os.path.join(self.dir_treeforge, 'logs'))},
            'gene_trees': {'dir': Path(os.path.join(self.dir_treeforge, 'gene_trees'))}
        }

    def _iteration_dirs(self, iteration):
        """Make sure that the directory for an iteration exists."""
        iter_key = f'iter_{iteration}'
        if iter_key not in self.master_dict:
            self.master_dict[iter_key] = {
                'mafft': {'dir': Path(os.path.join(self.dir_treeforge, f'iter_{iteration}', 'mafft'))},
                'tree': {
                    'dir'      : Path(os.path.join(self.dir_treeforge, f'iter_{iteration}', 'tree')),
                    'trimmed'  : Path(os.path.join(self.dir_treeforge, f'iter_{iteration}', 'tree', 'trimmed'))
                }
            }
            self.master_dict[iter_key]['mafft']['dir'].mkdir(parents=True, exist_ok=True)
            self.master_dict[iter_key]['tree']['dir'].mkdir(parents=True, exist_ok=True)
            self.master_dict[iter_key]['tree']['trimmed'].mkdir(parents=True, exist_ok=True)

    def _load_previous(self) -> None:
        """Load data from previous run.
        This will require dev flags to be set.
        I'm not sure that the save/load functionality works correctly but it's a start.
        Suppose i'm really the only one that needs to know about this."""
        skip_mapping = {
            'b': 'blast',
            'm': 'mcl',
            'a': 'mafft',
            't': 'tree',
            'p': 'prune',
            'w': 'write_tree',
            'r': 'prank',
            's': 'super'
        }
        if self.dir_logs.exists():
            file = Path(os.path.join(self.dir_logs, 'previous_run.json'))
            if not file.exists():
                return
            with open(file, 'r') as f:
                prev_data = json.load(f)
                
            for flag, step in skip_mapping.items():
                skip_value = flag in self.skip
                setattr(self, f'{step}_skip', skip_value)
                
                if skip_value and step in ['mafft', 'tree']:
                    for key in prev_data.keys():
                        if key.startswith('iter_'):
                            iter_num = int(key.split('_')[1])
                            if step in prev_data[key]:
                                self._iteration_dirs(iter_num)
                                self.master_dict[key][step] = prev_data[key][step]
                elif skip_value and step in prev_data:
                    self.master_dict[step] = prev_data[step]

    def _make_dirs(self):
        """Create all the directories."""
        def create_dirs(d: dict):
            for key, value in d.items():
                if isinstance(value, dict):
                    create_dirs(value)
                elif isinstance(value, Path) and key == 'dir':
                    value.mkdir(parents=True, exist_ok=True)

        create_dirs(self.master_dict)

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

    def _prev_print(self, step_name):
        """Print metrics from previous run.
        More dev stuff."""
        
        def print_nested_dict(d, exclude_keys=None):
            """Print nested dict as sorted list, excluding certain keys."""
            if exclude_keys is None:
                exclude_keys = set()
            
            result = []
            for key, value in sorted(d.items()):
                if key in exclude_keys:
                    continue
                if isinstance(value, dict):
                    for nested_key, nested_value in sorted(value.items()):
                        if isinstance(nested_value, (str, int, float, bool)):
                            result.append((nested_key, nested_value))
                        elif isinstance(nested_value, list):
                            result.append((nested_key, len(nested_value)))
                elif isinstance(value, (str, int, float, bool)):
                    if key not in self.printClass.key_translate:
                        continue
                    result.append((key, value))
                elif isinstance(value, list):
                    if key not in self.printClass.key_translate:
                        continue
                    result.append((key, len(value)))
            return result
        
        def print_tree_metrics(metrics):
            """Print Tree step metrics in same format as regular run."""
            if not isinstance(metrics, dict):
                return
            
            if 'trim_tips' in metrics:
                self.printout('metric', 'Trim Tips')
                trim_data = metrics['trim_tips']
                summary_metrics = print_nested_dict(trim_data, exclude_keys={'file_details'})
                if summary_metrics:
                    self.printout('metric', summary_metrics)
                if 'file_details' in trim_data and isinstance(trim_data['file_details'], list):
                    for detail in trim_data['file_details']:
                        if isinstance(detail, dict):
                            detail_metrics = print_nested_dict(detail)
                            if detail_metrics:
                                self.printout('metric', detail_metrics)
            
            if 'mask_tips' in metrics:
                self.printout('metric', 'Mask Tips')
                mask_data = metrics['mask_tips']
                summary_metrics = print_nested_dict(mask_data, exclude_keys={'file_details'})
                if summary_metrics:
                    self.printout('metric', summary_metrics)
                if 'file_details' in mask_data and isinstance(mask_data['file_details'], list):
                    for detail in mask_data['file_details']:
                        if isinstance(detail, dict):
                            detail_metrics = print_nested_dict(detail)
                            if detail_metrics:
                                self.printout('metric', detail_metrics)
            
            if 'cut_branches' in metrics:
                self.printout('metric', 'Cut Branches')
                cut_data = metrics['cut_branches']
                summary_metrics = print_nested_dict(cut_data, exclude_keys={'file_details'})
                if summary_metrics:
                    self.printout('metric', summary_metrics)
                if 'file_details' in cut_data and isinstance(cut_data['file_details'], list):
                    for detail in cut_data['file_details']:
                        if isinstance(detail, dict):
                            detail_metrics = print_nested_dict(detail)
                            if detail_metrics:
                                self.printout('metric', detail_metrics)
            
            if 'write_tree' in metrics:
                self.printout('metric', 'Write Tree')
                write_metrics = print_nested_dict(metrics['write_tree'])
                if write_metrics:
                    self.printout('metric', write_metrics)
        
        if step_name.lower() == 'tree':
            for key in self.master_dict.keys():
                if key.startswith('iter_') and 'tree' in self.master_dict[key]:
                    metrics = self.master_dict[key]['tree']
                    print_tree_metrics(metrics)
        elif step_name.lower() == 'mafft':
            for key in self.master_dict.keys():
                if key.startswith('iter_') and 'mafft' in self.master_dict[key]:
                    metrics = self.master_dict[key]['mafft']
                    if 'mafft' in metrics:
                        mafft_metrics = print_nested_dict(metrics['mafft'])
                        if mafft_metrics:
                            self.printout('metric', mafft_metrics)
        else:
            if step_name.lower() in self.master_dict:
                metrics = self.master_dict[step_name.lower()]
                if step_name.lower() == 'blast':
                    blast_metrics = print_nested_dict(metrics, exclude_keys={'files'})
                    if blast_metrics:
                        self.printout('metric', blast_metrics)
                elif step_name.lower() == 'mcl':
                    mcl_metrics = print_nested_dict(metrics)
                    if mcl_metrics:
                        self.printout('metric', mcl_metrics)
                elif step_name.lower() == 'prune':
                    if 'prune' in metrics:
                        prune_metrics = print_nested_dict(metrics['prune'])
                        if prune_metrics:
                            self.printout('metric', prune_metrics)
                elif step_name.lower() == 'prank':
                    prank_metrics = print_nested_dict(metrics)
                    if prank_metrics:
                        self.printout('metric', prank_metrics)
                elif step_name.lower() == 'astral':
                    astral_metrics = print_nested_dict(metrics)
                    if astral_metrics:
                        self.printout('metric', astral_metrics)

    def _execute(self, step_name, step_method, skip_flag):
        """Check that a step isnt being skipped before running it.
        Even more dev stuff."""
        if not getattr(self, skip_flag, False):
            return step_method()
        self.printout('subtitle', f'{step_name} Skipped')
        self._prev_print(step_name)
        return {'status': 'skipped'}

    def _clutter_check(self):
        if self.save and self.clutter:
            self.printout('error', 'Cannot Save with Clutter flag set')
            self.printout('error', 'Retaining files as safety measure')
            return
        if self.clutter:
            self.printout('metric', 'Performing final cleanup')
            
            if self.dir_treeforge.exists():
                for item in self.dir_treeforge.iterdir():
                    if item.is_dir() and item.name not in {'logs', 'gene_trees'}:
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
        """Move the fai files to the TreeForge/fai directory.
        Makes things look cleaner."""
        files_fai = glob(os.path.join(self.dir_base, '*.fai'))
        dir_fai = os.path.join(self.dir_treeforge, 'fai')
        os.makedirs(dir_fai, exist_ok=True)
        for fai in files_fai:
            fai_dest = os.path.join(dir_fai, os.path.basename(fai))
            shutil.move(fai, fai_dest)
        if os.path.exists(os.path.join(self.dir_base, 'phyx.logfile')):
            os.remove(os.path.join(self.dir_base, 'phyx.logfile'))

    def _finish_pipeline(self, gene_tree_count=None):
        """Finish the pipeline and cleanup files."""
        self.printout('title', 'TreeForge Complete')
        self._move_fai()
        final_cleanup(
            root_dir=self.dir_treeforge,
            logs_dir=self.dir_logs,
            printout=self.printout,
            protected_files=None,
            protected_dirs=None,
            protected_exts=None
        )
        self._save_csv()
        
        if gene_tree_count is None:
            gene_tree_count = len(glob(os.path.join(self.dir_treeforge, 'gene_trees', '*.tre')))
        
        out_dict = {'Species Tree'   : os.path.join(self.dir_treeforge, 'SpeciesTree.tre'),
                    'Metrics CSV'    : os.path.join(self.dir_treeforge, 'summary.csv'),
                    'Gene Trees'     : os.path.join(self.dir_treeforge, 'gene_trees'),
                    'Gene Tree Count': gene_tree_count,
                    }
        if self.output_super_matrix:
            out_dict['SuperMatrix'] = os.path.join(self.dir_treeforge, 'SuperMatrix.tre')
        self.printout('final', out_dict)

    def _save_csv(self):
        """Save the metrics to a csv file.
        Really needs some work. A lot of metrics here that arent useful, and probably 
        a lot of good metrics I'm not keeping track of."""
        def iter_dct(dct):
            result = {}
            for k,v in dct.items():
                if v is None:
                    result[k] = ''
                elif isinstance(v, list):
                    result[k] = len(v)
                elif isinstance(v, str):
                    result[k] = v
                elif isinstance(v, dict):
                    result.update(iter_dct(v))
                else:
                    result[k] = v
            return result
                    
        flattened_dct = {k: {} for k in self.master_dict.keys()}
        for k,v in self.master_dict.items():
            flattened = iter_dct(v)
            flattened_dct[k] = flattened
        
        csv_file = Path(os.path.join(self.dir_treeforge, f'summary.csv'))
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            header_row = []
            value_row = []
            
            root_key_order = ['blast', 'mcl']
            existing_iterations = []
            for key in flattened_dct.keys():
                if key.startswith('iter_'):
                    existing_iterations.append(key)
            existing_iterations.sort(key=lambda x: int(x.split('_')[1]))
            root_key_order.extend(existing_iterations)
            root_key_order.extend(['prune', 'prank', 'astral'])
            
            for root_key in root_key_order:
                if root_key in flattened_dct:
                    nested_dict = flattened_dct[root_key]
                    header_row.append(root_key)
                    value_row.append('')
                    
                    for metric_key, metric_value in nested_dict.items():
                        if metric_key not in self.printClass.key_translate:
                            continue
                        else:
                            metric_key = self.printClass.key_translate[metric_key]
                        header_row.append(metric_key)
                        value_row.append(str(metric_value))
            
            writer.writerow(header_row)
            writer.writerow(value_row)
        
        self.printout('metric', f'Metrics saved to {csv_file}')
        return flattened_dct

    def _contig_renamer(self):
        """Rename the contigs/transcripts/scaffolds.
        This should take care of any naming schemes or weird characters.
        Prevents treeforge and the stuff it wraps from crashing"""
        contig_renamer = Contig(self.dir_base, self.master_dict['base']['temp_fasta'])
        self.renamed_map = contig_renamer.name_mapping
        contig_renamer.rename_contigs(self.files_fasta)
        return contig_renamer.return_dict
    
    def run(self):
        self._contig_renamer()
        self._execute('BLAST', self.blast, 'blast_skip')
        mcl_result = self._execute('MCL', self.mcl, 'mcl_skip')
        self._cleanup_intermediate_files('blast')
        
        if not self.mcl_skip and mcl_result and mcl_result.get('mcl_count', 0) == 0:
            self.printout('error', 'No clusters produced by MCL')
            self.printout('error', 'Try adjusting the MCL parameters')
            # self._finish_pipeline(gene_tree_count=0)
            return
        
        continue_iterations = True
        
        for i in range(self.iterations):
            if not continue_iterations:
                self.printout('info', f'Skipping remaining iterations (iter_{i} to iter_{self.iterations-1})')
                break
                
            self._execute('MAFFT', self.mafft, 'mafft_skip')
            
            if (f'iter_{self.current_iter}' in self.master_dict and 
                'mafft' in self.master_dict[f'iter_{self.current_iter}'] and
                'mafft' in self.master_dict[f'iter_{self.current_iter}']['mafft'] and
                self.master_dict[f'iter_{self.current_iter}']['mafft']['mafft'].get('status') == 'no_files_to_process'):
                self.printout('info', f'No files to process in iteration {int(self.current_iter) + 1}, stopping iterations early')
                continue_iterations = False
            else:
                tree_result = self._execute('Tree', self.tree, 'tree_skip')
                self._cleanup_intermediate_files('mafft')
                self._cleanup_intermediate_files('mcl')
                
                if tree_result and 'write_tree' in tree_result:
                    modified_count = tree_result['write_tree'].get('modified_files_count', 0)
                    self.printout('info', f'Iteration {self.current_iter}: {modified_count} files modified')
                    
                    if modified_count == 0:
                        self.printout('info', f'No files modified in iteration {int(self.current_iter) + 1}')
                        continue_iterations = False
            
            if self.current_iter > 0:
                prev_iter = self.current_iter - 1
                self._cleanup_intermediate_files('mafft', iter_override=prev_iter)
                self._cleanup_intermediate_files('tree', iter_override=prev_iter)
            
            self.current_iter += 1 if self.current_iter < self.iterations - 1 else 0
            
        self._execute('Prune', self.prune, 'prune_skip')
        for i in range(self.iterations):
            self.current_iter = i
            if f'iter_{i}' in self.master_dict:
                self._cleanup_intermediate_files('tree')
                self._cleanup_preserved_files('prune')
            
        self._execute('PRANK', self.prank, 'prank_skip')
        self._cleanup_intermediate_files('prank')
        self._cleanup_preserved_files('prank')
        self._cleanup_preserved_files('mcl')
        self._execute('Astral', self.astral, 'astral_skip')
        self._cleanup_preserved_files('astral')
        self._finish_pipeline()

    @bilge_crew()
    def blast(self):
        self.printout('subtitle', 'BLAST')
        renamed_files = [Path(f) for f in glob(os.path.join(self.dir_treeforge, 'temp_fasta', '*')) if f.endswith(tuple(['.fasta', '.fa', '.fas', '.fna']))]
        blast = Blast(self.dir_base,                                                # Input Base Directory
                      self.dir_treeforge,                                           # TreeForge Directory
                      renamed_files,                                                # List of renamed FASTA files
                      self.master_dict['blast']['files']['concatenated'],           # Concatenated FASTA file
                      self.master_dict['blast']['files']['raw_blast'],              # Raw BLAST file
                      self.threads,                                                 # Number of threads
                      self.log,                                                     # Log level
                      self.hc,                                                      # Highlight color
                      self.bc,                                                      # Background color
                      self.blast_evalue,                                            # BLAST E-value threshold
                      self.blast_max_targets)                                       # BLAST max target sequences
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
                  self.min_seq_length)                                                  # Minimum sequence length
        mcl.run()
        return mcl.return_dict

    @bilge_crew()
    def mafft(self):
        self.printout('subtitle', f'MAFFT - {self.current_iter + 1}')
        self._iteration_dirs(self.current_iter)
        mafft = Mafft(self.dir_base,                                                    # Input Base Directory
                      self.master_dict[f'iter_{self.current_iter}']['mafft']['dir'],    # MAFFT Directory
                      self.threads,                                                     # Number of threads
                      self.log,                                                         # Log level
                      self.hc,                                                          # Highlight color
                      self.bc,                                                          # Background color
                      self.mafft_maxiter,                                               # MAFFT max iterations
                      self.pxclsq_threshold,                                            # pxclsq probability threshold
                      self.thread_divisor)                                              # Thread division factor
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
                    clutter=self.clutter)                                               # Pass clutter flag
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
                      self.prune_min_tree_leaves)                                       # Min tree leaves for pruning
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
                      self.bootstrap_replicates)                                        # IQ-TREE bootstrap replicates
        prank.run()
        return prank.return_dict

    @bilge_crew()
    def astral(self):
        self.printout('subtitle', 'Astral')
        astral = Astral(self.master_dict['prank']['dir'],                               # PRANK Directory
                        self.master_dict['super']['dir'],                               # Supertree Directory
                        self.dir_treeforge,                                             # TreeForge Directory
                        self.renamed_map,                                               # Renamed Map
                        self.output_super_matrix,                                       # Output supermatrix
                        self.super_bootstrap,                                           # Supermatrix bootstrap replicates
                        self.threads,                                                   # Number of threads
                        self.log,                                                       # Log level
                        self.hc,                                                        # Highlight color
                        self.bc)                                                        # Background color
        astral.run()
        return astral.return_dict
