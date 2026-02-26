"""
                  ╭──────────────────╮
          ╭───────┴──────╮           │   
      ╭───┴───╮       ╭──┴──╮     ╭──┴──╮
    ╭─┴─╮   ╭─┴─╮   ╭─┴─╮   │   ╭─┴─╮   │
    T   r   e   e   -   F   o   r   g   e

                                 EC BRETZ
"""

import os
import sys
import argparse
import re
from glob                import glob
from pathlib             import Path
from datetime            import datetime
from core.utils.logo     import print_logo, print_help, print_args, print_hcluster_warning
from core.datahub        import DataHub
from core.utils.deps     import Deps
from core.utils.printout import PrintOut
from core.utils.config   import ConfigManager
from core.utils.constants import FASTA_EXTENSIONS


MAJOR = 0
MINOR = 4
PATCH = 1

VERSION = f"{MAJOR}.{MINOR}.{PATCH}"

class TreeForge:
    def __init__(self):
        self.early_log_lines = []
        self._log_writer = lambda text: self.early_log_lines.append(text)
        nocolor_init = '--nocolor' in sys.argv

        self.highlight_color, self.background_color = print_logo(VERSION, log_func=self._log_writer, nocolor=nocolor_init)
        self.printClass                             = PrintOut('', self.highlight_color, self.background_color)
        self.printout                               = self.printClass.printout
        self.config_manager                         = ConfigManager(self.highlight_color, self.background_color)

        self.original_defaults                      = self.config_manager.get_defaults_dict().copy()
        self.original_defaults['config']            = None
        self.original_defaults['config_create']     = False
        self.original_defaults['log']               = 3
        self.original_defaults['subprocess_logs']   = False
        self.original_defaults['help']              = False
        self.original_defaults['version']           = False
        self.original_defaults['nocolor']          = False

    def parser(self, config_values=None):
        if len(sys.argv) == 1:
            print_help(self.highlight_color, self.original_defaults)
            sys.exit()
        parser = argparse.ArgumentParser(description="TreeForge", add_help=False)

        self.defaults_dict                    = self.config_manager.get_defaults_dict()
        self.defaults_dict['config']          = None
        self.defaults_dict['config_create']   = False
        self.defaults_dict['log']             = 3
        self.defaults_dict['subprocess_logs'] = False
        self.defaults_dict['help']            = False
        self.defaults_dict['version']         = False
        self.defaults_dict['nocolor']        = False

        if config_values:
            for key, value in config_values.items():
                if key in self.defaults_dict:
                    self.defaults_dict[key] = value
        


        def error(message):
            print_help(self.highlight_color, self.defaults_dict)
            message       = message.replace('unknown option:', '')
            message       = message.strip()
            chunks        = []
            current_chunk = ""
            words         = message.split()
            for word in words:
                if len(current_chunk) + len(word) + 1 <= 58:
                    if current_chunk:
                        current_chunk += " " + word
                    else:
                        current_chunk = word
                else:
                    if current_chunk:
                        chunks.append(current_chunk)
                    current_chunk = word
            if current_chunk:
                chunks.append(current_chunk)
            for chunk in chunks:
                self.printout('error', chunk)
            self.printout('error', 'Try --help for more information.')
            sys.exit(2)
        parser.error = error

        # Basic

        parser.add_argument("--input-dir",                  "-d",   type=str,    help="Directory of FASTA files",                   default=self.defaults_dict['input_dir'])
        parser.add_argument("--iter",                       "-i",   type=int,    help="Number of iterations",                       default=self.defaults_dict['iter'])
        parser.add_argument("--threads",                    "-t",   type=int,    help="Number of threads",                          default=self.defaults_dict['threads'])
        parser.add_argument("--clutter",                    "-c",                help="Remove Intermediate Files",                  default=self.defaults_dict['clutter'],                  action="store_true")
        parser.add_argument("--output-dir",                 "-o",   type=str,    help="Output directory",                           default=self.defaults_dict['output_dir'])
        parser.add_argument("--subprocess-logs",            "-sl",               help="Save subprocess stdout/stderr to log files", default=self.defaults_dict['subprocess_logs'],         action="store_true")

        # BLAST
        parser.add_argument("--blast-evalue",               "-be",  type=float,  help="BLAST E-value threshold",                    default=self.defaults_dict['blast_evalue'])
        parser.add_argument("--blast-max-targets",          "-bm",  type=int,    help="BLAST max target sequences",                 default=self.defaults_dict['blast_max_targets'])
        
        # MCL
        parser.add_argument("--mcl-hit-frac-cutoff",        "-mf",  type=float,  help="Hit fraction cutoff for MCL",                default=self.defaults_dict['mcl_hit_frac_cutoff'])
        parser.add_argument("--mcl-minimum-taxa",           "-mt",  type=int,    help="Minimum taxa for MCL clustering",            default=self.defaults_dict['mcl_minimum_taxa'])
        parser.add_argument("--mcl-inflation",              "-mi",  type=float,  help="MCL inflation parameter",                    default=self.defaults_dict['mcl_inflation'])
        parser.add_argument("--mcl-perfect-identity",       "-mp",  type=float,  help="Perfect identity threshold",                 default=self.defaults_dict['mcl_perfect_identity'])
        parser.add_argument("--mcl-coverage-threshold",     "-mc",  type=float,  help="Coverage threshold for identical sequences", default=self.defaults_dict['mcl_coverage_threshold'])
        parser.add_argument("--mcl-min-seq-length",         "-ml",  type=int,    help="Minimum sequence length",                    default=self.defaults_dict['mcl_min_seq_length'])
        
        # MAFFT    
        parser.add_argument("--mafft-maxiter",              "-mm",  type=int,    help="MAFFT max iterations",                       default=self.defaults_dict['mafft_maxiter'])
        parser.add_argument("--mafft-pxclsq-threshold",     "-mx",  type=float,  help="pxclsq probability threshold (MAFFT)",       default=self.defaults_dict['mafft_pxclsq_threshold'])
        parser.add_argument("--mafft-thread-divisor",       "-md",  type=int,    help="Thread division factor",                     default=self.defaults_dict['mafft_thread_divisor'])
        
        # Tree
        parser.add_argument("--tree-start-from-prev",       "-tg",               help="Reiterative iqtree start tree",              default=self.defaults_dict.get('tree_start_from_prev', False), action="store_true")
        parser.add_argument("--tree-relative-cutoff",       "-tr",  type=float,  help="Relative cutoff for trimming tips",          default=self.defaults_dict['tree_relative_cutoff'])
        parser.add_argument("--tree-absolute-cutoff",       "-ta",  type=float,  help="Absolute cutoff for trimming tips",          default=self.defaults_dict['tree_absolute_cutoff'])
        parser.add_argument("--tree-branch-cutoff",         "-tb",  type=float,  help="Branch cutoff for cutting branches",         default=self.defaults_dict['tree_branch_cutoff'])
        parser.add_argument("--tree-mask-paralogs",         "-tm",  type=str,    help="Mask paraphyletic tips",                     default=self.defaults_dict['tree_mask_paralogs'],       choices=['y', 'n'])
        parser.add_argument("--tree-outlier-ratio",         "-to",  type=float,  help="Ratio threshold for outlier detection",      default=self.defaults_dict['tree_outlier_ratio'])
        parser.add_argument("--tree-max-trim-iterations",   "-ti",  type=int,    help="Maximum trimming iterations",                default=self.defaults_dict['tree_max_trim_iterations'])
        parser.add_argument("--tree-min-subtree-taxa",      "-ts",  type=int,    help="Minimum taxa for valid subtrees",            default=self.defaults_dict['tree_min_subtree_taxa'])
        parser.add_argument("--tree-min-leaves",            "-tl",  type=int,    help="Minimum leaves for valid tree",              default=self.defaults_dict['tree_min_leaves']) 
        
        # Prune
        parser.add_argument("--prune-orthocutoff",          "-po",  type=int,    help="Ortholog Minimum Taxa Cutoff",               default=self.defaults_dict['prune_orthocutoff'])
        parser.add_argument("--prune-relative-cutoff",      "-pr",  type=float,  help="Relative tip cutoff for pruning",            default=self.defaults_dict['prune_relative_cutoff'])
        parser.add_argument("--prune-absolute-cutoff",      "-pa",  type=float,  help="Absolute tip cutoff for pruning",            default=self.defaults_dict['prune_absolute_cutoff'])
        parser.add_argument("--prune-outlier-ratio",        "-por", type=float,  help="Outlier ratio for pruning",                  default=self.defaults_dict['prune_outlier_ratio'])
        parser.add_argument("--prune-max-trim-iterations",  "-pmt", type=int,    help="Max trim iterations for pruning",            default=self.defaults_dict['prune_max_trim_iterations'])
        parser.add_argument("--prune-min-tree-leaves",      "-pml", type=int,    help="Min tree leaves for pruning",                default=self.defaults_dict['prune_min_tree_leaves'])
        
        # PRANK
        parser.add_argument("--prank-seqtype",              "-ps",  type=str,   help="Sequence type for PRANK",                     default=self.defaults_dict['prank_seqtype'],            choices=['dna', 'aa']) #maybe some day
        parser.add_argument("--prank-pxclsq-threshold",     "-pp",  type=float, help="pxclsq probability threshold",                default=self.defaults_dict['prank_pxclsq_threshold'])
        parser.add_argument("--prank-bootstrap",            "-pb",  type=int,   help="IQ-TREE bootstrap replicates",                default=self.defaults_dict['prank_bootstrap'])
        
        # Super
        parser.add_argument("--super-bootstrap",            "-sb",  type=int,   help="Supermatrix bootstrap replicates",            default=self.defaults_dict['super_bootstrap'])
        parser.add_argument("--bes-support",                "-bs",  type=float, help="Molecular distance support",                  default=self.defaults_dict['bes_support'])
        parser.add_argument("--super-matrix",               "-sm",              help="Output supermatrix",                          default=self.defaults_dict['super_matrix'],            action="store_true")
        
        # HCluster
        parser.add_argument("--hcluster-enabled",           "-hc",              help="Run Hierarchical Clustering",                 default=self.defaults_dict['hcluster_enabled'],        action="store_true")
        parser.add_argument("--hcluster-id",                "-hi",  type=float, help="HCluster id",                                 default=self.defaults_dict['hcluster_id'])
        parser.add_argument("--hcluster-iddef",             "-hid", type=int,   help="HCluster iddef",                              default=self.defaults_dict['hcluster_iddef'])
        parser.add_argument("--hcluster-tree",              "-hg",  type=str,   help="Path to custom guide tree",                   default=self.defaults_dict['hcluster_tree'])
        parser.add_argument("--hcluster-tool",              "-ht",  type=str,   help="Clustering tool (vsearch or mmseqs2)",        default=self.defaults_dict['hcluster_tool'],           choices=['vsearch', 'mmseqs2'])
        parser.add_argument("--hcluster-use-busco",         "-hub",             help="Use BUSCO-filtered files for clustering",     default=self.defaults_dict['hcluster_use_busco'],      action="store_true")
        
        # BUSCO
        parser.add_argument("--busco-evalue",               "-bce", type=float, help="BUSCO BLAST E-value threshold",               default=self.defaults_dict['busco_evalue'])
        parser.add_argument("--busco-max-targets",          "-bct", type=int,   help="BUSCO BLAST max target sequences",            default=self.defaults_dict['busco_max_targets'])
        parser.add_argument("--busco-coverage-threshold",   "-bcc", type=float, help="BUSCO coverage threshold",                    default=self.defaults_dict['busco_coverage_threshold'])
        # Configuration
        parser.add_argument("--config",                             type=str,   help="Path to configuration file",                  default=None,                                           nargs="?", const=True)
        parser.add_argument("--config-create",                      type=str,   help="Create a configuration template",             default=False,                                          nargs="?", const="config.yaml")
        parser.add_argument("--config-save",                        type=str,   help="Save current arguments to config file",       default=False,                                          nargs="?", const="config.yaml")
        
        # Standard
        parser.add_argument("--version",                    "-v",               help="Print version",                                                                                       action="store_true",)
        parser.add_argument("--nocolor",                                        help="Disable colored terminal output",                                                             default=self.defaults_dict['nocolor'],                    action="store_true")
        parser.add_argument("--log",                        "-l",   type=int,   help=argparse.SUPPRESS,                             default=self.defaults_dict['log'],                      choices=[0, 1, 2, 3, 4])
        parser.add_argument("--help",                       "-h",               help=argparse.SUPPRESS,                             default=self.defaults_dict['help'],                     action="store_true")
        
        return parser.parse_args()
    
    def _file_exists(self, dir):
        if os.path.exists(dir):
            files = []
            for ext in FASTA_EXTENSIONS:
                files.extend(glob(os.path.join(dir, f'*{ext}')))
            if not files:
                self.printout('error', 'No FASTA files found in directory')
                sys.exit(1)
            return files
        else:
            self.printout('error', 'Directory does not exist')
            sys.exit(1)

    def run(self):
        args = self.parser()

        nocolor = getattr(args, 'nocolor', False)
        if nocolor:
            self.highlight_color = ''
            self.background_color = ''
            self.printClass.set_nocolor(True)
            self.config_manager.set_nocolor(True)

        if args.help:
            print_help(self.highlight_color, self.original_defaults, nocolor=nocolor)
            sys.exit()
        if args.version:
            print(f"TreeForge v{VERSION}")
            sys.exit()
        
        if getattr(args, 'config_create', False):
            config_filename = args.config_create if args.config_create != "config.yaml" else "config.yaml"
            if not config_filename.endswith('.yaml'):
                config_filename += '.yaml'
            config_path = Path(args.output_dir) / config_filename if args.output_dir else Path(config_filename)
            self.config_manager.create_config(config_path)
            sys.exit(0)

        if args.output_dir:
            Path(args.output_dir).mkdir(parents=True, exist_ok=True)

        config_values = None
        if args.config is not None:
            if args.config is True:
                self.printout('error', 'No configuration file path provided with --config')
                self.printout('error', 'Usage: --config /path/to/config.yaml')
                self.printout('error', 'Or use --config-create to create a new config file')
                sys.exit(1)
            else:
                config_path = Path(args.config)
                if config_path.exists():
                    config_values = self.config_manager.load_config(config_path)
                else:
                    self.printout('error', f'Config file not found: {config_path}')
                    sys.exit(1)
        
        if config_values:
            args = self.parser(config_values)
            nocolor = getattr(args, 'nocolor', False)
            if nocolor:
                self.highlight_color = ''
                self.background_color = ''
                self.printClass.set_nocolor(True)
                self.config_manager.set_nocolor(True)

        # you can load a config with --config config.yaml and then save a copy with updated values with --config-save config2.yaml
        config_save_used = getattr(args, 'config_save', False)
        if config_save_used:
            config_dict = {}
            for key, value in args.__dict__.items():
                if key in self.config_manager.get_defaults_dict():
                    config_dict[key] = value

            config_filename = args.config_save if args.config_save != "config.yaml" else "config.yaml"
            if not config_filename.endswith('.yaml'):
                config_filename += '.yaml'
            config_path = Path(args.output_dir) / config_filename if args.output_dir else Path(config_filename)
            self.config_manager.save_config(config_dict, config_path)
            out_name = str(config_path.absolute()) if len(str(config_path.absolute())) < 32 else '...' + str(config_path.absolute())[-29:]
            self.printout('info', f"Config saved to {out_name}")
        
        passed_args = {k: v for k, v in args.__dict__.items()
                       if k in self.original_defaults and v != self.original_defaults[k] and v != '-'}
        
        args.highlight_color  = self.highlight_color
        args.background_color = self.background_color

        print_args(args, args.highlight_color, passed_args, log_func=self._log_writer, nocolor=nocolor)

        # Print warning if hierarchical clustering is enabled
        # Still very much under construction and minimally tested
        if args.hcluster_enabled:
            print_hcluster_warning(args.highlight_color, log_func=self._log_writer, nocolor=nocolor)

        deps = Deps(args.log, args.highlight_color, args.background_color, nocolor=nocolor)
        deps.check_deps()
        self._file_exists(args.input_dir)
        dataHub = DataHub(args)
        
        if hasattr(dataHub, 'printClass') and hasattr(dataHub.printClass, 'log_handle') and dataHub.printClass.log_handle:
            for line in self.early_log_lines:
                dataHub.printClass.log_handle.write(line + '\n')
            dataHub.printClass.log_handle.flush()
        
        dataHub.run()

if __name__ == "__main__":
    # Test args for debugging
    # sys.argv = [
    #     'treeforge.py',
    #     '-d', '/home/eric/treeTest2',
    #     '-i', '1',
    #     '-t', '1',
    #     '-hf', '0.3',
    #     '-mt', '10',
    #     '-po', '20',
    #     '-ps', 'aa',
    #     '-l', '3'
    # ]
    try:
        main = TreeForge()
        main.run()
    except KeyboardInterrupt:
        print("\n\nProcess interrupted by user.")
        os._exit(130)