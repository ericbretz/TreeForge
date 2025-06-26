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
from core.utils.logo import print_logo, print_help, print_args
from core.datahub import DataHub
from core.utils.deps import Deps
from core.utils.printout import PrintOut


'''
Skip flags:
b: Skip BLAST
m: Skip MCL
a: Skip MAFFT
t: Skip Tree
p: Skip Prune
w: Skip Write
r: Skip Prank
s: Skip Astral
'''

MAJOR = 0
MINOR = 2
PATCH = 0

VERSION = f"{MAJOR}.{MINOR}.{PATCH}"

class TreeForge:
    def __init__(self):
        self.highlight_color, self.background_color = print_logo(VERSION)
        self.printClass = PrintOut('', self.highlight_color, self.background_color)
        self.printout = self.printClass.printout

    def parser(self):
        if len(sys.argv) == 1:
            print_help(self.highlight_color)
            sys.exit()
        parser = argparse.ArgumentParser(description="TreeForge", add_help=False)

        self.defaults_dict = {
            'dir'                     : os.getcwd(),
            'iter'                    : 3,
            'threads'                 : 2,
            'clutter'                 : False,
            'blast_evalue'            : 10.0,
            'blast_max_targets'       : 1000,
            'mcl_hit_frac_cutoff'     : 0.3,
            'mcl_minimum_taxa'        : 10,
            'mcl_inflation'           : 1.4,
            'mcl_perfect_identity'    : 100.0,
            'mcl_coverage_threshold'  : 0.9,
            'mcl_min_seq_length'      : 300,
            'mafft_maxiter'           : 1000,
            'mafft_pxclsq_threshold'  : 0.1,
            'mafft_thread_divisor'    : 4,
            'tree_relative_cutoff'    : 0.2,
            'tree_absolute_cutoff'    : 0.3,
            'tree_branch_cutoff'      : 0.02,
            'tree_mask_paralogs'      : 'n',
            'tree_outlier_ratio'      : 20.0,
            'tree_max_trim_iterations': 10,
            'tree_min_subtree_taxa'   : 4,
            'tree_min_leaves'         : 4,
            'prune_orthocutoff'       : 20,
            'prune_relative_cutoff'   : 0.2,
            'prune_absolute_cutoff'   : 0.3,
            'prank_seqtype'           : 'aa',
            'prank_pxclsq_threshold'  : 0.3,
            'prank_bootstrap'         : 1000,
            'super_bootstrap'         : 1000,
            'log'                     : 3,
            'help'                    : False,
            'version'                 : False,
            'clutter'                 : False,
            'SAVE'                    : False,
            'SKIP'                    : '-',
        }

        def error(message):
            print_help(self.highlight_color)
            message = message.replace('ambiguous option:', '')

            # message = 'debug sladja 23123123 qw asd qwe ads eww ddd d f rasd qwwqec sasd e342m debug debug debug'

            message = message.strip()
            
            chunks = []
            current_chunk = ""
            words = message.split()
            
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
        
        # Basic Parameters
        parser.add_argument("--dir",                        "-d",   type=str,    help="Directory of FASTA files",           default=self.defaults_dict['dir']                      )
        parser.add_argument("--iter",                       "-i",   type=int,    help="Number of iterations",               default=self.defaults_dict['iter']                                )
        parser.add_argument("--threads",                    "-t",   type=int,    help="Number of threads",                  default=self.defaults_dict['threads']                                )
        parser.add_argument("--clutter",                    "-c",                help="Remove Intermediate Files",          action="store_true", default=self.defaults_dict['clutter']                      )

        # BLAST Stage Parameters
        parser.add_argument("--blast-evalue",               "-be",  type=float,  help="BLAST E-value threshold",            default=self.defaults_dict['blast_evalue']                             )
        parser.add_argument("--blast-max-targets",          "-bm",  type=int,    help="BLAST max target sequences",         default=self.defaults_dict['blast_max_targets']                             )
        
        # MCL Stage Parameters
        parser.add_argument("--mcl-hit-frac-cutoff",        "-hf",  type=float,  help="Hit fraction cutoff for MCL",        default=self.defaults_dict['mcl_hit_frac_cutoff']                              )
        parser.add_argument("--mcl-minimum-taxa",           "-mt",  type=int,    help="Minimum taxa for MCL clustering",    default=self.defaults_dict['mcl_minimum_taxa']                               )
        parser.add_argument("--mcl-inflation",              "-mi",  type=float,  help="MCL inflation parameter",            default=self.defaults_dict['mcl_inflation']                              )
        parser.add_argument("--mcl-perfect-identity",       "-mp",  type=float,  help="Perfect identity threshold",         default=self.defaults_dict['mcl_perfect_identity']                            )
        parser.add_argument("--mcl-coverage-threshold",     "-mc",  type=float,  help="Coverage threshold for identical sequences", default=self.defaults_dict['mcl_coverage_threshold']                      )
        parser.add_argument("--mcl-min-seq-length",         "-ml",  type=int,    help="Minimum sequence length",            default=self.defaults_dict['mcl_min_seq_length']                              )
        
        # MAFFT Stage Parameters    
        parser.add_argument("--mafft-maxiter",              "-mm",  type=int,    help="MAFFT max iterations",               default=self.defaults_dict['mafft_maxiter']                             )
        parser.add_argument("--mafft-pxclsq-threshold",     "-mx",  type=float,  help="pxclsq probability threshold (MAFFT)",default=self.defaults_dict['mafft_pxclsq_threshold']                             )
        parser.add_argument("--mafft-thread-divisor",       "-md",  type=int,    help="Thread division factor",             default=self.defaults_dict['mafft_thread_divisor']                                )
        
        # Tree Stage Parameters
        parser.add_argument("--tree-relative-cutoff",       "-tr",  type=float,  help="Relative cutoff for trimming tips",  default=self.defaults_dict['tree_relative_cutoff']                              )
        parser.add_argument("--tree-absolute-cutoff",       "-ta",  type=float,  help="Absolute cutoff for trimming tips",  default=self.defaults_dict['tree_absolute_cutoff']                              )
        parser.add_argument("--tree-branch-cutoff",         "-tb",  type=float,  help="Branch cutoff for cutting branches", default=self.defaults_dict['tree_branch_cutoff']                             )
        parser.add_argument("--tree-mask-paralogs",         "-tm",  type=str,    help="Mask paraphyletic tips",             default=self.defaults_dict['tree_mask_paralogs'], choices=['y', 'n']          )
        parser.add_argument("--tree-outlier-ratio",         "-to",  type=float,  help="Ratio threshold for outlier detection", default=self.defaults_dict['tree_outlier_ratio']                          )
        parser.add_argument("--tree-max-trim-iterations",   "-ti",  type=int,    help="Maximum trimming iterations",        default=self.defaults_dict['tree_max_trim_iterations']                               )
        parser.add_argument("--tree-min-subtree-taxa",      "-ts",  type=int,    help="Minimum taxa for valid subtrees",    default=self.defaults_dict['tree_min_subtree_taxa']                                )
        parser.add_argument("--tree-min-leaves",            "-tl",  type=int,    help="Minimum leaves for valid tree",      default=self.defaults_dict['tree_min_leaves']                                ) 
        
        # Prune Stage Parameters
        parser.add_argument("--prune-orthocutoff",          "-po",  type=int,  help="Ortholog Minimum Taxa Cutoff",       default=self.defaults_dict['prune_orthocutoff']                               )
        parser.add_argument("--prune-relative-cutoff",      "-pr",  type=float,  help="Relative tip cutoff for pruning",    default=self.defaults_dict['prune_relative_cutoff']                              )
        parser.add_argument("--prune-absolute-cutoff",      "-pa",  type=float,  help="Absolute tip cutoff for pruning",    default=self.defaults_dict['prune_absolute_cutoff']                              )
        
        # PRANK Stage Parameters
        parser.add_argument("--prank-seqtype",              "-ps",  type=str,   help="Sequence type for PRANK",             default=self.defaults_dict['prank_seqtype'], choices=['dna', 'aa']       )
        parser.add_argument("--prank-pxclsq-threshold",     "-pp",  type=float, help="pxclsq probability threshold (PRANK)",default=self.defaults_dict['prank_pxclsq_threshold']                               )
        parser.add_argument("--prank-bootstrap", "-pb",  type=int,   help="IQ-TREE bootstrap replicates",        default=self.defaults_dict['prank_bootstrap']                              )
        
        # Super Stage Parameters
        parser.add_argument("--super-bootstrap",            "-sb",  type=int,   help="Supermatrix bootstrap replicates",    default=self.defaults_dict['super_bootstrap']                              )
        
        # Standard Arguments
        parser.add_argument("--version",                    "-v",               help="Print version",                       action="store_true"                       )
        parser.add_argument("--log",                        "-l",   type=int,   help=argparse.SUPPRESS,                     default=self.defaults_dict['log'],         choices=[0, 1, 2, 3, 4])
        parser.add_argument("--help",                       "-h",               help=argparse.SUPPRESS,                     action="store_true", default=self.defaults_dict['help']                       )
        parser.add_argument("--SKIP",                               type=str,   help=argparse.SUPPRESS,                     default=self.defaults_dict['SKIP']                               )
        parser.add_argument("--SAVE",                                           help=argparse.SUPPRESS,                     action="store_true", default=self.defaults_dict['SAVE']        )
        
        return parser.parse_args()
    
    def run(self):
        args                  = self.parser()
        passed_args           = {k: v for k, v in args.__dict__.items() if v != self.defaults_dict[k] and v != '-'}
        args.highlight_color  = self.highlight_color
        args.background_color = self.background_color
        if args.help:
            print_help(args.highlight_color)
            sys.exit()
        if args.version:
            print(f"TreeForge v{VERSION}")
            sys.exit()
        print_args(args, args.highlight_color, passed_args)
        deps = Deps(args.log, args.highlight_color, args.background_color)
        deps.check_deps()
        dataHub = DataHub(args)
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
    main = TreeForge()
    main.run()