import os
import sys
import argparse
from core.utils.logo import print_logo, print_help, print_args
from core.datahub import DataHub
from core.utils.deps import Deps


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

VERSION = "0.1.1"

class TreeForge:
    def __init__(self):
        self.highlight_color, self.background_color = print_logo(VERSION)

    def parser(self):
        if len(sys.argv) == 1:
            print_help(self.highlight_color)
            sys.exit()
        parser = argparse.ArgumentParser(description="TreeForge", add_help=False)
        parser.add_argument("--dir",             "-d", type=str,    help="Directory of FASTA files",         default=os.getcwd()                       )
        parser.add_argument("--iter",            "-i", type=int,    help="Number of iterations",             default=3                                 )
        parser.add_argument("--threads",         "-t", type=int,    help="Number of threads",                default=2                                 )
        parser.add_argument("--hit-frac-cutoff", "-f", type=float,  help="Hit fraction cutoff",              default=0.3                               )
        parser.add_argument("--minimum-taxa",    "-m", type=int,    help="Minimum taxa for MCL clustering",  default=10                                )
        parser.add_argument("--orthocutoff",     "-o", type=float,  help="Ortholog Minimum Taxa Cutoff",     default=20                                )
        parser.add_argument("--seqtype",         "-s", type=str,    help="Sequence type",                    default='aa',      choices=['dna', 'aa']  )
        parser.add_argument("--clutter",         "-c",              help="Remove Intermediate Files",        action="store_true"                       )
        parser.add_argument("--version",         "-v",              help="Print version",                    action="store_true"                       )
        parser.add_argument("--log",             "-l", type=int,    help=argparse.SUPPRESS,                  default=3,         choices=[0, 1, 2, 3, 4])
        parser.add_argument("--help",            "-h",              help=argparse.SUPPRESS,                  action="store_true"                       )
        parser.add_argument("--SKIP",                  type=str,    help=argparse.SUPPRESS,                  default='-'                               )
        parser.add_argument("--SAVE",                               help=argparse.SUPPRESS,                  action="store_true", default=False        )
        return parser.parse_args()
    
    def run(self):
        args                  = self.parser()
        args.highlight_color  = self.highlight_color
        args.background_color = self.background_color
        if args.help:
            print_help(args.highlight_color)
            sys.exit()
        if args.version:
            print(f"TreeForge v{VERSION}")
            sys.exit()
        print_args(args, args.highlight_color)
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
    #     '-f', '0.3',
    #     '-m', '10',
    #     '-o', '20',
    #     '-s', 'aa',
    #     '-l', '3'
    # ]
    main = TreeForge()
    main.run()