import os
import shutil
from pathlib import Path
from core.utils.printout import PrintOut
from core.treeutils.newick import parse
from core.treeutils.utils import get_front_labels, get_name

class Prune:
    def __init__(self,
                 dir_base,
                 dir_prune,
                 dir_trimmed,
                 dir_ortho,
                 orthocutoff,
                 threads,
                 log,
                 hc,
                 bc):
        
        self.dir_base    = dir_base
        self.dir_prune   = dir_prune
        self.dir_trimmed = dir_trimmed
        self.dir_ortho   = dir_ortho
        self.orthocutoff = orthocutoff
        self.threads     = threads
        self.log         = log
        self.hc          = hc
        self.bc          = bc

        self.subtree_files       = [Path(os.path.join(dir_trimmed, f)) for f in os.listdir(dir_trimmed) if f.endswith('.subtree')]
        self.relative_tip_cutoff = 0.2
        self.absolute_tip_cutoff = 0.3
        self.minimum_taxa        = orthocutoff
        self.threads             = threads
        self.outpath             = Path(os.path.join(self.dir_prune, '1to1ortho'))
        self.ortho1to1_files     = []
        self.outpath.mkdir(parents=True, exist_ok=True)

        self.printClass          = PrintOut(log, hc, bc)
        self.printout            = self.printClass.printout
        self.return_dict         = {}

    def run(self):
        self.printout('metric', 'Finding 1-to-1 orthologs')
        for st in self.subtree_files:
            with open(st, 'r') as infile:
                intree = parse(infile.readline())
            curroot = intree
            pp_trees = []
            if self.get_front_score(curroot) >= self.minimum_taxa:
                self.ortho1to1_files.append(st)

        self.write_ortho1to1_files()
        self.return_dict['prune'] = {
            'ortho1to1'      : self.ortho1to1_files,
            'ortho1to1_dir'  : self.outpath,
            'ortho1to1_count': len(self.ortho1to1_files),
        }
        return self.return_dict
    
    def write_ortho1to1_files(self):
        self.printout('metric', 'Writing Trees')
        for st in self.ortho1to1_files:
            outfile = os.path.join(self.outpath, Path(st).stem + '_1to1ortho.tre')
            shutil.copy(st, outfile)
            name = outfile if len(outfile) < 35 else  '...' + outfile[-35:]

    def get_front_score(self, node):
        front_labels = get_front_labels(node)
        num_labels = len(front_labels)
        num_taxa = len(set([get_name(i) for i in front_labels]))
        if num_taxa == num_labels:
            return num_taxa
        return -1