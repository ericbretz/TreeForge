import os
import sys
import subprocess
import ete3
from core.utils.printout import PrintOut

class Astral:
    def __init__(self,
                 dir_prank,
                 dir_super,
                 dir_treeforge,
                 threads,
                 log,
                 hcolor,
                 bcolor):
        self.dir_prank     = dir_prank
        self.dir_super     = dir_super
        self.dir_treeforge = dir_treeforge
        self.threads       = threads
        self.log           = log
        self.hcolor        = hcolor
        self.bcolor        = bcolor

        self.cln_files       = ' '.join([os.path.join(self.dir_prank, f) for f in os.listdir(self.dir_prank) if f.endswith('.treefile')])
        self.concat_tree     = os.path.join(self.dir_super, 'concat.tree')
        self.finaltree       = os.path.join(self.dir_treeforge, 'FinalTree.tree')
        self.return_dict     = {'concat_tree': self.concat_tree, 'finaltree': self.finaltree, 'num_trees': len(self.cln_files.split())}

        self.printClass = PrintOut(log, hcolor, bcolor)
        self.printout   = self.printClass.printout

    def run(self):
        self.tree_concat()
        self.astral()
        return self.return_dict
    
    def tree_concat(self):
        """
        Concatenate gene trees.
        """
        self.printout('metric', 'Concatenating gene trees')
        with open(self.concat_tree, 'w') as outfile:
            for file in self.cln_files.split():
                try:
                    tree = ete3.Tree(file)
                    for node in tree.traverse():
                        if node.is_leaf():
                            if '.' in node.name:
                                node.name = node.name.split('.')[0]
                    outfile.write(tree.write(format=1) + '\n')
                except Exception as e:
                    self.printout('error', f'Failed to process tree file {file}: {str(e)}')
                    continue

    def astral(self):
        """
        Run ASTRAL.
        """
        self.printout('metric', 'ASTRAL')
        cmd         = f'astral -i {self.concat_tree} -o {self.finaltree}'
        return_code = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if return_code.returncode != 0:
            self.printout('error', 'ASTRAL failed')
            sys.exit(1)