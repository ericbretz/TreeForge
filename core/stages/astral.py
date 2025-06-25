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
                 renamed_map,
                 threads,
                 log,
                 hcolor,
                 bcolor):
        self.dir_prank     = dir_prank
        self.dir_super     = dir_super
        self.dir_treeforge = dir_treeforge
        self.renamed_map   = renamed_map
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
        Concatenate gene trees and rename leaves to original file names.
        """
        self.printout('metric', 'Concatenating gene trees')
        unique_keys = {k.split('@')[0] for k in self.renamed_map.keys()}
        unqique_leafs = {}
        with open(self.concat_tree, 'w') as outfile:
            for file in self.cln_files.split():
                try:
                    tree = ete3.Tree(file)
                    for node in tree.traverse():
                        if node.is_leaf():
                            key = str(node.name).split('.')[0]
                            if key not in unqique_leafs:
                                unqique_leafs[key] = node.name
                            matching_keys = [k for k in self.renamed_map.keys() if k.startswith(key)]
                            if matching_keys:
                                original_name, source_file = self.renamed_map[matching_keys[0]]
                                source_name = source_file.split('.')[0]
                                node.name = source_name
                            else:
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