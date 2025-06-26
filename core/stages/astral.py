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
        self.concat_tree     = os.path.join(self.dir_super, 'concat.tre')
        self.finaltree       = os.path.join(self.dir_treeforge, 'FinalTree.tre')
        self.dir_gene_trees  = os.path.join(self.dir_treeforge, 'gene_trees')
        self.return_dict     = {'concat_tree': self.concat_tree, 'gene_trees': self.dir_gene_trees, 'finaltree': self.finaltree, 'num_trees': len(self.cln_files.split())}

        self.printClass = PrintOut(log, hcolor, bcolor)
        self.printout   = self.printClass.printout

    def run(self):
        self.tree_concat()
        self.astral()
        return self.return_dict
    
    def tree_concat(self):
        """
        Concatenate gene trees and rename leaves to original file names.
        Saves both a directory of gene trees and a single concatenated tree.
        """
        self.printout('metric', 'Concatenating gene trees')
        unique_keys = {k.split('@')[0] for k in self.renamed_map.keys()}
        unqique_leafs = {}
        os.makedirs(self.dir_gene_trees, exist_ok=True)
        with open(self.concat_tree, 'w') as outfile:
            for file in self.cln_files.split():
                try:
                    tree = ete3.Tree(file)
                    tree_gene = tree.copy(method='deepcopy')
                    for node in tree.traverse():
                        if node.is_leaf():
                            key = str(node.name).split('_')[0]
                            if key not in unqique_leafs:
                                unqique_leafs[key] = node.name
                            matching_keys = [k for k in self.renamed_map.keys() if k.startswith(key)]
                            if matching_keys:
                                original_name, source_file = self.renamed_map[matching_keys[0]]
                                source_name = source_file.split('.')[0]
                                node.name = source_name
                            else:
                                node.name = node.name.split('.')[0]
                    for node in tree_gene.traverse():
                        if node.is_leaf():
                            node.name = self.renamed_map[node.name.replace('_', '@')][0]
                    tree_gene.write(format=1, outfile=os.path.join(self.dir_gene_trees, file.split('/')[-1].split('_')[0] + '.tre'))
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