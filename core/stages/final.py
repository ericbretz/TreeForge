import os
import sys
import ete3
import shutil
import subprocess
from core.utils.printout import PrintOut
from core.utils.bes import BES

class Astral:
    def __init__(self,
                 dir_prank,
                 dir_super,
                 dir_treeforge,
                 renamed_map,
                 output_super_matrix,
                 super_bootstrap,
                 threads,
                 log,
                 hcolor,
                 bcolor,
                 bes_support):
        self.dir_prank           = dir_prank
        self.dir_super           = dir_super
        self.dir_treeforge       = dir_treeforge
        self.bes_support         = bes_support
        self.renamed_map         = renamed_map
        self.output_super_matrix = output_super_matrix
        self.super_bootstrap     = super_bootstrap
        self.threads             = threads
        self.log                 = log
        self.hcolor              = hcolor
        self.bcolor              = bcolor

        self.cln_files           = ' '.join([os.path.join(self.dir_prank, f) for f in os.listdir(self.dir_prank) if f.endswith('.treefile')])
        self.pxcat_files         = ' '.join([os.path.join(self.dir_prank, f) for f in os.listdir(self.dir_prank) if f.endswith('-cln')])
        self.concat_tree         = os.path.join(self.dir_super,     'concat.tre')
        self.finaltree           = os.path.join(self.dir_treeforge, 'SpeciesTree.coalescent.tre')
        self.supermatrix         = os.path.join(self.dir_super,     'SuperMatrix.treefile')
        self.genetree            = os.path.join(self.dir_treeforge, 'SuperMatrix.tre')
        self.dir_gene_trees      = os.path.join(self.dir_treeforge, 'gene_trees')
        self.s_matrix_dir        = os.path.join(self.dir_super,     'super.matrix')
        self.s_model_dir         = os.path.join(self.dir_super,     'super.model')
        self.return_dict         = {'concat_tree': self.concat_tree, 
                                    'gene_trees' : self.dir_gene_trees,
                                    'finaltree'  : self.finaltree,
                                    'num_trees'  : len(self.cln_files.split())}

        self.printClass = PrintOut(log, hcolor, bcolor)
        self.printout   = self.printClass.printout

    def run(self):
        self.file_check()
        self.tree_concat()
        self.astral()
        if self.output_super_matrix:
            self.pxcat()
            self.iqtree()
        self.bes()
        return self.return_dict
    
    def file_check(self):
        if len(self.cln_files.split()) == 0:
            self.printout('error', 'No gene trees found')
            self.printout('error', 'Try adjusting parameters')
            sys.exit(1)
        if len(self.pxcat_files.split()) == 0:
            self.printout('error', 'No PXCAT files found')
            self.printout('error', 'Try adjusting parameters')
            sys.exit(1)
    
    def tree_concat(self):
        """
        Concatenate gene trees and rename leaves to original file names.
        Saves both a directory of gene trees and a concatenated tree.
        """
        self.printout('metric', 'Concatenating gene trees')
        os.makedirs(self.dir_gene_trees, exist_ok=True)
        with open(self.concat_tree, 'w') as outfile:
            for file in self.cln_files.split():
                try:
                    tree = ete3.Tree(file)
                    tree_gene = tree.copy(method='deepcopy')
                    for node in tree.traverse():
                        if node.is_leaf():
                            lookup_key = node.name.replace('_', '@')
                            if lookup_key in self.renamed_map:
                                _, source_file = self.renamed_map[lookup_key]
                                node.name = source_file.split('.')[0]
                            else:
                                node.name = node.name.split('.')[0]
                    for node in tree_gene.traverse():
                        if node.is_leaf():
                            lookup_key = node.name.replace('_', '@')
                            if lookup_key in self.renamed_map:
                                _, source_file = self.renamed_map[lookup_key]
                                node.name = source_file.split('.')[0]
                            else:
                                node.name = node.name.split('.')[0]
                    tree_gene.write(format=1, outfile=os.path.join(self.dir_gene_trees, file.split('/')[-1].split('_')[0] + '.tre'))
                    outfile.write(tree.write(format=1) + '\n')
                except Exception as e:
                    self.printout('error', f'Failed to process tree file {file}: {str(e)}')
                    continue

    def astral(self):
        """
        Run ASTRAL.
        """
        self.printout('metric', 'Astral')
        cmd         = f'astral -i {self.concat_tree} -o {self.finaltree}'
        return_code = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if return_code.returncode != 0:
            self.printout('error', 'Astral failed')
            self.printout('error', return_code.stderr.decode('utf-8'))
            sys.exit(1)

    def pxcat(self):
        self.printout('metric', 'PXCAT')
        cmd         = f'pxcat -s {self.pxcat_files} -p {self.s_model_dir} -o {self.s_matrix_dir}'
        return_code = subprocess.run(cmd, shell=True)
        if return_code.returncode != 0:
            self.printout('error', 'PXCAT failed')
            self.printout('error', return_code.stderr.decode('utf-8'))
            sys.exit(1)
        if os.path.exists("phyx.logfile"):
            os.remove("phyx.logfile")

    def iqtree(self):
        self.printout('metric', 'SuperMatrix')
        file_dir    = os.path.join(self.dir_super, 'SuperMatrix')
        cmd         = f'iqtree2 -s {self.s_matrix_dir} -spp {self.s_model_dir} -nt {self.threads} -bb {self.super_bootstrap} -m GTR+G -pre {file_dir} -redo'
        return_code = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if return_code.returncode != 0:
            self.printout('error', 'IQ-TREE SuperMatrix failed')
            self.printout('error', return_code.stderr.decode('utf-8'))
            sys.exit(1)
        if os.path.exists(self.supermatrix):
            shutil.move(self.supermatrix, self.genetree)

    def bes(self):
        self.printout('metric', 'Convert Coalescent Distance to Molecular Distance')
        bes_runner = BES(printout=self.printout)
        bes_runner.run(self.finaltree, self.concat_tree, self.bes_support, self.dir_treeforge)