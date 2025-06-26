import os
import sys
import subprocess
import ete3
from core.utils.printout import PrintOut

class SuperMatrix:
    def __init__(self,
                 dir_prank,
                 dir_super,
                 dir_tf,
                 threads,
                 log,
                 hcolor,
                 bcolor,
                 super_bootstrap):
        self.dir_prank = dir_prank
        self.dir_super = dir_super
        self.dir_tf    = dir_tf
        self.threads   = threads
        self.log       = log
        self.hcolor    = hcolor
        self.bcolor    = bcolor
        self.super_bootstrap = super_bootstrap

        self.intree          = os.path.join(self.dir_super, 'PrankAligned.treefile')
        self.outtree         = os.path.join(self.dir_super, 'SuperMatrix.tre')
        self.intree_renamed  = os.path.join(self.dir_super, 'PrankAlignedRenamed.tre')
        self.outtree_renamed = os.path.join(self.dir_super, 'SuperMatrixRenamed.tre')
        self.finaltree       = os.path.join(self.dir_tf,    'FinalTree.tre')
        self.s_matrix_dir    = os.path.join(self.dir_super, 'super.matrix')
        self.s_model_dir     = os.path.join(self.dir_super, 'super.model')
        self.cln_files       = ' '.join([os.path.join(self.dir_prank, f) for f in os.listdir(self.dir_prank) if f.endswith('.fas-cln')])
        self.return_dict     = {}

        self.printClass = PrintOut(log, hcolor, bcolor)
        self.printout   = self.printClass.printout

    def run(self):
        self.pxcat()
        self.iqtree()
        self.astral()
        return self.return_dict

    def pxcat(self):
        self.printout('metric', 'PXCAT')
        cmd         = f'pxcat -s {self.cln_files} -p {self.s_model_dir} -o {self.s_matrix_dir}'
        return_code = subprocess.run(cmd, shell=True)
        if return_code.returncode != 0:
            self.printout('error', 'PXCAT failed')
            self.printout('error', return_code.stderr.decode('utf-8'))
            sys.exit(1)

    def iqtree(self):
        self.printout('metric', 'IQ-TREE')
        file_dir    = os.path.join(self.dir_super, 'PrankAligned')
        cmd         = f'iqtree2 -s {self.s_matrix_dir} -spp {self.s_model_dir} -nt {self.threads} -bb {self.super_bootstrap} -m GTR+G -pre {file_dir} -redo'
        return_code = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if return_code.returncode != 0:
            self.printout('error', 'IQ-TREE failed')
            self.printout('error', return_code.stderr.decode('utf-8'))
            sys.exit(1)

    def astral(self):
        self.printout('metric', 'ASTRAL')
        cmd         = f'astral -i {self.intree} -o {self.finaltree}'
        return_code = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if return_code.returncode != 0:
            self.printout('error', 'ASTRAL failed')
            self.printout('error', return_code.stderr.decode('utf-8'))
            sys.exit(1)