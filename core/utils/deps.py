import os
import sys
import shutil
from core.utils.printout import PrintOut

class Deps:
    def __init__(self, log_level, hcolor, bcolor):
        self.hcolor     = hcolor
        self.bcolor     = bcolor
        self.printClass = PrintOut(log_level, hcolor, bcolor)
        self.printout   = self.printClass.printout

        self.deps = {
            'mafft'  : 'mafft',
            'mcl'    : 'mcl',
            'prank'  : 'prank',
            'iqtree2': 'iqtree2',
            'pxcat'  : 'pxcat',
            'pxclsq' : 'pxclsq',
            'blast'  : 'blastn',
        }

    def check_deps(self):
        self.printout('title', 'Checking dependencies')
        missing_deps = []
        
        for dep,cmd in self.deps.items():
            if not shutil.which(cmd):
                missing_deps.append(dep)
            else:
                self.printout('metric', {dep: 'found'})
        
        if missing_deps:
            self.printout('error', f'Missing dependencies: {", ".join(missing_deps)}')
            sys.exit(1)

if __name__ == '__main__':
    deps = Deps(1, 'green', 'blue')
    deps.check_deps()