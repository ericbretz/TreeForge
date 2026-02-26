import os
import sys
import shutil
from pathlib import Path
from core.utils.printout import PrintOut

class Deps:
    def __init__(self, log_level, hcolor, bcolor, nocolor=False):
        self.hcolor     = '' if nocolor else hcolor
        self.bcolor     = '' if nocolor else bcolor
        self.printClass = PrintOut(log_level, self.hcolor, self.bcolor)
        if nocolor:
            self.printClass.set_nocolor(True)
        self.printout   = self.printClass.printout

        self.deps = {
            'mafft'      : 'mafft',
            'mcl'        : 'mcl',
            'prank'      : 'prank',
            'iqtree2'    : 'iqtree2',
            # 'pxcat'      : 'pxcat',  # Replaced with concat.py
            'pxclsq'     : 'pxclsq',
            'blast'      : 'blastn',
            'blastx'     : 'blastx',
            'makeblastdb': 'makeblastdb',
            'astral'     : 'astral',
            'vsearch'    : 'vsearch',
            'diamond'    : 'diamond',
            'mmseqs'     : 'mmseqs',
        }

    def check_deps(self):
        self.printout('title', 'Checking dependencies')
        missing_deps = []
        
        for dep,cmd in self.deps.items():
            if not shutil.which(cmd):
                missing_deps.append(dep)
            else:
                self.printout('metric', {dep: 'found'})
        
        self.check_busco_database()
        
        if missing_deps:
            self.printout('error', f'Missing dependencies: {", ".join(missing_deps)}')
            sys.exit(1)
    
    def check_busco_database(self):
        """Check if BUSCO database file exists"""
        busco_db_path = Path(__file__).resolve().parent / 'euk_db.faa'
        
        if busco_db_path.exists():
            self.printout('metric', {'busco_database': 'found'})
        else:
            self.printout('warning', f'BUSCO database not found at {busco_db_path}')
            self.printout('warning', 'BUSCO functionality will not be available')
            self.printout('info', 'To enable BUSCO, download euk_db.faa and place it in core/utils/')

if __name__ == '__main__':
    deps = Deps(1, 'green', 'blue')
    deps.check_deps()