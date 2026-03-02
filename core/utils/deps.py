import os
import sys
import shutil
from pathlib import Path
from core.utils.printout import PrintOut

class Deps:
    def __init__(self, log_level, hcolor, bcolor, nocolor=False, astral_jar=None, hcluster_tool=None):
        self.hcolor       = '' if nocolor else hcolor
        self.bcolor       = '' if nocolor else bcolor
        self.printClass   = PrintOut(log_level, self.hcolor, self.bcolor)
        if nocolor:
            self.printClass.set_nocolor(True)
        self.printout     = self.printClass.printout
        self.astral_jar   = astral_jar
        self.hcluster_tool = hcluster_tool

        self.deps = {
            'mafft'      : 'mafft',
            'mcl'        : 'mcl',
            'prank'      : 'prank',
            'iqtree2'    : 'iqtree2',
            # 'pxcat'      : 'pxcat',  # Replaced with concat.py
            'pxclsq'     : 'pxclsq',
            # 'blast'      : 'blastn', # replaced this with mmseqs2 for all-by-all
            # 'blastx'     : 'blastx',
            # 'makeblastdb': 'makeblastdb',
            # 'diamond'    : 'diamond', # replaced with mmseqs
            'mmseqs'     : 'mmseqs',
        }

    def check_deps(self):
        self.printout('title', 'Checking dependencies')
        missing_deps = []
        
        for dep, cmd in self.deps.items():
            if not shutil.which(cmd):
                missing_deps.append(dep)
            else:
                self.printout('metric', {dep: 'found'})

        if self.hcluster_tool == 'vsearch':
            if not shutil.which('vsearch'):
                missing_deps.append('vsearch')
            else:
                self.printout('metric', {'vsearch': 'found'})

        self._check_astral(missing_deps)
        if self.hcluster_tool is not None:
            self.check_busco_database()
        
        if missing_deps:
            self.printout('error', f'Missing dependencies: {", ".join(missing_deps)}')
            sys.exit(1)

    def _check_astral(self, missing_deps: list) -> None:
        if self.astral_jar:
            jar_path = Path(self.astral_jar)
            if jar_path.is_file() and shutil.which('java'):
                self.printout('metric', {'astral': 'found'})
                self.printout('metric', {'java': 'found'})
                return
            if not jar_path.is_file():
                self.printout('warning', f'ASTRAL not found: {self.astral_jar}')
            elif not shutil.which('java'):
                self.printout('warning', 'java not found')
        if shutil.which('astral'):
            self.printout('metric', {'astral': 'found'})
            self.astral_jar = None
        else:
            missing_deps.append('astral')
    
    def check_busco_database(self):
        """Check if BUSCO database file exists"""
        busco_db_path = Path(__file__).resolve().parent / 'euk_db.faa'
        
        if busco_db_path.exists():
            self.printout('metric', {'busco db': 'found'})
        else:
            self.printout('warning', f'BUSCO database not found at {busco_db_path}')
            self.printout('warning', 'BUSCO functionality will not be available')
            sys.exit(1)
            # self.printout('info', 'To enable BUSCO, download euk_db.faa and place it in core/utils/')

if __name__ == '__main__':
    deps = Deps(1, 'green', 'blue')
    deps.check_deps()