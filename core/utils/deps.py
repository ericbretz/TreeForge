import os
import sys
import shutil
import importlib
from core.utils.printout import PrintOut

class Deps:
    def __init__(self, log_level, hcolor, bcolor):
        self.hcolor    = hcolor
        self.bcolor    = bcolor
        self.printClass = PrintOut(log_level, hcolor, bcolor)
        self.printout   = self.printClass.printout

        self.deps = {
            'mafft'  : 'mafft',
            'mcl'    : 'mcl',
            'blast'  : 'blastn',
            'prank'  : 'prank',
            'iqtree2': 'iqtree2',
            'pxcat'  : 'pxcat',
            'pxclsq' : 'pxclsq',
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
        self.check_pip()

    def check_pip(self):
        requirements_file = os.path.join(os.path.dirname(__file__), '..', '..', 'requirements.txt')
        self.printout('title', 'Checking pip dependencies')
        
        if not os.path.exists(requirements_file):
            self.printout('error', f'Requirements file not found: {requirements_file}')
            sys.exit(1)
        
        missing_deps = []
        
        # Map package names to their actual import names
        package_import_map = {
            'biopython': 'Bio'
        }
        
        with open(requirements_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    package_name = line.split('>=')[0].split('<=')[0].split('==')[0].split('!=')[0].split('~=')[0].strip()
                    
                    # Use the import name if it's different from the package name
                    import_name = package_import_map.get(package_name, package_name)
                    
                    try:
                        importlib.import_module(import_name)
                        self.printout('metric', {package_name: 'found'})
                    except ImportError:
                        missing_deps.append(package_name)
        
        if missing_deps:
            self.printout('error', f'Missing pip dependencies: {", ".join(missing_deps)}')
            sys.exit(1)

if __name__ == '__main__':
    deps = Deps(1, 'green', 'blue')
    deps.check_deps()
    deps.check_pip()