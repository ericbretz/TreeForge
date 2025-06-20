import os
import sys
import subprocess
from typing import Dict, List
from core.utils.printout import PrintOut
from concurrent.futures import ThreadPoolExecutor, as_completed

"""
PARALLEL FLOW:

ThreadPoolExecutor
    │
    ├─ Worker 1: .fa → MAFFT → .aln → pxclsq → .cln → iqtree2 → .iqtree
    │
    ├─ Worker 2: .fa → MAFFT → .aln → pxclsq → .cln → iqtree2 → .iqtree
    │
    ├─ Worker 3: .fa → MAFFT → .aln → pxclsq → .cln → iqtree2 → .iqtree
    │
    └─ Worker N: .fa → MAFFT → .aln → pxclsq → .cln → iqtree2 → .iqtree
"""

class Mafft:
    def __init__(self,
                 dir_base,
                 dir_iteration,
                 threads,
                 log,
                 hc,
                 bc):

        self.dir_base           = dir_base
        self.dir_iteration      = dir_iteration
        self.threads            = threads
        self.log                = log
        self.hc                 = hc
        self.bc                 = bc
        self.printClass         = PrintOut(log, hc, bc)
        self.printout           = self.printClass.printout

        self.fasta_files        = [os.path.join(self.dir_iteration, f) for f in os.listdir(self.dir_iteration) if f.endswith('.fa')]
        self.aln_files          = []
        self.cln_files          = []
        self.iqtree_files       = []
        self.total_time         = 0
        self.return_dict        = {}

    def run(self):
        self.max_workers = min(self.threads // 4, len(self.fasta_files))
        self.threads_per_job = max(4, self.threads // max(1, self.max_workers))
        results = self.mafft_thread()
        self.return_mafft(results)
        return self.return_dict

    def mafft(self, fasta_file: str) -> Dict[str, List[str]]:

        aln_file = fasta_file.replace('.fa', '.aln')
        self.aln_files.append(aln_file)
        cmd = f'mafft --auto --maxiterate 1000 --thread {self.threads_per_job} {fasta_file}'
        mafft = subprocess.run(cmd, stdout=open(aln_file, 'w'), stderr=subprocess.PIPE, shell=True)
        if mafft.returncode != 0:
            self.printout('error', 'Mafft failed')
            self.printout('error', mafft.stderr.decode('utf-8'))
            sys.exit(1)
        else:
            pass

        cln_file     = self.pxclsq(aln_file)
        iqtree_file  = self.iqtree(cln_file)
        phyx_logfile = os.path.join(self.dir_base, "phyx.logfile")
        if os.path.exists(phyx_logfile): 
            os.remove(phyx_logfile)

        return_dict = {
            'aln_files': [aln_file],
            'cln_files': [cln_file],
            'iqtree_files': [iqtree_file]
        }

        return return_dict
        
    def mafft_thread(self) -> Dict[str, List[str]]:
        combined_results = {
            'aln_files': [],
            'cln_files': [],
            'iqtree_files': []
        }
        
        total_files = len(self.fasta_files)
        completed_files = 0
        self.printout('progress', f"MAFFT progress: {completed_files:0{len(str(total_files))}d}/{total_files} files completed")
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {executor.submit(self.mafft, fasta_file): fasta_file for fasta_file in self.fasta_files}
            
            for future in as_completed(futures):
                result = future.result()
                completed_files += 1
                fasta_file = futures[future]
                for key in combined_results:
                    combined_results[key].extend(result[key])
                self.printout('progress', f"MAFFT progress: {completed_files:0{len(str(total_files))}d}/{total_files} files completed")
        
        print()
        return combined_results
    
    def pxclsq(self, aln_file):
        cln_file = aln_file.replace('.aln', '.cln')
        self.cln_files.append(cln_file)
        cmd = f'pxclsq -s {aln_file} -p 0.1 -o {cln_file}'
        pxclsq = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if pxclsq.returncode != 0:
            self.printout('error', 'Pxclsq failed')
            self.printout('error', pxclsq.stderr.decode('utf-8'))
            sys.exit(1)
        else:
            pass
        return cln_file

    def iqtree(self, cln_file):
        iqtree_file = cln_file.replace('.cln', '.iqtree')
        self.iqtree_files.append(iqtree_file)
        cmd = f'iqtree2 -s {cln_file} -nt {self.threads_per_job} -m GTR+G -redo -fast'
        iqtree = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if iqtree.returncode not in [0,2]:
            self.printout('error', 'IQ-TREE failed')
            self.printout('error', iqtree.returncode)
            self.printout('error', iqtree.stderr.decode('utf-8'))
            return iqtree_file
        else:
            pass
        return iqtree_file
    
    def return_mafft(self, results):
        self.return_dict['mafft'] = {
            'aln_files'   : results['aln_files'],
            'cln_files'   : results['cln_files'],
            'iqtree_files': results['iqtree_files'],
        }
        return