import os
import sys
import subprocess
from pathlib             import Path
from typing              import Dict, List, Optional
from core.stages.base_stage import BaseStage
from core.utils.sublogger import run_stage_subprocess, run_logged_subprocess
from concurrent.futures  import ThreadPoolExecutor, as_completed

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

class Mafft(BaseStage):
    def __init__(self,
                 dir_base,
                 dir_iteration,
                 threads,
                 log,
                 hc,
                 bc,
                 mafft_maxiter,
                 pxclsq_threshold,
                 thread_divisor,
                 prev_trimmed_tree_dir: Optional[str] = None,
                 use_prev_trees_as_start: bool = False,
                 fast_mode: bool = False,
                 subprocess_dir=None,
                 shared_printClass=None):

        super().__init__(log, hc, bc, threads, subprocess_dir, shared_printClass)
        self.dir_base           = Path(dir_base)
        self.dir_iteration      = Path(dir_iteration)
        self.mafft_maxiter      = mafft_maxiter
        self.pxclsq_threshold   = pxclsq_threshold
        self.thread_divisor     = thread_divisor
        self.fast_mode          = fast_mode  # For hierarchical clustering guide trees only

        self.prev_trimmed_tree_dir = Path(prev_trimmed_tree_dir) if prev_trimmed_tree_dir else None
        self.use_prev_trees_as_start = use_prev_trees_as_start

        self.fasta_files        = [str(self.dir_iteration / f) for f in os.listdir(self.dir_iteration) if f.endswith('.fa')]
        self.aln_files          = []
        self.cln_files          = []
        self.iqtree_files       = []
        self.total_time         = 0

    def run(self):
        """
        Run MAFFT.
        """
        self.max_workers     = max(1, min(self.threads // self.thread_divisor, len(self.fasta_files)))
        self.threads_per_job = max(self.thread_divisor, self.threads // max(1, self.max_workers))
        
        if not self.fasta_files:
            self.printout('info', 'No .fa files found to process')
            self.return_mafft({
                'aln_files'      : [],
                'cln_files'      : [],
                'iqtree_files'   : [],
                'files_processed': 0,
                'status'         : 'no_files_to_process'
            })
            return self.return_dict
        
        if os.path.exists('phyx.logfile'):
            os.remove('phyx.logfile')

        results = self.mafft_thread()
        self.return_mafft(results)
        return self.return_dict

    def mafft(self, fasta_file: str) -> Dict[str, List[str]]:
        """
        Run MAFFT.
        """
        aln_file = fasta_file.replace('.fa', '.aln')
        self.aln_files.append(aln_file)
        
        # Use faster settings for hierarchical clustering guide trees
        if self.fast_mode:
            cmd = f'mafft --retree 2 --maxiterate 0 --thread {self.threads_per_job} {fasta_file}'
        else:
            cmd = f'mafft --auto --maxiterate {self.mafft_maxiter} --thread {self.threads_per_job} {fasta_file}'
        
        try:
            if self.subprocess_dir:
                # For logged subprocess, capture output and write to file
                result = run_logged_subprocess(cmd, self.subprocess_dir, f'mafft_{Path(fasta_file).stem}', shell=True, check=True)
                with open(aln_file, 'wb') as f:
                    f.write(result.stdout)
            else:
                mafft = subprocess.run(cmd, stdout=open(aln_file, 'w'), stderr=subprocess.PIPE, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            self.printout('error', 'Mafft failed')
            if e.stderr:
                self.printout('error', e.stderr.decode('utf-8'))
            sys.exit(1)

        cln_file     = self.pxclsq(aln_file)
        iqtree_file  = self.iqtree(cln_file)

        return_dict = {
            'aln_files'   : [aln_file],
            'cln_files'   : [cln_file],
            'iqtree_files': [iqtree_file]
        }

        return return_dict
        
    def mafft_thread(self) -> Dict[str, List[str]]:
        """
        Run MAFFT in parallel.
        """
        combined_results = {
            'aln_files'   : [],
            'cln_files'   : [],
            'iqtree_files': []
        }
        
        total_files          = len(self.fasta_files)
        completed_files      = 0
        completed_percentage = f'{completed_files/total_files*100:.1f}%'
        self.printout('progress', f"MAFFT progress: {completed_files:0{len(str(total_files))}d}/{total_files} ({completed_percentage}) completed")
        try:
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                futures = {executor.submit(self.mafft, fasta_file): fasta_file for fasta_file in self.fasta_files}
                
                for future in as_completed(futures):
                    result           = future.result()
                    completed_files += 1
                    for key in combined_results:
                        combined_results[key].extend(result[key])
                    completed_percentage = f'{completed_files/total_files*100:.1f}%'
                    self.printout('progress', f"MAFFT progress: {completed_files:0{len(str(total_files))}d}/{total_files} ({completed_percentage}) completed")
        except KeyboardInterrupt:
            self.printout('error', 'Process interrupted by user')
            raise
        
        print()
        return combined_results
    
    def pxclsq(self, aln_file):
        """
        Run PXCLSQ.
        """
        cln_file = aln_file.replace('.aln', '.cln')
        self.cln_files.append(cln_file)
        cmd = f'pxclsq -s {aln_file} -p {self.pxclsq_threshold} -o {cln_file}'
        
        try:
            if self.subprocess_dir:
                run_logged_subprocess(cmd, self.subprocess_dir, f'pxclsq_{Path(aln_file).stem}', shell=True, check=True)
            else:
                subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            self.printout('error', 'Pxclsq failed')
            if e.stderr:
                self.printout('error', e.stderr.decode('utf-8'))
            sys.exit(1)
        
        return cln_file

    def iqtree(self, cln_file):
        """
        Run IQ-TREE.
        """
        iqtree_file = cln_file.replace('.cln', '.iqtree')
        self.iqtree_files.append(iqtree_file)

        def run_iqtree(cmd: str, log_suffix: str = '') -> subprocess.CompletedProcess:
            if self.subprocess_dir:
                return run_logged_subprocess(cmd, self.subprocess_dir, f'iqtree_{Path(cln_file).stem}{log_suffix}', shell=True, check=False)
            else:
                return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

        base_cmd = f'iqtree2 -s {cln_file} -nt {self.threads_per_job} -m GTR+G -redo -fast'

        starting_tree_used = False
        cmd_to_run = base_cmd

        if self.use_prev_trees_as_start and self.prev_trimmed_tree_dir is not None:
            cluster_id = Path(cln_file).stem
            prev_tree = self.prev_trimmed_tree_dir / f"{cluster_id}.subtree"

            if prev_tree.exists() and prev_tree.is_file() and prev_tree.stat().st_size > 0:
                cmd_to_run = f"{base_cmd} -t {prev_tree}"
                starting_tree_used = True
            else:
                self.printout('info', f'No previous tree found for {cluster_id}, running IQ-TREE without starting tree')

        iqtree = run_iqtree(cmd_to_run, '_with_start_tree' if starting_tree_used else '')

        if iqtree.returncode not in [0, 2] and starting_tree_used:
            self.printout('warning', f'IQ-TREE failed for {Path(cln_file).name} with starting tree, retrying without starting tree')
            iqtree = run_iqtree(base_cmd)

        if iqtree.returncode not in [0, 2]:
            self.printout('error', 'IQ-TREE failed')
            self.printout('error', iqtree.returncode)
            self.printout('error', iqtree.stderr.decode('utf-8'))
            return iqtree_file

        return iqtree_file
    
    def return_mafft(self, results):
        """
        Return MAFFT results.
        """
        self.return_dict['mafft'] = {
            'aln_files'   : results['aln_files'],
            'cln_files'   : results['cln_files'],
            'iqtree_files': results['iqtree_files'],
        }
        return