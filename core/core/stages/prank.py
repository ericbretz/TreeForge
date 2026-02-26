import os
import sys
import subprocess
from pathlib                import Path
from threading              import Lock
from core.treeutils.newick  import parse
from core.stages.base_stage import BaseStage
from core.utils.sublogger   import run_stage_subprocess, run_logged_subprocess
from core.utils.fasta_io    import read_fasta_sequences, write_single_sequence
from queue                  import Queue, Empty
from concurrent.futures     import ThreadPoolExecutor, as_completed


"""
Run prank, pxcslq, iqtree2. This will dynamically adjust the number of workers based on the number of threads.
it also reallocates workers after each stage to make things move faster.
"""

class Prank(BaseStage):
    def __init__(self,
                 dir_ortho1to1,
                 dir_prank,
                 infile_ending,
                 dna_aa,
                 concat_fasta,
                 threads,
                 log,
                 hc,
                 bc,
                 prank_pxclsq_threshold,
                 bootstrap_replicates,
                 subprocess_dir=None,
                 shared_printClass=None):
        super().__init__(log, hc, bc, threads, subprocess_dir, shared_printClass)
        self.dir_ortho1to1          = Path(dir_ortho1to1)
        self.dir_prank              = Path(dir_prank)
        self.infile_ending          = infile_ending
        self.dna_aa                 = dna_aa
        self.concat_fasta           = concat_fasta
        self.prank_pxclsq_threshold = prank_pxclsq_threshold
        self.bootstrap_replicates   = bootstrap_replicates

        self.prank_workers          = max(1, int(threads * 0.4))
        self.pxclsq_workers         = max(1, int(threads * 0.2))
        self.iqtree_workers         = max(1, int(threads * 0.4))
        self.iqtree_threads_per_job = 2

    def run(self):
        """
        Run PRANK.
        """
        self.write_fasta_from_tree()
        self.process_alignments()
        if os.path.exists('phyx.logfile'):
            os.remove('phyx.logfile')
        return self.return_dict

    def write_fasta_from_tree(self):
        """
        Write FASTA from tree.
        """
        self.printout('metric', 'Writing FASTA from tree')
        seq_mapping = {}
        for tree in self.dir_ortho1to1.glob('*.tre'):
            with open(tree, 'r') as infile:
                intree = parse(infile.readline())
            leaf_count = 0
            for node in intree.iternodes():
                if node.is_leaf():
                    leaf_count += 1
                    parts       = node.name.split('_')
                    if len(parts) >= 2:
                        if len(parts) == 2:
                            fasta_id = f"{parts[0]}@{parts[1]}"
                        else:
                            fasta_id = f"{parts[0]}_{parts[1]}@{'_'.join(parts[2:])}"
                        seq_mapping[fasta_id] = (tree.stem, node.name)
        created_files = set()
        
        sequences = read_fasta_sequences(self.concat_fasta)
        for fasta_id, sequence in sequences.items():
            if fasta_id in seq_mapping:
                tree_name, orig_name = seq_mapping[fasta_id]
                outfile_path = self.dir_prank / f"{tree_name}.fa"
                mode = 'w' if str(outfile_path) not in created_files else 'a'
                write_single_sequence(orig_name, sequence, outfile_path, mode=mode)
                if mode == 'w':
                    created_files.add(str(outfile_path))
        
        self.printout('metric', {'fasta_from_tree': len(created_files)})

    def run_prank_single(self, fa_file):
        """
        Run PRANK
        Single Threaded
        """
        out_base = fa_file.with_suffix('')
        cmd = f"prank -d={fa_file} -o={out_base}.aln"
        
        run_stage_subprocess(cmd, f'PRANK {fa_file.stem}', self.subprocess_dir, self.printout)
            
        aln_file = Path(f"{out_base}.aln.best.fas")
        new_name = aln_file.with_name(f"{aln_file.stem.replace('.aln.best', '')}.fas")
        aln_file.rename(new_name)
        return new_name

    def run_pxclsq_single(self, fas_file):
        """
        Run PXCLSQ.
        Single Threaded
        """
        out_file = fas_file.with_name(f"{fas_file.name}-cln")
        cmd = f"pxclsq -s {fas_file} -p {self.prank_pxclsq_threshold} -o {out_file}"
        
        run_stage_subprocess(cmd, f'PRANK pxclsq {fas_file.stem}', self.subprocess_dir, self.printout)
        return out_file

    def run_iqtree_single(self, cln_file):
        """
        Run IQ-TREE.
        Single Threaded
        """
        if self.bootstrap_replicates >= 1000:
            cmd = f"iqtree2 -s {cln_file} -nt {self.iqtree_threads_per_job} -bb {self.bootstrap_replicates} -redo -m GTR+G"
        else:
            cmd = f"iqtree2 -s {cln_file} -nt {self.iqtree_threads_per_job} -redo -m GTR+G"
        
        try:
            if self.subprocess_dir:
                result = run_logged_subprocess(cmd, self.subprocess_dir, f'prank_iqtree_{cln_file.stem}', shell=True, check=False)
                if result.returncode not in [0, 2]:
                    raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)
            else:
                iqtree = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                if iqtree.returncode not in [0, 2]:
                    raise subprocess.CalledProcessError(iqtree.returncode, cmd, iqtree.stdout, iqtree.stderr)
        except subprocess.CalledProcessError as e:
            self.printout('error', 'IQ-TREE2 failed')
            if e.stderr:
                self.printout('error', e.stderr.decode('utf-8'))
            sys.exit(1)
        return cln_file.with_suffix('.treefile')

    def run_optimized(self, files):
        """
        Run PRANK, PXCLSQ, IQ-TREE in parallel.
        """
        results          = []
        total_files      = len(files)
        completed_prank  = 0
        completed_pxclsq = 0
        completed_iqtree = 0
        prank_queue      = Queue()
        pxclsq_queue     = Queue()
        iqtree_queue     = Queue()
        
        progress_lock = Lock()
        
        def update_progress(stage, increment=1):
            nonlocal completed_prank, completed_pxclsq, completed_iqtree
            with progress_lock:
                if stage == 'prank':
                    completed_prank += increment
                    percentage       = f'{completed_prank/total_files*100:.1f}%'
                    self.printout('progress', f"Prank: {completed_prank}/{total_files} ({percentage})")
                elif stage == 'pxclsq':
                    completed_pxclsq += increment
                    percentage        = f'{completed_pxclsq/total_files*100:.1f}%'
                    self.printout('progress', f"Pxclsq: {completed_pxclsq}/{total_files} ({percentage})")
                elif stage == 'iqtree':
                    completed_iqtree += increment
                    percentage        = f'{completed_iqtree/total_files*100:.1f}%'
                    self.printout('progress', f"Iqtree2: {completed_iqtree}/{total_files} ({percentage})")

        for file in files:
            prank_queue.put(file)

        def prank_worker():
            while True:
                try:
                    fa_file = prank_queue.get(timeout=1)
                    if fa_file is None:
                        break
                    
                    fas_result = self.run_prank_single(fa_file)
                    pxclsq_queue.put(fas_result)
                    update_progress('prank')
                    
                except Empty:
                    break

        def pxclsq_worker():
            while True:
                try:
                    fas_file = pxclsq_queue.get(timeout=1)
                    if fas_file is None:
                        break
                    
                    cln_result = self.run_pxclsq_single(fas_file)
                    iqtree_queue.put(cln_result)
                    update_progress('pxclsq')
                    
                except Empty:
                    break

        def iqtree_worker():
            while True:
                try:
                    cln_file = iqtree_queue.get(timeout=1)
                    if cln_file is None:
                        break
                    
                    tree_result = self.run_iqtree_single(cln_file)
                    results.append(tree_result)
                    update_progress('iqtree')
                    
                except Empty:
                    break

        try:
            with ThreadPoolExecutor(max_workers=self.prank_workers + self.pxclsq_workers + self.iqtree_workers) as executor:
                prank_futures = [executor.submit(prank_worker) for _ in range(self.prank_workers)]
                
                for future in as_completed(prank_futures):
                    try:
                        future.result()
                    except Exception as e:
                        self.printout('error', f'Prank worker failed: {e}')
                        sys.exit(1)
                
                for _ in range(self.prank_workers):
                    prank_queue.put(None)
                
                pxclsq_futures = [executor.submit(pxclsq_worker) for _ in range(self.pxclsq_workers)]
                
                for future in as_completed(pxclsq_futures):
                    try:
                        future.result()
                    except Exception as e:
                        self.printout('error', f'Pxclsq worker failed: {e}')
                        sys.exit(1)
                
                iqtree_futures = [executor.submit(iqtree_worker) for _ in range(self.iqtree_workers)]
                
                for future in as_completed(iqtree_futures):
                    try:
                        future.result()
                    except Exception as e:
                        self.printout('error', f'Iqtree worker failed: {e}')
                        sys.exit(1)
        except KeyboardInterrupt:
            self.printout('error', 'Process interrupted by user')
            raise

        self.printout('metric', {'prank_trees_generated': len(results)})
        return results

    def _reallocate_resources_after_prank(self):
        """
        Reallocate resources after PRANK.
        Calculate new PXCLSQ and IQ-TREE workers based on available threads.
        """
        available_threads   = self.threads - self.prank_workers
        new_pxclsq_workers  = max(1, int(available_threads * 0.6))
        new_iqtree_workers  = max(1, available_threads - new_pxclsq_workers)
        self.pxclsq_workers = new_pxclsq_workers
        self.iqtree_workers = new_iqtree_workers

    def _reallocate_resources_after_pxclsq(self):
        available_threads   = self.threads - self.prank_workers - self.pxclsq_workers
        new_iqtree_workers  = max(1, available_threads // 2)
        self.iqtree_workers = new_iqtree_workers

    def process_alignments(self):
        fa_files = list(self.dir_prank.glob('*.fa'))        
        if not fa_files:
            self.printout('error', 'No FASTA files found for processing')
            return
        
        trees = self.run_optimized(fa_files)
        
        self.return_dict.update({
            'alignments': list(self.dir_prank.glob('*.fas')),
            'cleaned_alignments': list(self.dir_prank.glob('*-cln')),
            'trees': trees
        })