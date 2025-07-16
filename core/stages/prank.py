import os
import sys
import pysam
import subprocess
from pathlib                import Path
from threading              import Lock
from core.treeutils.newick  import parse
from core.utils.printout    import PrintOut
from queue                  import Queue, Empty
from concurrent.futures     import ThreadPoolExecutor, as_completed


"""
Run prank, pxcslq, iqtree2. This will dynamically adjust the number of workers based on the number of threads.
it also reallocates workers after each stage to make things move faster.
"""

class Prank:
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
                 bootstrap_replicates):
        self.dir_ortho1to1          = Path(dir_ortho1to1)
        self.dir_prank              = Path(dir_prank)
        self.infile_ending          = infile_ending
        self.dna_aa                 = dna_aa
        self.concat_fasta           = concat_fasta
        self.threads                = threads
        self.log                    = log
        self.hc                     = hc
        self.bc                     = bc
        self.prank_pxclsq_threshold = prank_pxclsq_threshold
        self.bootstrap_replicates   = bootstrap_replicates
        self.return_dict            = {}

        self.prank_workers          = max(1, int(threads * 0.4))
        self.pxclsq_workers         = max(1, int(threads * 0.2))
        self.iqtree_workers         = max(1, int(threads * 0.4))
        self.iqtree_threads_per_job = 2

        self.printClass             = PrintOut(log, hc, bc)
        self.printout               = self.printClass.printout

    def run(self):
        self.write_fasta_from_tree()
        self.process_alignments()
        return self.return_dict

    def write_fasta_from_tree(self):
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
        
        with pysam.FastxFile(self.concat_fasta) as fasta:
            for entry in fasta:
                if entry.name in seq_mapping:
                    tree_name, orig_name = seq_mapping[entry.name]
                    outfile_path         = os.path.join(self.dir_prank, f"{tree_name}.fa")
                    mode = 'w' if outfile_path not in created_files else 'a'
                    with open(outfile_path, mode) as outfile:
                        outfile.write(f">{orig_name}\n{entry.sequence}\n")
                    if mode == 'w':
                        created_files.add(outfile_path)
        
        self.printout('metric', f'Created {len(created_files)} FASTA files from tree data')

    def run_prank_single(self, fa_file):
        out_base = fa_file.with_suffix('')
        cmd      = f"prank -d={fa_file} -o={out_base}.aln"
        prank    = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        
        if prank.returncode != 0:
            self.printout('error', 'PRANK failed')
            self.printout('error', prank.stderr.decode('utf-8'))
            sys.exit(1)
            
        aln_file = Path(f"{out_base}.aln.best.fas")
        new_name = aln_file.with_name(f"{aln_file.stem.replace('.aln.best', '')}.fas")
        aln_file.rename(new_name)
        return new_name

    def run_pxclsq_single(self, fas_file):
        out_file = fas_file.with_name(f"{fas_file.name}-cln")
        cmd      = f"pxclsq -s {fas_file} -p {self.prank_pxclsq_threshold} -o {out_file}"
        pxclsq   = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        
        if pxclsq.returncode != 0:
            self.printout('error', 'Pxclsq failed')
            sys.exit(1)
        return out_file

    def run_iqtree_single(self, cln_file):
        cmd    = f"iqtree2 -s {cln_file} -nt {self.iqtree_threads_per_job} -bb {self.bootstrap_replicates} -redo -m GTR+G"
        iqtree = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        
        if iqtree.returncode not in [0,2]:
            self.printout('error', 'IQ-TREE2 failed')
            self.printout('error', iqtree.stderr.decode('utf-8'))
            sys.exit(1)
        return cln_file.with_suffix('.treefile')

    def run_optimized(self, files):
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

        self.printout('metric', f'Generated {len(results)} trees')
        return results

    def _reallocate_resources_after_prank(self):
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