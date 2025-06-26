import os
import sys
from pathlib import Path
import pysam
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from concurrent.futures import wait
from core.treeutils.newick import parse
from core.utils.printout import PrintOut

"""
PARALLEL FLOW:

ThreadPoolExecutor
    │
    ├─ Worker 1: .fa → PRANK → .fas → pxclsq → -cln → iqtree2 → .treefile
    │
    ├─ Worker 2: .fa → PRANK → .fas → pxclsq → -cln → iqtree2 → .treefile
    │
    ├─ Worker 3: .fa → PRANK → .fas → pxclsq → -cln → iqtree2 → .treefile
    │
    └─ Worker N: .fa → PRANK → .fas → pxclsq → -cln → iqtree2 → .treefile
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
        self.dir_ortho1to1 = Path(dir_ortho1to1)
        self.dir_prank     = Path(dir_prank)
        self.infile_ending = infile_ending
        self.dna_aa        = dna_aa
        self.concat_fasta  = concat_fasta
        self.threads       = threads
        self.log           = log
        self.hc            = hc
        self.bc            = bc
        self.prank_pxclsq_threshold = prank_pxclsq_threshold
        self.bootstrap_replicates    = bootstrap_replicates
        self.return_dict   = {}

        self.printClass    = PrintOut(log, hc, bc)
        self.printout      = self.printClass.printout

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
                    parts = node.name.split('_')
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
                    outfile_path = os.path.join(self.dir_prank, f"{tree_name}.fa")
                    mode = 'w' if outfile_path not in created_files else 'a'
                    # if '_' in orig_name:
                    #     name_parts = orig_name.rsplit('_', 1)
                    #     orig_name = name_parts[0]
                    # elif '@' in orig_name:
                    #     orig_name = orig_name.split('@')[0]
                    with open(outfile_path, mode) as outfile:
                        outfile.write(f">{orig_name}\n{entry.sequence}\n")
                    if mode == 'w':
                        created_files.add(outfile_path)

    def run_parallel(self, func, files, max_workers, step_name):
        results         = []
        total_files     = len(files)
        completed_files = 0
        
        self.printout('metric', step_name)
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(func, f): f for f in files}
            
            for future in as_completed(futures):
                result           = future.result()
                completed_files += 1
                input_file       = futures[future]
                cluster_name     = input_file.stem
                self.printout('metric', {'cluster': cluster_name.split("_")[0], 'completed': completed_files, 'total': total_files})
                results.append(result)
        return results

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
        threads_per_job = max(1, self.threads // 4)
        cmd             = f"iqtree2 -s {cln_file} -nt {threads_per_job} -bb {self.bootstrap_replicates} -redo -m GTR+G"
        iqtree          = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        
        if iqtree.returncode not in [0,2]:
            self.printout('error', 'IQ-TREE failed')
            self.printout('error', iqtree.stderr.decode('utf-8'))
            sys.exit(1)
        return cln_file.with_suffix('.treefile')

    def run_staged_parallel(self, files, max_workers):
        results = []
        total_files = len(files)
        completed_prank = 0
        completed_pxclsq = 0
        completed_iqtree = 0
        
        total_steps = total_files * 3
        completed_steps = 0

        self.printout('progress', f"PRANK: {completed_steps:0{len(str(total_steps))}d}/{total_steps} (P:{completed_prank:0{len(str(total_files))}d}/{total_files} px:{completed_pxclsq:0{len(str(total_files))}d}/{total_files} IQ:{completed_iqtree:0{len(str(total_files))}d}/{total_files})")
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            prank_futures = {
                executor.submit(self.run_prank_single, f): f 
                for f in files
            }
            
            running_pxclsq = set()
            running_iqtree = set()
            pxclsq_futures = {}
            iqtree_futures = {}
            
            while prank_futures or pxclsq_futures or iqtree_futures:
                done, _ = wait(
                    list(prank_futures.keys()) + 
                    list(pxclsq_futures.keys()) + 
                    list(iqtree_futures.keys()),
                    return_when='FIRST_COMPLETED'
                )
                
                for future in done:
                    if future in prank_futures:
                        input_file = prank_futures[future]
                        cluster_name = input_file.stem.split("_")[0]
                        prank_result = future.result()
                        completed_prank += 1
                        completed_steps += 1
                        
                        self.printout('progress', f"{completed_steps}/{total_steps} (P:{completed_prank}/{total_files} px:{completed_pxclsq}/{total_files} IQ:{completed_iqtree}/{total_files})")
                        
                        if len(running_pxclsq) < max_workers:
                            pxclsq_future = executor.submit(self.run_pxclsq_single, prank_result)
                            pxclsq_futures[pxclsq_future] = prank_result
                            running_pxclsq.add(prank_result)
                        
                        del prank_futures[future]
                        
                    elif future in pxclsq_futures:
                        input_file = pxclsq_futures[future]
                        cluster_name = input_file.stem.split("_")[0]
                        pxclsq_result = future.result()
                        completed_pxclsq += 1
                        completed_steps += 1
                        
                        self.printout('progress', f"{completed_steps}/{total_steps} (P:{completed_prank}/{total_files} px:{completed_pxclsq}/{total_files} IQ:{completed_iqtree}/{total_files})")
                        
                        running_pxclsq.remove(input_file)
                        
                        if len(running_iqtree) < max(1, self.threads // 4):
                            iqtree_future = executor.submit(self.run_iqtree_single, pxclsq_result)
                            iqtree_futures[iqtree_future] = pxclsq_result
                            running_iqtree.add(pxclsq_result)
                        
                        del pxclsq_futures[future]
                        
                    elif future in iqtree_futures:
                        input_file = iqtree_futures[future]
                        cluster_name = input_file.stem.split("_")[0]
                        tree_result = future.result()
                        completed_iqtree += 1
                        completed_steps += 1
                        
                        self.printout('progress', f"{completed_steps}/{total_steps} (P:{completed_prank}/{total_files} px:{completed_pxclsq}/{total_files} IQ:{completed_iqtree}/{total_files})")
                        
                        running_iqtree.remove(input_file)
                        results.append(tree_result)
                        del iqtree_futures[future]
                
                while len(running_pxclsq) < max_workers and completed_prank > len(pxclsq_futures) + completed_pxclsq:
                    for prank_result in [r for r in results if r not in running_pxclsq and r not in running_iqtree]:
                        pxclsq_future = executor.submit(self.run_pxclsq_single, prank_result)
                        pxclsq_futures[pxclsq_future] = prank_result
                        running_pxclsq.add(prank_result)
                        break
                
                iqtree_max = max(1, self.threads // 4)
                while len(running_iqtree) < iqtree_max and completed_pxclsq > len(iqtree_futures) + completed_iqtree:
                    for pxclsq_result in [r for r in results if r not in running_iqtree]:
                        iqtree_future = executor.submit(self.run_iqtree_single, pxclsq_result)
                        iqtree_futures[iqtree_future] = pxclsq_result
                        running_iqtree.add(pxclsq_result)
                        break
        return results

    def process_alignments(self):
        fa_files    = list(self.dir_prank.glob('*.fa'))
        max_workers = max(1, self.threads - 2)
        trees       = self.run_staged_parallel(fa_files, max_workers)
        
        self.return_dict.update({
            'alignments': list(self.dir_prank.glob('*.fas')),
            'cleaned_alignments': list(self.dir_prank.glob('*-cln')),
            'trees': trees
        })