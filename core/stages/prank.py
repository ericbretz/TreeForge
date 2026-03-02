import os
import sys
import subprocess
from pathlib                import Path
from core.treeutils.newick  import parse
from core.stages.base_stage import BaseStage
from core.utils.sublogger   import run_stage_subprocess, run_logged_subprocess
from core.utils.fasta_io    import read_fasta_sequences, write_single_sequence
from concurrent.futures     import ThreadPoolExecutor, as_completed

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
                 subprocess_dir    = None,
                 shared_printClass = None):
        super().__init__(log, hc, bc, threads, subprocess_dir, shared_printClass)
        self.dir_ortho1to1          = Path(dir_ortho1to1)
        self.dir_prank              = Path(dir_prank)
        self.infile_ending          = infile_ending
        self.dna_aa                 = dna_aa
        self.concat_fasta           = concat_fasta
        self.prank_pxclsq_threshold = prank_pxclsq_threshold
        self.bootstrap_replicates   = bootstrap_replicates

        self.iqtree_threads_per_job = 2
        self.prank_workers          = threads
        self.pxclsq_workers         = threads
        self.iqtree_workers         = max(1, threads // self.iqtree_threads_per_job)

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
            for node in intree.iternodes():
                if node.is_leaf():
                    parts = node.name.split('_')
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
        seqtype = '-protein' if self.dna_aa == 'aa' else '-DNA'
        out_base = fa_file.with_suffix('')
        cmd = f"prank -d={fa_file} -o={out_base}.aln {seqtype}"
        
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
        iqtree_model = 'LG+G' if self.dna_aa == 'aa' else 'GTR+G'
        if self.bootstrap_replicates >= 1000:
            cmd = f"iqtree2 -s {cln_file} -nt {self.iqtree_threads_per_job} -bb {self.bootstrap_replicates} -redo -m {iqtree_model}"
        else:
            cmd = f"iqtree2 -s {cln_file} -nt {self.iqtree_threads_per_job} -redo -m {iqtree_model}"
        
        try:
            if self.subprocess_dir:
                result = run_logged_subprocess(cmd, self.subprocess_dir, f'prank_iqtree_{cln_file.stem}', shell=True, check=False)
            else:
                result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if result.returncode not in [0, 2]:
                raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)
        except subprocess.CalledProcessError as e:
            self.printout('error', 'IQ-TREE2 failed')
            if e.stderr:
                self.printout('error', e.stderr.decode('utf-8'))
            sys.exit(1)
        return cln_file.with_suffix('.treefile')

    def run_optimized(self, files):
        """
        Run PRANK, PXCLSQ, IQ-TREE sequentially
        """
        def _run_phase(fn, items, n_workers, label):
            total     = len(items)
            completed = 0
            results   = []
            with ThreadPoolExecutor(max_workers=n_workers) as executor:
                futures = {executor.submit(fn, item): item for item in items}
                for future in as_completed(futures):
                    try:
                        results.append(future.result())
                        completed += 1
                        self.printout('progress', f"{label}: {completed}/{total} ({completed/total*100:.1f}%)")
                    except Exception as e:
                        self.printout('error', f'{label} worker failed: {e}')
                        sys.exit(1)
            return results

        try:
            fas_files = _run_phase(self.run_prank_single,  files,     self.prank_workers,  'Prank')
            cln_files = _run_phase(self.run_pxclsq_single, fas_files, self.pxclsq_workers, 'Pxclsq')
            trees     = _run_phase(self.run_iqtree_single, cln_files, self.iqtree_workers,  'Iqtree')
        except KeyboardInterrupt:
            self.printout('error', 'Process interrupted by user')
            raise

        self.printout('metric', {'prank_trees_generated': len(trees)})
        return trees

    def process_alignments(self):
        fa_files = list(self.dir_prank.glob('*.fa'))        
        if not fa_files:
            self.printout('error', 'No FASTA files found for processing')
            return
        
        trees = self.run_optimized(fa_files)
        
        self.return_dict.update({
            'alignments'        : list(self.dir_prank.glob('*.fas')),
            'cleaned_alignments': list(self.dir_prank.glob('*-cln')),
            'trees'             : trees
        })