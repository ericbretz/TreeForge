import sys
import subprocess
from pathlib import Path
from core.stages.base_stage import BaseStage
from core.utils.sublogger import run_stage_subprocess
from core.utils.fasta_io import concatenate_fasta_files

class Blast(BaseStage):
    def __init__(self, dir_base, dir_treeforge, files_fasta, files_concatenated, files_raw_blast, threads, log, hc, bc, blast_evalue, blast_max_targets, subprocess_dir=None, shared_printClass=None):
        super().__init__(log, hc, bc, threads, subprocess_dir, shared_printClass)
        self.dir_base           = Path(dir_base)
        self.dir_treeforge      = Path(dir_treeforge)
        self.files_fasta        = files_fasta
        self.files_concatenated = files_concatenated
        self.files_raw_blast    = files_raw_blast
        self.blast_evalue       = blast_evalue
        self.blast_max_targets  = blast_max_targets
        self.contig_count       = 0
        self.return_dict        = {'contig_count': self.contig_count}

    def run(self):
        self.fasta_concatenation()
        self.make_blast_db()
        self.all_by_all()
        return self.return_dict

    def fasta_concatenation(self):
        '''
        Concatenate all FASTA files into one file.
        '''
        self.printout('metric', 'FASTA Concatenation')
        self.contig_count = concatenate_fasta_files(self.files_fasta, self.files_concatenated)
        self.return_dict['contig_count'] = self.contig_count
        return
    
    def make_blast_db(self):
        '''
        Build BLAST database from concatenated FASTA.
        '''
        self.printout('metric', 'Build Database')
        cmd = f'makeblastdb -in {self.files_concatenated} -parse_seqids -dbtype nucl -out {self.files_concatenated}'
        run_stage_subprocess(cmd, 'BLAST makeblastdb', self.subprocess_dir, self.printout)
        return
        
    def all_by_all(self):
        '''
        Run BLAST all-by-all search.
        '''
        self.printout('metric', 'All-by-All Comparison')
        cmd = f'blastn -db {self.files_concatenated} -query {self.files_concatenated} -evalue {self.blast_evalue} -num_threads {self.threads} -max_target_seqs {self.blast_max_targets} -out {self.files_raw_blast} -outfmt "6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore"'
        run_stage_subprocess(cmd, 'BLAST blastn', self.subprocess_dir, self.printout)
        return