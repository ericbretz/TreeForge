import sys
import subprocess
import pysam
from core.utils.printout import PrintOut

class Blast:
    def __init__(self, dir_base, dir_treeforge, files_fasta, files_concatenated, files_raw_blast, threads, log, hc, bc):
        self.dir_base           = dir_base
        self.dir_treeforge      = dir_treeforge
        self.files_fasta        = files_fasta
        self.files_concatenated = files_concatenated
        self.files_raw_blast    = files_raw_blast
        self.threads            = threads
        self.contig_count       = 0
        self.return_dict = {'contig_count': self.contig_count}
        self.printClass = PrintOut(log, hc, bc)
        self.printout = self.printClass.printout

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
        lines = []
        for fasta in self.files_fasta:
            fasta = pysam.FastaFile(fasta)
            reference = fasta.references
            for i, name in enumerate(reference):
                seq = fasta.fetch(name)
                lines.append(f">{name}\n{seq}")
                self.contig_count += 1
        self.return_dict['contig_count'] = self.contig_count
        with open(self.files_concatenated, 'w') as f:
            f.write("\n".join(lines))
        return
    
    def make_blast_db(self):
        '''
        Build BLAST database from concatenated FASTA.
        '''
        self.printout('metric', 'Build Database')
        cmd   = f'makeblastdb -in {self.files_concatenated} -parse_seqids -dbtype nucl -out {self.files_concatenated}'
        blast = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if blast.returncode != 0:
            self.printout('error', 'BLAST makeblastdb failed')
            sys.exit(1)
        else:
            return
        
    def all_by_all(self):
        '''
        Run BLAST all-by-all search.
        '''
        self.printout('metric', 'All-by-All Comparison')
        cmd   = f'blastn -db {self.files_concatenated} -query {self.files_concatenated} -evalue 10 -num_threads {self.threads} -max_target_seqs 1000 -out {self.files_raw_blast} -outfmt "6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore"'
        blast = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if blast.returncode != 0:
            self.printout('error', 'BLAST blastn failed')
            sys.exit(1)
        else:
            return