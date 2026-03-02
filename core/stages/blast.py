import sys
from pathlib import Path
from core.stages.base_stage import BaseStage
from core.utils.sublogger import run_stage_subprocess
from core.utils.fasta_io import concatenate_fasta_files

_SEQTYPE_TO_MMSEQS = {'aa': 'prot', 'nuc': 'nucl'}
_MMSEQS_TO_SEQTYPE = {'prot': 'aa',  'nucl': 'nuc'}

class Blast(BaseStage):
    def __init__(self, 
        dir_base, 
        dir_treeforge, 
        files_fasta, 
        files_concatenated, 
        files_raw_blast, 
        threads, 
        log, 
        hc, 
        bc, 
        blast_evalue, 
        blast_max_targets, 
        seqtype           = None,
        subprocess_dir    = None,
        shared_printClass = None):
        super().__init__(log, hc, bc, threads, subprocess_dir, shared_printClass)
        self.dir_base           = Path(dir_base)
        self.dir_treeforge      = Path(dir_treeforge)
        self.files_fasta        = files_fasta
        self.files_concatenated = Path(files_concatenated)
        self.files_raw_blast    = files_raw_blast
        self.blast_evalue       = blast_evalue
        self.blast_max_targets  = blast_max_targets
        self.seqtype_override   = seqtype
        self.seq_type           = None
        self.contig_count       = 0
        self.return_dict        = {'contig_count': self.contig_count, 'seqtype': None}

    def run(self):
        self.fasta_concatenation()
        if self.seqtype_override is not None:
            self.seq_type = _SEQTYPE_TO_MMSEQS[self.seqtype_override]
        else:
            self.detect_seq_type()
        self.return_dict['seqtype'] = _MMSEQS_TO_SEQTYPE[self.seq_type]
        self.printout('metric', [('seqtype', self.return_dict['seqtype'])])
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

    def detect_seq_type(self):
        '''
        Detect whether input sequences are protein or nucleotide.
        '''
        PROTEIN_ONLY    = set('EFILPQZ')
        NUCLEOTIDE_ONLY = set('U')
        seq_chars       = []
        max_sample      = 50000

        with open(self.files_concatenated, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('>'):
                    continue
                seq_chars.extend(c.upper() for c in line if c.isalpha())
                if len(seq_chars) >= max_sample:
                    break

        if not seq_chars:
            self.printout('error', 'FASTA file does not contain sequence data.')
            sys.exit(1)

        sample = set(seq_chars[:max_sample])
        if sample & PROTEIN_ONLY:
            self.seq_type = 'prot'
        elif sample & NUCLEOTIDE_ONLY:
            self.seq_type = 'nucl'
        else:
            self.seq_type = 'nucl'

    def all_by_all(self):
        '''
        Run all-by-all comparison
        '''
        self.printout('metric', 'All-by-All Comparison')
        search_type = '1' if self.seq_type == 'prot' else '3'
        tmp_dir     = self.files_concatenated.parent / 'mmseqs_tmp'
        cmd = [
            'mmseqs', 'easy-search',
            str(self.files_concatenated),
            str(self.files_concatenated),
            str(self.files_raw_blast),
            str(tmp_dir),
            '--search-type',   search_type,
            '--threads',       str(self.threads),
            '-e',              str(self.blast_evalue),
            '--max-seqs',      str(self.blast_max_targets),
            '--format-output', 'query,qlen,target,tlen,qframe,pident,nident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits',
            '--format-mode',   '0',
        ]
        run_stage_subprocess(cmd, 'MMseqs2 easy-search', self.subprocess_dir, self.printout)
        return
