import os
import sys
import pysam
import subprocess
import numpy as np
from pathlib             import Path
from typing              import Dict, Tuple
from core.utils.printout import PrintOut

class MCL:
    def __init__(self,
                 dir_base,
                 dir_treeforge,
                 dir_mcl,
                 dir_iteration,
                 files_raw_blast,
                 files_concatenated,
                 hit_frac_cutoff,
                 minimum_taxa,
                 threads,
                 log,
                 hc,
                 bc,
                 mcl_inflation,
                 perfect_identity,
                 coverage_threshold,
                 min_seq_length):
        self.printClass         = PrintOut(log, hc, bc)
        self.printout           = self.printClass.printout
        self.dir_base           = dir_base
        self.dir_treeforge      = dir_treeforge
        self.dir_mcl            = dir_mcl
        self.dir_iteration      = dir_iteration
        self.files_raw_blast    = Path(files_raw_blast)
        self.files_concatenated = Path(files_concatenated)
        self.hit_frac_cutoff    = hit_frac_cutoff
        self.minimum_taxa       = minimum_taxa
        self.threads            = threads
        self.log                = log
        self.hc                 = hc
        self.bc                 = bc
        self.mcl_inflation      = mcl_inflation
        self.perfect_identity   = perfect_identity
        self.coverage_threshold = coverage_threshold
        self.min_seq_length     = min_seq_length
    
        self.blast_mcl_out      = Path(os.path.join(self.dir_mcl, self.files_raw_blast.with_suffix(f'.hit-frac{self.hit_frac_cutoff}.minusLogEvalue').name))
        self.mcl_out            = Path(os.path.join(self.dir_mcl, self.files_raw_blast.with_suffix(f'.hit-frac{self.hit_frac_cutoff}_I{self.mcl_inflation}_e5').name))
        self.ident_out          = Path(os.path.join(self.dir_mcl, self.files_raw_blast.with_suffix(f'.ident').name))
        self.mcl_count          = 0

        self.return_dict        = {
                                    'blast_mcl_out' : str(self.blast_mcl_out),
                                    'mcl_out'       : str(self.mcl_out),
                                    'ident_out'     : str(self.ident_out),
                                    'mcl_count'     : self.mcl_count
                                }

    def run(self):
        self.blast_to_mcl()
        self.mcl()
        self.mcl_to_fasta()
        self.return_mcl()
        return self.return_dict
    
    def get_taxon_name(self, taxon: str) -> str:
        """Get taxon name from sequence ID."""
        return taxon.split("@")[0]

    def get_minus_log_evalue(self, raw_value: str) -> float:
        """Get minus log evalue from raw value."""
        try:
            value = float(raw_value)
            if value <= 0:
                return 180.0
            result = -np.log10(value)
            if result < 0.01:
                return 0.01
            return result
        except (ValueError, TypeError):
            return 180.0

    def blast_to_mcl(self):
        """Convert BLAST output to MCL input."""
        self.printout('metric', 'BLAST to MCL')
        outname = self.blast_mcl_out
        current_hits: Dict[Tuple[str, str], Tuple[float, float, float, float, float]] = {}
        max_minus_log_evalue: Dict[Tuple[str, str], float] = {}

        with open(self.files_raw_blast) as infile, \
             open(outname, 'w') as outfile, \
             open(self.ident_out, 'w') as ident_file:
            
            for line in infile:
                if len(line) < 3:
                    continue
                    
                spls = line.strip().split("\t")
                query, hit = spls[0], spls[2]
                
                if query == hit:
                    continue
                    
                query_taxon, hit_taxon = self.get_taxon_name(query), self.get_taxon_name(hit)
                
                pident = float(spls[5])
                qlen, qstart, qend = map(float, (spls[1], spls[10], spls[11]))
                slen, sstart, send = map(float, (spls[3], spls[12], spls[13]))
                
                if pident == self.perfect_identity and query_taxon != hit_taxon:
                    if (abs(qend - qstart) / qlen >= self.coverage_threshold or abs(send - sstart) / slen >= self.coverage_threshold) and \
                       qlen >= self.min_seq_length and slen >= self.min_seq_length:
                        ident_file.write(line)
                
                hit_key = (query, hit)
                minus_log_evalue = self.get_minus_log_evalue(spls[14])
                
                if hit_key in current_hits:
                    prev_qstart, prev_qend, prev_sstart, prev_send, prev_evalue = current_hits[hit_key]
                    current_hits[hit_key] = (
                        min(qstart, prev_qstart),
                        max(qend, prev_qend),
                        min(sstart, prev_sstart),
                        max(send, prev_send),
                        max(minus_log_evalue, prev_evalue)
                    )
                else:
                    current_hits[hit_key] = (qstart, qend, sstart, send, minus_log_evalue)
                    max_minus_log_evalue[hit_key] = minus_log_evalue
            
            for (query, hit), (qstart, qend, sstart, send, evalue) in current_hits.items():
                perc_qrange = (qend - qstart + 1) / qlen
                perc_srange = (send - sstart + 1) / slen
                
                if perc_qrange >= self.hit_frac_cutoff and perc_srange >= self.hit_frac_cutoff:
                    outfile.write(f"{query}\t{hit}\t{max_minus_log_evalue[(query, hit)]}\n")
        
        return str(outname)
    
    def mcl(self):
        """Run MCL."""
        self.printout('metric', 'MCL')
        cmd = f"mcl {self.blast_mcl_out} --abc -te {self.threads} -tf 'gq(5)' -I {self.mcl_inflation} -o {self.mcl_out}"
        mcl = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if mcl.returncode != 0:
            self.printout('error', 'MCL failed')
            sys.exit(1)
        else:
            return
        
    def mcl_to_fasta(self):
        """Convert MCL output to FASTA."""
        self.printout('metric', 'MCL to fasta')
        clusterDICT = {}
        with open(self.mcl_out, 'r') as infile:
            for line in infile:
                if len(line) < 3: 
                    continue
                spls = line.strip().split('\t')
                if len(set(i.split("@")[0] for i in spls)) >= self.minimum_taxa:
                    self.mcl_count += 1
                    clusterID = str(self.mcl_count)
                    for seqID in spls:
                        clusterDICT[seqID] = clusterID
        with pysam.FastxFile(self.files_concatenated) as fasta:
            for entry in fasta:
                seqid = entry.name
                if seqid in clusterDICT:
                    clusterID = clusterDICT[seqid]
                    filename = os.path.join(self.dir_iteration, f"cluster{clusterID}.fa")
                    with open(filename, "a") as outfile:
                        outfile.write(f">{seqid}\n{entry.sequence}\n")
                        
    def return_mcl(self):
        """Return MCL results."""
        self.return_dict = {
            'blast_mcl_out' : str(self.blast_mcl_out),
            'mcl_out'       : str(self.mcl_out),
            'ident_out'     : str(self.ident_out),
            'mcl_count'     : self.mcl_count
        }
        return