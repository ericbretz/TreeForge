import os
import csv
import tempfile
import numpy as np
from pathlib                import Path
from typing                 import Dict, Tuple
from core.stages.base_stage import BaseStage
from core.utils.sublogger   import run_stage_subprocess
from core.utils.fasta_io    import read_fasta_sequences, write_single_sequence

class MCL(BaseStage):
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
                 min_seq_length,
                 subprocess_dir=None,
                 shared_printClass=None):
        super().__init__(log, hc, bc, threads, subprocess_dir, shared_printClass)
        self.dir_base           = Path(dir_base)
        self.dir_treeforge      = Path(dir_treeforge)
        self.dir_mcl            = Path(dir_mcl)
        self.dir_iteration      = Path(dir_iteration)
        self.files_raw_blast    = Path(files_raw_blast)
        self.files_concatenated = Path(files_concatenated)
        self.hit_frac_cutoff    = hit_frac_cutoff
        self.minimum_taxa       = minimum_taxa
        self.mcl_inflation      = mcl_inflation
        self.perfect_identity   = perfect_identity
        self.coverage_threshold = coverage_threshold
        self.min_seq_length     = min_seq_length
    
        self.blast_mcl_out      = self.dir_mcl / self.files_raw_blast.with_suffix(f'.hit-frac{self.hit_frac_cutoff}.minusLogEvalue').name
        self.mcl_out            = self.dir_mcl / self.files_raw_blast.with_suffix(f'.hit-frac{self.hit_frac_cutoff}_I{self.mcl_inflation}_e5').name
        self.ident_out          = self.dir_mcl / self.files_raw_blast.with_suffix(f'.ident').name
        self.mcl_count          = 0

        self.return_dict        = {
                                    'blast_mcl_out' : str(self.blast_mcl_out),
                                    'mcl_out'       : str(self.mcl_out),
                                    'ident_out'     : str(self.ident_out),
                                    'mcl_count'     : self.mcl_count
                                }

    def run(self):
        """
        Run MCL.
        """
        self.blast_to_mcl()
        self.mcl()
        self.mcl_to_fasta()
        self.return_mcl()
        return self.return_dict
    
    def get_taxon_name(self, taxon: str) -> str:
        """
        Get taxon name from sequence ID.
        Split sequence ID by '@' and return the first part
        """
        return taxon.split("@")[0]

    def get_minus_log_evalue(self, raw_value: str) -> float:
        """
        Get minus log evalue from raw value.
        """
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
        
        FLUSH_THRESHOLD = 100000
        current_hits: Dict[Tuple[str, str], Tuple[float, float, float, float, float, float, float]] = {}
        temp_files      = []
        
        with open(self.files_raw_blast) as infile, \
             open(self.ident_out, 'w') as ident_file:
            
            reader = csv.reader(infile, delimiter='\t')
            
            for spls in reader:
                if len(spls) < 15:
                    continue
                
                query, hit = spls[0], spls[2]
                
                if query == hit:
                    continue
                
                query_taxon, hit_taxon = self.get_taxon_name(query), self.get_taxon_name(hit)
                
                pident = float(spls[5])
                qlen, qstart, qend = float(spls[1]), float(spls[10]), float(spls[11])
                slen, sstart, send = float(spls[3]), float(spls[12]), float(spls[13])
                
                if pident == self.perfect_identity and query_taxon != hit_taxon:
                    if (abs(qend - qstart) / qlen >= self.coverage_threshold or 
                        abs(send - sstart) / slen >= self.coverage_threshold) and \
                       qlen >= self.min_seq_length and slen >= self.min_seq_length:
                        ident_file.write('\t'.join(spls) + '\n')
                
                hit_key = (query, hit)
                minus_log_evalue = self.get_minus_log_evalue(spls[14])
                
                if hit_key in current_hits:
                    prev_qstart, prev_qend, prev_sstart, prev_send, prev_evalue, stored_qlen, stored_slen = current_hits[hit_key]
                    current_hits[hit_key] = (
                        min(qstart, prev_qstart),
                        max(qend, prev_qend),
                        min(sstart, prev_sstart),
                        max(send, prev_send),
                        max(minus_log_evalue, prev_evalue),
                        stored_qlen,
                        stored_slen
                    )
                else:
                    current_hits[hit_key] = (qstart, qend, sstart, send, minus_log_evalue, qlen, slen)
                
                if len(current_hits) >= FLUSH_THRESHOLD:
                    temp_file = self._flush_hits_to_temp(current_hits)
                    temp_files.append(temp_file)
                    current_hits.clear()
        
        if current_hits:
            temp_file = self._flush_hits_to_temp(current_hits)
            temp_files.append(temp_file)
        
        self._merge_temp_files(temp_files, outname)
        
        for temp_file in temp_files:
            try:
                os.unlink(temp_file)
            except:
                pass
        
        return str(outname)
    
    def _flush_hits_to_temp(self, hits_dict: Dict) -> str:
        """Flush current hits to a temporary file."""
        temp_fd, temp_path = tempfile.mkstemp(suffix='.mcl_temp', dir=self.dir_mcl)
        
        with os.fdopen(temp_fd, 'w') as f:
            for (query, hit), (qstart, qend, sstart, send, evalue, qlen, slen) in hits_dict.items():
                perc_qrange = (qend - qstart + 1) / qlen
                perc_srange = (send - sstart + 1) / slen
                
                if perc_qrange >= self.hit_frac_cutoff and perc_srange >= self.hit_frac_cutoff:
                    f.write(f"{query}\t{hit}\t{evalue}\n")
        
        return temp_path
    
    def _merge_temp_files(self, temp_files: list, output_file: str):
        """Merge temporary files into final output, handling duplicates."""
        if not temp_files:
            open(output_file, 'w').close()
            return
        
        if len(temp_files) == 1:
            os.rename(temp_files[0], output_file)
            return
        
        final_hits: Dict[Tuple[str, str], float] = {}
        
        for temp_file in temp_files:
            with open(temp_file, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    if len(row) >= 3:
                        query, hit, evalue = row[0], row[1], float(row[2])
                        key = (query, hit)
                        if key in final_hits:
                            final_hits[key] = max(final_hits[key], evalue)
                        else:
                            final_hits[key] = evalue
        
        with open(output_file, 'w') as f:
            for (query, hit), evalue in final_hits.items():
                f.write(f"{query}\t{hit}\t{evalue}\n")
        
        return str(output_file)
    
    def mcl(self):
        """Run MCL."""
        self.printout('metric', 'MCL')
        cmd = f"mcl {self.blast_mcl_out} --abc -te {self.threads} -tf 'gq(5)' -I {self.mcl_inflation} -o {self.mcl_out}"
        run_stage_subprocess(cmd, 'MCL', self.subprocess_dir, self.printout)
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
        
        sequences = read_fasta_sequences(self.files_concatenated)
        for seqid, sequence in sequences.items():
            if seqid in clusterDICT:
                clusterID = clusterDICT[seqid]
                filename = self.dir_iteration / f"cluster{clusterID}.fa"
                write_single_sequence(seqid, sequence, filename, mode='a')
                        
    def return_mcl(self):
        """Return MCL results."""
        self.return_dict = {
            'blast_mcl_out'  : str(self.blast_mcl_out),
            'mcl_out'        : str(self.mcl_out),
            'ident_out'      : str(self.ident_out),
            'mcl_count'      : self.mcl_count,
            'mcl_inflation'  : self.mcl_inflation,
            'mcl_hit_frac'   : self.hit_frac_cutoff,
            'mcl_eval_power' : 5,
        }
        return