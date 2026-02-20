import os
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple
from collections import defaultdict
from core.stages.base_stage import BaseStage
from core.utils.fasta_io import read_fasta_sequences, write_fasta_sequences

class GeneCluster(BaseStage):
    def __init__(self,
                 dir_base,
                 dir_treeforge,
                 dir_genecluster,
                 files_raw_blast,
                 files_concatenated,
                 minimum_taxa,
                 threads,
                 log,
                 hc,
                 bc,
                 shared_printClass=None):
        super().__init__(log, hc, bc, threads, shared_printClass=shared_printClass)
        self.dir_base           = Path(dir_base)
        self.dir_treeforge      = Path(dir_treeforge)
        self.dir_genecluster    = Path(dir_genecluster)
        self.files_raw_blast    = Path(files_raw_blast)
        self.files_concatenated = Path(files_concatenated)
        self.minimum_taxa       = minimum_taxa
        
        self.gene_count         = 0
        self.gene_fasta_files   = []
        self.return_dict        = {
            'gene_count': self.gene_count,
            'gene_files': self.gene_fasta_files
        }
    
    def run(self):
        """
        Run gene-based clustering for HCluster mode.
        """

        # self.printout('metric', 'Gene-based clustering for HCluster mode')
        
        self.dir_genecluster.mkdir(parents=True, exist_ok=True)
        
        blast_data    = self.parse_blast_results()
        gene_groups   = self.group_by_gene(blast_data)
        gene_clusters = self.select_best_per_species(gene_groups)
        
        self.write_gene_fastas(gene_clusters)
        
        self.return_dict['gene_count'] = self.gene_count
        self.return_dict['gene_files'] = self.gene_fasta_files
        
        self.printout('metric', [('gene_clusters_created', self.gene_count)])
        
        return self.return_dict
    
    def parse_blast_results(self) -> List[Dict]:
        """
        Parse BLAST results.
        """

        self.printout('metric', 'Parsing BLAST results')
        
        blast_hits = []
        
        if not self.files_raw_blast.exists():
            self.printout('error', f'BLAST results file not found: {self.files_raw_blast}')
            return blast_hits
        
        with open(self.files_raw_blast, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if len(line.strip()) < 3:
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 16:
                    self.printout('warning', f'Malformed BLAST line {line_num}: insufficient columns')
                    continue
                
                try:
                    query   = parts[0]
                    qlen    = int(parts[1])
                    subject = parts[2]
                    slen    = int(parts[3])
                    pident  = float(parts[5])
                    length  = int(parts[7])
                    
                    if query == subject:
                        continue
                    
                    query_species, query_gene = self._parse_sequence_id(query)
                    subject_species, subject_gene = self._parse_sequence_id(subject)
                    
                    if not query_gene or not subject_gene:
                        continue
                    
                    blast_hits.append({
                        'query'          : query,
                        'query_species'  : query_species,
                        'query_gene'     : query_gene,
                        'qlen'           : qlen,
                        'subject'        : subject,
                        'subject_species': subject_species,
                        'subject_gene'   : subject_gene,
                        'slen'           : slen,
                        'pident'         : pident,
                        'length'         : length
                    })
                
                except (ValueError, IndexError) as e:
                    self.printout('warning', f'Error parsing BLAST line {line_num}: {str(e)}')
                    continue
        
        self.printout('metric', [('blast_hits_parsed', len(blast_hits))])
        return blast_hits
    
    def _parse_sequence_id(self, seq_id: str) -> Tuple[str, str]:
        """
        Parse sequence ID.
        Split sequence ID by '@' and return the first and second parts.
        """
        try:
            parts = seq_id.split('@')
            if len(parts) == 2:
                return parts[0], parts[1]
        except Exception:
            pass
        return '', ''
    
    def group_by_gene(self, blast_hits: List[Dict]) -> Dict[str, List[Dict]]:
        """
        Group hits by gene ID.
        Return a dictionary with gene IDs as keys and lists of hits as values.
        """
        self.printout('metric', 'Grouping hits by gene')
        
        gene_groups = defaultdict(list)
        
        for hit in blast_hits:
            gene_id = hit['query_gene']
            gene_groups[gene_id].append(hit)
        
        self.printout('metric', [('unique_genes_found', len(gene_groups))])
        return dict(gene_groups)
    
    def select_best_per_species(self, gene_groups: Dict[str, List[Dict]]) -> Dict[str, Dict[str, str]]:
        """
        Select best match per species for each gene.
        Return a dictionary with gene IDs as keys and dictionaries of sequence IDs as values.
        """

        # self.printout('metric', 'Selecting best match per species for each gene')
        
        gene_clusters = {}
        
        for gene_id, hits in gene_groups.items():
            species_hits = defaultdict(list)
            
            if hits:
                query_species = hits[0]['query_species']
                query_seq = hits[0]['query']
                species_hits[query_species].append({
                    'seq_id': query_seq,
                    'pident': 100.0,
                    'length': hits[0]['qlen']
                })
            
            for hit in hits:
                subject_species = hit['subject_species']
                species_hits[subject_species].append({
                    'seq_id': hit['subject'],
                    'pident': hit['pident'],
                    'length': hit['slen']
                })
            
            best_matches = {}
            for species, matches in species_hits.items():
                sorted_matches = sorted(
                    matches,
                    key=lambda x: (x['pident'], x['length']),
                    reverse=True
                )
                best_match = sorted_matches[0]
                best_matches[best_match['seq_id']] = best_match['seq_id']
            
            if len(best_matches) >= self.minimum_taxa:
                gene_clusters[gene_id] = best_matches
        
        self.printout('metric', [('genes_selected', len(gene_clusters)), ('minimum_taxa', self.minimum_taxa)])
        return gene_clusters
    
    def write_gene_fastas(self, gene_clusters: Dict[str, Dict[str, str]]):
        """
        Write gene FASTA files for each gene cluster.
        """
        self.printout('metric', 'Writing gene FASTA files')
        
        try:
            sequences = read_fasta_sequences(self.files_concatenated)
        except Exception as e:
            self.printout('error', f'Error reading concatenated FASTA: {str(e)}')
            return
        
        for gene_id, seq_ids in gene_clusters.items():
            output_file = self.dir_genecluster / f'{gene_id}.fa'
            
            try:
                gene_sequences = {}
                for seq_id in seq_ids.keys():
                    if seq_id in sequences:
                        gene_sequences[seq_id] = sequences[seq_id]
                    else:
                        self.printout('warning', f'Sequence {seq_id} not found in concatenated FASTA')
                
                write_fasta_sequences(gene_sequences, output_file)
                self.gene_count += 1
                self.gene_fasta_files.append(str(output_file))
            
            except Exception as e:
                self.printout('error', f'Error writing gene file {output_file}: {str(e)}')
                continue
        
        self.printout('metric', [('gene_fasta_files_written', self.gene_count)])
