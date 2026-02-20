import sys
import os
import pysam
import subprocess
import glob
import shutil
from pathlib import Path
from core.stages.base_stage import BaseStage
from core.utils.sublogger import run_stage_subprocess
from core.stages.hcluster import HCluster

class Busco(BaseStage):
    def __init__(self,
                 dir_base,
                 dir_treeforge,
                 dir_busco,
                 dir_iteration,
                 files_fasta,
                 threads,
                 log,
                 hc,
                 bc,
                 busco_evalue,
                 busco_max_targets,
                 busco_coverage_threshold,
                 hcluster_id,
                 hcluster_iddef,
                 species_tree_file=None,
                 id_to_stem=None,
                 subprocess_dir=None,
                 shared_printClass=None):
        super().__init__(log, hc, bc, threads, subprocess_dir, shared_printClass)
        self.dir_base           = Path(dir_base)
        self.dir_treeforge      = Path(dir_treeforge)
        self.dir_busco_root     = Path(dir_busco)
        self.dir_busco          = Path(dir_busco) / 'busco'
        self.dir_iteration      = Path(dir_iteration)
        self.files_fasta        = files_fasta
        self.busco_evalue       = busco_evalue
        self.busco_max_targets  = busco_max_targets
        self.busco_coverage_threshold = busco_coverage_threshold
        self.hcluster_id        = hcluster_id
        self.hcluster_iddef     = hcluster_iddef
        self.species_tree_file  = species_tree_file
        self.id_to_stem         = id_to_stem if id_to_stem is not None else {}
        
        self.busco_count        = 0
        self.return_dict        = {
                                    'busco_count'      : self.busco_count,
                                    'busco_output_dir' : str(self.dir_busco)
                                }

    def run(self):
        """
        Run BUSCO extraction.
        """
        self.extract_busco_genes()
        return self.return_dict
    
    def run_with_clustering(self):
        """
        Extract BUSCO genes and run hierarchical clustering.
        """
        self.extract_busco_genes()
        self.run_hierarchical_clustering()
        return self.return_dict

    def extract_busco_genes(self):
        """
        Extract BUSCO genes using DIAMOND.
        """
        self.printout('metric', 'Extracting BUSCO genes using DIAMOND')
        
        if not self.files_fasta:
            self.printout('error', 'No input FASTA files provided')
            return
        
        missing_files = [f for f in self.files_fasta if not os.path.exists(f)]
        if missing_files:
            self.printout('error', f'Input FASTA files not found: {missing_files}')
            return
        
        current_file_dir     = Path(__file__).resolve().parent
        treeforge_script_dir = current_file_dir.parent.parent
        self.busco_db_path   = treeforge_script_dir / 'core' / 'utils' / 'euk_db.faa'
        
        if not self.busco_db_path.exists():
            self.printout('error', f'BUSCO database not found at {self.busco_db_path}')
            self.printout('error', 'Please ensure the BUSCO database file exists in core/utils/')
            return
        
        try:
            with open(self.busco_db_path, 'r') as f:
                first_line = f.readline().strip()
                if not first_line.startswith('>'):
                    self.printout('error', f'Invalid BUSCO database format: {self.busco_db_path}')
                    return
        except Exception as e:
            self.printout('error', f'Error reading BUSCO database: {str(e)}')
            return
        
        self.busco_db_dir = self.dir_busco / 'busco_db'
        
        self.dir_busco.mkdir(parents=True, exist_ok=True)
        
        self.make_busco_blast_db()
        
        self.busco_species_files = []
        total_files              = len(self.files_fasta)
        for idx, species_fasta in enumerate(self.files_fasta, 1):
            try:
                species_busco_file = self.process_species_fasta(species_fasta, idx, total_files)
                if species_busco_file:
                    self.busco_species_files.append(species_busco_file)
            except Exception as e:
                self.printout('error', f'Error processing {species_fasta}: {str(e)}')
                continue
        
        self.printout('progress', f'Processing complete: {len(self.busco_species_files)}/{total_files} files processed')
        print()
        
        self.busco_count                  = len(self.busco_species_files)
        self.return_dict['busco_count']   = self.busco_count
        self.return_dict['busco_files']   = self.busco_species_files
    
    def process_species_fasta(self, species_fasta, current=None, total=None):
        """
        Process a species FASTA file to extract BUSCO genes.
        """
        species_name = os.path.splitext(os.path.basename(species_fasta))[0]
        if current is not None and total is not None:
            self.printout('progress', f'Processing {species_name} ({current}/{total})')
        else:
            self.printout('progress', f'Processing {species_name}')
        
        species_blast_results = self.dir_busco / f'{species_name}_busco_blast.txt'
        
        try:
            self.busco_blast_search_species(species_fasta, species_blast_results)
        except Exception as e:
            self.printout('error', f'BLAST search failed for {species_name}: {str(e)}')
            return None
        
        species_busco_fasta = self.dir_busco / f'{species_name}_busco.fa'
        try:
            self.create_species_busco_fasta(species_fasta, species_blast_results, species_busco_fasta)
        except Exception as e:
            self.printout('error', f'Failed to create BUSCO FASTA for {species_name}: {str(e)}')
            return None
        
        if os.path.exists(species_busco_fasta) and os.path.getsize(species_busco_fasta) > 0:
            return species_busco_fasta
        else:
            self.printout('warning', f'No BUSCO hits found for {species_name}')
            if os.path.exists(species_busco_fasta):
                os.remove(species_busco_fasta)
            return None
    
    def busco_blast_search_species(self, species_fasta, output_file):
        """
        Run BUSCO DIAMOND blastx search.
        """
        db_path = self.busco_db_dir / "busco_db"
        cmd = f'diamond blastx --query {species_fasta} --db {db_path}.dmnd --out {output_file} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --evalue {self.busco_evalue} --max-target-seqs {self.busco_max_targets} --threads {self.threads}'
        
        try:
            result = run_stage_subprocess(cmd, f'BUSCO DIAMOND {Path(species_fasta).stem}', self.subprocess_dir, self.printout, check=False)
            if result.returncode != 0:
                self.printout('error', f'DIAMOND blastx failed for {species_fasta}')
                if result.stderr:
                    stderr_text = result.stderr if isinstance(result.stderr, str) else result.stderr.decode('utf-8')
                    self.printout('error', stderr_text)
                raise RuntimeError('DIAMOND search failed')
            else:
                if os.path.exists(output_file) and os.path.getsize(output_file) == 0:
                    self.printout('warning', f'No DIAMOND results for {os.path.basename(species_fasta)}')
        except Exception as e:
            self.printout('error', f'Error running DIAMOND search for {species_fasta}: {str(e)}')
            raise
    
    def create_species_busco_fasta(self, species_fasta, blast_results, output_fasta):
        """
        Create a species BUSCO FASTA file from DIAMOND blastx results.
        """
        busco_hit_ids = set()
        
        try:
            with open(blast_results, 'r') as infile:
                for line_num, line in enumerate(infile, 1):
                    if len(line.strip()) < 3:
                        continue
                        
                    parts = line.strip().split('\t')
                    if len(parts) >= 12:
                        query_id = parts[0]
                        evalue = float(parts[10])
                        
                        if evalue <= self.busco_evalue:
                            busco_hit_ids.add(query_id)
                    else:
                        self.printout('warning', f'Malformed BLAST line {line_num}: {line.strip()}')
        except Exception as e:
            self.printout('error', f'Error parsing BLAST results: {str(e)}')
            raise
        
        if not busco_hit_ids:
            self.printout('warning', f'No BUSCO hits found in {blast_results}')
            return
        
        try:
            with open(output_fasta, 'w') as outfile:
                with pysam.FastxFile(species_fasta) as fasta:
                    for entry in fasta:
                        if entry.name in busco_hit_ids:
                            outfile.write(f">{entry.name}\n{entry.sequence}\n")
        except Exception as e:
            self.printout('error', f'Error creating BUSCO FASTA: {str(e)}')
            raise
    
    
    def make_busco_blast_db(self):
        """
        Build BUSCO DIAMOND database.
        """
        self.printout('metric', 'Build BUSCO Database')
        self.busco_db_dir.mkdir(parents=True, exist_ok=True)
        
        db_path = self.busco_db_dir / "busco_db"
        cmd = f'diamond makedb --in {self.busco_db_path} --db {db_path}'
        
        try:
            result = run_stage_subprocess(cmd, 'BUSCO DIAMOND makedb', self.subprocess_dir, self.printout, check=False)
            if result.returncode != 0:
                self.printout('error', 'DIAMOND makedb failed')
                if result.stderr:
                    stderr_text = result.stderr if isinstance(result.stderr, str) else result.stderr.decode('utf-8')
                    self.printout('error', stderr_text)
                raise RuntimeError('DIAMOND database creation failed')
            else:
                return
        except Exception as e:
            self.printout('error', f'Error creating DIAMOND database: {str(e)}')
            raise
    
    
    
    def run_hierarchical_clustering(self):
        """
        Run hierarchical clustering.
        """
        self.printout('metric', 'Running Hierarchical Clustering')
        
        if not self.species_tree_file:
            self.printout('error', 'No species tree provided for hierarchical clustering')
            return
        
        if not hasattr(self, 'busco_species_files') or not self.busco_species_files:
            self.printout('error', 'No BUSCO species files to process')
            return
        
        species_file_map = {}
        for f in self.busco_species_files:
            fname = os.path.basename(f)
            species_name = fname.split('.')[0]

            if species_name not in species_file_map:
                species_file_map[species_name] = []
            species_file_map[species_name].append(str(f))

        hcluster = HCluster(
            self.dir_base,
            self.dir_treeforge,
            self.dir_busco,
            self.busco_species_files,
            self.species_tree_file,
            self.threads,
            self.log,
            self.hc,
            self.bc,
            self.hcluster_id,
            self.hcluster_iddef,
            species_file_map=species_file_map,
            id_to_stem=self.id_to_stem
        )
        
        hcluster.run()
        
        self.return_dict.update(hcluster.return_dict)
        self.return_dict['busco_files_processed'] = len(self.busco_species_files)
    
    
    

