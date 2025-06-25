import os
from Bio import SeqIO
import string
import random

"""
This renames the transcripts in fasta files so that treeforge doesnt lock up or crash. including whatever it wraps.
Its a basic map of original names to new names.

Later the map will be used to rename species trees after the source file names.

file: file1.fasta
    contigA@scaffold01238 -> AAAAAA@00000001 -> file1
    contigB_scaffold01239 -> AAAAAA@00000002 -> file1
    con3a91231jafold01240 -> AAAAAA@00000003 -> file1

file: file2.fasta
    contigD@scaffold01241 -> AAAAAA@00000004 -> file2
    contigE@scaffold01242 -> AAAAAA@00000005 -> file2
    contigF@scaffold01243 -> AAAAAA@00000006 -> file2

file: file3.fasta
    contigG@scaffold01244 -> AAAAAA@00000007 -> file3
    contigH@scaffold01245 -> AAAAAA@00000008 -> file3
    contigI@scaffold01246 -> AAAAAA@00000009 -> file3

Final Tree:
    ┌── file1
   ┌┤
───┤└── file3
   │
   └─── file2
"""

class Contig:
    def __init__(self, dir_base=".", dir_temp_fasta="."):
        self.dir_base        = dir_base
        self.dir_temp_fasta  = dir_temp_fasta
        self.name_mapping    = {}
        self.reverse_mapping = {}
        self.current_prefix  = None
        self.counter         = 0
        self.return_dict     = {}
        
        os.makedirs(self.dir_temp_fasta, exist_ok=True)

    def _generate_prefix(self):
        """Creates a random 5-character prefix using uppercase letters."""
        return ''.join(random.choices(string.ascii_uppercase, k=5))

    def _generate_new_name(self, prefix):
        """Creates a new simplified name in the format PREFIX@COUNTER."""
        new_name = f"{prefix}@{self.counter:08d}"
        self.counter += 1
        return new_name

    def rename_contigs(self, fasta_files):
        """
        Renames contigs in fasta files with simplified names.
        This should take care of any formatting or weird characters people use for transcript names
        """
        renamed_files = []
        total_contigs = 0
        
        for fasta_file in fasta_files:
            if not os.path.exists(fasta_file):
                pass
                continue
                
            prefix = self._generate_prefix()
            self.current_prefix = prefix
            
            base_name   = os.path.basename(fasta_file)
            name, ext   = os.path.splitext(base_name)
            output_file = os.path.join(self.dir_temp_fasta, f"{name}_renamed{ext}")
            
            file_contigs = self._process_fasta_file(fasta_file, output_file, prefix)
            renamed_files.append(output_file)
            total_contigs += file_contigs
            
            self.counter = 0
        
        self.return_dict = {
            'renamed_files'  : renamed_files,
            'total_contigs'  : total_contigs,
            'files_processed': len(renamed_files),
            'name_mappings'  : len(self.name_mapping)
        }
        return renamed_files

    def _process_fasta_file(self, input_file, output_file, prefix):
        """Processes a single fasta file and renames the contigs."""
        source_file  = os.path.basename(input_file)
        contig_count = 0
        
        with open(output_file, 'w') as out_handle:
            for record in SeqIO.parse(input_file, "fasta"):
                original_name = record.id
                new_name      = self._generate_new_name(prefix)
                
                self.name_mapping[new_name]         = (original_name, source_file)
                self.reverse_mapping[original_name] = (new_name, source_file)
                
                record.id          = new_name
                record.description = ""
                SeqIO.write(record, out_handle, "fasta")
                contig_count += 1
        
        return contig_count

    def get_original_name(self, new_name):
        """
        Gets the original contig name from the simplified name.
        """
        mapping = self.name_mapping.get(new_name)
        return mapping[0] if mapping else None

    def get_source_file(self, new_name):
        """
        Gets the source file name for a simplified contig name.
        """
        mapping = self.name_mapping.get(new_name)
        return mapping[1] if mapping else None

    def get_original_info(self, new_name):
        """
        Gets both original name and source file for a simplified contig name.
        """
        return self.name_mapping.get(new_name, (None, None))

    def get_new_name(self, original_name):
        """
        Gets the simplified contig name from the original name.
        """
        mapping = self.reverse_mapping.get(original_name)
        return mapping[0] if mapping else None

    def get_new_info(self, original_name):
        """
        Gets both simplified name and source file for an original contig name.
        """
        return self.reverse_mapping.get(original_name, (None, None))

    def save_mapping(self, output_file):
        """
        Saves the name mapping to a file for later use.
        """
        with open(output_file, 'w') as f:
            f.write("New_Name\tOriginal_Name\tSource_File\n")
            for new_name, (original_name, source_file) in self.name_mapping.items():
                f.write(f"{new_name}\t{original_name}\t{source_file}\n")

    def load_mapping(self, mapping_file):
        """
        Loads a name mapping from a file.
        """
        self.name_mapping.clear()
        self.reverse_mapping.clear()
        
        with open(mapping_file, 'r') as f:
            next(f)
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    new_name, original_name, source_file = parts[:3]
                    self.name_mapping[new_name] = (original_name, source_file)
                    self.reverse_mapping[original_name] = (new_name, source_file)

    def revert_fasta_file(self, renamed_file, output_file):
        """
        Reverts a renamed fasta file back to original names.
        """
        with open(output_file, 'w') as out_handle:
            for record in SeqIO.parse(renamed_file, "fasta"):
                new_name = record.id
                original_name = self.get_original_name(new_name)
                
                if original_name:
                    record.id = original_name
                else:
                    pass
                SeqIO.write(record, out_handle, "fasta")

    def ID_iteration(self):
        """Legacy method kept for compatibility."""
        pass