import os
import string
import random
import pysam
from pathlib import Path
from typing import Tuple, Optional, Dict

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


class BiDirectionalMapping:
    def __init__(self):
        self._forward: Dict[str, Tuple[str, str]] = {}  # new_name -> (original_name, source_file)
        self._reverse: Dict[str, Tuple[str, str]] = {}  # original_name -> (new_name, source_file)
    
    def add(self, new_name: str, original_name: str, source_file: str) -> None:
        """Add a bidirectional mapping entry."""
        self._forward[new_name]      = (original_name, source_file)
        self._reverse[original_name] = (new_name, source_file)
    
    def get_by_new_name(self, new_name: str) -> Optional[Tuple[str, str]]:
        """Get (original_name, source_file) by new_name."""
        return self._forward.get(new_name)
    
    def get_by_original_name(self, original_name: str) -> Optional[Tuple[str, str]]:
        """Get (new_name, source_file) by original_name."""
        return self._reverse.get(original_name)
    
    def clear(self) -> None:
        """Clear all mappings."""
        self._forward.clear()
        self._reverse.clear()
    
    def __len__(self) -> int:
        """Return the number of mappings."""
        return len(self._forward)
    
    def items_forward(self):
        """Iterate over forward mappings."""
        return self._forward.items()
    
    def items_reverse(self):
        """Iterate over reverse mappings."""
        return self._reverse.items()


class Contig:
    def __init__(self, dir_base=".", dir_temp_fasta="."):
        self.dir_base        = Path(dir_base)
        self.dir_temp_fasta  = Path(dir_temp_fasta)
        self.mapping         = BiDirectionalMapping()
        self.current_prefix  = None
        self.counter         = 0
        self.return_dict     = {}
        
        self.dir_temp_fasta.mkdir(parents=True, exist_ok=True)
    
    @property
    def name_mapping(self) -> Dict[str, Tuple[str, str]]:
        return self.mapping._forward
    
    @property
    def reverse_mapping(self) -> Dict[str, Tuple[str, str]]:
        return self.mapping._reverse

    def _generate_prefix(self):
        return ''.join(random.choices(string.ascii_uppercase, k=5))

    def _generate_new_name(self, prefix):
        new_name = f"{prefix}@{self.counter:08d}"
        self.counter += 1
        return new_name

    def rename_contigs(self, fasta_files):
        renamed_files = []
        total_contigs = 0
        
        for fasta_file in fasta_files:
            if not os.path.exists(fasta_file):
                pass
                continue
                
            prefix = self._generate_prefix()
            self.current_prefix = prefix
            
            fasta_path   = Path(fasta_file)
            name         = fasta_path.stem
            ext          = fasta_path.suffix
            output_file  = self.dir_temp_fasta / f"{name}_renamed{ext}"

            file_contigs = self._process_fasta_file(fasta_file, output_file, prefix)
            renamed_files.append(output_file)
            total_contigs += file_contigs
            
            self.counter = 0
        
        self.return_dict = {
            'renamed_files'  : renamed_files,
            'total_contigs'  : total_contigs,
            'files_processed': len(renamed_files),
            'name_mappings'  : len(self.mapping)
        }
        return renamed_files

    def _process_fasta_file(self, input_file, output_file, prefix):
        source_file  = os.path.basename(input_file)
        contig_count = 0
        
        with open(output_file, 'w') as out_handle:
            with pysam.FastxFile(str(input_file)) as fasta:
                for entry in fasta:
                    original_name = entry.name
                    new_name      = self._generate_new_name(prefix)
                    
                    self.mapping.add(new_name, original_name, source_file)
                    
                    out_handle.write(f">{new_name}\n{entry.sequence}\n")
                    contig_count += 1
        
        return contig_count

    def get_original_name(self, new_name: str) -> Optional[str]:
        """Get original name from new name."""
        result = self.mapping.get_by_new_name(new_name)
        return result[0] if result else None

    def get_source_file(self, new_name: str) -> Optional[str]:
        """Get source file from new name."""
        result = self.mapping.get_by_new_name(new_name)
        return result[1] if result else None

    def get_original_info(self, new_name: str) -> Tuple[Optional[str], Optional[str]]:
        """Get (original_name, source_file) from new name."""
        return self.mapping.get_by_new_name(new_name) or (None, None)

    def get_new_name(self, original_name: str) -> Optional[str]:
        """Get new name from original name."""
        result = self.mapping.get_by_original_name(original_name)
        return result[0] if result else None

    def get_new_info(self, original_name: str) -> Tuple[Optional[str], Optional[str]]:
        """Get (new_name, source_file) from original name."""
        return self.mapping.get_by_original_name(original_name) or (None, None)

    def save_mapping(self, output_file: str) -> None:
        """Save name mapping to file."""
        with open(output_file, 'w') as f:
            f.write("New_Name\tOriginal_Name\tSource_File\n")
            for new_name, (original_name, source_file) in self.mapping.items_forward():
                f.write(f"{new_name}\t{original_name}\t{source_file}\n")

    def load_mapping(self, mapping_file: str) -> None:
        """Load name mapping from file."""
        self.mapping.clear()
        
        with open(mapping_file, 'r') as f:
            next(f)  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    new_name, original_name, source_file = parts[:3]
                    self.mapping.add(new_name, original_name, source_file)

    def revert_fasta_file(self, renamed_file, output_file):
        with open(output_file, 'w') as out_handle:
            with pysam.FastxFile(str(renamed_file)) as fasta:
                for entry in fasta:
                    new_name      = entry.name
                    original_name = self.get_original_name(new_name)
                    
                    if original_name:
                        out_handle.write(f">{original_name}\n{entry.sequence}\n")
                    else:
                        out_handle.write(f">{new_name}\n{entry.sequence}\n")