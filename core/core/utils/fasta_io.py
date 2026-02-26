import pysam
from pathlib import Path
from typing import Dict, List, Union, Set


def read_fasta_sequences(fasta_file: Union[str, Path]) -> Dict[str, str]:
    sequences = {}
    with pysam.FastxFile(str(fasta_file)) as fasta:
        for entry in fasta:
            sequences[entry.name] = entry.sequence
    return sequences

def write_fasta_sequences(sequences: Dict[str, str], output_file: Union[str, Path], mode: str = 'w') -> None:
    with open(str(output_file), mode) as f:
        for seq_id, sequence in sequences.items():
            f.write(f">{seq_id}\n{sequence}\n")

def write_single_sequence(seq_id: str, sequence: str, output_file: Union[str, Path], mode: str = 'w') -> None:
    with open(str(output_file), mode) as f:
        f.write(f">{seq_id}\n{sequence}\n")

def filter_fasta_by_ids(
    input_file : Union[str, Path],
    output_file: Union[str, Path],
    keep_ids   : Set[str],
    mode       : str = 'w'
) -> int:
    count = 0
    with pysam.FastxFile(str(input_file)) as fasta_in:
        with open(str(output_file), mode) as fasta_out:
            for entry in fasta_in:
                if entry.name in keep_ids:
                    fasta_out.write(f">{entry.name}\n{entry.sequence}\n")
                    count += 1
    return count

def concatenate_fasta_files(
    input_files: List[Union[str, Path]],
    output_file: Union[str, Path]
) -> int:
    total_count = 0
    lines       = []
    for input_file in input_files:
        fasta      = pysam.FastaFile(str(input_file))
        references = fasta.references
        for name in references:
            seq = fasta.fetch(name)
            lines.append(f">{name}\n{seq}")
            total_count += 1
    with open(str(output_file), 'w') as f:
        f.write("\n".join(lines))
    return total_count

def get_fasta_sequence_count(fasta_file: Union[str, Path]) -> int:
    count = 0
    with pysam.FastxFile(str(fasta_file)) as fasta:
        for _ in fasta:
            count += 1
    return count

def get_fasta_ids(fasta_file: Union[str, Path]) -> List[str]:
    ids = []
    with pysam.FastxFile(str(fasta_file)) as fasta:
        for entry in fasta:
            ids.append(entry.name)
    return ids

def split_fasta_by_mapping(
    input_file        : Union[str, Path],
    output_dir        : Union[str, Path],
    id_to_file_mapping: Dict[str, str]
) -> Dict[str, int]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    file_handles = {}
    file_counts = {}
    
    try:
        with pysam.FastxFile(str(input_file)) as fasta:
            for entry in fasta:
                if entry.name in id_to_file_mapping:
                    filename = id_to_file_mapping[entry.name]
                    
                    if filename not in file_handles:
                        filepath = output_dir / filename
                        file_handles[filename] = open(filepath, 'w')
                        file_counts[filename] = 0
                    
                    file_handles[filename].write(f">{entry.name}\n{entry.sequence}\n")
                    file_counts[filename] += 1
    finally:
        for handle in file_handles.values():
            handle.close()
    
    return file_counts
