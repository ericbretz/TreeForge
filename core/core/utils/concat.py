import sys
from Bio import SeqIO


class SequenceConcatenator:
    _DNA_CHARS = frozenset('ACGTRYMKSWHBVDN')

    def __init__(self, uppercase=False, printout=None):
        self.uppercase      = uppercase
        self.printout       = printout
        self.taxa           = {}
        self.partition_info = []

    @staticmethod
    def _parse_species_name(full_id):
        if '@' in full_id:
            parts = full_id.split('@')
            if len(parts) == 2 and parts[0].isupper() and parts[1].replace('_', '').isdigit():
                return parts[0]
        elif '_' in full_id:
            parts = full_id.split('_', 1)
            if len(parts) == 2 and parts[0].isupper() and parts[1][0].isdigit():
                return parts[0]
        return full_id

    def add_alignment(self, filepath):
        sequences  = {}
        seq_length = None
        seq_type   = 'DNA'

        try:
            records = list(SeqIO.parse(filepath, 'fasta'))
        except Exception as e:
            self.printout('error', f"Error reading {filepath}: {e}")
            sys.exit(1)

        if not records:
            self.printout('error', f"No sequences found in {filepath}")
            sys.exit(1)

        lengths = [len(rec.seq) for rec in records]
        if len(set(lengths)) > 1:
            self.printout('error', f"Sequences in {filepath} are not aligned")
            sys.exit(1)

        seq_length = lengths[0]

        for rec in records:
            species_name = self._parse_species_name(rec.id)

            if species_name in sequences:
                continue

            seq = str(rec.seq)
            sequences[species_name] = seq

            if seq_type == 'DNA':
                seq_upper = seq.upper().replace('-', '').replace('N', '')
                if seq_upper and not set(seq_upper).issubset(self._DNA_CHARS):
                    seq_type = 'AA'

        file_index = len(self.partition_info)
        new_taxa   = [t for t in sequences if t not in self.taxa]

        for taxon in list(self.taxa.keys()) + new_taxa:
            if taxon not in self.taxa:
                self.taxa[taxon] = []

            while len(self.taxa[taxon]) < file_index:
                _, prev_start, prev_end, _ = self.partition_info[len(self.taxa[taxon])]
                prev_length = prev_end - prev_start + 1
                self.taxa[taxon].append(('-' * prev_length, len(self.taxa[taxon])))

            if taxon in sequences:
                seq = sequences[taxon]
                if self.uppercase:
                    seq = seq.upper()
                self.taxa[taxon].append((seq, file_index))
            else:
                self.taxa[taxon].append(('-' * seq_length, file_index))

        if self.partition_info:
            start = self.partition_info[-1][2] + 1
        else:
            start = 1
        end = start + seq_length - 1

        self.partition_info.append((filepath, start, end, seq_type))

    def write_concatenated(self, output_path):
        with open(output_path, 'w') as f:
            for taxon_id, seq_parts in self.taxa.items():
                while len(seq_parts) < len(self.partition_info):
                    missing_idx = len(seq_parts)
                    _, start, end, _ = self.partition_info[missing_idx]
                    seq_parts.append(('-' * (end - start + 1), missing_idx))

                f.write(f">{taxon_id}\n")
                for seq, _ in seq_parts:
                    f.write(seq)
                f.write('\n')

    def write_partition_file(self, output_path):
        with open(output_path, 'w') as f:
            for i, (_, start, end, seq_type) in enumerate(self.partition_info, 1):
                f.write(f"{seq_type}, p{i} = {start}-{end}\n")

    def concatenate_files(self, file_list_path, output_matrix, output_partition=None):
        with open(file_list_path, 'r') as f:
            files = [line.strip() for line in f if line.strip()]

        if not files:
            self.printout('error', f"No files found in {file_list_path}")
            sys.exit(1)

        self.printout('metric', [('alignments', len(files))])

        total = len(files)
        for i, filepath in enumerate(files, 1):
            self.add_alignment(filepath)
            self.printout('progress', f"Concatenating: {i:0{len(str(total))}d}/{total} alignments processed")
        print()

        self.printout('metric', [('concat_matrix', output_matrix)])
        self.write_concatenated(output_matrix)

        if output_partition:
            self.printout('metric', [('partition_file', output_partition)])
            self.write_partition_file(output_partition)

        total_length = self.partition_info[-1][2] if self.partition_info else 0
        num_taxa     = len(self.taxa)

        return num_taxa, total_length
