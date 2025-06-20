from typing import List, Iterator
from dataclasses import dataclass
from Bio.Seq import Seq

offset = 33

@dataclass
class Sequence:
	name: str = ""
	seq: str = ""
	qualstr: str = ""
	qualarr: List[int] = None

	def __post_init__(self):
		if self.qualarr is None:
			self.qualarr = []

	def set_qualstr(self, qual: str) -> None:
		self.qualstr = qual
		if not self.qualarr:
			self.qualarr = [ord(j) - offset for j in self.qualstr]
			for quality_score in self.qualarr:
				assert 0 < quality_score <= 41, "Change the offset in seq.py\nqual"

	def get_fastq(self) -> str:
		return f"@{self.name}\n{self.seq}\n+\n{self.qualstr}\n"

	def get_fasta(self) -> str:
		return f">{self.name}\n{self.seq}\n"

	def rev_comp(self) -> None:
		seq_obj = Seq(self.seq)
		self.seq = str(seq_obj.reverse_complement())

def fastq_generator(infile) -> Iterator[Sequence]:
	if isinstance(infile, str):
		with open(infile, 'r') as f:
			yield from _fastq_generator_internal(f)
	else:
		yield from _fastq_generator_internal(infile)

def _fastq_generator_internal(infile) -> Iterator[Sequence]:
	line = infile.readline()
	while line:
		if line.startswith('@'):
			name = line[1:].strip()
			seq = infile.readline().strip()
			infile.readline().strip()
			qual = infile.readline().strip()
			tseq = Sequence(name=name, seq=seq)
			tseq.set_qualstr(qual)
			yield tseq
		line = infile.readline()

def read_fasta_file(filename: str) -> List[Sequence]:
	seqlist = []
	templab = ""
	tempseq = ""
	first = True
	
	with open(filename, "r") as fl:
		for line in fl:
			if line.startswith(">"):
				if not first:
					seqlist.append(Sequence(templab, tempseq))
				templab = line.strip()[1:]
				tempseq = ""
				first = False
			else:
				tempseq += line.strip()
	
	seqlist.append(Sequence(templab, tempseq))
	return seqlist