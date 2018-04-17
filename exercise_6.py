#1. Define a Protein Class

class Protein(object):
	"""This class represent proteins with their identifier and sequence
		And several defined methods"""
	def __init__(self, identifier, sequence):

		self.identifier = identifier
		self.sequence = sequence
	
	def get_identifier(self):
		"""Return the identifier of the protein"""
		return self.identifier

	def get_sequence(self):
		"""Return the sequence the protein"""
		return self.sequence

	def get_mw(self):
		"""Return the molecular weight of the protein"""
		aminoacid_mw = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 
		'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
		
		molecular_weight = sum(aminoacid_mw[aa] for aa in self.sequence)
		return molecular_weight

	def has_subsequence(self, protein_subsequence):
		"""Return if the protein sequence has the entered subsequence"""
		if protein_subsequence.upper() in self.sequence:
			return True
		else:
			return False

	def get_length(self):
		"""Return the length of the protein"""
		return len(self.sequence)


#2. Fasta_iterator generator function that returns Protein objects

def FASTA_iterator( fasta_filename ):
	""" A generator function that reads a fasta file, and each iteration returns a Protein object """

	fi = open( fasta_filename, "r" )
	sequence = ""
	identifier = ""

	for line in fi:
		line = line.strip()
		if line.startswith(">"):
			if sequence != "":
				yield Protein(identifier, sequence)
				sequence = ""
			
			identifier = line[1:]
			
		else:
			sequence += line

	if sequence != "":
		yield Protein(identifier, sequence)

