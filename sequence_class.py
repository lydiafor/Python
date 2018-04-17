import sequence_data

class Sequence(object):
	"""
	Generic class to define a BioPolymer.
	A biopolymer is composed by an identifier and a sequence of monomers
	"""

	alphabet = set()	# Class attribute to define the possibles monomers
	mw = {}			# Class attribute to define the monomer molecular weights

	def __init__(self, identifier, sequence):

		self.__identifier = identifier

		for letter in sequence:
			if letter not in self.alphabet:
				raise IncorrectSequenceLetter(letter, self.__class__.__name__)
				 

		self.__sequence = sequence
		self.__mw = None

	def get_identifier(self):
		"""
		Getter method to obtain the identifier of the sequence
		"""
		return self.__identifier

	def get_sequence(self):
		"""
		Getter method to get the sequence string
		"""
		return self.__sequence

	def get_mw(self):
		"""
		Calculate the molecular weight of the Sequence instance
		as the sum of the molecular weight of all the monomers
		Return a float number
		"""
		if self.__mw is None:
			self.__mw = sum( self.mw[letter] for letter in self.__sequence )
		return self.__mw

	def has_subsequence(self, sequence_obj ):
		"""
		Check if the sequence of sequence_obj is contained in 
		the Sequence object
		"""
		return sequence_obj.get_sequence() in self.__sequence

	def __len__(self):
		return len(self.__sequence)

	def __lt__(self, other):
		return self.get_mw() < other.get_mw()

	def __eq__(self, other):
		return self.get_sequence() == other.get_sequence()

	def __getitem__(self, key):
		return self.__sequence[key]

	def __add__(self, other):
		if type(self) == type(other):
			return type(self)(identifier = "%s+%s" %(self.get_identifier(),
                                                                     other.get_identifier()),
                                              sequence = self.__sequence + other.get_sequence())
		else:
			raise TypeError("It is not possible to concatenate different types of sequences")

	def __contains__(self, string):
		return string in self.__sequence

	def __hash__(self):
		return hash((self.__identifier, self.__sequence))
 
#Sequences subclasses
class ProteinSequence(Sequence):
	"""
	Protein Sequence object. The monomers of ProteinSequence are
	the aminoacids
	"""
	
	#Override all specific Class attributes for ProteinSequence
	alphabet = set(sequence_data.protein_letters)
	mw = sequence_data.protein_weights


class NucleotideSequence(Sequence):
	"""
	Nucleotide Sequence Object. The monomers of NucleotideSequence
	are not defined, as they can be ribonucleotides or deoxyribonucleotides
	"""

	complement = {}
	start_codons = set()
	end_codons = set()

	translation_dict =  {}
	
	def translate(self):
		started = False
		translated_sequence = ""
		sequence = self.get_sequence()
		for i in range(0,len(sequence),3):
			codon = sequence[i:i+3]
			if started is False:
				if codon in self.start_codons:
					started = True
					translated_sequence = self.translation_dict[codon]
			elif codon in self.stop_codons:
				break
			else:
				translated_sequence+=self.translation_dict[codon]
		return ProteinSequence( identifier = self.get_identifier()+"_translated",
					sequence = translated_sequence )


class DNASequence(NucleotideSequence):

	alphabet = set(sequence_data.dna_letters)
	mw = sequence_data.dna_weights

	complement = sequence_data.dna_complement
	stop_codons = set(sequence_data.dna_stop_codons)
	start_codons = set(sequence_data.dna_start_codons)

	translation_dict = sequence_data.dna_table
	
	def transcribe(self):
		return RNASequence( identifier = self.get_identifier()+"_transcribed", sequence = self.get_sequence().replace("T","U") )


class RNASequence(NucleotideSequence):
	alphabet = set(sequence_data.rna_letters)
	mw = sequence_data.rna_weights
	complement = sequence_data.rna_complement

	stop_codons = set(['UAA', 'UAG', 'UGA'])
	start_codons = set(['UUG', 'CUG', 'AUG'])

	translation_dict = sequence_data.rna_table

	def reverse_transcribe(self):
                return DNASequence( identifier = self.get_identifier()+"_transcribed", sequence = self.get_sequence().replace("U","T") )

#Value Error exception subclass
class IncorrectSequenceLetter(Exception):
	""" A ValueError subclass to handle with problems when adding a new sequence instances """
	def __init__(self, letter, class_name):
		self.letter = letter
		self.class_name = class_name

	def __str__(self):
		return "The sequence item %s is not found in the alphabet of class %s\n" %(self.letter, self.class_name)


