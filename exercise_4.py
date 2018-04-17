#Fasta iterator function

def FASTA_iterator( fasta_filename ):
	""" A generator function that reads a fasta file, and each iteration returns a tupple (identifier, sequence) """

	fi = open( fasta_filename, "r" )
	sequence = ""
	identifier = ""
	for line in fi:
		line = line.strip()
		if line.startswith(">"):
			if sequence != "":
				yield (identifier, sequence)
				sequence = ""
			
			identifier = line[1:]
			
		else:
			sequence += line

	if sequence != "":
		yield(identifier, sequence)

#1. Maximum length

def get_max_sequence_length_from_FASTA_file ( fasta_filename ):
	""" Given a multiline Fasta file, return the length of the sequence with the maximum length. """
	return max([ len(seq[1]) for seq in FASTA_iterator(fasta_filename) ])
	 
#2. Minimum length

def get_min_sequence_length_from_FASTA_file ( fasta_filename ):
	""" Given a multiline Fasta file, return the length of the sequence with the minimum length. """
	return min([ len(seq[1]) for seq in FASTA_iterator(fasta_filename) ])

#3. Tuples of sequences with maximum length

def get_longest_sequence_from_FASTA_file( fasta_filename ):
	"""Given a Fasta file, it returns a list of tuples (identifier, sequence)
	corresponding to sequence(s) with maximum length, sorted by identifier """
	max_len = get_max_sequence_length_from_FASTA_file(fasta_filename)
	longest_seqs = [(identifier, sequence) for identifier, sequence in FASTA_iterator( fasta_filename) if len(sequence) == max_len ]
	return sorted(longest_seqs, key=lambda x: x[0])

#4. Tuple of sequences with minimum length

def get_shortest_sequence_from_FASTA_file( fasta_filename ): 
	"""Given a Fasta file, it returns a list of tuples (identifier, sequence)
	corresponding to sequence(s) with minimum length, sorted by identifier """
	min_len = get_min_sequence_length_from_FASTA_file(fasta_filename)
	shortest_seqs = [(identifier, sequence) for identifier, sequence in FASTA_iterator( fasta_filename) if len(sequence) == min_len ]
	return sorted(shortest_seqs, key=lambda x: x[0])

#5. Molecular weights of all proteins

def get_molecular_weights (fasta_filename):
	""" Given a protein Fasta file, returns a dictionary with molecular weights of all proteins.
	Keys would be protein identifiers and associated values with float type """
	
	aminoacid_mw = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 
	'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}

	molecular_weight = {}
	for ids, seq in FASTA_iterator(fasta_filename):
		molecular_weight[ids] =  sum(aminoacid_mw[x] for x in seq)

	return molecular_weight

#6. Highest molecular weight

def get_sequence_with_max_molecular_weight( fasta_filename ):
	""" Given a protein FASTA file, returns a tupple (identifier, sequence) of the proteins	
	with Highest molecular weight. If there are more than one proteins with same maximum 
	molecular weight, just return the first one """

	id_list = [identifier for identifier, sequence in FASTA_iterator(fasta_filename)]
	molecular_weight = get_molecular_weights(fasta_filename)
	heaviest_seqs = [(identifier, sequence) for identifier, sequence in FASTA_iterator(fasta_filename) if molecular_weight[identifier] == max(molecular_weight.values())]
	return sorted(heaviest_seqs, key= lambda x: id_list.index(x[0]))[0]

