#1. Fasta generator function
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

#2. Fasta dictionary generator function

def compare_fasta_file_identifies( fasta_filename_list ):
	"""given a list of Fasta files, this function returns a dictionary that contains the following keys: 
		intersection: set with common identifiers found in all files
		union: set with all indetifiers (unique)
		frequency: dictionary with all identifiers as keys and number of files in 
					which it appears as values (int)
		Specific: dictionary with name of inputs files as keys and a set with
					specific identifies as values."""

	output_dict = { "intersection" : set(),
			"union" : set(),
			"frequency" : dict(),
			"specific" : dict()
			}

	ids_prot = {}

	#Store identifiers into a dicctionary
	for fasta in fasta_filename_list:
		ids_prot[fasta] = set()
		for iden, seq in FASTA_iterator(fasta):
			ids_prot[fasta].add(iden.upper()) #not difference betwen lower and upper cases

	#Get the intersection, union, frequency and specific.
	n = 0
	for fasta, ids in ids_prot.items():
		n += 1
		if n == 1:
			intersect = set(ids)
			union = set(ids)
	
		else:
			intersect.intersection_update(ids)
			union.update(ids)

		for element in ids:
			if element not in output_dict["frequency"]:
				output_dict["frequency"][element] = 1
			else:
				output_dict["frequency"][element] += 1

		specific = set(ids)
		for fasta2, ids2 in ids_prot.items():
			if fasta2 == fasta:
				continue
			else:
				specific.difference_update(ids2)

		output_dict["specific"][fasta] = specific

	output_dict["intersection"] = intersect
	output_dict["union"] = union

	return output_dict