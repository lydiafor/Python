#1. Extract data proteins in a file
def count_sequences_by_residue_threshold(filename, residue, threshold=0.03):
	""" given a multi-Fasta file, returns the total number of proteins having 
	a relative frequency higher or equal than a given threshold for a given residue """
	fd = open(filename, "r")
	proteins = {}
	total_prot = 0

	for line in fd:
		line = line.strip()
		if line.startswith(">"):
			prot_id = line[1:]
			proteins[prot_id] = ""
	
		else:
			proteins[prot_id] += line

	for seq in proteins.values():
		residue_count = seq.count(residue)
		length = len(seq)
		frequency = residue_count / length
		if frequency >= threshold:
				total_prot += 1

	fd.close()
	return total_prot
	
print(count_sequences_by_residue_threshold("prot.fasta", "L"))			

#2. Extract protein first and final aa and their absolute frequency

def print_sequences_tails(filename, output_filename, first_n=10, last_m=10):
	""" given a protein multi-Fasta file, return an outfile with the 
	protein identifier first and final n aa and their absulute frequency """

	fd = open(filename,"r")
	fo = open(output_filename, "w")
	proteins = {}
	prot_n = 0

	for line in fd:
		line = line.strip()
		if line.startswith(">"):
			prot_id = line[1:]
			prot_n += 1
			proteins[prot_id] = ""
		
		else:
			proteins[prot_id] += line

	fo.write("The file " + filename + " contains " + str(prot_n) + 
				" proteins. Here we show the code of the protein, the first " +
				str(first_n) + " aminoacids of each protein and the last " + 
				str(last_m) + " aminoacids\n")

	for id in proteins:
		seq = proteins[id]
		first_aa = seq[:first_n]
		last_aa = seq[-last_m:]
		total_aa = first_aa + last_aa
		output = ""
		for aa in total_aa:
			if aa not in output:
				output += aa + ":"
				aa_count = seq.count(aa)
				output += str(aa_count) + ","

		fo.write(id + "\t" + first_aa + "\t" + last_aa + "\t" + output + "\n")

	fd.close()
	fo.close()

