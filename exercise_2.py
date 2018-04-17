#1. Number of subsequences contained in multi-line fasta file

def calculate_aminoacid_frequencies(fasta_filename, subsequences_filename, number_of_repetitions, output_filename):
	""" This functions takes as inputs a multi-line protein Fasta file, and a subssequences file,
        calculates the proportion of proteins in FASTA cotaining at least N-times each of the sub-sequences, 
	te save it in an output file separated by tabulators and ordered by the proportion values"""

	ff = open(fasta_filename, "r")
	fs = open(subsequences_filename, "r")
	fo = open(output_filename, "w")

	seqs = {}
	count_subseqs = {}
	prot_n = 0
	subseq_n = 0

	#Store fasta information into a dictionary, count number of proteins
	for line in ff:
		line = line.strip()
		if line.startswith(">"):
			prot_id = line[1:]
			prot_n += 1
			seqs[prot_id] = ""
		else:
			seqs[prot_id] += line

	#Store subsequences file information into a list and count number of subsequences
	for line in fs:
		line = line.strip()
		count_subseqs.setdefault(line, 0)
		subseq_n += 1

	#Write the number of proteins and subsequences in output file
	fo.write("#Number of proteins:\t\t%3s\n#Number of subsequences:\t%3s\n#subsequence proportions:\n" 
		%(prot_n, subseq_n))
	
	#iterate both data structures and check the number of times a protein contains a motif
	for subseq in count_subseqs:
		for seq in seqs.values():
			count = seq.count(subseq)
			if count >= number_of_repetitions:
				count_subseqs[subseq] += 1

	#Write the results and computhe the relative frequencie in sorted values.
	for subseq, count in sorted(count_subseqs.items(), key=lambda x: x[1], reverse=True):
		rel_freq = count / prot_n
		fo.write("%s\t%10d\t%.4f\n" %(subseq, count, rel_freq))

	ff.close()
	fs.close()
	fo.close()

calculate_aminoacid_frequencies("prot.fasta", "subseq.txt", 1, "test.txt")

