#Import modules
from sequence import *
import re
import sys
import os
import argparse
import gzip
import random

#Define command-line arguments parser

parser = argparse.ArgumentParser(description="""This program read the input DNA Fasta file(s) 
											and calculate the length and molecular weight of their corresponding proteins.""")

parser.add_argument('-o', '--output-file',
					dest = "outfile",
					action = "store",
					default = None,
					help = "Output File")

parser.add_argument('-i', '--input',
					dest = "infile",
					action = "store",
					default = None,
					help = "Input File")

parser.add_argument('-v', '--vervose',
					dest = "verbose",
					action = "store_true",
					default = False,
					help = "Print log in stderr")

parser.add_argument('-p', '--pattern',
					dest = "pattern",
					action = "store",
					default = False,
					help = "Pritn only sequences with match to that pattern")

parser.add_argument('-r', '--random',
					dest = "number",
					action = "store",
					type = int,
					default = False,
					help = "Integer specifying the number of sequences to show")

options = parser.parse_args()

##Function definitions

#Fasta iterator

def FASTA_iterator(filename, sequence_class):

	if filename.endswith(".gz"):
		fd = gzip.open(filename, 'rt')
	else:
		fd = open(filename)

	sequence = ""
	for line in fd:
		if line[0]==">":
			if len(sequence)>0:
				try:
					yield sequence_class(identifier, sequence)
				except IncorrectSequenceLetter as e:
					sys.stderr.write(str(e))

			identifier = line[1:].strip()
			sequence = ""
		else:
			sequence+=line.strip()
	
	if len(sequence)>0:
		try:
			yield sequence_class(identifier, sequence)
		
		except IncorrectSequenceLetter as e:
			sys.stderr.write(str(e))

		finally:
			fd.close()


#Get the correspond fasta files of input

def get_input_file(input):
	""" Handling with different kind of input: only fasta or gunzip fasta files 
	for a given path or for the current directory """
	fasp = re.compile('.fa$|.fasta$|.fa.gz$|.fasta.gz$')
	path = input

	if path == None:
		path = os.getcwd()

	if os.path.isdir(path):
		fasta_files = [f for f in os.listdir(path) if fasp.search(f) is not None]
		os.chdir(path)
	else:
		fasta_files = [input]

	return fasta_files

#Open different kind of output files

def get_output_file(output):
	""" Handling with different kind of output: None, normal file or gunzip file """

	if output == None:
		outfile = sys.stdout
	else:
		if ".gz" in output:
			outfile = gzip.open(output, "wt")
		else:
			outfile = open(output, "w")

	return outfile

#Function to obtain the Protein sequences of each DNA sequence of FASTA file, its length and molecular weight from a list of fasta file

def get_protein_len_mw(fasta_list, pattern=False):
	"""Return a list with the identifier of Proteins, sequence length and molecular weight 
	from a list of fasta files. If there is a pattern return only the sequences that match with that pattern
	Also return the number of fasta files, total sequences and protein sequences"""

	total_prot = []
	n = 0
	m = 0

	for fasta in fasta_list:
		n += 1
		for dna in FASTA_iterator(fasta, DNASequence):
			m += 1
			prot = dna.translate()
			len_prot = len(prot)
			mw_prot = prot.get_mw()

		#Check if there is a pattern and output only those sequences which contains that pattern
			if pattern != False:
				p = re.compile(pattern)
				if p.search(prot.get_sequence()) is not None:
					total_prot.append(tuple([prot.get_identifier(), len_prot, mw_prot]))
				else:
					continue
			else:
				total_prot.append(tuple([prot.get_identifier(), len_prot, mw_prot]))

	l = len(total_prot)
	return total_prot, n, m, l 

#Function to get N random protein sequences

def get_n_proteins(prot_list, n):
	"""Return N random proteins from the proteins lists """
	if n > len(prot_list):
		print("""Random number is higher than actual selected proteins. Returning all proteins""")
		random_prot = prot_list
	else:
		random_prot = random.sample(prot_list, n)
	return random_prot


#Function to get the progress log of the program

def get_progres_log(n_files, fasta_files, n_sequences, n_prots, n_random, outfile):
	"""Return a message in standard error with the progression log of the program """
	if outfile == None:
		outfile = "STDOUT"

	sys.stderr.write("%d FASTA files found:\n" %(n_files))
	for fasta in fasta_files:
		sys.stderr.write("%s\n" %(fasta))

	sys.stderr.write("""\n%d DNA sequences found.\n
%d Protein sequences matched regular expression.\n
%d Protein sequences selected randomly\n
Printing output to %s.\n
Program finished correctly.\n\n"""
					%(n_sequences, n_prots, n_random, outfile))

#Function to print data into outfile

def print_results(prot_list, outfile):
	"""Printing the results to selected outfile"""
	for data in sorted(prot_list, key=lambda x: x[1], reverse=True):
		data = "%s\t%d\t%.3f\n" %(data[0], data[1], data[2])
		outfile.write(data)


#Execution Control

if __name__ == "__main__":
	fasta_files = get_input_file(options.infile)
	outfile = get_output_file(options.outfile)

	info_prots = get_protein_len_mw(fasta_files, options.pattern)
	prots = info_prots[0]
	n_files = info_prots[1]
	n_sequences = info_prots[2]
	n_proteins = info_prots[3]

	if options.number:
		prots = get_n_proteins(prots, options.number)
		n_random = len(prots)
	else:
		n_random = 0

	if options.verbose:
		get_progres_log(n_files, fasta_files, n_sequences, n_proteins, n_random, options.outfile)

	print_results(prots, outfile)