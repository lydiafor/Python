# Python

In this repository you can find the exercises for the python subject of the Master in Bioinformatics for the Health Science (UPF)

**Exercise 1:** 

1. Given a multi-line protein FASTA file, return an intweger corresponding to the total number of protein sequences having a relative frequency higher or equal that a given threshold for a given residue.

2. Given a protein FASTA file, save on a file the protein identifier, first N-aminoacids, the last M-aminoacids and the absolute frequency in the protein of each of the N/M-aminoacids. The first three fields must be separated by a tabulator, the absolute frequency of the residues must have the format residue (frequency and must be separated by comma). One protein by line. The first line must have a summatory of the information.

  fil: *prot.fasta / dna.fasta*

**Exercise 2:** 

Create a script that, given a multi-line protein FASTA file and a "sub-sequences" file, calculates the proportion of proteins in the FASTA file containing at least N-times each of the sub-sequences. Saving it in an output file with a specified format and ordered by the proportion value.
  
  file: *prot.fasta, subseq.txt*
  
**Exercise 3:**

1. Build a generator function that reads a Fasta file, and returns a tupple with (identifier, sequence) in each iteration

2. Given a list of FASTA files, create a function that returns a dictionary with the following keys and their associated values:
  
   - *intersection:* a set iwth the common identifiers found in all the files.
   - *union*: set with all identifiers (uniques) found in all files.
   - *frequency:* dictionary with all the identifiers as keys and the number of files in which it appear as values.
   - *specific:* dictionary with the name of input files as keys and a set with the specific identifiers as values.
   
   file: *sample_fasta1.fa, sample_fasta2.fa, sample_fasta3.fa*
   
**Exercise 4:**

Given a multiline FASTA file, compute:
  1. Length of the sequence with maximum length
  2. Length of the sequence with minimum length
  3. List of tuples (identifier, sequence) corresponding to the sequences with maximum length. Sorted by identifier
  4. List of tuples (identifier, sequence) corresponding to the sequences with minimum length. Sorted by identifier
  5. For proteins, returns a dictionary with the molecular weights of all proteins in file. 
  6. For proteins, a tuple (identifier, sequence) of the protein with the highest molecular weight
  
  file: *prot.fasta*
  
**Exercise 5**

Create a python script that calculates the mean of the minimum distance between any two residues pairs found in the same chain of a PDB. The script, when executed by command line, should output in standard output the mean distance for each chain (with 4 decimal positions). The python script should use a single argument corresponding to the PDB file path to use. This command line argument is optional. If the PDB file path is not defined, read the PDB file from standard input.
  
  file: *1e9h.pdb*

**Exercise 6**

1. Define a new class named Protein, with the identifier and sequence as attributes, and with the methods: get_identifier, get_sequence, get_mw, has_subsequence, get_length.

2. Modify the Fasta generator function to return protein objects instead of tuples.

  file: *prot.fasta*

**Exercise 7**

1. Read the input DNA Fasta file and calculate the length and molecular weight of their corresponding proteins.
Command line parameters:

    -o --output-file: output file. If not defined, output is printed to standard output. If
defined, it prints the output to this file. If output file endswith “.gz”, output must be
gzipped.

    -i --input: if this argument is defined and it corresponds to a file, use it as input. If it
corresponds to a directory, use all “.fa”, “.fasta”, “.fa.gz” and “.fasta.gz” files in the defined directory. If this argument is not defined, use all “.fa” and “.fasta” files in the current directory.

    -v --verbose: if this argument is defined, the progression log is printed to standard error. If it is not used, the log is not showed.

    -p --pattern: string with a regular expression. If defined, output only the sequences having the given regular expression in the ProteinSequence.

    -r --random: integer defining the number of sequences to be printed in the output. If it not set, the complete output must be printed. If it is defined, a random selection of the defined size has to be printed.
 
Output should be sorted by sequence length, from longest to shortest.

2. Define and implement the following classes: Sequence, ProteinSequence, NucleotideSequence, DNASequence and RNASequence.
Conditions:
- alphabet must be a class attribute that specifies the possible alphabet of the sequence

- when creating a ne Sequence instance, it must check that the sequence is correct by cheking in the alphabet. If not, raise an exception with the incorrect letter.

- Required data for attributes are in sequence_dictionary.py

Methods: get_identifier, get_sequence, get_mw, has_subsequence, translate, transcribe, reverse_transcribe, len(sequence), sequence1 == sequence2, sequence1 != sequence2, sequence + sequence, sequence[i], in operator, compating sequences according to their molecular weights, hash the classses.

3. Createa new ValueError exception subclass names *IncorrectSequenceLetter*, when a letter is not found in the alphabet and class name of the sequence. Modify the Fasta generator function to skip sequences having incorrect letters (print a meesage erros in the standard error and continue with the next sequence).

  files: *sequence_class.py, sequence_data.py, dna.fasta*
