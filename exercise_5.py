#1. Calculate mean of the minimum distance between any two residues pairs.

import sys
import math

def calculate_pdb_chain_mean_minimum_distances(pdb_file_path):
	"""Calculate the mean of the minimum distance between any two residues pairs found in the same chain of a pdb
		The script works with command-line arguments, given the pdb_path """

	pdb = open(pdb_file_path)
	dic = dict()
	output = dict()

	#Store the positions of the residues of the aminoacids into a dictionary
	for line in pdb:
		if line.startswith("ATOM"):
			#Definition of attributes
			elements = line.split()
			chain = elements[4]
			n_aa = elements[5]
			x_pos = float(elements[6])
			y_pos = float(elements[7])
			z_pos = float(elements[8])
			pos_set = (x_pos, y_pos, z_pos)

			if chain not in dic:
				dic[chain] = {}

			if n_aa not in dic[chain].keys():
				dic[chain][n_aa] = [pos_set]
			else:
				dic[chain][n_aa].append(pos_set)
	
	#Iterate the dictionary to calculate minimum distances
	for chain in dic.keys():
		values = list(dic[chain].keys())
		n_tot = len(values)
		min_dist_tot = []
		n = 0
		
		#For each residue position
		while n < n_tot - 1:
			coords = dic[chain][values[n]]

			#Iterate again residues positions to get all comparisons
			for m in range(n+1, n_tot): 
				coords2 = dic[chain][values[m]]
				res_dist = []

				#For each coordenate in position n calculate the distances with each coordinate in position m
				for x,y,z in coords:
					distances = [ ( (x - x2)**2 + (y - y2)**2 + (z - z2)**2 ) for x2,y2,z2 in coords2 ]
					res_dist.extend((distances))

				#Get the minimum distance for each n,m pair of residues
				min_dist = math.sqrt(min(res_dist))
				min_dist_tot.append(min_dist)
		
			n += 1
		
		#Compute the mean od distances for each chain
		mean_min_dist = sum(min_dist_tot) / len(min_dist_tot)

		output[chain] = mean_min_dist

	return output

#Execution control
if __name__ == "__main__":

	try:
		file = sys.argv[1]
	except:
		file = input("Enter a path of a pdb: ")

	dic = calculate_pdb_chain_mean_minimum_distances(file)

	for chain, mean_dist in dic.items():
		print("%s: %.4f" %(chain, mean_dist))