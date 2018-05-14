#!/usr/bin/env python 
__version__= '0.0.0'

#function to get xyz coordinates out of *.log or *.xyz files
def get_xyz():

	return
	
#function to reverse order of IRC1
def reverse_IRC():

	return

#function that centers xyz structure by center of mass or centroid, takes structure and center of mass or centroid position and returns centered structure
def center_xyz(structure, center):
	
	return structure

#function that finds center of mass
def find_centerofmass():

	return

#function that finds centroid
def find_centroid():

	return

	
#function that calculates RMSD between two xyz structures
def RMSD():
	
	return

	

def main():
	
	import argparse
	import sys
	
	"""
	Assume same atom count and numbering!
	Assume right direction of IRCs (make function that reverses IRCs!) 
	This means the last point of IRC1 is the reference
	This means the first point of IRC2 is used for calculating the rotational matrix!
	"""
	parser =  argparse.ArgumentParser(usage='%(prog)s [options] IRC_1 IRC_2',  description='description')
	
	
	args = parser.parse_args()
	
	#get structures from IRC1
	
	#get structures from IRC2
	
	# Translate all structures according to either center of mass or centroid
	
	# grep last structure of IRC1
	# grep first structure of IRC2
	
	# Calculate rotational matrix
	
	# Apply rotation to all structures in IRC2
	
	# Print one file with translated IRC1 and translated and rotated IRC2 structures
	
	
	return



if __name__ == "__main__":
    main()
