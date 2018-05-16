#!/usr/bin/env python 
__version__= '0.0.0'

import numpy as np
import copy

#Periodic table

periodic_table = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
    "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
    "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]

#Atomic mass according to CIAAW http://www.ciaaw.org/ for most common isotope, checked up to O

atomic_mass = [0.0000,1.0078,4.0026,7.0160,9.0122,11.009,12.000,14.003,15.995,18.998,20.180,22.990,"Mg",26.982,27.977,30.974,"S","Cl",39.948,39.098,40.078, 44.956,47.867,50.941,51.996,54.938,55.845,58.933,58.693,63.546,65.38,69.723,72.630,74.922,"Se","Br","Kr","Rb","Sr","Y","Zr",
    "Nb","Mo","Tc","Ru",102.91,106.42,107.87,"Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
    "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]

	
#get mass of molecule from atom list and also returns a list with masses of each atom for mass weighted centering
def get_mass(atoms):
	mass_list =[]
	for item in atoms:
		mass_list.append(atomic_mass[periodic_table.index(item)])
	print(mass_list)
	mass=sum(mass_list)
	return mass, mass_list


#function to get xyz coordinates and atom list out of *.log or *.xyz files
def get_xyz():

	return

#checks if a string contains only an integer
def isInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False	
	
#gets xyz structures from xyz file and provides a list of numpy arrays with xyz data of each structure, number of atoms and an atom list 
def xyz_from_xyz(file):
	n_atoms = 0
	atoms = []
	structures = []
	input = open(file, 'r')
	#search for number of atoms
	for line in input:
		if isInt(line.strip()):
			n_atoms=int(line)
			break	
	else: #exits if no line with number of atoms was found
		exit('No xyz coordinates found in file: ' + file)
			
	#skip one line
	input.readline()
	
	# now there should be n_atoms lines of coordinates WHAT IF NOT???
	for i in range(n_atoms):
		l=input.readline().split()
		
		atoms.append(l[0]) #get atom symbol and append to atom list
		
		coords=[float(x) for x in l[1:]] #convert line to list of floats
		coords=np.array([coords]) #create array with coords
		try: #try append, doesn't work if XYZ doesn't exist yet
			XYZ=np.concatenate((XYZ,coords), axis=0)
		except NameError:
			XYZ=coords
			
			
	structures.append(XYZ) #append first structure to structures list
	del XYZ #get rid of that for the next structure
	
	
	#now search for more structures
	
	
		
	
	return structures, n_atoms, atoms
	
#function to reverse order of IRC, takes list of numpy arrays and reverses them
def reverse_IRC(list):
	list.reverse()
	return list

#function that centers xyz structure by center of mass or centroid, takes structure and mode of centering and returns centered structure
def center_xyz(structure, atoms, mode):
	
	#first calculate translation vector V by calculating center and reversing the vector
	if mode == 'c':
		V=-1*find_centroid(structure)
	elif mode == 'm':
		V=-1*find_centerofmass(structure, atoms)
	else:
		exit('something went wrong')
	
	#now add vector to each atom
	structure=structure+V
	return structure

#function that finds center of mass, input is geometry P and atom list
def find_centerofmass(P, atoms):
	S=copy.deepcopy(P)
	#first get mass of molecule and mass list
	mass,mass_list=get_mass(atoms)
	#now multiply each line by it's mass
	for i in range(len(mass_list)):
		S[i,:] *=mass_list[i]	
	C=S.mean(axis=0)
	C[:]=C[:]/(mass/len(mass_list))
	return C

#function that finds centroid from numpy array, returns numpy array with position of centroid
def find_centroid(P):
	C=P.mean(axis=0)
	return C

	
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
	
	
	
	Arguments:
	
    mass weight vs centroid
	include H or not
	reverse IRC1 
	reverse IRC2
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
