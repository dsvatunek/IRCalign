#!/usr/bin/env python 
__version__= '0.1.0'

import numpy as np
import copy

#Periodic table

periodic_table = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
    "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
    "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]

#Atomic mass according to CIAAW http://www.ciaaw.org/ for most common isotope, checked up to Rb

atomic_mass = [0.0000,1.0078,4.0026,7.0160,9.0122,11.009,12.000,14.003,15.995,18.998,19.992,22.990,23.985,26.982,27.977,30.974,31.972,34.969,39.962,38.963,39.962,44.956,47.948,50.944,51.941,54.938,55.935,58.933,57.935,62.930,63.929,68.926,73.921,74.922,79.917,78.918,83.911,"Rb","Sr","Y","Zr",
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


#function to get xyz coordinates and atom list out of *.log or *.xyz files, right now only xyz
def get_xyz(file):
	structures,n_atoms,atoms=xyz_from_xyz(file)
	return structures, n_atoms, atoms

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
		exit('Error: No xyz coordinates found in file: ' + file)
			
	#skip one line
	input.readline()
	
	# now there should be n_atoms lines of coordinates WHAT IF NOT???
	for i in range(n_atoms):
		l=input.readline().split()
		
		if l[0] in periodic_table:
			atoms.append(l[0]) #get atom symbol and append to atom list
		else:
			exit('Error: something is wrong with the first structure in file: '+file)
		coords=[float(x) for x in l[1:]] #convert line to list of floats
		coords=np.array([coords]) #create array with coords
		try: #try append, doesn't work if XYZ doesn't exist yet
			XYZ=np.concatenate((XYZ,coords), axis=0)
		except NameError:
			XYZ=coords
				
	structures.append(XYZ) #append first structure to structures list
	del XYZ #get rid of that for the next structure
		
	#now search for more structures
	
	for line in input:
		#start extracting if atom number line is found
		try:
			if int(line.strip()) == n_atoms:
				#read one line to skip title
				input.readline()
				
				# now there should be n_atoms lines of coordinates WHAT IF NOT???
				for i in range(n_atoms):
					l=input.readline().split()
					coords=[float(x) for x in l[1:]]
					coords=np.array([coords])
					try: #try append, doesn't work if XYZ doesn't exist yet
						XYZ=np.concatenate((XYZ,coords), axis=0)
					except NameError:
						XYZ=coords
				structures.append(XYZ)
				del XYZ
		except ValueError:
			pass
				
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

	
#function that calculates RMSD between two xyz structures, takes structures as (m*3) numpy array
def calc_RMSD(A, B):
	RMSD=0.0
	for i in range(3):
		for x in range(len(A)):
			RMSD += (A[x,i]-B[x,i])**2
	return np.sqrt(RMSD/len(A))

	#takes structure A and reference structure R, calculates optimal rotation (R) using SVD for alignment based on Kabsch algorithm!
def find_rotation(A, R):
	# Computation of the covariance matrix
	C = np.dot(np.transpose(A), R)
	
	u, s, vh = np.linalg.svd(C)
	
	#assure right handed coordinate system)
	if (np.linalg.det(u) * np.linalg.det(vh)) < 0.0:
		s[-1] = -s[-1]
		u[:, -1] = -u[:, -1]
	
	#calculate rotation R
	rot = np.dot(u, vh)
		
	return rot
	
# rotates structure, takes strcuture A and rotation R, returns rotated structure B
def rotate(A, rot):
	B=np.dot(A, rot)
	return B

# does partial procrustes analysis, takes structure A and reference R and returns RMSD
def procrustes(A, R):
	rot=find_rotation(A, R)
	#calculates rotated structure B
	B=rotate(A, rot)
	RMSD=calc_RMSD(B,R)
	return RMSD
	
#appends structure to file
def print_xyz(A, atoms, file):
	
	file.write(str(len(atoms))+ "\n\n")
	for i in range(len(A)):
		file.write("{0:2s} {1:15.12f} {2:15.12f} {3:15.12f}\n".format(atoms[i], A[i, 0], A[i, 1], A[i, 2]))
		
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
	#parser =  argparse.ArgumentParser(usage='%(prog)s [options] IRC_1 IRC_2',  description='description')
	#args = parser.parse_args()
	
	
	#HARDCODE FOR TESTING
	
	center_mode='c'
	irc1='irc.xyz'
	irc2='test.xyz'
	
	
	
	#get structures from IRC1
	structures_irc1,n_atoms_irc1,atoms_irc1=get_xyz(irc1)
	
	#get structures from IRC2
	structures_irc2,n_atoms_irc2,atoms_irc2=get_xyz(irc2)
	
	#check if IRCs are compatible
	if atoms_irc1 == atoms_irc2:
		pass
	else:
		exit('Error: Structures don\'t match between IRCs')
	
	# Translate all structures according to either center of mass or centroid
	for i in range(len(structures_irc1)):
		structures_irc1[i]=center_xyz(structures_irc1[i], atoms_irc1, center_mode)
	
	for i in range(len(structures_irc2)):
		structures_irc2[i]=center_xyz(structures_irc2[i], atoms_irc1, center_mode)
		
	""" 
	find out orientation of IRCs automatically by
	checking all 4 different possibilities	
	and reversing IRCs if necessary
	"""
	
	rmsd1=procrustes(structures_irc1[0], structures_irc2[0])
	print('rmsd1=' +  str(rmsd1))
	rmsd2=procrustes(structures_irc1[-1], structures_irc2[0])
	print('rmsd2=' +  str(rmsd2))
	rmsd3=procrustes(structures_irc1[0], structures_irc2[-1])
	print('rmsd3=' +  str(rmsd3))
	rmsd4=procrustes(structures_irc1[-1], structures_irc2[-1])
	print('rmsd4=' +  str(rmsd4))
	
	if rmsd1 < rmsd2 and rmsd1 < rmsd3 and rmsd1 < rmsd4:
		structures_irc1=reverse(structures_irc1)		
		print('Reversed IRC1')
	elif rmsd2 < rmsd1 and rmsd2 < rmsd3 and rmsd2 < rmsd4:
		pass		
	elif rmsd3 < rmsd1 and rmsd2 < rmsd2 and rmsd3 < rmsd4:
		structures_irc1=reverse(structures_irc1)
		structures_irc2=reverse(structures_irc2)
		print('Reversed IRC1 and IRC2')		
	elif rmsd4 < rmsd1 and rmsd2 < rmsd2 and rmsd4 < rmsd3:
		structures_irc2=reverse(structures_irc2)
		print('Reversed IRC2')			
	else:
		print('Warning: Was not able to find correct orientation of IRCs, assuming there are correctly aligned!')
	
	# Calculate rotational matrix
	rot=find_rotation(structures_irc2[0],structures_irc1[-1])
	
	# Apply rotation to all structures in IRC2
	for i in range(len(structures_irc2)):
		structures_irc2[i]=rotate(structures_irc2[i],rot)
	
	# Print one file with translated IRC1 and translated and rotated IRC2 structures
	output = open(irc1 + '_' + irc2 + '.xyz','w')

	#print IRC1
	for i in range(len(structures_irc1)):
		print_xyz(structures_irc1[i],atoms_irc1, output)
	#print IRC2
	for i in range(len(structures_irc2)):
		print_xyz(structures_irc2[i],atoms_irc2, output)
	#close file
	output.close
		
	
	return



if __name__ == "__main__":
    main()
