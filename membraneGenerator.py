#! /opt/rh/python27/root/usr/bin/python2.7

'''
This is a program that creates a membrane file and inserts a protein structure into it if given by the user.
This script need python 2.7 to work
If you have different versions of python installed, the following can be run to use python2.7:
scl enable python27 bash

Dependencies:   1. Requires vmd, chimera and R installed in bin
		2. hacked amberlipid2charmm.py script that converts to amber and recognizes POP residue (instead of POPC)
		3. rpy2 package installed on python
		4. bio3d package installed on R

'''

#Libraries needed

import sys,os
import numpy as np
from Bio.PDB import *
import xpdb
import argparse
from random import uniform
from random import randint
from math import pi
from subprocess import call
from rpy2 import *
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

#global variable intialization
pdbp = 1 
io = 1
bad_res = {}
counters = {}
counters['Good'] = 0
counters['Bad'] = 0

class ResSelect(Select):
    def accept_residue(self, res):
        try:
            bad_res[(res.id,res.resname)]
            counters['Good']+=1
            return False
        except KeyError:
            counters['Bad']+=1
            return True

class OppSelect(Select):
    def accept_residue(self, res):
        try:
            bad_res[(res.id,res.resname)]
            return True
        except KeyError:
            return False

def AlignProtein(ref,sample):
	'''This function uses a refernce protein provided by the user to align the new protein'''

	#write a script "chimera_align.py" to help alignment of the proteins
	with open("chimera_align.py","wb") as of:
		of.write('import os\n')
		of.write('import sys\n')
		of.write('from chimera import runCommand as rc\n')
		of.write('ref = sys.argv[-2]\n')
		of.write('sample = sys.argv[-1]\n')
		of.write('rc("open " + ref)\n')
		of.write('rc("open " + sample)\n')
		of.write('rc("mm #0 #1")\n')
		of.write('rc("write #1 "+ "chimera_aligned_"+sample)\n')
		of.write('rc("close all")\n')
		of.write('rc("stop now")\n')
		of.write('print "Chimera output written"\n')	
	
	if os.path.exists("chimera_align.py"):
		#call the script in command line to run the above written script
		call(['chimera','--nogui','chimera_align.py',ref,sample])
		#remove the script
		call(['rm','chimera_align.py'])
	else:
		print "Something went wrong. Exiting..."
		sys.exit()

def ReformatProtein(protein):
	'''This function reformats the protein into amber readable format'''

	bio3d = importr('bio3d')
	read = robjects.r['read.pdb']
	pdb = read(protein)

	convert = robjects.r['convert.pdb']
	new = convert(pdb,type="amber")

	write = robjects.r['write.pdb']
	write(new,file="amber_formatted"+protein)

	print 'Protein converted into amber readable format...'

def SetUp():
	'''This function sets up readers and parsers'''

	global pdbp,io
	pdbp = PDBParser(QUIET = True, PERMISSIVE=True, structure_builder=xpdb.SloppyStructureBuilder())
	io = xpdb.SloppyPDBIO()

def CreateMembrane(size,membrane_file):
	'''This function creates the membrane of size provided by the user'''

	print "Creating Membrane..."
	FNULL = open(os.devnull,'w')

	#Write a script "create_membrane.tcl"
	with open("create_membrane.tcl","wb") as of:
		of.write("set m_size [lindex $argv 1]\n")
		of.write("set b_max [expr $m_size/2 + 2]\n")
		of.write("set b_min [expr -1 * $b_max]\n")
		of.write("package require membrane\n")
		of.write("package require solvate\n")
		of.write("package require autoionize\n")
		of.write('membrane -l POPC -x $m_size -y $m_size -o [lindex $argv 0] -top "c36"\n')
		of.write('set sel [atomselect top "name P and resname POPC"]\n')
		of.write("set xmin [lindex [measure minmax $sel] 0 0]\n")
		of.write("set ymin [lindex [measure minmax $sel] 0 1]\n")
		of.write("set xmax [lindex [measure minmax $sel] 1 0]\n")
		of.write("set ymax [lindex [measure minmax $sel] 1 1]\n")
		of.write("$sel delete\n")
		of.write("# Solvate box\n")
		of.write('set b_array [list [list $xmin $ymin -40] [list $xmax $ymax 40]]\n')
		of.write('solvate [lindex $argv 0].psf [lindex $argv 0].pdb -b 1.5 -minmax $b_array -o [lindex $argv 0]\n')
		of.write('# Identify average z position of 21st and 31st Oxygen in POPC residue\n')
		of.write('set sel [atomselect top "name O21 O31 and resname POPC and z < 0"]\n')
		of.write('set zmin [lindex [measure center $sel] 2]\n')
		of.write('$sel delete\n')
		of.write('set sel [atomselect top "name O21 O31 and resname POPC and z > 0"]\n')
		of.write('set zmax [lindex [measure center $sel] 2]\n')
		of.write('$sel delete\n')
		of.write('# Set all betas zero in preparation of identify\n')
		of.write('set all [atomselect top all]\n')
		of.write('$all set beta 0\n')
		of.write('$all delete\n')
		of.write('# Identify water molecules in membrane\n')
		of.write('set wat [atomselect top "water and noh and z < $zmax and z > $zmin"]\n')
		of.write('$wat set beta 1\n')
		of.write('$wat delete\n')
		of.write('set del [atomselect top "beta 1"]\n')
		of.write('foreach seg [$del get segname] res [$del get resid] {\n')
		of.write('    set small_sel [atomselect top "segname $seg and resid $res"]\n')
		of.write('    $small_sel set beta 1\n')
		of.write('    $small_sel delete\n')
		of.write('}\n')
		of.write('$del delete\n')
		of.write('package require psfgen\n')
		of.write('readpsf [lindex $argv 0].psf\n')
		of.write('coordpdb [lindex $argv 0].pdb\n')
		of.write('set del [atomselect top "beta 1"]\n')
		of.write('foreach seg [$del get segname] res [$del get resid] {\n')
		of.write('  delatom $seg $res\n')
		of.write('}\n')
		of.write('$del delete\n')
		of.write('writepsf [lindex $argv 0].psf\n')
		of.write('writepdb [lindex $argv 0].pdb\n')
		of.write('autoionize -psf [lindex $argv 0].psf -pdb [lindex $argv 0].pdb -sc 0.15 -cation POT -o [lindex $argv 0]\n')
		of.write('set sel [atomselect top "resname POPC"]\n')
		of.write('$sel set chain L\n')
		of.write('foreach res [lsort -unique [$sel get residue]] {\n')
		of.write('  incr count\n')
		of.write('  set temp [atomselect top "residue $res"]\n')
		of.write('  $temp set resid $count\n')
		of.write('  $temp delete\n')
		of.write('}\n')
		of.write('$sel delete\n')
		of.write('set out [atomselect top "beta 0"]\n')
		of.write('$out writepdb [lindex $argv 0].pdb\n')
		of.write('$out delete\n')
		of.write('quit\n')

	#Call the above written script to run
	call(['vmd','-dispdev','text','-e','create_membrane.tcl', '-args', membrane_file, str(size)], stdout=FNULL)

	#Remove the script
	call(['rm','create_membrane.tcl'])
	print "Membrane created..."
	print "Correcting residue names for membrane file created..."
	SetUp()
	try:
		membrane = pdbp.get_structure('Membrane', membrane_file+'.pdb')
	except IOError:
		print "Unable to read in file"

	#Rename water molecules to match amber specifications
	for residue in membrane[0]['W']:
		residue.resname = 'WAT'
		for atom in residue:
			if atom.id == 'OH2':
				atom.fullname = ' O'
	#Rename ions to match amber specifications
	for residue in membrane[0]['I']:
		if residue.resname == 'CLA':
			residue.resname = 'Cl-'
			for atom in residue:
				atom.fullname = 'Cl-'    
		if residue.resname == 'POT':
			residue.resname = 'K+'
			for atom in residue:
				atom.fullname = 'K+'
	io.set_structure(membrane)
	io.save(membrane_file + '.pdb')

	#Create amber version of membrane file
	call(['python', '/home/vvenkat6/ChemoreceptorProject/Scripts/charmmlipid2amber.py','-c','/home/vvenkat6/ChemoreceptorProject/Scripts/charmmlipid2amber/charmmlipid2amber.csv','-i',membrane_file+'.pdb','-o',membrane_file+'.amber.pdb']) #Hard-coded path

	print "Amber format membrane generated..."

	#Removing temporary files:
	call(['rm',membrane_file+'.pdb'])
	call(['rm','amber_formattedchimera_aligned_sample'])
	call(['rm','chimera_aligned_sample'])

def MergePdb(membrane_file,sample):
	'''This function merges the membrane with the protein file'''

	SetUp()
	global bad_res, counters

	try:
		membrane = pdbp.get_structure('Membrane', membrane_file+'.amber.pdb')
	except IOError:
		print "Unable to read in amber membrane file"
		sys.exit()

	try:
		protein = pdbp.get_structure('Protein','amber_formattedchimera_aligned_'+sample)
	except IOError:
		print "Unable to read in amber protein file"
		sys.exit()

	print "Merging the protein and membrane file and removing conflicting molecules.."

	#Prepare membrane file for searching
	atomlist = Selection.unfold_entities(membrane, 'A')
	ns = NeighborSearch(atomlist)

	#Remove all waters that are within 4 angstroms of the protein
	for atom in Selection.unfold_entities(protein, 'A'):
		a,b,c = atom.coord
		neighbors = ns.search(np.array([a,b,c]), 3.5, level = 'R')
		for residue in neighbors:
			membrane[0][residue.get_parent().id][residue.id].full_id=True
			bad_res[(residue.id,residue.resname)]=True

	print "Water molecules within 4 angstroms of the protein have been removed"

	try:
		membrane[0].add(protein[0].get_list()[0])
	except:
		print "Not able to write"

	io.set_structure(membrane)
	io.save('mem_pro.pdb', ResSelect())
	call(['rm',membrane_file+'.amber.pdb'])
	call(['rm',membrane_file+'.log'])
	call(['rm',membrane_file+'.psf'])

def CheckIons(bad_res,membrane):
	'''Function to remove ions and ensure neutrality'''
	count_c = 0
	count_k = 0
	for key in bad_res:
		if key[1]=='Cl-':
                	count_c+=1
		elif key[1]==' K+':
			count_k+=1
	if count_c > count_k:
		num_delete = count_c - count_k
		while num_delete>0:
			rand_index = randint(0,len(membrane[0]['I'].get_list())-1)
                        cand_res = membrane[0]['I'].get_list()[rand_index]
			try:
                                bad_res[(cand_res.id,cand_res.resname)]
                        except KeyError:
                                if cand_res.resname == ' K+':
                                        bad_res[(cand_res.id,cand_res.resname)]=True
                                        num_delete -= 1
	if count_k > count_c:
                num_delete = count_k-count_c
                while num_delete>0:
                        rand_index = randint(0,len(membrane[0]['I'].get_list())-1)
                        cand_res = membrane[0]['I'].get_list()[rand_index]
                        try:
                                bad_res[(cand_res.id,cand_res.resname)]
                        except KeyError:
                                if cand_res.resname == 'Cl-':
                                        bad_res[(cand_res.id,cand_res.resname)]=True
                                        num_delete -= 1
	return(membrane)

	
def AddTer(output,input):
	'''This functions adds appropriate TER cards to the protein'''
	f = open(input,'r')
	o = open(output,'w')

	print "Adding appropriate TER cards"
	current = "1"
	counter = 0
	for line in f:
        	line = line.strip()
       		word = line.split()
        	if word[0] == 'TER':
                	counter = 1
        	if len(word) > 1 and counter == 0:
                	if word[3] == "PC" or word[3] == "PA" or word[3] == "OL":
                        	next = word[5]
                        	if current != next and word[3] == "PA":
                                	o.write( "TER\n")
                                	current = next
                        	o.write(line+'\n')
        	if counter == 1:
                	o.write(line+'\n')
	o.close()
	f.close()
	call(['rm',input])
	
def CheckLipid(output,input):
	'''This function makes sure that the hanging head/tails that are not attached are removed'''
	
	print "Performing final check..."
	f = open(input,'r')
	o = open(output,'w')

	count = 0
	counter = 0
	resnum_current_set = set()
	resnum_set = set()
	residue_set = set()

	for line in f:
        	line = line.strip()
        	word = line.split()

        	if word[0] == 'TER' and counter == 0:
                	if len(residue_set) < 3:
                        	for x in resnum_current_set:
                                	if x in resnum_set:
                                        	resnum_set.remove(x)
                	resnum_current_set = set()
                	residue_set = set()

        	if word[0]!='TER' and counter == 0:             
                	residue_set.add(word[3])
                	resnum_current_set.add(word[5])
                	resnum_set.add(word[5])
                	if word[3] != 'PA' and word[3] != 'PC' and word[3] != 'OL':
                        	counter = 1
                        	break
	f.close()
	f = open(input,'r')
	prev_word = 'ATOM'
	counter = 0
	for line in f:
        	line = line.strip()
        	word = line.split()
        	if word[0] == 'TER':
                	if counter == 0:
                        	o.write(line+'\n')
                        	counter = 1
        	elif word[0] == 'END':
                	o.write(line)
        	elif word[5] in resnum_set:
                	o.write(line+'\n')
                	counter = 0
        	prev_word = word[0]
	o.close()
	f.close() 
	call(['rm',input])

def AddLigands(ligand_list,num_ligands,membrane_file,membrane_size):
	'''This function adds ligands to the membranes'''
	print "Adding Ligands to the membrane"
	SetUp()
	global bad_res, counters

	min_x = -1 * membrane_size/2
	min_y = -1 * membrane_size/2
	z = membrane_size/2

	ligand_spacing = membrane_size/num_ligands
	
	for i in range(1,len(ligand_list)+1):
		print "Inserting ligand number : " + str(i)
		i = i - 1
		try:
                	membrane = pdbp.get_structure('Membrane', membrane_file+'.amber.pdb')
		except IOError:
                	print "Unable to read in amber membrane file"
                	sys.exit()

		#prepare mebrane file for searching
		atomlist = Selection.unfold_entities(membrane, 'A')
        	ns = NeighborSearch(atomlist)

		for x in range(0,num_ligands):
            		for y in range(0,num_ligands):
                		#Calculate displacement vector
               			displace = Vector(min_x+x * ligand_spacing + ligand_spacing/2 , min_y+y * ligand_spacing + ligand_spacing/2, pow(-1,x+y)*3.0/4*z )
		
                		#Read in ligand
                		try:
                    			ligand = pdbp.get_structure('Ligand', ligand_list[i])
                		except IOError:
                    			print "Unable to read in ligand file"
					sys.exit()

                		#Find centroid of ligand
                		centroid = Vector(0,0,0)
                		counter = 0
                		for atom in ligand[0].get_list()[0].get_list()[0]:
                    			centroid = centroid + atom.get_vector()
                    			counter += 1
                		centroid = centroid/counter

                		#Recenter ligand at zero
                		for atom in ligand[0].get_list()[0].get_list()[0]:
                    			atom.set_coord(atom.get_vector()-centroid)

                		#Get rotation matrix from random random unit vector and angle
                		theta = uniform(0,2*pi)
                		ux = uniform(-1,1)
                		uy = uniform(-1,1)
                		uz = uniform(-1,1)
                		norm = pow(pow(ux,2)+pow(uy,2)+pow(uz,2),.5)
                		u = Vector(ux/norm,uy/norm,uz/norm)
                		rot_m = rotaxis(theta,u)

				#Rotate ligand
                		for atom in ligand[0].get_list()[0].get_list()[0]:
                    			atom.set_coord(atom.get_vector().left_multiply(rot_m))

                		#Add displacement to ligand
                		for atom in ligand[0].get_list()[0].get_list()[0]:
                    			atom.set_coord(atom.get_vector()+displace)

                		#Remove all waters that are within 4 angstroms of the ligands
                		for atom in Selection.unfold_entities(ligand, 'A'):
                    			a,b,c = atom.coord
                    			neighbors = ns.search(np.array([a,b,c]), 3.5, level = 'R')
                    			for residue in neighbors:
                        			membrane[0][residue.get_parent().id][residue.id].full_id=True
                        			bad_res[(residue.id,residue.resname)]=True

                		#Add ligand to pdb file
                		ligand[0].get_list()[0].id = 'l' #Change chain name to letter a-whatever
                		#Change residue number
                		temp_id = ligand[0].get_list()[0].get_list()[0].get_id()
                		ligand[0].get_list()[0].get_list()[0].id = (temp_id[0],x*num_ligands+y+1,temp_id[2])
                		ligand[0].get_list()[0].id
                		try:
                    			membrane[0].add(ligand[0].get_list()[0])
                		except:
                    			membrane[0]['l'].add(ligand[0].get_list()[0].get_list()[0])
			
			#CheckIons(bad_res,membrane)	

			#write out membrane with ligands
			io.set_structure(membrane)
			out_file = membrane_file + '_' + ligand_list[i]+'.pdb'
			io.save(out_file, ResSelect())
			
			#Adding TER card to the files
			AddTer("TER"+out_file,out_file)
			#CheckLipid("Final"+out_file,"TER"+out_file)			
			call(['mv',"TER"+out_file,out_file])

def main():
	parser = argparse.ArgumentParser(description='Program to create a membrane file of given length and insert a protein if the file is provided')
	parser.add_argument('-s','--size',dest='size',help='Enter the membrane size in A (int)',required=True)
	parser.add_argument('-p','--protein',dest='sample',help='Enter the location of the pdb file that needs to be inserted',required=False,default='Na')
	parser.add_argument('-r','--reference',dest='ref',help='Enter the location of the reference pdb file that can be used as a reference to align the protein',required=False,default='Na')
	parser.add_argument('-l','--ligand',nargs='+',dest='ligand',help='Enter the location of the ligand file to be inserted',default='Na')
	parser.add_argument('-n','--number',dest='lig_num',help='Enter the number of ligands per x to insert into the membrane. Use this mode if the ligand coordinates are not known.',default=1)
	parser.add_argument('-o','--output',dest='output',help='Enter the output file name',required=False,default='output.pdb')
	args = parser.parse_args()

	if(os.path.exists(args.sample) and (args.ref == 'Na' or not os.path.exists(args.ref))):
		print "The reference file is missing or unable to be opened"
		parser.print_help()
		sys.exit()

	if(not os.path.exists(args.sample) and os.path.exists(args.ref)):
		print "The pdb file is missing or unable to be opened"
		parser.print_help()
		sys.exit()

	if(args.ligand != 'Na' and int(args.lig_num) > 0 and (args.ref != 'Na' or args.sample != 'Na')):
		print "When picking multiple ligands, proteins cannot be inserted by this program."
		print "Please try again"
		sys.exit()


	membrane_file = 'POPC_'+str(args.size)+'A'

	#Function to create membrane
	CreateMembrane(args.size,membrane_file)

	#Function to insert only n ligands in random locations inside the membrane
	if(args.ligand != 'Na' and int(args.lig_num) > 0):
		#Function to add ligands
		AddLigands(args.ligand,int(args.lig_num),membrane_file,int(args.size))	

		print "The membrane files with ligands have been saved in the current directory"
		sys.exit()
	
	#Function to insert protein into file if protein and reference pdb are given
	if(args.sample != 'Na' and args.ref != 'Na'):
		#Function to align new protein to reference protein
		AlignProtein(args.ref,args.sample)

		#Function to reformat protein pdb into amber readable format
		ReformatProtein("chimera_aligned_"+args.sample)
		
		#Function to merge the protein and membrane file
		MergePdb(membrane_file,args.sample)

		#Function to insert TER cards appropriately
		AddTer("temp_output","mem_pro.pdb")

		#Function to make sure the heads and tails are all attached
		CheckLipid(args.output,temp_output)
		
		print "The program ran successfully"
		print "The output is located in: " + args.output

        else:
                print "Program did not detect protein to be inserted!"
		print "Your file has been created!"
		sys.exit()

P
