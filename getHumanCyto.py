# We have a data file containing expression quantities and 
# localization parameters for different yeast genes.  We want to
# only take lines that have cytoplasm in the location and that have
# an expression quantity measured.

import sys

import shutil
import os
import time
import datetime
import math
import urllib
from array import array
import re

def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False


inNameSL = 'HPRD_FLAT_FILES_041310/humanProteomeSubcellularLocalization.txt'
inNameSQ = 'HPRD_FLAT_FILES_041310/humanProteomeSequence.txt'

OutName = 'humanCyto.fasta'
    
# read in input sequences and process one at a time, then write the output to the output file
infileSL = open(inNameSL,'r')
infileSQ = open(inNameSQ,'r')
outfile = open(OutName,'w')


yid = ''
ydict = {}
ensList = []

goCyto = '(GO:0005737)'
pid = '00000'

pidSub = '00000'
pidEns = ''

matlabBufSize = 4095


for line in infileSQ:
	
	line = line.strip('\n')
	
	
	# See if header or sequence
	if line[0] == '>':
		components = re.split('\|',line)
		
		pidSub=components[1]
		pidEns=components[2]
		
	else:
		
		if pidEns in ydict:
			print pidEns
		ydict[pidEns]=line
		
	


for line in infileSL:
	
	line = line.strip('\n');
	components = re.split('\t',line)
	
	# Shortcut for important items
	pidNum = components[0]
	pidSubNum = components[1]
	pidEns2 = components[2]
	pidLoc = components[9]
	
	# Check to see if it's a cytosolic protein
	searchgo = re.search(goCyto,pidLoc)
	if searchgo != None:
		
		#If no subnumber, assign to be #_1
		# if pidSubNum == '-':
			# pidSubNum = pidNum + '_1'
			
		# If it hasnt already been added, add it to file, and make sure it wont overflow matlab buffer size
		if pidEns2 not in ensList:
			
			if len(ydict[pidEns2]) <= matlabBufSize:
				
				outfile.write('>'+pidEns2+'\n')
				outfile.write(ydict[pidEns2])
				outfile.write('\n')
				ensList.append(pidEns2)
				
	
	
# close files at the end
infileSL.close()
infileSQ.close()
outfile.close()
