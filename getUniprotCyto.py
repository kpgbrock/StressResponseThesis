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


inNameSL = 'uniprotSlimeMold.txt'
OutName = 'slimeMoldCyto.fasta'
    
# read in input sequences and process one at a time, then write the output to the output file
infileSL = open(inNameSL,'r')
outfile = open(OutName,'w')



ydict = {}

matlabBufSize = 4095
acID = ''
protCyto = False
onSeq = False
seq = ''

for line in infileSL:
	
	line = line.strip('\n')
		
	# See what type of line it is
	lheader = line[0:2]
	
	# if lheader == '//':
		# newProt = True
	
	# Store accession ID if available
	if lheader == 'AC':
		components = re.split('\W+',line)
		acID = components[1]
		
	# If it's a line that would include subcellular localization, search for Cytoplasm
	elif lheader == 'CC':
		isCyto = re.search('Cytoplasm',line)
		if isCyto != None:
			protCyto = True
	# See if we're about to start feeding the sequence on the next line	
	elif lheader == 'SQ':
		onSeq = True
		seq = ''
	# See if we've come to the end of the protein
	elif lheader == '//':
		onSeq = False
		if protCyto:
			ydict[acID] = seq
		protCyto = False
		acID = ''
	elif onSeq:
		seq = ''.join([seq,''.join(re.split('\s+',line))])
		
	


for key in ydict:
	if len(ydict[key]) <= matlabBufSize:
		
		outfile.write('>'+key+'\n')
		outfile.write(ydict[key])
		outfile.write('\n')
	else:
		print key
				
	
	
# close files at the end
infileSL.close()
outfile.close()
