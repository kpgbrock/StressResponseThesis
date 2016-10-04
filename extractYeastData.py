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


InName = 'yeastGFPfull.txt'
OutName = 'yeastGFPfiltered.txt'
    
# read in input sequences and process one at a time, then write the output to the output file
infile = open(InName,'r')
outfile = open(OutName,'w')

isFirst = True

# Go through each EC ID
for line in infile:

	if isFirst:
		outfile.write(line)
		isFirst = False
	else:	
		# Take out newlines for formatting purposes
		line = line.strip('\n')
		components = re.split('\t',line)
		
		# Check to make sure valid compound id
		check = re.search('cytoplasm',components[8])
		if check != None:
			if is_number(components[6]):
				outfile.write(line+'\n')
   	 
    
# close files at the end
infile.close()
outfile.close()
