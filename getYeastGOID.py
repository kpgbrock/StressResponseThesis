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


InName = 'yeastGOFullData.txt'
OutNameF = 'yeastGOFunction.txt'
OutNameP = 'yeastGOProcess.txt'
    
# read in input sequences and process one at a time, then write the output to the output file
infile = open(InName,'r')
outfileF = open(OutNameF,'w')
outfileP = open(OutNameP,'w')

# isFirst = True

# Go through each EC ID

yid = ''
ydictF = {}
ydictP = {}

isCyto = False

for line in infile:
	
	line = line.strip('\n');
	components = re.split('\t',line)
	
	# Shortcut for important items
	yidTemp = components[0]
	yidCFP = components[3]
	yidCyto = components[4]
	yidGO = components[5]
	
	# Check to see if we've started on a new protein
	if yidTemp != yid:
		yid = yidTemp
		isCyto = False
	
	# Check to see if it's a cytosolic protein
	if yidCyto == 'cytoplasm':
		isCyto = True
	
	if isCyto and (yidCFP == 'F'):
		x = re.search('GO:(\d*)',yidGO)
		if x is None:
			print line
			print yidGO
			continue
		goid = x.group(1)
		
		if goid in ydictF:
			ydictF[goid].append(yid)
		else:
			ydictF[goid] = [yid]
	
	if isCyto and (yidCFP == 'P'):
		x = re.search('GO:(\d*)',yidGO)
		if x is None:
			print 'hello'
			print yidGO
			continue
		goid = x.group(1)
		if goid in ydictP:
			ydictP[goid].append(yid)
		else:
			ydictP[goid] = [yid]	
	
	# if isFirst:
		# outfile.write(line)
		# isFirst = False
	# else:	
		#Take out newlines for formatting purposes
		# line = line.strip('\n')
		# components = re.split('\t',line)
		
		#Check to make sure valid compound id
		# check = re.search('cytoplasm',components[8])
		# if check != None:
			# if is_number(components[6]):
				# outfile.write(line+'\n')
   	 
for k in ydictF.keys():
	outfileF.write(k+'\t')
	for x in ydictF[k]:
		outfileF.write(','+'\''+x+'\''+'\t')
	outfileF.write('\n')

for k in ydictP.keys():
	outfileP.write(k+'\t')
	for x in ydictP[k]:
		outfileP.write(','+'\''+x+'\''+'\t')
	outfileP.write('\n')
	
# close files at the end
infile.close()
outfileF.close()
outfileP.close()
