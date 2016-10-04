# CpdIdLookup.py
# Link compound ID's with compound names through Kegg screenscraper
# Used to help simplify metabolite naming system

import sys

import shutil
import os
import time
import datetime
import math
import urllib
from array import array
import re

def searchForID(searchterm,line,RelAddr):

	pmidAddr = RelAddr + "".join(line[0:len(line)])
	print pmidAddr
   	address = urllib.urlopen(pmidAddr)
	
    # Find and store names
   	l = address.readlines();
   	res = ""
   	for i in range(len(l)):
		l[i] = l[i].strip('\n')
		res = res + "".join(l[i])
	print res
   	m = re.search(searchterm,res)
	return m


# Beginning of url to access KEGG compound entry
RelAddr = 'http://www.ncbi.nlm.nih.gov/nuccore/'

InName = 'C:/Users/Kelly/Desktop/purkinjeGeneBank.txt'
OutName = 'C:/Users/Kelly/Desktop/purkinjeList.txt'
    
# read in input sequences and process one at a time, then write the output to the output file
infile = open(InName,'r')
outfile = open(OutName,'w')

# Go through each EC ID
for line in infile:

    # Take out newlines for formatting purposes
    line = line.strip('\n')
    mRNA = line;
	
    # Check to see if in RefSeq format
    check = re.search('NM.*',line)
	
	# If not, attempt to get mRNA RefSeq
    if check == None:
		mRNA = searchForID('reference mRNA sequence \((?P<name>.*)\)',line,RelAddr)
		
    if mRNA != None:
		print mRNA
		protein = searchForID('\((?P<name>NP_.*)\)',mRNA,RelAddr)
		outfile.write(protein + '\n')
   	 
    
# close files at the end
infile.close()
outfile.close()
address.close()
