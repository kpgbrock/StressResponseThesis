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

# Beginning of url to access KEGG compound entry
RelAddr = 'http://www.genome.jp/dbget-bin/www_bget?cpd:'

InName = 'C:/Users/Kelly/Desktop/LabWork/rawCPD.txt'
OutName = 'C:/Users/Kelly/Desktop/LabWork/metNames.txt'
    
# read in input sequences and process one at a time, then write the output to the output file
infile = open(InName,'r')
outfile = open(OutName,'w')

# Go through each EC ID
for line in infile:

    # Take out newlines for formatting purposes
    line = line.strip('\n')
    
    # Check to make sure valid compound id
    check = re.search('C.*',line)
    if check == None:
   	 outfile.write(line + '\n')
    else:
   	 # Open KEGG compound entry
   	 pmidAddr = RelAddr + "".join(line[1:len(line)-1])
   	 address = urllib.urlopen(pmidAddr)
    
    
   	 # Find and store names
   	 l = address.readlines();
   	 res = ""
   	 for i in range(len(l)):
   		 l[i] = l[i].strip('\n')
   		 res = res + "".join(l[i])
   	 m = re.search('Name.*?: solid">(?P<name>.*?)(;|(</td></tr>))',res)
   	 outfile.write(m.group('name') + '\n')
   	 
    
# close files at the end
infile.close()
outfile.close()
address.close()
