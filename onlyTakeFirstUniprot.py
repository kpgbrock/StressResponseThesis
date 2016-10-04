
import sys

import shutil
import os
import time
import datetime
import math
import urllib
from array import array
import re

inName = 'C:/Users/Kelly/Desktop/testNP.txt'
outName = 'C:/Users/Kelly/Desktop/test.txt'
infile = open(inName,'r')
outfile = open(outName,'w')

# Go through each EC ID
for line in infile:

    # Take out newlines for formatting purposes
	line = line.strip('\n')
	if line == '\n':
		outfile.write('\n')
	elif line == 'Data not found':
		outfile.write('\n')
	else:
		line = re.split('\|',line)
		line = line[0]
		outfile.write(line + '\n')
	
infile.close()
outfile.close()