# Kelly Brock
# England Rotation
# Extracts cytoplasm-based 

import sys
import re

def KD(seq):
	kd = 0;
	for i in range(len(seq)):
		AA = seq[i]
		if AA == 'I':
			kd = kd + 4.5
		elif AA == 'X':
			kd = kd - 0.4
		elif AA == 'V':
			kd = kd + 4.2
		elif AA == 'L':
			kd = kd + 3.8
		elif AA == 'F':
			kd = kd + 2.8
		elif AA == 'C':
			kd = kd + 2.5
		elif AA == 'M':
			kd = kd + 1.9
		elif AA == 'A':
			kd = kd + 1.8
		elif AA == 'G':
			kd = kd - 0.4
		elif AA == 'T':
			kd = kd - 0.7
		elif AA == 'W':
			kd = kd - 0.9 
		elif AA == 'S':
			kd = kd - 0.8
		elif AA == 'Y':
			kd = kd - 1.3
		elif AA == 'P':
			kd = kd - 1.6
		elif AA == 'H':
			kd = kd - 3.2
		elif AA == 'E':
			kd = kd - 3.5
		elif AA == 'Q':
			kd = kd - 3.5
		elif AA == 'D':
			kd = kd - 3.5
		elif AA == 'N':
			kd = kd - 3.5
		elif AA == 'K':
			kd = kd - 3.9
		elif AA == 'R':
			kd = kd - 4.5
	return kd*1.0/len(seq)
	
def readFasta(filename):
    """reads in a FASTA sequence"""

    stream = open(filename)
    seqdict = {}
    seq = []
    yorf = ''
	
    for line in stream:
		if line.startswith(">"):
			
			if seq != []:
				seqdict[yorf] = "".join(seq)
			seq = []
			yorf = re.findall('Y\w+',line)
			if yorf == []:
				yorf = 'N/A'
			else:
				yorf = yorf[0]
			continue
			
		seq.append(line.rstrip())
		
    
    seqdict[yorf] = "".join(seq)
    
    stream.close()
    return seqdict
	
def readYGFP(seqdict,infile,outfile):
	seqmin = 150
	seqmax = 200
	kdseqmin = -1
	kdseqmax = 0
	stream = open(infile,'r')
	streamout = open(outfile,'w')
	
	locsummary = 8
	orf = 1
	for line in stream:
		row = re.split('\t',line)
		if len(row) <= locsummary:
			continue
		
		if re.findall('cytoplasm',row[locsummary]) != []:
			if row[orf] in seqdict:
				#print len(seqdict[row[orf]])
				#print KD(seqdict[row[orf]])
				if (len(seqdict[row[orf]]) >= seqmin) & (len(seqdict[row[orf]]) <= seqmax) & (KD(seqdict[row[orf]]) > kdseqmin) & (KD(seqdict[row[orf]]) < kdseqmax):
					streamout.write('>' + row[orf] + '\n')
					streamout.write(seqdict[row[orf]] + '\n')
	stream.close()
	streamout.close()
	
def main():
	fastafile = 'orf_trans.txt'
	locfile = 'yeastgfp.txt'
	outfile = 'yeastCytoplasmFiltered2.txt'
	seqdict = readFasta(fastafile)
	readYGFP(seqdict,locfile,outfile)

main()