#!/usr/bin/env python
#author: mhnrchs
#partly inspired from dany's alf.genomes2dna.py

import sys
from optparse import OptionParser
import re
from Bio import SeqIO

FASTA_HEADER_PAT = re.compile('^(G\d+_SE\d+), sequence type: (.*), locus: (-?\d+)$') # from dany's script

def setBreakpoint(bps, cls, tailClass):
	if bps[cls] == -1: return # breakpoint was already evident, nothing to do
	if bps[cls] == 0:	# class not yet seen, set class of neighbor at tail
		bps[cls] = tailClass
	else:
		if bps[cls] != tailClass:
			bps[cls] = -1 # breakpoint evident, set it

if __name__ == '__main__':
	usage = 'usage: %prog <ALF SIMULATED DNA files> <homology matrix>\n\
Call with several files as arguments or append files into one.\n\
Be sure to input files in correct order though.'
	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()

	if len(args) < 1:
		parser.print_help()
		exit(1)


	## parse sequence input
	seq_name = ""
	seqs = [] # seqs [i][j] = [seq_name, locus, length, positionInFile, class] for gene j in genome i
	locus = 0
	for arg in args[:-1]: # iterate over input files except last one
		posInFile = 0
		atoms = []
		for rec in SeqIO.parse(open(arg), 'fasta'): # iterate over records
			posInFile += 1
			desc = FASTA_HEADER_PAT.match(rec.description)
			seq_name = re.search('_(.*), s', rec.description).group(1)
			locus = int(desc.group(3))
			atoms.append ([seq_name, locus,len(rec)-1,posInFile,0])
		seqs.append(atoms)
	# sort by locus for each sequence
	for i in range (len(seqs)):
		seqs[i] = sorted(seqs[i], key = lambda x: abs(x[1]))
	
	
	## parse homologies
	f = open(args[-1])
	line = f.readlines()[1]
	line = line[6:-2]
	homologies = eval(line)	# homologies[a][b][c] = list of genes in genome b homologous to gene c in genome a
	nrGenomes = len(homologies)
	currentClass = 1
	for i in range(nrGenomes):
		atoms = seqs[i]
		for (name, locus, length, posInFile, cls) in atoms:
			if (cls > 0): continue # was already classified
			seqs[i][abs(locus)-1][4] = currentClass
			for j in range(nrGenomes):
				for homLocus in homologies[i][j][abs(locus)-1]:
					seqs[j][homLocus-1][4] = currentClass
			currentClass += 1
			
	for i in range (100):
		for j in range (5):
			print seqs[j][i][4],
		print


	## find evident breakpoints
	breakpoints = [0] * (currentClass) # if breakpoints[0][10] == -1, there is a breakpoint before members of class 10
	# idea: for every segment, let c be its class. check the neighbor class. If two predeccors don't match, there is a breakpoint.
	seqNr = 0
	for atoms in seqs:
		for j in range(1,len(atoms)-1):
			locus = atoms[j][0]
			nextLocus = atoms[j+1][0]
			cls = atoms[j][4]
			neighborClass = sys.maxint 	# for segments at sequence start/end
			if locus > 0: # check right neighbor
				neighborClass = atoms[j-1][4]
				setBreakpoint(breakpoints, cls, neighborClass)
			else: # check left neighbor
				neighborClass = atoms[j+1][4]
				setBreakpoint(breakpoints, cls, neighborClass)
		seqNr += 1
		
	print breakpoints

	## build result by joining segments with no breakpoint between them
	result = []
	count = 0
	strand = ''
	seqNr = 0
	for atoms in seqs:
		curPos = 0
		endPos = -1 # this needs to be -1 for init due to the +1 in else below
		lastPos = atoms[0][3]
		lastName = atoms[0][0]
		for (name, locus, length, posInFile, cls) in atoms:
			lastAtom = atoms[abs(locus)-2] # -1 for the one before, another -1 because indexing starts at 0
			strand = '+' if lastAtom[1] >= 0 else '-'
			if breakpoints[cls] == -1: # make a new segment if a breakpoint is evident
				count += 1
				result.append((name, count, locus-1, strand, curPos, endPos, lastAtom[4]))
				curPos = endPos + 1
				endPos = curPos
			else: # go on to next locus
				endPos += 1
			endPos += length
			lastPos = posInFile
		# append last segment
		count += 1
		lastLocus = atoms[-1][1]
		lastClass = atoms[-1][4]
		result.append((lastName, count, lastLocus, '+' if lastLocus >= 0 else '-', curPos, endPos, lastClass))
		seqNr += 1
	
	## print result
	for (name,count,locus,strand,start,end, cls) in result:
		print (name + '\t' +		# sequence name
			str(count) + '\t' +		# atom number
			str(cls) + '\t' +		# class nr
			strand + '\t' +			# strand
			str(start) + '\t' + 	# atom start
			str(end))				# atom end
