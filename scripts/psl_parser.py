import logging

class Psl_data:
	"""
	PSL file format
	http://www.ensembl.org/info/website/upload/psl.html
	"""
	
	## Initialise with default values
	def __init__(self, matches= -1, misMatches = -1, repMatches = -1, nCount = -1, qNumInsert = -1, qBaseInsert = -1, tNumInsert = -1, tBaseInsert = -1, strand = "none", qName = None, qSize = -1, qStart = -1, qEnd = -1, tName = None, tSize = -1, tStart = -1, tEnd = -1, blockCount = 0 , blockSizes = [], qStarts = [], tStarts = []):
		self.matches = matches
		self.misMatches = misMatches
		self.repMatches = repMatches
		self.nCount = nCount
		self.qNumInsert = qNumInsert
		self.qBaseInsert = qBaseInsert
		self.tNumInsert = tNumInsert
		self.tBaseInsert = tBaseInsert
		self.strand = strand
		self.qName = qName
		self.qSize = qSize
		self.qStart = qStart
		self.qEnd = qEnd
		self.tName = tName
		self.tSize = tSize
		self.tStart = tStart
		self.tEnd = tEnd
		self.blockCount = blockCount
		self.blockSizes = blockSizes
		self.qStarts = qStarts
		self.tStarts = tStarts

## Initialise with full list(single psl line) of string represented values
	def __init__(self, parameterlist):
		self.matches = int(parameterlist[0])
		self.misMatches = int(parameterlist[1])
		self.repMatches = int(parameterlist[2])
		self.nCount = int(parameterlist[3])
		self.qNumInsert = int(parameterlist[4])
		self.qBaseInsert = int(parameterlist[5])
		self.tNumInsert = int(parameterlist[6])
		self.tBaseInsert = int(parameterlist[7])
		self.strand = parameterlist[8]
		self.qName = parameterlist[9]
		self.qSize = int(parameterlist[10])
		self.qStart = int(parameterlist[11])
		self.qEnd = int(parameterlist[12])
		self.tName = parameterlist[13]
		self.tSize = int(parameterlist[14])
		self.tStart = int(parameterlist[15])
		self.tEnd = int(parameterlist[16])
		self.blockCount = int(parameterlist[17])
		## using generators for potential big number of blocks
		self.blockSizes = (int(n) for n in parameterlist[18].split(",")[:-1])
		self.qStarts = (int(n) for n in parameterlist[19].split(",")[:-1])
		self.tStarts = (int(n) for n in parameterlist[20].split(",")[:-1])


## string representation of the class
	def __str__(self):
		return "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17} {18} {19} {20}".format(self.matches, self.misMatches, self.repMatches, self.nCount, self.qNumInsert, self.qBaseInsert, self.tNumInsert, self.tBaseInsert, self.strand, self.qName, self.qSize, self.qStart, self.qEnd, self.tName, self.tSize, self.tStart, self.tEnd, self.blockCount, list(self.blockSizes), list(self.qStarts), list(self.tStarts))		
		

## parse psl-file, returns a list
def parse_psl(filename):
	data = []
	with open(filename, "r") as input:
			## read line by line and parse non comments to Psl-data
			for line in input.readlines():
				if ("##" in line):
					continue
				## append to entry list, line needs to be split for the constructor
				data.append(Psl_data(line.split()))
	return data
		
## generator function for Psl-entries from the specified file
def gen_psl_content(filename):
	with open(filename, "r") as input:
		for line in input.readlines():
			## ignore comments
			if ("##" in line):
				continue
			## generator expression
			yield Psl_data(line.split())


if __name__ == "__main__":
	import argparse
	import sys
	(_, filename) = sys.argv
	## some testing expressions
	"""test = parse_psl(filename)
	for entry in test:
		print(entry)"""
	gen_test = gen_psl_content(filename)
	"""for entry in gen_test:
		print(entry)"""
	for _ in range(5):
		print(gen_test.next())
