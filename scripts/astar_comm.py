from operator import itemgetter
from collections import deque
import networkx as nx
import pdb
import time
 
#Input datatype is .psl file with variable number of alignments 
input = open("seq.fa.psl")
#input = open("drosPSLtest.psl")
content = input.readlines()
input.close()

array = [] 
# minimal required matches in alignment. maybe dominik can help here
minMatches = 100
# array of names and position of the aligned sequences in the alignment
sequences = []
#array of raw breakpoints
breakpoints = []
# start of the sequences of one alignment used for extracting the blocks
seq1Start = 0
seq2Start = 0
# boolean which is true if array already contains sequence in question
containsSeq = False
# array of breakpoints after mapping
foundBr = []		
#set of breakpoints used while mapping	
activeBr = set()
# graph with weighted bpś
comGraph = []
#graph with bpś that have a weight over the newSeqCost
newGraph = []
# minimal Segment length
minSegLen = 400
# cost of opening a new segment
newSeqCost = 40
# array of the edges
edges = []
#list of list of edges which contain a problem for the algorithm
edgeList = []

#function returns the start/relative position of the given sequence in the composed sequence
def getSeqPos(sequences, index):
	if sequences:
		return int(sequences[index-1][2]) + int(sequences[index-1][3])
	else:
		return 0	

#appending alignments to the array which have a higher or equal number of matches than 100	
for line in content:
	splitLine = line.split()
	if ("##" not in splitLine[0] and int(splitLine[0])>=minMatches):
		array.append(splitLine)
				
#adding first sequence to sequences	
if not sequences:
	sequences.append([array[0][9],0,array[0][10],0])	
		
#adding sequencename, position, lenght and relative sequence position	
for x in range (0, len(array)):
	for y in range(0, len(sequences)):
		if sequences[y][0] == array[x][9]:
			Seq1Start = sequences[y][3]
			containsSeq = True
			break
	if not containsSeq:
		sequences.append([array[x][9],len(sequences),array[x][10],getSeqPos(sequences, len(sequences))])
		seq1Start = getSeqPos(sequences, len(sequences))
	containsSeq = False
	for y in range(0, len(sequences)):
		if sequences[y][0] == array[x][13]:
			seq2Start = sequences[y][3]
			containsSeq = True
			break
	if not containsSeq:
		sequences.append([array[x][13],len(sequences),array[x][14],getSeqPos(sequences, len(sequences))])
		seq2Start = getSeqPos(sequences, len(sequences))
		
	#splitting the array of blocks as i treat every block as its own alignment	
	blocks1 = array[x][19][:-1].split(",")
	for s in range(0,len(blocks1)):
		blocks1[s] = int(blocks1[s])+seq1Start 
	s=0
	blocks2 = array[x][20][:-1].split(",")
	for s in range(0,len(blocks2)):
		blocks2[s] = int(blocks2[s])+seq2Start 
	s=0
	blocksLen = array[x][18][:-1].split(",")

	#list of breakpoints which contains: brStart, brEnd, 2BrStart, 2BrEnd
	if(blocks1[0]>blocks2[0]):
		for c in range(0, len(blocks2)):
			if int(blocksLen[c]) > minSegLen:
				breakpoints.append([blocks2[c], blocks2[c]+int(blocksLen[c]), blocks1[c], blocks1[c]+int(blocksLen[c])])
			#breakpoints.append([int(array[x][12])+seq1Start, int(array[x][11])+seq1Start, int(array[x][16])+seq2Start, int(array[x][15])+seq2Start, blocks2, blocks1])
	else:
		for c in range(0, len(blocks1)):
			# block have to be above minSeqLen to be considered
			if int(blocksLen[c]) > minSegLen:
				breakpoints.append([blocks1[c], blocks1[c]+int(blocksLen[c]), blocks2[c], blocks2[c]+int(blocksLen[c])])	
		#breakpoints.append([int(array[x][11])+seq1Start, int(array[x][12])+seq1Start, int(array[x][15])+seq2Start, int(array[x][16])+seq2Start, blocks1, blocks2])

breakpoints.sort(key=lambda x: x[0])

#mapping breakpoints
x=0
for x in range (len(breakpoints)):
	activeBr.add(breakpoints[x][0])
	activeBr.add(breakpoints[x][1])
	activeBr.add(breakpoints[x][2])
	activeBr.add(breakpoints[x][3])
	foundBr.append(breakpoints[x][0])
	foundBr.append(breakpoints[x][1])
	foundBr.append(breakpoints[x][2])
	foundBr.append(breakpoints[x][3])
x=0	
print(len(activeBr))

#actual mapping of the breakpoints which takes super long. dont know how to make this faster
#bei minSeglen=200 383671 gemappte Brs (1h) ausgehend von 6651 zu anfang
while activeBr:
	br = activeBr.pop()	
	for x in range(len(breakpoints)):
			if br > breakpoints[x][0] and br < breakpoints[x][1]:
				if breakpoints[x][2]+br-breakpoints[x][0] not in foundBr:
					activeBr.add(breakpoints[x][2]+br-breakpoints[x][0])
				foundBr.append(breakpoints[x][2]+br-breakpoints[x][0])

			
foundBr.sort()	




comGraph.append([foundBr[0],1,False])
x=0
#giving weight to the bpś dependent on how many other bpś are close. 
for x in range(len(foundBr)):
	if 	foundBr[x] == foundBr[x-1]:
		comGraph[-1][1] += 20	
	else:
		if foundBr[x-1]+20>foundBr[x]:
			comGraph[-1][1] += foundBr[x-1]+20-foundBr[x]
		comGraph.append([foundBr[x],1,False])	
x=0		
#adding brś to the new Graph if their weight is bigger than the cost for opening a new segment
for x in range (len(comGraph)):
	if comGraph[x][1]>=newSeqCost:
		newGraph.append(comGraph[x])
				
#return the cost of an edge from x to y
def getEdgeCost(x,y):
	cost = newSeqCost
	for i in range (x+1,y):
			cost+=newGraph[i][1]	
	return cost

#returning clostest node before x that could have an edge to x
def getClosest(x):
	y=x-1
	while (1):
		if newGraph[x][0]-newGraph[y][0] > minSegLen:
			break
		y=y-1
	return y	

#returns next node
def getNext(x):
	y=0
	if x != len(newGraph):
		y= x+1
	else:
		y=x
	return y	
#returns node before x no matter how far 
def getClose(x):
	y=0
	if x != 0:
		y=x-1
	return y		
		
	
x=0
#creating the edges and also seperating the edges in lists which are not trivial, so the algorithm can be applied independently
for x in range (0, len(newGraph)):
	y=x+1
	#condition for a new problem
	if x != len(newGraph)-1 and newGraph[x][0]>newGraph[0][0]+minSegLen and int(newGraph[getNext(x)][0]) - int(newGraph[getClose(x)][0]) > minSegLen:
		edges.append([newGraph[getClosest(x)][0], newGraph[x][0], getEdgeCost(getClosest(x),x)])
		edgeList.append(edges)
		edges = []
	#creating new edge	
	while (y < len(newGraph) and newGraph[y][0]<newGraph[x][0]+3*minSegLen):	
		if newGraph[y][0]>=newGraph[x][0]+minSegLen:
			newGraph[y][2] = True
			edges.append([newGraph[x][0],newGraph[y][0],getEdgeCost(x,y)])
		y+=1	

			
print("edges done")		

#approximates the distance to the target		
def heuristic (v1, v2):
	x = (int)(v2-v1)/(3*minSegLen)
	return x*newSeqCost

#creation on the graph and start of time measurement	
diGraph = nx.DiGraph()
start = time.time()
o=0
for edges in edgeList:
	#if the range of edges is shorter than the segment lenght, we merge the array to one edge
	if  edges [-1][1] - edges[0][0] < minSegLen:
		mergedEdge = [(int)((edges [-1][1] + edges[0][0])/2), (int)((edges [-1][1] + edges[0][0])/2), 0]
		#this is just to count how many non trivial problem we have
		if len(edges) >=3:
			o+=1
		edges = []
		edges.append(mergedEdge)
	diGraph.add_weighted_edges_from(edges)
	segmentation = nx.astar_path(diGraph, edges[0][0], edges[-1][1],heuristic)
	segCost = nx.astar_path_length(diGraph, edges[0][0], edges[-1][1])
	diGraph.clear()
	print(segmentation)
	print(segCost)
	
end = time.time()
print(end-start)	
#print(sequences)	
#print(mappedBr)
#print(graph)
#print(comGraph)
#print(newGraph)
#print(edges)
#print(breakpoints)		
#print(foundBr)
print(o)
