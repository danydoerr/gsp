#!/usr/bin/env python
#author: dany,mhnrchs

from sys import stdout, stderr, exit
from optparse import OptionParser
from os.path import basename
from copy import deepcopy
import re


import networkx as nx
from Bio import SeqIO

FASTA_HEADER_PAT = re.compile('^G\d+_(SE\d+), sequence type: (.*), locus: (-?\d+)$') # from dany's script
FILE_PAT = re.compile('SE(\d+)')


def readSequences(files):
    seqs = [] # seqs [i][j] = [seq_name, locus, length, positionInFile, class] for gene j in genome i
    for f in files: # iterate over input files (except last one, it's the homology matrix)
        atoms = []
        pos = 1
        for rec in SeqIO.parse(open(f), 'fasta'): # iterate over records
            desc = FASTA_HEADER_PAT.match(rec.description)
            seq_name = desc.group(1)
            locus = int(desc.group(3))
            atoms.append([seq_name, locus, len(rec), pos, None])
            pos += 1
        # sort atoms by locus
        atoms.sort(key = lambda x: abs(x[1]))
        seqs.append(atoms)
    return seqs


def readHomologies(data):
    f = open(data)
    line = f.readlines()[1]
    line = line[6:-2]
    # homologies[a][b][c] = list of genes in genome b homologous to gene c in genome a
    return eval(line)


def assignClasses(seqs, homologies):

    G = nx.Graph()
    for i in xrange(len(homologies)):
        for j in xrange(i, len(homologies)):
            for x in xrange(len(homologies[i][j])):
                for y in homologies[i][j][x]:
                    G.add_edge((i, x), (j, y-1))

    f = -1
    cls2pos = list()
    for C in nx.connected_components(G):
        f += 1
        cls2pos.append(list(C))
        for Gx, x in C:
            seqs[Gx][x][4] = f

    return cls2pos

def mergeAtoms(seqs, cls2pos):
    # create copy of seqs
    seqs = deepcopy(seqs)

    for c in xrange(len(cls2pos)):
        pos = cls2pos[c]
        right_neighbors = set()
        i = 0
        while len(right_neighbors) < 2 and i < len(pos):
            Gx, x = pos[i]
            if seqs[Gx][x][1] > 0:
                if x == len(seqs[Gx])-1:
                    right_neighbors.add(None)
                else:
                    right_neighbors.add((seqs[Gx][x+1][4], seqs[Gx][x+1][1] < 0
                        and -1 or 1))
            else:
                if x < 1:
                    right_neighbors.add(None)
                else:
                    right_neighbors.add((seqs[Gx][x-1][4], seqs[Gx][x-1][1] < 0
                        and -1 or 1))
            i += 1

        if len(right_neighbors) == 1 and next(iter(right_neighbors)) != None:
            # merge atoms into one
            print >> stderr, '++ merge neighboring families %s and %s' %(c,
                    next(iter(right_neighbors))[0])
            for Gx, x in pos:
                if seqs[Gx][x][1] < 0:
                    new_a = seqs[Gx][x-1]
                    new_a[2] += seqs[Gx][x][2]
                    if type(new_a[3]) == int:
                        new_a[3] = [new_a[3]]
                    if type(seqs[Gx][x][3]) == int:
                        seqs[Gx][x][3] = [seqs[Gx][x][3]]
                    new_a[3] = seqs[Gx][x][3] + new_a[3]
                    y = x
                else:
                    new_a = seqs[Gx][x]
                    new_a[2] += seqs[Gx][x+1][2]
                    if type(new_a[3]) == int:
                        new_a[3] = [new_a[3]]
                    if type(seqs[Gx][x+1][3]) == int:
                        seqs[Gx][x+1][3] = [seqs[Gx][x+1][3]]
                    new_a[3].extend(seqs[Gx][x+1][3])
                    y = x+1
                old_c = seqs[Gx][y][4]
                while y < len(seqs[Gx]) and seqs[Gx][y][4] == old_c:
                    seqs[Gx][y] = new_a
                    y += 1
    return seqs


def printSeqs(seqs, out):
    print >> out, '\t'.join(('GENOME', 'START', 'END', 'STRAND', 'FAMILY', 'POS_IN_FASTA'))
    for atoms in seqs:
        start = 0
        for i in xrange(len(atoms)):
            if not i or atoms[i-1] != atoms[i]:
                Gx, locus, length, pos, c = atoms[i]
                if type(pos) != int:
                    pos = ','.join(map(str, pos))
                print >> out, '\t'.join(map(str, (Gx, start+1, start+length,
                    locus < 0 and '-' or '+', c, pos)))
                start += length

if __name__ == '__main__':
    usage = 'usage: %prog <ALF SIMULATED DNA FILE 1> .. <ALF SIMULATED DNA ' + \
            'FILE N> <homology matrix>'
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()

    if len(args) < 2:
        parser.print_help()
        exit(1)

    ## parse sequence input
    files = sorted(args[:-1], key=lambda x:
            int(FILE_PAT.match(basename(x)).group(1)))
    seqs = readSequences(files)

    ## parse homologies
    homologies = readHomologies(args[-1])

    ## group atoms into classes according to homology assignment
    cls2pos = assignClasses(seqs, homologies)

    mergedSeqs = mergeAtoms(seqs, cls2pos)

    ## print result
    printSeqs(mergedSeqs, stdout)
    
