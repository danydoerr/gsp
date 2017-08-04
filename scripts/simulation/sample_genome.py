#!/usr/bin/env python

from sys import stdout,stderr,exit,argv
from optparse import OptionParser
from os.path import basename
from random import gammavariate as gamma, choice, randint
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import logging

CODONS = ['GCT', 'GCC', 'GCA', 'GCG', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', \
        'AGG', 'AAT', 'AAC', 'GAT', 'GAC', 'TGT', 'TGC', 'CAA', 'CAG', 'GAA', \
        'GAG', 'GGT', 'GGC', 'GGA', 'GGG', 'CAT', 'CAC', 'ATT', 'ATC', 'ATA', \
        'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG', 'AAA', 'AAG', 'ATG', 'TTT', \
        'TTC', 'CCT', 'CCC', 'CCA', 'CCG', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', \
        'AGC', 'ACT', 'ACC', 'ACA', 'ACG', 'TGG', 'TAT', 'TAC', 'GTT', 'GTC', \
        'GTA', 'GTG']
CODON_START = ['ATG']
CODON_STOP = [ 'TAA', 'TGA', 'TAG']

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

DEFAULT_MIN_LENGTH = 100
DEFAULT_CHROMOSOME_NO = 6

if __name__ == '__main__':

    usage = 'usage: %prog [options] <#CONSERVED SEGMENTS>'
    parser = OptionParser(usage=usage)

    parser.add_option('-c', '--chromosomes', dest='chrs',type=int,
            default=DEFAULT_CHROMOSOME_NO, metavar='INT',
            help='Number of chromosomes. [default=%default]')

    parser.add_option('-l', '--minlength', dest='min', type=int,
            default=DEFAULT_MIN_LENGTH, metavar='INT',
            help='Miminimum length of a conserved segment. [default=%default]')

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        exit(1)

    
    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.ERROR)
    ch.setFormatter(logging.Formatter('!! %(message)s'))
    
    cf = logging.FileHandler('%s.log' %(basename(argv[0]).rsplit('.', 1)[0]), mode='w')
    cf.setLevel(logging.INFO)
    cf.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t++ %(message)s'))

    LOG.addHandler(cf)
    LOG.addHandler(ch)

    #
    # main 
    #

    random_dna   = lambda x: ''.join(choice('ACGT') for _ in xrange(x))
    random_codon = lambda x: ''.join(choice(CODONS) for _ in xrange(x))

    n = int(args[0])
    marker_len = [int(options.min + gamma(0.7, 600)) for _ in xrange(n)]
    marker_len = [max(m/3 * 3, 6) for m in marker_len]
    inter_marker_len = map(lambda x: x >= options.min and x or 0,
            [int(gamma(0.7, 5000)) for _ in xrange(n+options.chrs)])
    chromosome_len = [gamma(7.5, 1) for _ in xrange(options.chrs)]
    s = sum(chromosome_len)
    chromosome_len = [int(round(n*x/s)) for x in chromosome_len]
    chromosome_len[randint(0, options.chrs-1)] += n-sum(chromosome_len)

    recs = list()
    chr_id = i = m = p = 0
    while m < n:
        if chromosome_len[chr_id] and inter_marker_len[i]:
            recs.append(SeqRecord(Seq(random_dna(inter_marker_len[i]),
                generic_dna), id='intermarker_%s|%s:%s|chromosome|%s|strand|+' %(i+1,
                    p+1, p+inter_marker_len[i], chr_id), description=''))
            p += inter_marker_len[i] 
        i += 1
        
        for _ in xrange(chromosome_len[chr_id]):
            recs.append(SeqRecord(Seq(choice(CODON_START) + \
                    random_codon(marker_len[m]/3) + choice(CODON_STOP),
                generic_dna), id='marker_%s|%s:%s|chromosome|%s|strand|+' %(m+1,
                        p+1, p+marker_len[m], chr_id), description=''))
            p += marker_len[m]
            m += 1
            if inter_marker_len[i]:
                recs.append(SeqRecord(Seq(random_dna(inter_marker_len[i]),
                    generic_dna), id='intermarker_%s|%s:%s|chromosome|%s|strand|+' %(i+1,
                        p+1, p+inter_marker_len[i], chr_id), description=''))
                p += inter_marker_len[i] 
            i += 1
        chr_id += 1

    SeqIO.write(recs, stdout, 'fasta')

