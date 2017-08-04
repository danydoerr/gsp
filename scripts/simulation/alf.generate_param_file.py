#!/usr/bin/env python

from sys import stdout,stderr,exit,argv
from optparse import OptionParser
from os.path import basename, join
import subprocess
import logging

DEFAULT_PAM = 10
DEFAULT_NSPECIES = 10

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


if __name__ == '__main__':

    usage = 'usage: %prog [options] <SIMULATED ROOT GENOME>'
    parser = OptionParser(usage=usage)

    parser.add_option('-p', '--pam', dest='pam',type=int, default=DEFAULT_PAM,
            metavar='INT', help='evolutionary distance [default=%default]')

    parser.add_option('-n', '--nspecies', dest='nspecies',type=int,
            default=DEFAULT_NSPECIES, metavar='INT',
            help='number of species in the sampled tree [default=%default]')

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

    seqTypes = list()

    for line in open(args[0]):
        if line.startswith('>'):
            if line.startswith('>marker'):
                seqTypes.append(1)
            else:
                seqTypes.append(2)
   
    dbName = '%s.db' %args[0].rsplit('.', 1)[0] 
    p = subprocess.Popen(['fasta2darwin', args[0], '-d', '-o',
        dbName], stdout=subprocess.PIPE)
    output, err = p.communicate()

    drw = """
mname := '%s_p%s_n%s':
unitIsPam := true:
realorganism := '%s'; # real genome db as first organism
treeType := 'BDTree':
birthRate := 0.01:
deathRate := 0.001:
mutRate := %s:        # PAM distance from origin to recent species (for random trees)
scaleTree := true: 
NSpecies := %s:

# parameters concerning the substitution models
substModels := [SubstitutionModel('CPAM'), SubstitutionModel('TN93', [0.8, 0.8, 0.8], [seq(0.25, 4)], true)];
indelModels := [IndelModel(0.000003, ZIPF, [1.821], 50), IndelModel(0.00003, ZIPF, [1.821], 50)];
rateVarModels := [RateVarModel(Gamma, 5, 0.01, 1), RateVarModel(Gamma, 3, 0.01, 5)];

seqTypes := [[1,1,1, 'marker'], [2,2,2, 'inter_marker']]:
seqTypeAssignments := [%s];

modelSwitchS := [[1.0, 0], [0, 1.0]]:
modelSwitchD := [[1.0, 0], [0, 1.0]]:

# parameters concerning gene duplication
geneDuplRate := 0.0001;
numberDupl := 5;
transDupl := 0.5;
fissionDupl := 0.1;
fusionDupl := 0.1;
P_pseudogene := 0.2;
ratefac_pseudogene := 0.9;
P_neofunc := 0.5;
ratefac_neofunc := 1.5;
life_neofunc := 10;
P_subfunc := 0.3;
ratefac_subfunc := 1.2;
life_subfunc := 10;

# parameters concerning gene loss
geneLossRate := 0.0001;
numberLoss := 5;

# parameters concerning genome rearrangement
invers := 0.004;
invSize := 200;
transloc := 0.002;
transSize := 100;
invtrans := 0.1;
numberFusion := 0;

# only output relevant data
simOutput := {'VP', 'Fasta'};
""" %(args[0].rsplit('.', 1)[0], options.pam, options.nspecies, dbName, options.pam, options.nspecies, ', '.join(map(str,
    seqTypes))) 

    print >> stdout, drw
