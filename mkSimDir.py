#!/usr/bin/env python
# mkSimDir creates a whole new simulation directory or updates
# an old one (needed on Windows). One needs one simulation
# directory per simultaneously running simulation. Since windows
# does not support symbolic links, mkSimDir copies the files
# necessary to simulate to the newly created directory. Linux,
# on the other hand, makes symbolic links, which makes it
# easier to update the program when there is a change. To update
# the program on windows, one needs to make mkSimDir on all
# simulation directories one has.
#
# Arg 1: Name (including path) of the new simulation directory.
#
# Author: Fredrik Edin, 2004
# Address: freedin@nada.kth.se
#
# Remember to update paths marked with # CHANGE # in accordance
# to the organization of your file system.

import sys, shutil, os, glob

# CHANGE (Do not write the ~ symbol): Location of STANDARDFILER catalog #
# Remember to use two backslashes when writing windows pathdelimiters
home = "C:\\Documents and Settings\\freedin.admin\\Desktop\\Neuron\\STANDARDFILER\\"
neuronhome = 'C:\\nrn56\\'
dest = sys.argv[1]

# If the directory already exists, just copy the
# files into the directory and do not make a new one.
if not os.path.exists( dest ):
	os.mkdir( dest )


# Copy all .mod files. They don't work with symbolic links
modfiles = glob.glob( home+'*.mod' )
for name in modfiles:
	shutil.copy( name, dest )

# Compile .mod files
os.chdir( dest )
if sys.platform == 'win32':
	os.system( 'sh ' + neuronhome + 'lib\\mknrndll.sh ' + neuronhome )
else:
	os.system( 'nrnivmodl' )

# Link the rest of the files (UNIX/Linux) or copy them
files = ['MultiModuleWMNet.hoc','Net.hoc','Results.hoc','View.hoc','ECellIAF.hoc','ECell2NMDA.hoc','findSims.m','getTree.m','getTreeSpec.m','ICellIAF.hoc','ICell2NMDA.hoc','MyRandom.hoc','Simulator.hoc','showSimulation.m','showSERIES.m','SimLoop.py','totCurr.txt','RandomSeed.py','doneSoFar.txt','rC1v.txt','rC2v.txt','rNMDA.txt','rNMDA2.txt','rAMPA.txt','PQ.hoc']
for name in files:
	if sys.platform == 'win32':
		shutil.copy( home+name, name )
	else:
		os.symlink( home+name, name )

# Create a DATA directory for storage of simulation resutls
os.mkdir( 'DATA' )
os.chdir( '..' )

