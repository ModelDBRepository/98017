#!/usr/bin/env python 

# Comment:
# SimLoop.py program loops over simulations. The following steps are performed:  
# 1. One writes all parameters needed for a series of simulations in a file
#    whose file name contains the word SERIES. Let's call this file the SERIES
#    file. The SERIES file can be saved anywhere in the file system.
# 2. One enters a simulation catalog, i.e. a catalog where SimLoop and the
#    simulation program exists.
# 3. One runs SimLoop with the SERIES file as argument. SimLoop copies the
#    SERIES file to a file named SERIES in the simulation directory.
# 4. SimLoop updates the counter within the original SERIES file (but not 
#    SERIES in the simulation directory) to the next position,
#    so that other SimLoop instances can start new simulations from the same
#    SERIES file.
# 5. The simulation program reads parameters from SERIES. The simulation is 
#    performed.
# 6. A new directory with a name based on the current date and time is created
#    to save the results of a simulation. All files needed to run the 
#    simulation, as well as the results files are copied into the new 
#    directory.
# 7. SimLoop writes the name of the directory where results are saved onto 
#    the SERIES file for future reference.
# 8. Old files (files ending with ~) are removed.
# 9. Steps 3-8 are repeated until all simulations have been performed.
#10. The contents of SERIES is copied to the original SERIES file.

# Author: Fredrik Edin, 2004
# Address: freedin@nada.kth.se



import sys, os, shutil, glob, string
from time import time, localtime, strftime, gmtime
from string import split

# Function which reads rows in the SERIES file until it comes to a row 
# not beginning with #.
def skipComment():
	pos = f.tell()
	#print pos
	i = 1
	s = f.readline()
	#print str( i ) + ' ' + s
	while s and s[0] == '#':
		pos = f.tell()
		s = f.readline()
		i = i + 1
		#print str( i ) + ' ' + str(pos) + ' ' + s
	f.seek( pos )

# Initiation
finished = 0
iter = 1

# Step 3 above
name = sys.argv[1]
if sys.platform != 'win32':
	pathdelimiter = '/'
	home = '/'
	pathstring = '/DATA/'
else:
	pathdelimiter = '\\'
	home = 'C:\\'
	pathstring = '\\DATA\\'

# Transform name into name with absolute path
if name[0] != '~' and string.find( name, home ):
	name = os.getcwd()+pathdelimiter+name
else:
	name = 'C:\\Documents and Settings\\Fredrik Edin\\Mina dokument\\Neuron\\STANDARDFILER\\' + name[1:] # CHANGE #

if name == os.getcwd() + pathdelimiter + 'SERIES.txt':
	print 'ERROR: The SERIES file should be something other than the'
	print 'file SERIES in the simulation directory'
	sys.exit()	


# Determine file name, stripped of path
pathpos = string.rfind( name, pathdelimiter )
if pathpos >= 0:
	onlyname = name[pathpos+1:]
else:
	onlyname = name



# Loop with steps 3-8 above
while finished == 0:

	# Step 3 above
	# copy file to SERIES (which MultiModuleWMNetXP.hoc reads) without
	# counter updated to next position
	shutil.copy( name, 'SERIES.txt' )

	# Step 4 above
	f = open( name, 'r+' )

	# read comments at head of file
	skipComment()
	
	# Read counter position
	s = f.readline()
	currpos = f.tell()
	curr = int( f.readline() )
        
	# skip to current simulation (the one to be simulated by the 
	# simulation program)
	rg = range( curr )
	skipComment()
	for i in rg:
		f.readline()
		skipComment()

	# read first row of current simulation to learn the number
	# of rows of the simulation blocks
	firstrowpos = f.tell()
	firstrow = f.readline()
	restoffile = f.read()
	f.close()

	# Update counter
	num = int( split(firstrow)[2] )
	newcurr = curr + num
	strcurr = '0000' + str( newcurr )
	strcurr = strcurr[len( strcurr )-5:] + '\n'
	f = open( name, 'r+')
	f.seek( currpos )
	f.write( strcurr )
	f.close()
	print "Simulation #: " + str( iter )
	starttime = time()
	# End of step 4

	# Step 5 above
	if sys.platform == 'win32':
		os.system( 'nrniv MultiModuleWMNetXP.hoc' )
		#print "No simulation"
	elif not string.find( sys.platform, 'linux' ):  # 0 = found on place 0
		os.system( 'i686/special MultiModuleWMNet.hoc' )
	elif not string.find( sys.platform, 'sun' ):    # 0 = found on place 0
		os.system('nice sparc/special MultiModuleWMNet.hoc')

	dt = gmtime( time() - starttime )
	if dt[3] != 0:
	    tstr = "%H:%M:%S"
	else:
	    tstr = "%M:%S"

	print "Simulation finished. Elapsed time: " + strftime( tstr,dt )

	# Step 6 above
	now=os.getcwd() + pathstring + strftime( "%y%m%d.%H%M.%S",localtime( time() ) )

	os.mkdir( now )
	shutil.copy( 'SERIES.txt', now+pathdelimiter+onlyname )
	files_to_copy = ['SimLoopXP.py','ECell.hoc','ICell.hoc','MultiModuleWMNetXP.hoc','MyRandom.hoc','Net.hoc','Results.hoc','Simulator.hoc','View.hoc','rAMPA.txt','rNMDA.txt','RandomSeed.py'] + glob.glob( '*.mod' ) + glob.glob( 'LFP*.txt' )
	for file in files_to_copy:
		shutil.copy( file, now )
	
	files_to_move = ['C.txt','CKey.txt','Q.txt','QKey.txt','Connections.txt','Params.txt','ParamKey.txt','Parameters.txt','APs.txt','SERIES.txt'] + glob.glob( 'CURRENT*.txt' ) + glob.glob( 'X*' ) + glob.glob( 'SEED' )
	for file in files_to_move:
		shutil.copy( file, now )
		#os.remove( file )

	
	# Steps 7 above.
	# skip to current simulation (the one just simulated by the 
	# simulation program) and
	f = open( name, 'r+')
	skipComment()
	s = f.readline()
	s = f.readline()
	rg = range( curr )
	skipComment()
	for i in rg:
		s = f.readline()
		skipComment()
	firstrowpos = f.tell()
	firstrow = f.readline()
	restoffile = f.read()
	f.seek( firstrowpos )

	# a) write path to directory containing results onto SERIES
	firstrow = firstrow[:len(firstrow)-1]+'\t'+now+'\n'
	f.write( firstrow )
	f.write( restoffile )

	# b) see whether there are any simulations left
	f.seek( firstrowpos )
	f.readline()
	skipComment()
	rg = range( num-1 )
	for i in rg:
		f.readline()
		skipComment()

	s = f.read()
	if len(s) < 2:
		finished = 1
	f.close()
	# End of step 7

	# Remove old files
	# Step 8 above.
	if os.path.exists( name+'~' ):
		os.remove( name+'~' )
	filestoremove = glob.glob( '*~' )
	for fil in filestoremove:
		os.remove( fil )
	iter = iter + 1

# Old stuff not used at the moment
# save series name in AllSeries file
#fpath = 'C:\Neuron\mixed_synapses_temp\multimodulenet'
#f = open( fpath + 'AllSeries', 'a')
#f.write( name + '\n' )
#f.close()
#os.remove( fpath + '*~' )

