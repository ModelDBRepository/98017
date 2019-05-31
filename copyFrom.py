#!/usr/bin/env python
# This program moves files from a remote computer. If you have a subdirectory
# tree which is a copy of the subdirectory tree at dayhoff, files are
# tomatically moved to the subdirectory tree on your computer, and all
# names of search paths will be updated. Only UNIX.
#
# Arg 1: The name of the remote host
# Arg 2: The name of the series located in the frS catalog below
# Arg 3: The name relative to the frF catalog of the catalog in which the
#        simulation files are located
#(Arg 4: The location where the series file is to be placed relative to
#        the toS catalog, unless the path starts with '/' or '~', in
#        which case the path is absolute or relative to the home catalog.)
#
# Change path names according to your own path organization. See # CHANGE #
# below.
# 
# Author: Fredrik Edin, 2004
# Address: freedin@nada.kth.se



import sys, os
from string import replace, rfind

# print usage
if len( sys.argv ) < 3:
	print "Usage:"
	print "Arg 1: Name of remote host"
	print "Arg 2: Location of remote SERIES"
	print "Arg 3: Location of remote simulation catalog"
	print "Arg 4: Location of SERIES in new file system"
	sys.exit()
	

# Determines operating system
if sys.platform == 'win32':
	pathdelimiter = '\\'
else:
 	pathdelimiter = '/'

# This copies a series file from a remote host and saves the file name 
# of the series file in scpstr
remotehost = sys.argv[1]
serienamn = sys.argv[2]

# frS: Location of series file in remote host relative to home catalog
frS = 'Neuron/SIM_SERIER/' + serienamn                     # CHANGE #

# frF: Location of simulation catalog in remote host relative to home catalog
frF = 'Neuron/MYSIMS/' + sys.argv[3] + '/DATA/*'           # CHANGE #

# toS: Where you want to move series
toS = '/disk0/NOBACKUP/SIM_SERIER'                         # CHANGE #
if len( sys.argv ) > 4:
	tmp = sys.argv[4]
	if tmp[0] == '/' or tmp[0] == '~':
		toS = tmp
	else:
		toS = toS + '/' + tmp

# toF: Where you want to move simulation catalog
toF = '/disk0/NOBACKUP/' + sys.argv[3] + '/DATA/'      # CHANGE #

# copy files from remote host
scpstr = 'scp -ru ' + remotehost + ':' + frS + ' ' + toS
os.system( scpstr )
scpstr = 'scp -ru ' + remotehost + ':' + frF + ' ' + toF
os.system( scpstr )
thisdir = os.getcwd()
os.chdir( toF )
filestomove = glob.glob( 'DATA' + pathdelimiter + '*' )
for name in filestomove:
	shutil.move( name, '.' )
os.rmtree( 'DATA' )
os.chdir( thisdir )

# Change path names in files
fr = '/home/freedin/Neuron/MYSIMS/'                    # CHANGE #
to = '/disk0/NOBACKUP/'                                # CHANGE #
if rfind( toS, '/' ) != len( toS ) - 1:
	toS = toS + '/'
f = open(toS + serienamn, 'r')
g = open('CHSERIES', 'w')
s = f.readline()
while s:
	#ind = rfind( s, '\t' )
	#s1 = s[:ind]
	#s2 = s[ind+1:]
	#s2n = replace( s2, fr, to )
	#g.write( s1 )
	#g.write( s2n )
	str = replace( s, fr, to )
	g.write( str )	

	# copy files in the series to host computer
	#scpstr = 'scp -r
	s = f.readline()
f.close()
g.close()
shutil.move( 'CHSERIES', toS+serienamn )
	
# save series name in AllSeries file
# Old stuff not used right now
#f = open('/afs/nada.kth.se/home/o/u1sxc4xo/Neuron/Program/STANDARDFILER/AllSeries', 'a')                             # CHANGE #
#f.write( toS + serienamn + '\n' )
#f.close()
