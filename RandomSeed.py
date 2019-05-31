#!/usr/bin/env python 

import sys, os, shutil, glob, string
from time import time, localtime, strftime, gmtime
from string import split

f = open( "SEED.txt", 'w' )
a = list(str(long(time()*1000)))
a.reverse()
a = ''.join(a[:9])
f.write(a)
f.close()


