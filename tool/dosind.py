#!/usr/bin/env python
import glob, sys
print [ max ( ([ int( x [ x.index('.') + 1 : ] ) for x in glob.glob(str(i) + '/doshistory/weight.*') ] + [-1]) ) for i in xrange( 8 if len(sys.argv) < 2 else int(sys.argv[1]) ) ]
