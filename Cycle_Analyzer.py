#!/usr/bin/env python

import os

thelist=[1431]+list(range(677,684))
for j in thelist:
	print("\nhistograms_Run%05d.root\n" % (j))
	os.system("python3 reconstruction.py configFile.txt -r %d" % (j))


