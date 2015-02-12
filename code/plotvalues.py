# -*- coding: cp1252 -*-
from __future__ import division
#----------------------------#
import numpy as N
import pylab as P
import sys, os
import os.path as path
#----------------------------#
def read_file_values(f = None):
	if f is None:
		return
	table = {}
	while True:
		line = f.readline()
		if len(line) == 0:
		    break        
		result = line.split("=")
		if len(result) != 2:
		    continue
		name,weight = result
		name = name.strip()
		#name = name.lower()

		weight = weight.strip()
		if '(' in weight or ')' in weight:
                        continue
		if '[' in weight or ']' in weight:
                        continue
		if not table.has_key(name):
		    table[name] = []
		table[name].insert(len(table[name])-1,float(weight))
	return table
#----------------------------#
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "USAGE:\tpython "+path.basename(sys.argv[0])+" <filepath> <valuename>"
		exit(1)
	path = sys.argv[1]
	valuename = sys.argv[2].lower()
	f = open(path)
	table = read_file_values(f)
	order = [valuename]
	for name in order:
		P.plot(N.array(table[name]), linewidth = 2)
		print table[name]
	P.legend(["valuename"])
	P.show()
