#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import argv


InFile = open (argv[1], 'r')
OutFile = open(argv[2], 'w')

LineNumber = 0
for line in InFile:
      Line = line.strip('\n')
      
      column = line.split('\t')
      
      start = 0
      
      OutputColumn = ("{}\t{}\t{}\t{}\t".format(column[0],start,column[1],column[2]))
      
      OutFile.write(OutputColumn + "\n")
      
      LineNumber = LineNumber + 1

      
InFile.close()
OutFile.close()


import os
os.system('cat chromInfo.bed | sort -k 1,1 -k2,2n| uniq > chromSort.bed')
os.system('cat dm3_exons.txt | sort -k 1,1 -k2,2n| uniq > dm3_exonsSort.bed')
