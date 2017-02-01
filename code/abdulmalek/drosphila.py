#!/usr/bin/env python


InputName= 'chromInfo.txt'

Input = open(InputName,'r')

LineNumber = 0
OutFileName =  'InputName.bed'

OutFile = open(OutFileName, 'w')


for line in Input:
      Line = line.strip('\n')
      
      column = line.split('\t')
      
      start = 0
      
      OutputColumn = ("{}\t{}\t{}\t{}\t".format(column[0],start,column[1],column[2]))
      OutFile.write(OutputColumn + "\n")
      
      LineNumber = LineNumber + 1
      
Input.close()
OutFile.close
