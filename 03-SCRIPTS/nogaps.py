#!/usr/bin/env python

# This program will:
#    1) open a file of FASTA formatted sequences
#    2) write two output files:
#       [filename]_wholeSeqs
#          sequences with a set number of gaps or fewer ('MaxGaps')
#       [filename]_withGaps
#          sequences with more than the maximum number of gaps

# To use this program:
# Enter the name of the file to be opened as an argument when running the program

# The file to be opened must be formatted:
#     1) with the name & sequence on separate lines
#    2) with the symbol '-' to indicate gaps in sequences
#     3) with no empty lines


import re
import sys

MaxGaps = 5 #Output file '_wholeSeqs' will only contain sequences with 5 or fewer gaps

InFileName = sys.argv[1]
OutFileName = sys.argv[1] + '_wholeSeqs'
OutFileName2 = sys.argv[1] + '_withGaps'


InFile = open(InFileName, 'r')
OutFile = open(OutFileName, 'a')
OutFile2 = open(OutFileName2, 'a')


for Line in InFile:
    Line=Line.strip('\n')
    if Line[0] == '>': #Do this only for the name line
           Name = Line

    else: #Do this for the sequence line
        Gaps = re.findall('\-',Line)     #Find the number of gaps in the sequence
        if len(Gaps) < MaxGaps + 1:        #If there are MaxGaps or fewer, write sequence and name to OutFile
            OutFile.write(Name + '\n')
            OutFile.write(Line + '\n')
        else:
            OutFile2.write(Name + '\n') #If there are more than MaxGaps, write name and sequence to OutFile2
            OutFile2.write(Line + '\n')


InFile.close()
OutFile.close()
OutFile2.close()
