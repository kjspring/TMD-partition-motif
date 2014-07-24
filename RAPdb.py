#!/usr/bin/python
########################################################
# Lipid Raft Associated Protein Database
# Author: Kevin Spring
# Date: 27 March 2013
# Language: Python
# Filename: RAPdb.py
########################################################
#
# Algorithm
# 1. Get a list of all lipid raft associated proteins in EBI
# 2. Extract only the ID numbers
# 3. For each ID number request data from UniProt
# 4. 


# libraries
from Bio import ExPASy
from Bio import SwissProt
import time
from Bio import SeqIO
from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from os.path import abspath
import re

# Load accession numbers and save as a list
source = "protein.list"
with open(source, 'r') as f:
    accessions = f.read().splitlines()

# iteratively go through the accession list and create a list of Bio.SwissProt.Record 
# objects containing the data
# http://www.biostars.org/p/17119/

proteinList = []
for x in accessions:
    handle = ExPASy.get_sprot_raw(x)
    try:
    	record = SwissProt.read(handle)
    except ValueException:
    	print "WARNING: Accession %s not found" % accession
    print "Downloading record complete"
    print len(proteinList)
    proteinList.append(record)
    print "Record added to list"
    time.sleep(1)

f.close()

--------------------
from Bio.SwissProt import KeyWList

dat = [] # initialize empty list
         # within each list have a protein
         # in that protein have a list of attributes
         # including ID, Number of transmembrane mentions

handle = open("keywlist.txt")
records = KeyWList.parse(handle)
for record in records:
    print(record['ID'])

# Go through each object in proteinList and save only those that have a single
# mention of TRANSMEM

noTMD = [] # initialize a list to keep track of how many TMD are in the protein
locList = [] # initialize list to save the transmembrane
tmdStart = []# initialize a list to save the start of the TMD
tmdEnd = [] # initialize a list to save the end of the TMD
entryName = [] # initialize a list to save the name of the sequence
index = [] # 
sequence = []

#http://stackoverflow.com/questions/2917372/how-to-search-a-list-of-tuples-in-python#

# Function Search for single-pass transmembrane
def transmembrane(proteinList, tmdStart, tmdEnd, sequence):
    for item in proteinList:
        index = [i for i, v in enumerate(item.features) if v[0] == "TRANSMEM"]
        noTMD = len(index)
        if noTMD == 1:
            tmdStart.append(item.features[index[0]][1]) #stores the start aa #
            tmdEnd.append(item.features[index[0]][2]) # stores the end aa #
            sequence.append(item.sequence)
            entryName.append(item.accessions[0])
             
transmembrane(proteinList, tmdStart, tmdEnd, sequence) # runs the function

locList.append(tmdStart)
locList.append(tmdEnd)

# Output the same as above, but extend the TMD sequence by 5 on each end.
outSPTMRAP = open("outSPTMRAP.fasta", "w")
i = 0
for e in sequence:
    if (locList[0][i] > 6):
        outSPTMRAP.write( ">sp|" + str(entryName[i]) + "|" + str(locList[0][i]-6) + "-" + str(locList[1][i]+5)  + "\n" + str(e[(locList[0][i]-6):(locList[1][i]+5)]) + "\n")
    else:
        outSPTMRAP.write( ">sp|" + str(entryName[i]) + "|0-" + str(locList[1][i]+5) + "\n" + str(e[0:(locList[1][i]+5)]) + "\n")
    i = i + 1

outSPTMRAP.close()

# Load the data

f_path = "outSPTMRAP.fasta"

# Open the data file and save it as a list
handle = open(f_path, "r")
tmdRecords = list(SeqIO.parse(handle, "fasta", IUPAC.extended_protein))
handle.close()
count = len(tmdRecords)

regex02 = "P[A-Z]{3,12}?G"
regex03 = "G[A-Z]{3,12}?P"

li01 = []

for i in range(0,count):
    search = re.search(regex02, str(tmdRecords[i].seq))
    if search != None:
        print tmdRecords[i]
        li01.append(tmdRecords[i])
    print i


li02 = []
for i in range(0,count):
    search = re.search(regex03, str(tmdRecords[i].seq))
    if search != None:
        print tmdRecords[i]
        li02.append(tmdRecords[i])
    print i

liAll = li01 + li02

# Remove duplicates
li = list(set(liAll))

# Output final list
output_handle = open("motifSeqFinal.fasta", "w")
SeqIO.write(li, output_handle, "fasta")
output_handle.close()

regex04 = "PC[A-Z]{2}G"

li_PCxxG = []
for i in range(0,count):
    search = re.search(regex04, str(tmdRecords[i].seq))
    if search != None:
        print tmdRecords[i]
        li_PCxxG.append(tmdRecords[i])
    print i
    
    
for i in liAll
    print i