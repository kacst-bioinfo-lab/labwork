#!/usr/bin/env python
# -*- coding: utf-8 -*-

def group(lst, n):# This is a generator function to read N items of a batch size n:
  for i in range(0, len(lst), n):
    val = lst[i:i+n]
    if len(val) == n:
      yield tuple(val)
    

def Spacer_e(string):# This is a function to sum the uppercase letters for character in string:
    return sum(c.isupper() for c in string)

   
def CNS_e(string):# This is a function to sum the lowercase letters for character in string:
    return sum(c.islower() for c in string)
       
output_file = open('nucleotide_counts.bed','w') # the output file
output_file.write('CA\tCB\tCR\tCB-CA\tCA-CR\tCB-CR\tSA\tSB\tSR\tSB-SA\tSA-SR\tSB-SR\n')#This is for the columns names


from Bio import SeqIO
all_seqs = []# empty variable to be used later to append all the seqs
for seq_record in SeqIO.parse("Taxo.fas", "fasta"):
    all_seqs.append(seq_record.seq.split()[0])

data = [line.strip().split() for line in open("log1.txt")]
for r, row in enumerate(group(data,3)):
        x=int(data[0][0])
        y=int(data[1][0])
        z=int(data[2][0])

for i, s in enumerate(group(all_seqs,3)):#group each three sequence in one batch of the size n & assigned it to batch s in i index
    A, B, R = s # This is to itertate inside each tupled batch & assign each item to its corresponding variable
    As = (A[x-1:y]).lower()
    Af = (A[y:z]).upper()
    As2 = (A[z:]).lower()
    Bs = (B[x-1:y]).lower()
    Bf = (B[y:z]).upper()
    Bs2 = (B[z:]).lower()
    Rs = (R[x-1:y]).lower()
    Rf = (R[y:z]).upper()
    Rs2 = (R[z:]).lower()
    CA = int(CNS_e(As) + CNS_e(As2))
    CB = int(CNS_e(Bs) + CNS_e(Bs2))
    CR = int(CNS_e(Rs) + CNS_e(Rs2))
    CB_CA = int(CB-CA)
    CA_CR = int(CA-CR)
    CB_CR = int(CB-CR)
    SA = int(Spacer_e(Af))
    SB = int(Spacer_e(Bf))
    SR = int(Spacer_e(Rf))
    SB_SA = int(SB-SA)
    SA_SR = int(SA-SR)
    SB_SR = int(SB-SR)
    #Printing the output into another file
    output_line = '%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%f\n' % \
    (CA, CB, CR, CB_CA, CA_CR, CB_CR, SA, SB, SR, SB_SA, SA_SR, SB_SR)
    output_file.write(output_line)
    if i is None:#Break condition I think I don't need it cuz it worked perfectley with out it.
       break
       output_file.close()
