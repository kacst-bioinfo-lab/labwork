#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''This is a program to Parse the Concerved elements [CNS] & Spacer of each sequence... where the CNS are represented by the lowercase letter,,, and the Spacers are represented by the uppercase letters. Here the sequences are grouped together and parsed for further computing'''



def group(lst, n):# This is a generator function to read N items of a batch size n:
  for i in range(0, len(lst), n):
    val = lst[i:i+n]
    if len(val) == n:
      yield tuple(val)
    

def Spacer_e(string):# This is a function to sum the uppercase letters for character in string:
    return sum(c.isupper() for c in string)


def CNS_e(string):# This is a function to sum the lowercase letters for character in string:
    return sum(c.islower() for c in string)
       
output_file = open('CNS_SPACER_Counts.bed','w') # the output file
output_file.write('CA\tCB\tCR\tCB-CA\tCA-CR\tCB-CR\tSA\tSB\tSR\tSB-SA\tSA-SR\tSB-SR\n')#This is for the columns names


from Bio import SeqIO
all_seqs = []# empty variable to be used later to append all the seqs
for seq_record in SeqIO.parse("2taxa_tree.fas", "fasta"):
    all_seqs.append(seq_record.seq.split()[0])

for i, s in enumerate(group(all_seqs,3)):#group each three sequence in one batch of the size n & assigned it to batch s in i index
    A, B, ROOT = s # This is to itertate inside each tupled batch & assign each item to its corresponding variable
    CA = CNS_e(A)# for line 29 to 40: just compute the sum of the spacer and CNS elements for each item
    CB = CNS_e(B)
    CR = CNS_e(ROOT)
    CB_CA = CB-CA
    CA_CR = CA-CR
    CB_CR = CB-CR
    SA = Spacer_e(A)
    SB = Spacer_e(B)
    SR = Spacer_e(ROOT)
    SB_SA = SB-SA
    SA_SR = SA-SR
    SB_SR = SB-SR
    #Printing the output into another file
    output_line = '%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%f\n' % \
    (CA, CB, CR, CB_CA, CA_CR, CB_CR, SA, SB, SR, SB_SA, SA_SR, SB_SR)
    output_file.write(output_line)
    if i is None:#Break condition I think I don't need it cuz it worked perfectley with out it.
        break
        output_file.close()
    
