#!/usr/bin/env python

def group(lst, n):# This is a generator function to read N items of a batch size n:
  for i in range(0, len(lst), n):
    val = lst[i:i+n]
    if len(val) == n:
      yield tuple(val)
    

def Spacer_e(string):# This is a function to sum the uppercase letters for character in string:
    return sum(c.isupper() for c in string)
  
def CNS_e(string):# This is a function to sum the lowercase letters for character in string:
    return sum(c.islower() for c in string)
       
output_file = open('nucleotide_counts_2.bed','w') # the output file
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
    Aslow = A[x:y]
    Afast = A[y:z]
    Aslow2 = A[z:]
    Bslow = B[x:y]
    Bfast = B[y:z]
    Bslow2 = B[z:]
    Rslow = R[x:y]
    Rfast = R[y:z]
    Rslow2 = R[z:]
    CA = CNS_e(Aslow) + CNS_e(Aslow2) + CNS_e(Afast)# for line 29 to 40: just compute the sum of the spacer and CNS elements for each item
    CB = CNS_e(Bslow) + CNS_e(Bslow2) + CNS_e(Bfast)
    CR = CNS_e(Rslow) + CNS_e(Rslow2) + CNS_e(Rfast)
    CB_CA = CB-CA
    CA_CR = CA-CR
    CB_CR = CB-CR
    SA = Spacer_e(Afast) + Spacer_e(Aslow) + Spacer_e(Aslow2)
    SB = Spacer_e(Bfast) + Spacer_e(Bslow) + Spacer_e(Bslow2)
    SR = Spacer_e(Rfast) + Spacer_e(Rslow) + Spacer_e(Rslow2)
    SB_SA = SB-SA
    SA_SR = SA-SR
    SB_SR = SB-SR
    #Printing the output into another file
    output_line = '%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%f\n' % \
        (CA, CB, CR, CB_CA, CA_CR, CB_CR, SA, SB, SR, SB_SA, SA_SR, SB_SR)
    output_file.write(output_line)
