###Python Code for Parsing FASTA file
A tremendous part of bioinformatics work requires dealing with the different types of file formats designed to
hold biological data. These data often are designed in a mannar which take multiple lines inside a file to descripte one object.
The most widely used format is **fasta** file, which contains different type of attributes descriping single object.In most cases bioinformaticians deals with fasta files that contain **multiple** sequences, making the task of parsing such files even harder.

Here is a python code to parse a **fasta** file that contains multiple sequences. The idea here is to extract the *sequenceRecord* for each record in the file, then group each N *sequenceRecord* together. what comes after is upto the objective of your study. For now we have a fasta file that contains a multiple sequence for a ROOT and its children represented by A and B. The sequences contain conserved elements and a spacers, represented by lowercase and uppercase letters, respectively.

