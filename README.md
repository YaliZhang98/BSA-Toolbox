# Sequence Analysis Toolbox

## Introduction
In this program, Flask, a popular framework in python is used to build a web application. This application is a toolbox that can get a uploaded fasta or fastq file or a pure nucleotide sequence submitted by users, and then performs different analysis and functions on nucleotide sequences.  

There are five different functions in this toolbox:  
1) GC content calculation. In this function, GC content of each sequence in the uploaded file will be calculated separately. In additon, the total GC content of all nucleotide sequences in file will be calculated as well.  
2) Sequence length detection. In this function, the maximum and minimum length of sequences in uploaded file will be found, and the specific sequence with the longest and shortest length also will be pointed at the same time. In addition, the average length of all sequences in file will be calculated.  
3) Reversed complement of nucleotide sequences. In this function, reversed complementary sequence of each nucleotide sequence in the uploaded file will be found respectively. All reversed complements will be output in a fasta format file.  
4) Translation of nucleotide sequence. In this function, all nucleotide sequence in the uploaded file will be translated to corresponding amino acid sequence separately. All amino acid sequences will be output in a fasta format file.  
5) GC content plot. In this function, all sequences in the upload file will be merged into one sequence, and the GC content changes will be ploted over this one sequence. The window size and step will be defined by user.  
  
## Usage






