# FastaSplitter
A simple Python3 application, capable of splitting fasta sequences into sequences of smaller lengths

USAGE: 

python3 Splitter.py {FastaFile} {Length} {Options} > {Output}

The program will take the Fasta File and split it into smaller sequences of the desired Length. 

Options:
  -o={X}; --overlap={X}   An overlap of {X} bases will be made between sequences.
  -c={Y}; --cov={Y}       The overlap will be calculated in order to generate the specified {Y} coverage.

If a sequence is not multiple of the given 'length', an additional sequence will be created using the last
'length' bases, making an overlap with the previous one.



