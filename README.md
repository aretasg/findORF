# find-open-reading-frames
A script to find open reading frames (ORF) in the nucleotide sequence in input as FASTA file. Searches all 6 reading frames of the sequence on assumption that 'ATG' is a start codon and 'TAG/TGA/TAA' is a stop codon.
For the script to work BioPython must be installed which can be do with "pip install biopython" command on the CLI.

find_ORF.py supports multiple CLI options supplied by the user:
* setting minimum nucleotide length for ORF detection (default = 300);
* translating ORFs into peptide sequences;
* selecting the translation table (see BioPython documentation);
* option to write an output to a FASTA file;
* option to include ORFs with ambiguous calls;
* option to search ORFs with a start codon specified by the user (default = 'ATG').

Example input and output files are provided (see foo files)

More information regarding the use and options will be added later and in the meantime use --help flag for information on all options and parameters.
