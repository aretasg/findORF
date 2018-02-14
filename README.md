# find-open-reading-frames
A script to find open reading frames (ORF) in the nucleotide sequence in an input as FASTA file. Searches all 6 reading frames of the sequence on the assumption that 'ATG' is a start codon and 'TAG/TGA/TAA' is a stop codon and should be suitable for reliable gene leads in prokaryotic genomes.

The script was written in and intended to be used with Python 3.6 although it should work with Python 2.7.
The script requires BioPython 1.7 to work which can be installed which can be do with "pip install biopython" command on the CLI.

To start using the script just simply run the controller.py script with find_ORF.py located in the same folder.
controller.py supports multiple CLI options supplied by the user:

* -i <file name> a flag to specify the filename and/or the directory (note that this is the only required option);
* -l <number> setting minimum nucleotide length for the ORF detection (default = 300);
* -p translating ORFs into peptide sequences;
* -t <number> selecting the translation table (see BioPython documentation: www.biopython.org);
* -o option to write an output to a FASTA file;
* -n option to include ORFs with ambiguous calls;
* -c <codon> option to search ORFs with a start codon specified by the user (default = 'ATG').

Example input and output files are provided (see foo files)
Example usage in CLI: "python controller.py -i seq_foo.fasta"

More information regarding the use and options will be added later and in the meantime use --help flag for information on all options and parameters.
