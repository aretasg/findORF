# findORF
A CLI tool to find open reading frames (findORF) in a DNA sequence as FASTA file. Searches all 6 reading frames of the sequence on the assumption that 'ATG' is a start codon and 'TAG/TGA/TAA' is a stop codon and should be suitable for reliable gene leads in prokaryotic/intron-free genomes.

**The script was written in and intended to be used with Python 3.6 although it should work with Python 2.7.**

To start using the script simply run the findORF.py script while fions.py and TABLE1.DAT are located in the same directory.
findORF.py supports multiple CLI options supplied by the user:

* -i <file name> a flag to specify the filename and/or the directory of a file (note that this is the only required option);
* -l <number> setting minimum nucleotide length for the ORF detection (default = 300);
* -p translating ORFs into peptide sequences;
* -o option to write an output to a FASTA file;
* -n option to include ORFs with ambiguous calls;
* -c <codon> option to search ORFs with a start codon specified by the user (default = 'ATG');
* -t <file name> selecting the translation table (default = TABLE1.DAT); You can replace TABLE1.DAT (GenBank Standard translation table) with a translation table of your preference that matches the default file format.

Example input and output files are provided (see foo files);
Example usage in CLI: "python findORF.py -i seq_foo.fasta -p -o"

Please use --help flag for more information on all options and parameters.
