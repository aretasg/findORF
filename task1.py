# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 15:11:52 2016

@author: kamsa
"""
import re
import argparse
import sys

''' Contents '''

''' Command line arguments parser.'''

parser = argparse.ArgumentParser(description='ORF finder: The program is intended to be used as a prokaryotic genome open reading frame (ORF) finder. ORFs are identified on the occurence of start and stop codons within the genome. The program allows for ORFs to be discovered in all 6 reading frames. By default the minimum length of an open ORF is set to 50 amino acids. This can be changed with --length argument. An output file will be returned containing the name and length of each ORF on one line, follwed by the amino acid sequence on the next, for each ORF. By default the naming convention of the output file follows the following format: (frame)_(min.length)_(inputfilename). This can be changed with the --output argument. By default the program finds ORFs in all 6 reading frames. The --frame argument can be used to specify an individual frame.')
parser.add_argument('input', nargs = '+',
                    help='The specified input file, must take the form of a txt file, where the first line containins the genome information, and the following lines contain the DNA seqeunce. E.g. .fasta')
parser.add_argument('-o', '--output', metavar='filename', 
                    help='Change the name of the output file. Remember to include the desired file extension eg. .fasta. The default file namin convention can be found in the description' )
parser.add_argument('-f', '--frame', 
                    help='The reading frame (RF) to be searched for ORFs can be specified. By default all reading frames are searched. One of the following arguments is can be used to search a single reading frame: f1 = forward RF1, f2 = forward RF2, f3 = forward RF3, r1 = Reverse RF1, r2 = Reverse RF2, r3 = Reverse RF3.'
                    ,metavar='frame')
parser.add_argument('-l', '--length', 
                    help = 'Allows the minimum length of the open reading frames to be specified in terms of amino acids. By default this value is set to 50. The argument can take the form of an integer between 1-39999.', 
                    type=int, choices=range(1, 39999), metavar='min_ORF_length')
args = parser.parse_args()

## The above code uses the argparse module in order to parse arguments to the
## command line. The input and frame argumensta are positional, whilst the 
## output and length arguments are optional.

''' Opening a file '''

for file in args.input:

    
    # Opening fasta file and forming a continous string.
    f = open(file) # Open fasta file and asign to variable f.
    next(f) # Skip the first line in the file, in this case skipping the header.
    ''' add if statement abobe for if header??'''
    dnaseq = '' # Empty string.
    
    ## Each line in the fasta file is concatenated to the previous to form a 
    ## continuous string of the DNA sequence. :-1 is used to remove \n from the 
    ## end of each line.
     
    for line in f: 
        dnaseq += line[:-1]   
    f.close() # Close file.
    fdnaseq = dnaseq.upper() # Converts all letters to upper case.
    
    
    
    
    
    ''' Finding the complementary sequence and reversing it. '''
    
    ## The .translate method is used to find the complemntary sequence. This allows
    ## every character in a string to be replaced by there key, whitch is 
    ## specified here in the dictionary oppbase. Finally the strand is reversed
    ## into the 5'-3' direction.
    
    oppbase = {ord('A'):'T', ord('C'):'G', ord('T'):'A', ord('G'):'C'} #http://stackoverflow.com/questions/18723493/python3-replace-using-dictionary
    trdnaseq = fdnaseq.translate(oppbase) 
    rdnaseq = trdnaseq[::-1] # reverse complemntary DNA into 5'-3' direction.
    
    
    
    
    
    ''' Identifying the 6 reading frames. '''
    
    ## The function bellow is used to generate a list of codons in all 6 reading
    ## frames for the forward and reverse strands. The function splits the string
    ## after every 3rd base and places the codon in a list. Each frame can be found
    ## by changing the starting point in the range (frame argument).
     
    def codonlist(seq, frame):
        codons = [] # Empty list                   # Move through seq in itterations
        for chr in range (frame - 1, len(seq), 3): # of 3, in the specified range.
            trinuc = seq[chr:chr+3] # Codon stored as variable.
            codons.append(trinuc)   # Variable appendid to list.
        return codons 
            
    fdnaseq1 = codonlist(fdnaseq, 1)  # 6 lists of codons are made for each reading
    fdnaseq2 = codonlist(fdnaseq, 2)  # frame.
    fdnaseq3 = codonlist(fdnaseq, 3)
    rdnaseq1 = codonlist(rdnaseq, 1)
    rdnaseq2 = codonlist(rdnaseq, 2)
    rdnaseq3 = codonlist(rdnaseq, 3)
    
    
    
    
    
    ''' Translating the reading frames. ''' 
    
    ## The dictionary bellow contains the 3 letter codons as values and the 
    ## corresponding amino acid letter as the key.
        
    aadict = {"AAA" : "K", "AAC" : "N", "AAG" : "K", "AAT" : "N",
    	"ACA" : "T", "ACC" : "T", "ACG" : "T", "ACT" : "T",
    	"AGA" : "R", "AGC" : "S", "AGG" : "R", "AGT" : "S",
    	"ATA" : "I", "ATC" : "I", "ATG" : "M", "ATT" : "I",
    	"CAA" : "Q", "CAC" : "H", "CAG" : "Q", "CAT" : "H",
    	"CCA" : "P", "CCC" : "P", "CCG" : "P", "CCT" : "P",
    	"CGA" : "R", "CGC" : "R", "CGG" : "R", "CGT" : "R",	
    	"CTA" : "L", "CTC" : "L", "CTG" : "L", "CTT" : "L",
    	"GAA" : "E", "GAC" : "D", "GAG" : "E", "GAT" : "D",
    	"GCA" : "A", "GCC" : "A", "GCG" : "A", "GCT" : "A",
    	"GGA" : "G", "GGC" : "G", "GGG" : "G", "GGT" : "G",
    	"GTA" : "V", "GTC" : "V", "GTG" : "V", "GTT" : "V",
    	"TAA" : "*", "TAC" : "Y", "TAG" : "*", "TAT" : "Y",
    	"TCA" : "S", "TCC" : "S", "TCG" : "S", "TCT" : "S",
    	"TGA" : "*", "TGC" : "C", "TGG" : "W", "TGT" : "C",
    	"TTA" : "L", "TTC" : "F", "TTG" : "L", "TTT" : "F" }
     
    ## The function bellow is used to translate the list of codons for each frame 
    ## generated previously. For each codon in a list the dictionary is used to 
    ## find the single letter key which corresponds to the codon with .get method.
    ## This key is then appeneded to a new list as the function itterates over each 
    ## codon from the input variable. The list of amino acid letters is then
    ## converted into a string and returned.
    
    def translator(listofcodons):
        aminoacidseq = []
        for codon in listofcodons:
            aminoacid = aadict.get(codon, 'x') #x = the value given to keys that are not found
            aminoacidseq.append(aminoacid) #forming list of amino acid letters
        aminoacidseqstring = ''.join(aminoacidseq) #Converts list to string
        return aminoacidseqstring
    
    faa1 = translator(fdnaseq1) # The translate function being used on the reading
    faa2 = translator(fdnaseq2) # frame strings generated previously.
    faa3 = translator(fdnaseq3)
    raa1 = translator(rdnaseq1)
    raa2 = translator(rdnaseq2)
    raa3 = translator(rdnaseq3)
    
    
    
    
    
    ''' Finding the ORF's. '''
    
    ## The re.findall function will be used to find all ORFs. The funcion searchs
    ## strings for a specified pattern and returns the results as a list. In this
    ## case the re function is being used to find sequences begginig with M, with 
    ## at least 50 amino acids in the middle (not incuding stop), and a stop (Z). 
     
    if args.length != True:
        minimumlength = 50
        regex = r'M[ARNDCEQGHILKMFPSTWYV]{' + str(minimumlength) + r',40000}\*'
    
    if args.length:
        minimumlength = args.length
        regex = r'M[ARNDCEQGHILKMFPSTWYV]{' + str(minimumlength) + r',40000}\*'
    
    ## Above command line argument is assigned to the variable minmumlenth. A regex
    ## has also been produced as a string that contains the minimumlength variable.
    ## This allows the minimum amino acid length to be set in the command line.
    
    forf1 = re.findall(regex, faa1)
    forf2 = re.findall(regex, faa2)
    forf3 = re.findall(regex, faa3)
    rorf1 = re.findall(regex, raa1)
    rorf2 = re.findall(regex, raa2)
    rorf3 = re.findall(regex, raa3)
    
    
    #if sys.argv[4] != int:
    #   print('The minimum amino acid length must be an integer, please refer to the manual for help.')
    
       
      
       
       
    ''' Outputting to file.'''
    # The bellow if statements are used to specify whitch frames will be written to file, whitch can be sepecified by the user. 
    if args.frame != True :  # Default all frames are written (see if statements further on)
        rframe = 'all'
    if args.frame:           # Only the user specified state will be written to file if argument is provided (see if statements further down)
        rframe = args.frame
     
    # The bellow if statements dictate the name of the file ddepending on whether an optional output file name is specified.    
    if args.output != True:  # Default - no output argument name is specified
        y = open('ORFs_{0}_{1}_{2}'.format(rframe, minimumlength, file), 'w') #open a writeable file. 
        x = 'ORFs_{0}_{1}_{2}'.format(rframe, minimumlength, file)
    if args.output:  # Output argument is specified.
        y = open('{0}'.format(args.output), 'w')
        x = '{0}'.format(args.output)
    # The xs here are used to help print an output summary coded later.    
    
        
    ## A function was produced in order to print the ORFs found to a fasta file.

    def writetofile(orf, orfnumber):
        i = 1
        for seq in orf :
            x = len(seq)  # length of open reading frame
            z = str(seq)
            organism = file.replace('.fasta','') # remove fasta from filename
            
            a = '\n'.join(z[f:f+60] for f in range(0, len(z), 60)) # This line uses the ORF string z, and devides it into lines of 60 for writting in Fasta format.
        
            y.write('> {0}_{1}_{2:03}, Frame: {1}, Length: {3}\n{4}\n'.format(
                    organism, orfnumber, i, x, a))
            i +=1
       
            
            
    summarytext = 'ORFs found in {0}, in {1} frame/s, with the minimum amino acid length: {2}.\nResults outputed to: {3}'.format(file, rframe, minimumlength, x)
    alllength = len(forf1) + len(forf2) + len(forf3) + len(rorf1) + len(rorf2) + len(rorf3)
            
    if rframe == 'all':     
        writetofile(forf1, 'F1') #The ORFs written to a fasta file in chronological 
        writetofile(forf2, 'F2') # order.
        writetofile(forf3, 'F3')
        writetofile(rorf1, 'R1')
        writetofile(rorf2, 'R2')
        writetofile(rorf3, 'R3')
        print('{0} {1}'.format(alllength, summarytext))
    elif rframe == 'f1':
        writetofile(forf1, 'F1')
        print('{0} {1}'.format(len(forf1), summarytext))
    elif rframe == 'f2':
        writetofile(forf2, 'F2')
        print('{0} {1}'.format(len(forf2), summarytext))
    elif rframe == 'f3':
        writetofile(forf3, 'F3')
        print('{0} {1}'.format(len(forf3), summarytext))
    elif rframe == 'r1':
        writetofile(rorf1, 'R1')
        print('{0} {1}'.format(len(rorf1), summarytext))
    elif rframe == 'r2':
        writetofile(rorf2, 'R2')
        print('{0} {1}'.format(len(rorf2), summarytext))
    elif rframe == 'r3':
        writetofile(rorf3, 'R3')
        print('{0} {1}'.format(len(rorf3), summarytext))
    else:
        print('Please select correct frame, refer to manual for usage by using the --help argument')
        sys.exit()
        
    
    y.close()
    
    
    
    
    
    
    


        

    
