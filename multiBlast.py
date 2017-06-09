# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import os,sys
import argparse as ap
import gzip
from Bio import SeqIO
import re


def getInputFiles(indir,inputFormat):
    
    for file in os.listdir(indir):
        
        if inputFormat == 'genbank':
            
            pass
        
        


def validateDir(directory):
    
    if os.path.isdir(directory):
        
        return True
    else:
        
        return False

def arguments():
    
    global indir
    global outdir
    global inputFormata
    
    desc = "A generalized tool for performing blasts on any number of genomes."
    
    parser = ap.ArgumentParser(description=desc)
    parser.add_argument('-i','--input',action='store',dest='indir')
    parser.add_argument('-o','--output',action='store',dest='outdir')
    parser.add_argument('-f','--format',action='store',dest='inputFormat')
    
    args = parser.parse_args()
    
    indir = args.indir
    outdir = args.outdir
    inputFormat = args.inputFormat
    
    #if not validateDir(indir):
    
        #print("Invalid Path to Input Directory")
        #parser.print_help()
        #sys.exit(1)
        
    if inputFormat in ['genbank','fasta']:
        
        if inputFormat == 'genbank':
            
            print("in Genebank")
    else:
    
        print("Invalid Input File Format")
        parser.print_help()
        sys.exit(1)
        
    return indir,outdir

def genbankToFasta(genbank_file):

    name = re.search(r'(G.+)_genomic.gbff.gz', genbank_file)
    
    fasta_file = name.group(1) + '.faa'
    
   
    g_file = gzip.open(genbank_file, 'rb')
    f_file = open(fasta_file, 'w')
    
    organism = ''
    
    for record in SeqIO.parse(g_file, 'genbank'):
        if record.features:
            for feature in record.features:
                if feature.type == 'source':
                    
                    if 'organism' in feature.qualifiers:
                        organism = str(feature.qualifiers['organism'][0])
                if feature.type == 'CDS':
                    
                    if 'protein_id' in feature.qualifiers and 'locus_tag' in feature.qualifiers and 'product' in feature.qualifiers and 'translation' in feature.qualifiers:
                        
                        f_file.write('>' + str(feature.qualifiers['protein_id'][0]) + '|')

                        f_file.write(str(feature.qualifiers['locus_tag'][0]) + '|')

                        f_file.write(str(feature.qualifiers['product'][0]) + '|[' + organism + ']')

                        translation = str(feature.qualifiers['translation'][0])
                        
                        for i in range(0, len(translation), 81):
                            translation = translation[:i] + '\n' + translation[i:]
                        f_file.write(translation + '\n')
    f_file.close()
    g_file.close()

if __name__ == "__main__":
    
    indir,outdir = arguments()

    