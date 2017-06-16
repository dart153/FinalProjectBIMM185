# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import os,sys
import subprocess as sp
import argparse as ap
import gzip
from Bio import SeqIO
import re

class Genome:
    
    def __init__(self,path,outdir):
        
        self.name = ''
        self.path = path
        self.blastdb = outdir+'/blastdb'
        
        
        #A list of the proteins
        proteins = []
        
        self.getName()
        self.makeBlastDB()
        
    def getName(self):
        
        filename = os.path.basename(self.path)
        
        self.name = re.match(r'([\w]+).',filename).group(1)
    
    def readFile(self):

        pass
    
    def makeBlastDB(self):
        
        
        self.blastdb = self.blastdb + '/{}'.format(self.name)
        
        if not os.path.exists(self.blastdb):
            os.makedirs(self.blastdb)
        
        command = 'gzcat {} | makeblastdb -dbtype prot -input_type fasta '.format(self.path)
        command += '-parse_seqids -hash_index -out ./{} -title {} -in -'.format(self.blastdb,self.name)
        
        os.system(command)
        
class Orthologs:
    
    def __init__(self,genome1,genome2):
        
        pass
    
    def runBlastP(self):
        
        pass
        
def getInputFiles(indir,inputFormat,outdir):
    
    global genomes
    
    for file in os.listdir(indir):
        
        #print(os.path.basename(file))
        
        if inputFormat == 'genbank' and os.path.isfile(file):
            
            infile = genbankToFasta(file,outdir)
            
            path = indir + '/' + infile
            
            genome = Genome(path,outdir)
            
            
            genomes.append(genome)
            
        else:
            
            path = indir + '/' + file
            genome = Genome(path,outdir)
            
            print(genome.name)
            
            genomes.append(genome)
        
        

'''
Validates the existence of the directory that is passed in as
an argument.

-- Used for testing the presence of the input directory containing the
 genomes for comparison
'''
def validateDir(directory):
    
    if os.path.isdir(directory):
        
        return True
    else:
        
        return False


def arguments():
    
    global indir
    global outdir
    global inputFormat
    global evalue
    global ident
    global scov
    
    desc = "A generalized tool for performing blasts on any number of genomes."
    
    parser = ap.ArgumentParser(description=desc)
    parser.add_argument('-i','--input',
                        action='store',
                        dest='indir',
                        help='The path to the input directory, containing the genomes to be compared.')
    parser.add_argument('-o','--output',
                        action='store',
                        dest='outdir',
                        help='The path to the output directory, in which the results will be stored.')
    parser.add_argument('-f','--format',
                        action='store',
                        dest='inputFormat',
                        help='The format of the input files. Can be either FASTA or Genbank')
    parser.add_argument('-e','--eval',
                        action='store',
                        dest='eval',
                        default=1.0,
                        type=float,
                        help='The maximum e-value to be considered for comparison')
    parser.add_argument('-id','--ident',
                        action='store',
                        dest='ident',
                        type=float,
                        default=0.25,
                        help='The minimum percent identity to be considered for homology.')
    parser.add_argument('-scov',
                        action='store',
                        dest='scov',
                        type=float,
                        default=0.8,
                        help='The minimum percentage of the smaller protein covered.')
    
    args = parser.parse_args()
    
    indir = args.indir
    outdir = args.outdir
    inputFormat = args.inputFormat
    evalue = args.eval
    ident = args.ident
    scov = args.scov
    
    #Check for input directory
    if not validateDir(indir):
    
        print("Invalid Path to Input Directory")
        parser.print_help()
        sys.exit(1)
    
    #Check for output directory
    if not outdir:
        
        print("No Output Directory Provided.")
        parser.print_help()
        sys.exit(1)
        

    #check if input format is valid
    if inputFormat in ['genbank','fasta']:
        
        getInputFiles(indir,inputFormat,outdir)
    else:
    
        print("Invalid Input File Format")
        parser.print_help()
        sys.exit(1)
        
    if ident < 0 or ident > 1:
        
        print("Invalid identity value.")
        parser.print_help()
        sys.exit(1)
        
        
        
    return indir,outdir



def genbankToFasta(genbank_file,outdir):

    name = re.search(r'(G.+)_genomic.gbff.gz', genbank_file)
    
    fasta_file = outdir+"/fastas/"+name.group(1) + '.faa'
    
   
    g_file = gzip.open(genbank_file, 'rb')
    f_file = open(fasta_file, 'w')
    
    organism = ''
    
    for record in SeqIO.parse(g_file, 'genbank'):
        if record.features:
            for feature in record.features:
                
                if feature.type == 'CDS':
 
                        f_file.write('>' + str(feature.qualifiers['protein_id'][0]) + '|')

                        f_file.write(str(feature.qualifiers['locus_tag'][0]) + '\n')

                        translation = str(feature.qualifiers['translation'][0])
                        
                        for i in range(0, len(translation), 81):
                            translation = translation[:i] + '\n' + translation[i:]
                        f_file.write(translation + '\n')
    f_file.close()
    g_file.close()
    
    return fasta_file

if __name__ == "__main__":
    
    #global variables
    global indir
    global outdir
    global evalue
    global ident
    global scov
    global genomes
    
    
    genomes = []
    
    indir,outdir = arguments()

    