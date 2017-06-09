# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import os,sys
import argparse as ap

def getInputFiles(indir):
    
    for file in os.listdir(indir):
        
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

if __name__ == "__main__":
    
    indir,outdir = arguments()

    