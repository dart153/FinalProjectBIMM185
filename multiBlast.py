#! /usr/bin/env python
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import os,sys
import numpy as np
import subprocess as sp
import argparse as ap
import gzip
import itertools
from Bio import SeqIO
import re

class Genome:

    def __init__(self,path,outdir):

        self.name = ''
        self.path = path
        self.blastdb = outdir+'/blastdb'


        #A list of the proteins
        self.proteins = []

        self.getName()
        self.readFile()
        self.makeBlastDB()

    def getName(self):

        filename = os.path.basename(self.path)

        self.name = re.match(r'([\w]+).',filename).group(1)

    def readFile(self):

        print(os.path.basename(self.path).split('.')[-1])

        if os.path.basename(self.path).split('.')[-1] == 'faa':

            for record in SeqIO.parse(self.path,'fasta'):

                self.proteins.append(record.id)
        else:

            with gzip.open(self.path) as f:

                for record in SeqIO.parse(f,'fasta'):

                    self.proteins.append(record.id)

    def makeBlastDB(self):

        self.blastdb = self.blastdb + '/{}'.format(self.name)

        if not os.path.exists(self.blastdb):
            os.makedirs(self.blastdb)

        command = ''

        if os.path.basename(self.path).split('.')[-1] == 'gz':

            command = 'gzcat {} | makeblastdb -dbtype prot -input_type fasta '.format(self.path)
            command += '-parse_seqids -hash_index -out ./{} -title {} -in -'.format(self.blastdb,self.name)


        else:

            command = 'makeblastdb -in {} -dbtype prot -input_type fasta '.format(self.path)
            command += '-parse_seqids -hash_index -out ./{} -title {}'.format(self.blastdb,self.name)

        os.system(command)


class Orthologs:

    def __init__(self,genome1,genome2,outdir,evalue,pident,scov):

        self.g1 = genome1
        self.g2 = genome2
        self.outdir = outdir +'/results'
        self.evalue = evalue
        self.pident = pident
        self.scov = scov

        #store the results

        self.columns = ['qseqid','sseqid','evalue','pident','length','slen']


        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        print('Evaluating {} and {} for orthologs'.format(self.g1.name,self.g2.name))

        #Run the blasts
        table1 = self.runBlastP(self.g1,self.g2)
        table2 = self.runBlastP(self.g2,self.g1)

        #Calculate orthologs
        self.BDBH(table1,table2)

    def runBlastP(self,g1,g2):

        results = '{}/{}_vs_{}.tsv'.format(self.outdir,g1.name,g2.name)

        command = ''
        if not os.path.exists(results):

            if os.path.basename(g1.path).split('.')[-1] == 'gz':

                command = 'gzcat {} | blastp -out {} '.format(g1.path,results)
                command += "-db {} -evalue {} -outfmt '6 qseqid sseqid evalue pident length slen' ".format(g2.blastdb,self.evalue)

            else:

                command = 'blastp -query {} -out {} '.format(g1.path,results)
                command += "-db {} -evalue {} -outfmt '6 qseqid sseqid evalue pident length slen'".format(g2.blastdb,self.evalue)

            os.system(command)

        #Read into table
        return self.loadToTable(results)

    def loadToTable(self,results):

        table = pd.read_table(results,header=None,sep='\t',names=self.columns)

        return self.processTable(table)

    def processTable(self,table):

        table['scov'] = np.NaN

        for index,row in table.iterrows():

            #print(row)
            sseqid= re.match(r'ref\|([\w\_\.]+)\|',str(row[1])).group(1) # I changed this because mine didn't have the ref part (not sure why)
            table.set_value(index,'sseqid',sseqid)

            scov= (float(row[4])/float(row[5]))*100
            table.set_value(index,'scov',scov)

        print(table.head())

        return table

    def BDBH(self,table1,table2):

        orthologs = open('{}/{}vs{}_ortho.txt'.format(self.outdir,self.g1.name,self.g2.name),'w+')

        orthologs.write('{}\t{}\n'.format('\t'.join(self.columns),'status'))

        queries = table1.qseqid.unique()

        for query in queries:

            #Subset and sort the list with only the particular query
            comparisons = self.subset(table1,query)

            if comparisons.shape[0] > 0:
                subject = comparisons['sseqid'].iloc[0]

                subjectHits = self.subset(table2,subject)

                if subjectHits.shape[0] > 0:
                    if query == subjectHits['sseqid'].iloc[0]:

                        orthologs.write('{}\t{}\n'.format(self.rowToString(subjectHits.iloc[0]),'BDBH'))
                    else:

                        self.topHits(subjectHits,orthologs)

        orthologs.close()

    def topHits(self,table,file):

        comparisons = table.loc[table['evalue'] == table['evalue'].iloc[0]]

        if comparisons.shape[0] > 0:
            for index,row in comparisons.iterrows():

                file.write('{}\t{}\n'.format(self.rowToString(row),'TOP'))


    def subset(self,table,query):

        #Subset the table based on each of the parameters
        comparisons = table.loc[table['qseqid'] == query]
        comparisons = comparisons.loc[comparisons['evalue'] <= self.evalue]
        comparisons = comparisons.loc[comparisons['pident'] >= self.pident]
        comparisons = comparisons.loc[comparisons['scov'] >= self.scov]

        comparisons = comparisons.sort_values(['evalue','pident','scov'], ascending=[True,False,False])

        return comparisons

    def rowToString(self,row):

        row = row.tolist()

        return '\t'.join([str(i) for i in row])


class Paralogs:

    def __init__(self, genome, outdir, evalue, ident, scov):

        self.g1 = genome
        self.outdir = outdir +'/results'
        self.evalue = evalue
        self.ident = ident
        self.scov = scov

        self.columns = ['qseqid','sseqid','evalue','ident','length','slen']

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        print('Evaluating {} for paralogs'.format(self.g1.name))

        #Run blast
        table = self.runBlastP(self.g1)
        print(table.head())
        print(table.describe())


    def runBlastP(self, g1):

        results = '{}/{}.tsv'.format(self.outdir,g1.name)

        command = ''
        if not os.path.exists(results):

            if os.path.basename(g1.path).split('.')[-1] == 'gz':

                command = 'gzcat {} | blastp -out {} '.format(g1.path,results)
                command += "-db {} -evalue {} -outfmt '6 qseqid sseqid evalue pident length slen' ".format(g1.blastdb,self.evalue)

            else:

                command = 'blastp -query {} -out {} '.format(g1.path,results)
                command += "-db {} -evalue {} -outfmt '6 qseqid sseqid evalue pident length slen'".format(g1.blastdb,self.evalue)

            os.system(command)

        return self.loadToTable(results)


        #load into pandas dataframe
        #get things evalue lower other 2 higher than cutoff
        #ignore self hits (subj & target are same)

    def loadToTable(self,results):

        table = pd.read_table(results,header=None,sep='\t',names=self.columns)

        return self.processTable(table)

    def processTable(self,table):

        table['scov'] = np.NaN

        for index,row in table.iterrows():

            sseqid= re.match(r'\|*([\w\_\.]+)\|',str(row[1])).group(1)
            table.set_value(index,'sseqid',sseqid)

            qseqid = re.match(r'\|*([\w\_\.]+)\|*',str(row[0])).group(1)
            table.set_value(index,'qseqid',sseqid)

            # ignore self hits



            scov= (float(row[4])/float(row[5]))*100
            table.set_value(index,'scov',scov)

        table = self.subset(table)

        return table


    def BDBH(self,table):

        #queries = self.table1.loc['qseqid'].unique()

        print(table.describe())


    def subset(self, table):
        comparisons=table
        comparisons = comparisons.loc[comparisons['evalue'] <= self.evalue]
        comparisons = comparisons.loc[comparisons['ident'] >= self.ident]
        comparisons = comparisons.loc[comparisons['scov'] >= self.scov]
        comparisons = comparisons.loc[comparisons['qseqid'] != comparisons['sseqid']]

        comparisons = comparisons.sort_values(['evalue','ident', 'scov'], ascending = [True, False, False])

        return comparisons


def getInputFiles(indir,inputFormat,outdir):

    global genomes

    for file in os.listdir(indir):

        if not os.path.isdir(file):

            if inputFormat == 'genbank':

                path = indir + '/' + file

                infile = genbankToFasta(path,outdir)


                genome = Genome(infile,outdir)


                genomes.append(genome)

            else:

                path = indir + '/' + file
                genome = Genome(path,outdir)

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
                        default=25.0,
                        help='The minimum percent identity to be considered for homology.')
    parser.add_argument('-scov',
                        action='store',
                        dest='scov',
                        type=float,
                        default=80.0,
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



    return indir,outdir



def genbankToFasta(genbank_file,outdir):

    name = re.search(r'(G.+)_genomic.gbff.gz', genbank_file)

    outdir = outdir + '/fastas'

    if not os.path.exists(outdir):
            os.makedirs(outdir)

    fasta_file = outdir+ '/'+name.group(1) + '.faa'


    g_file = gzip.open(genbank_file, 'rb')
    f_file = open(fasta_file, 'w')


    for record in SeqIO.parse(g_file, 'genbank'):

        if record.features:
            for feature in record.features:

                if feature.type == 'CDS':

                    if 'protein_id' in feature.qualifiers and 'locus_tag' in feature.qualifiers and 'translation' in feature.qualifiers:
                        f_file.write('>' + str(feature.qualifiers['protein_id'][0]) + '|')

                        f_file.write(str(feature.qualifiers['locus_tag'][0]))

                        translation = str(feature.qualifiers['translation'][0])

                        for i in range(0, len(translation), 81):
                            translation = translation[:i] + '\n' + translation[i:]
                        f_file.write(translation + '\n')
    f_file.close()
    g_file.close()

    return fasta_file

def findOrthologs():

    if len(genomes) > 1:
        comparisons = itertools.combinations(genomes,2)


        for comp in comparisons:

            ortholog = Orthologs(comp[0],comp[1],outdir,evalue,ident,scov)

def findParalogs():

    for genome in genomes:
        paralog = Paralogs(genome,outdir,evalue,ident,scov)

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
    findOrthologs()
    findParalogs()
