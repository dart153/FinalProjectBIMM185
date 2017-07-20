#! /usr/bin/env python
"""
multiBlast is a generalized tool that allows for multiple proteomes
in Genbank or FASTA format to be BLASTed and compared for both orthologs
and paralogs.
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


    def removeHomolog(self,id):

        if id in self.proteins:

            self.proteins.remove(id)


class Orthologs:

    def __init__(self,genome1,genome2,outdir,evalue,pident,scov):

        self.g1 = genome1
        self.g2 = genome2
        self.outdir = outdir +'/results'
        self.evalue = evalue
        self.pident = pident
        self.scov = scov

        self.ortho_ab = {}
        self.ortho_ba = {}


        #store the results
        self.columns = ['qseqid','sseqid','evalue','pident','length','slen']


        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        if not os.path.exists('{}/blast_results'.format(self.outdir)):

            os.makedirs('{}/blast_results'.format(self.outdir))

        print('Evaluating {} and {} for orthologs'.format(self.g1.name,self.g2.name))

        #Run the blasts

        #subject: g2;query: g1
        table1 = self.runBlastP(self.g1,self.g2)

        #subject: g1;query: g2
        table2 = self.runBlastP(self.g2,self.g1)

        #Calculate orthologs
        self.orthologTable = self.BDBH(table1,table2)

        #print(self.orthologTable)

        #print('ortho_ab: {}'.format(self.ortho_ab))
        #print('ortho_ba: {}'.format(self.ortho_ba))

    def runBlastP(self,g1,g2):

        results = '{}/blast_results/{}_vs_{}.tsv'.format(self.outdir,g1.name,g2.name)

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
            sseqid= re.match(r'ref\|([\w\_\.]+)\|',str(row[1])).group(1)
            table.set_value(index,'sseqid',sseqid)

            scov= (float(row[4])/float(row[5]))*100
            table.set_value(index,'scov',scov)

        return table

    def BDBH(self,table1,table2):

        orthologTable = pd.DataFrame(columns = self.columns)

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

                        #Proteins are not orphans
                        self.g1.removeHomolog(query)
                        self.g2.removeHomolog(subject)

                        if subject not in self.ortho_ab:

                            self.ortho_ab[subject] = []

                        self.ortho_ab[subject].append(query)

                        if query not in self.ortho_ba:

                            self.ortho_ba[query] = []

                        self.ortho_ba[query].append(subject)

                        #print(subjectHits.iloc[0])

                        #Add to table of orthologs
                        orthologTable = orthologTable.append(subjectHits.iloc[0],ignore_index=True)

                        #print(orthologTable)

                        #Write Pair to file
                        orthologs.write('{}\t{}\n'.format(rowToString(subjectHits.iloc[0]),'BDBH'))
                    else:

                        orthologTable = self.topHits(subjectHits,orthologs,orthologTable)

        orthologs.close()

        return orthologTable

    def topHits(self,table,file,orthologTable):

        comparisons = table.loc[table['evalue'] == table['evalue'].iloc[0]]

        if comparisons.shape[0] > 0:

            for index,row in comparisons.iterrows():

                subject = row['qseqid']
                query = row['sseqid']

                self.g2.removeHomolog(subject)
                self.g1.removeHomolog(query)

                if subject not in self.ortho_ab:

                    self.ortho_ab[subject] = []

                self.ortho_ab[subject].append(query)

                if query not in self.ortho_ba:

                    self.ortho_ba[query] = []

                self.ortho_ba[query].append(subject)

                orthologTable = orthologTable.append(row)

                file.write('{}\t{}\n'.format(rowToString(row),'TOP'))

        return orthologTable


    def subset(self,table,query):

        #Subset the table based on each of the parameters
        comparisons = table.loc[table['qseqid'] == query]
        comparisons = comparisons.loc[comparisons['evalue'] <= self.evalue]
        comparisons = comparisons.loc[comparisons['pident'] >= self.pident]
        comparisons = comparisons.loc[comparisons['scov'] >= self.scov]

        comparisons = comparisons.sort_values(['evalue','pident','scov'], ascending=[True,False,False])

        return comparisons



# paralogs class
class Paralogs:

    # constructor
    def __init__(self, genome, outdir, evalue, ident, scov):

        self.g1 = genome
        self.outdir = outdir +'/results'
        self.evalue = evalue
        self.ident = ident
        self.scov = scov
        self.pairs = []

        # columns for paralogs
        self.columns = ['qseqid','sseqid','evalue','ident','length','slen']

        # make output directory
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        if not os.path.exists('{}/blast_results'.format(self.outdir)):

            os.makedirs('{}/blast_results'.format(self.outdir))

        #Run blast
        table = self.runBlastP(self.g1)

        # write paralog data to a file
        self.writeParalogs(table)

    def runBlastP(self, g1):

        # output for blast data
        results = '{}/blast_results/{}.tsv'.format(self.outdir,g1.name)

        command = ''
        if not os.path.exists(results):

            # check if compressed
            if os.path.basename(g1.path).split('.')[-1] == 'gz':

                command = 'gzcat {} | blastp -out {} '.format(g1.path,results)
                command += "-db {} -evalue {} -outfmt '6 qseqid sseqid evalue pident length slen' ".format(g1.blastdb,self.evalue)

            else:

                command = 'blastp -query {} -out {} '.format(g1.path,results)
                command += "-db {} -evalue {} -outfmt '6 qseqid sseqid evalue pident length slen'".format(g1.blastdb,self.evalue)

            os.system(command)

        return self.loadToTable(results)


    def loadToTable(self,results):

        # load raw data into pandas table
        table = pd.read_table(results,header=None,sep='\t',names=self.columns)

        return self.processTable(table)

    def processTable(self,table):

        # get scov values
        table['scov'] = np.NaN

        for index,row in table.iterrows():

            scov = (float(row[4])/float(row[5]))*100
            table.set_value(index,'scov',scov)

            sseqid= re.match(r'ref\|([\w\_\.]+)\|',str(row[1])).group(1)
            table.set_value(index,'sseqid',sseqid)

        # get values in the range we want
        table = self.subset(table)

        return table


    def subset(self, table):
        comparisons = table
        comparisons = comparisons.loc[comparisons['evalue'] <= self.evalue]
        comparisons = comparisons.loc[comparisons['ident'] >= self.ident]
        comparisons = comparisons.loc[comparisons['scov'] >= self.scov]

        # sort
        comparisons = comparisons.sort_values(['evalue','ident', 'scov'], ascending = [True, False, False])

        return comparisons


    def writeParalogs(self, table):
        # write to tab separated output file
        output = open('{}/{}_para.tsv'.format(self.outdir,self.g1.name),'w')
        for index, row in table.iterrows():
            if row[0] != row[1]:
                if [row[1], row[0]] not in self.pairs:
                    if [row[0],row[1]] not in self.pairs:

                        self.pairs.append([row[0], row[1]])
                        output.write('{}\n'.format(rowToString(row)))

        output.close()

class Intersection:

    def __init__(self,genomes,orthologs):

        self.orthologs = orthologs
        self.ref = orthologs[0].g1
        self.genomes = genomes
        self.orthologDict = {}
        self.refComp = []
        self.common = set()

        self.loadGenomes(orthologs)
        self.findCommon()
        self.intersect(self.orthologs)

    def loadGenomes(self,orthologs):

        for ortholog in orthologs:

            #dictionary to find all the ortholog instances with particular genome
            if ortholog.g1.name not in self.orthologDict:

                self.orthologDict[ortholog.g1.name] = {}

            if ortholog.g2.name not in self.orthologDict:

                self.orthologDict[ortholog.g2.name] = {}

            self.orthologDict[ortholog.g1.name][ortholog.g2.name] = ortholog
            self.orthologDict[ortholog.g2.name][ortholog.g1.name] = ortholog


        #for key,value in self.orthologDict.iteritems():

            #print('{}: {}'.format(key.name,';'.join('{},{}'.format(x.g1.name,x.g2.name) for x in value)))

    def findCommon(self):

        #get only the ortholog objects that contain the reference genome
        self.refComp = self.orthologDict[self.ref.name]

        #create a set of the reference proteins from the first ortholog object
        orthologDF = self.refComp.values()[0].orthologTable

        #print(orthologDF)

        refList = orthologDF["sseqid"].tolist()
        self.common = set(refList)


    def intersect(self,orthologs):

        #flag to indicate whether gene is shared
        shared = 0


        for genome in self.genomes:

            if genome.name != self.ref.name:

                print('{}\t{}'.format(self.ref.name,genome.name))

                ortholog = self.orthologDict[self.ref.name][genome.name]

                possibleCommon = self.getRefOrthologs(self.ref,ortholog)

                self.common.intersection_update(possibleCommon)

                if not self.common:

                    print('Failed')

                    return None

        #Write the results to a file

        print('\n'.join(list(self.common)))
        '''
        #get a dictionary of sets containing the corresponding orthologs for the non-ref genomes
        orthologs = self.nonRefOrthologs(self.common)

        for i in range(len(self.genomes)-1):
            for j in range(i+1,len(self.genomes)):

                g1 = self.genomes[i]
                g2 = self.genomes[j]

                comparisons = self.orthologDict[g1.name][g2.name]

                #iterate row by row in the DataFrame
                for index,row in comparisons:

                    sseqid = row['sseqid']
                    qseqid = row['qseqid']
        '''


    def nonRefOrthologs(self,common):

        orthologs = {}

        #loop through genomes
        for genome in self.genomes:

            orthologs[genome] = set()
            #get the genes for which the refernce genes are orthologs
            ortholog = self.orthologDict[self.ref][genome.name]

            orthologs[genome] = self.getNonRefOrthologs(common,ortholog)

        return orthologs

    def getNonRefOrthologs(self,common,ortholog):

        nonRefOrthologs = set()

        orthologs = {}

        if self.ref == ortholog.g1:

            orthologs = ortholog.ortho_ab

        else:

            orthologs = ortholog.ortho_ba

        for gene in common:

            if gene in orthologs:

                nonRefOrthologs.update(orthologs[gene])

        return nonRefOrthologs

    def getRefOrthologs(self,genome,ortholog):

        #get the table of orthologs
        orthologTable = ortholog.orthologTable

        #prepare a set for the orthologs
        orthologs = set()

        if genome == ortholog.g1:

            orthologs = set(orthologTable["sseqid"].tolist())

        else:

            orthologs = set(orthologTable["qseqid"].tolist())

        return orthologs

        #for gene in self.common:

            #for genome in self.genomes:

                #ortholog = self.orthologDict[self.ref.name][genome.name]

                #if self.ref == ortholog.g1:

                    #if gene in ortholog.orthologTable['sseqid']:

                        #shared = 1
                        #refOrthologs[gene] = "null"
                    #else:

                        #shared = 0

                #if shared == 0:

                    #break

        #evaluate whether
        #for shared == 1:



def rowToString(row):

    row = row.tolist()

    return '\t'.join([str(i) for i in row])

'''
Determines the routine required depending on the input format.

--If genbank, then the file must first be converted to fasta
format
'''
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

'''
Reads the command line arguments and saves all
values for use by the program.
'''
def arguments():

    global indir
    global outdir
    global inputFormat
    global evalue
    global ident
    global scov

    desc = "A generalized tool for performing blasts on any number of genomes to determine homology."

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


'''
Converts Genbank formatted files into FASTA
formatted files for BLAST.
'''
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
'''
Creates an instance of the ortholog class for
each unique pair of genomes. Both genomes are BLASTed
against each other and orthologs are determined by
Bi-Directional Best Hit (BDBH) and Top Hit.
'''
def findOrthologs():

    if len(genomes) > 1:
        #comparisons = itertools.combinations(genomes,2)

        orthologs = []

        for i in range(len(genomes)-1):
            for j in range(i+1,len(genomes)):

                if genomes[i].name != genomes[j].name:

                    ortholog = Orthologs(genomes[i],genomes[j],outdir,evalue,ident,scov)

                    orthologs.append(ortholog)

        return genomes[0].name,orthologs

'''
Creates an instance of the paralog class for each
genome. Runs a self blast and searches for Paralogs
based on the user-inputed parameters.
'''
def findParalogs():

    for genome in genomes:
        paralog = Paralogs(genome,outdir,evalue,ident,scov)

'''
Compiles the orphan genes of each genome into a
file.
'''
def findOrphans():

    for genome in genomes:

        orphans = open('{}/results/{}_orph.txt'.format(outdir,genome.name),'w+')

        if inputFormat == 'fasta':
            orphans.write('{}'.format('\n'.join(genome.proteins)))
        else:
            for protein in genome.proteins:
                orphans.write('{}\n'.format('\t'.join(protein.split('|'))))

        orphans.close



'''
Main Method

--Runs the logic of the program through the
non class methods
'''
if __name__ == "__main__":

    #global variables
    global indir
    global outdir
    global inputFormat
    global evalue
    global ident
    global scov
    global genomes


    genomes = []

    indir,outdir = arguments()

    #Find orthologs
    ref, orthologs = findOrthologs()

    #Find paralogs
    findParalogs()

    #Find orphaned genes
    findOrphans()

    Intersection(genomes,orthologs)
