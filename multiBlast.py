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

'''
Class Genome is used to store all the relavent data required for each Genome
that is being compared. It also creates the BLAST database for each of the genomes.
'''
class Genome:

    '''
    Constructor for the genome class.

    Parameters:
        -path: path to the input file
        -outdir: path to the output directory
    '''
    def __init__(self,path,outdir):

        self.name = ''
        self.path = path
        self.blastdb = outdir+'/blastdb'

        #A list of the proteins
        self.proteins = []

        self.getName()
        self.readFile()
        self.makeBlastDB()

    '''
    Retrieves the genome name from the input file name
    '''
    def getName(self):

        filename = os.path.basename(self.path)

        self.name = re.match(r'([\w]+).',filename).group(1)

    '''
    Reads the input file and creates a list of all the proteins
    availible in the file.
    '''
    def readFile(self):

        if os.path.basename(self.path).split('.')[-1] == 'faa':

            for record in SeqIO.parse(self.path,'fasta'):

                self.proteins.append(record.id)
        else:

            with gzip.open(self.path) as f:

                for record in SeqIO.parse(f,'fasta'):

                    self.proteins.append(record.id)

    '''
    This function creates the BLAST database for the given genome file.
    '''
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


    '''
    This function removes the proteins that are found as orthologs in other genomes,so
    that at the end of the program, so called 'orphan' genes may be found.
    '''
    def removeHomolog(self,id):

        if id in self.proteins:

            self.proteins.remove(id)

'''
Class Orthologs is the class which runs blastp on every combination of
genomes in order to determine the orthologs.
'''
class Orthologs:

    '''
    Constructor for the Orthologs class.

    Parameters:
        -genome1: the first genome object
        -genome2: the second genome object
        -outdir: the path to the output directory
        -evalue: the maximum evalue to be considered for orthology
        -pident: the minimum percent identity to be considered for orthology
        -scov: the minimum subject coverage to be considered for orthology
    '''
    def __init__(self,genome1,genome2,outdir,evalue,pident,scov):

        self.g1 = genome1
        self.g2 = genome2
        self.outdir = outdir +'/results'
        self.evalue = evalue
        self.pident = pident
        self.scov = scov

        self.status = {}


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

    '''
    Runs blastp on the current combination of genomes.
    '''
    def runBlastP(self,g1,g2):

        #initialize the output file name for the blastp results
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

    '''
    Loads the output of blastp into PANDAS tables
    '''
    def loadToTable(self,results):

        table = pd.read_table(results,header=None,sep='\t',names=self.columns)

        return self.processTable(table)

    '''
    Process the table by adding the column 'scov' (subject coverage) and calculating
    the subject coverage. Also process the subject ids to not include 'ref' in
    the name.
    '''
    def processTable(self,table):

        table['scov'] = np.NaN

        for index,row in table.iterrows():

            sseqid = ''
            try:
                sseqid= re.match(r'ref\|([\w\_\.]+)\|',str(row[1])).group(1)
            except AttributeError:
                sseqid = row[1]
            table.set_value(index,'sseqid',sseqid)

            scov= (float(row[4])/float(row[5]))*100
            table.set_value(index,'scov',scov)

        return table

    '''
    The Bi-Directional Best Hit (BDBH) function of the Ortholog class. Is the
    highest level of confidence for orthology in this program.
    '''
    def BDBH(self,table1,table2):


        orthologTable = pd.DataFrame(columns = self.columns)

        orthologs = open('{}/{}_vs_{}_ortho.tsv'.format(self.outdir,self.g2.name,self.g1.name),'w+')

        orthologs.write('{}\t{}\n'.format('\t'.join(self.columns),'status(BDBH/TOP)'))

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



                        if subject not in self.status:

                            self.status[subject] = {}

                        self.status[subject][query] = 'BDBH'


                        #Add to table of orthologs
                        orthologTable = orthologTable.append(subjectHits.iloc[0],ignore_index=True)


                        #Write Pair to file
                        #orthologs.write('{}\t{}\n'.format(rowToString(subjectHits.iloc[0]),'BDBH'))
                    else:

                        orthologTable = self.topHits(subjectHits,orthologs,orthologTable)


        #print the table
        self.printTable(orthologTable,orthologs)

        orthologs.close()

        return orthologTable

    '''
    topHits is the function for the Top Hits method of orthology. It provides
    the lowest level of confidence for orthology.
    '''
    def topHits(self,table,file,orthologTable):

        comparisons = table.loc[table['evalue'] == table['evalue'].iloc[0]]

        if comparisons.shape[0] > 0:

            for index,row in comparisons.iterrows():

                query = row['sseqid']
                subject = row['qseqid']

                self.g1.removeHomolog(query)
                self.g2.removeHomolog(subject)



                #file.write('{}\t{}\n'.format(rowToString(row),'TOP'))

                if subject not in self.status:

                    orthologTable = orthologTable.append(row,ignore_index=True)
                    self.status[subject] = {}

                    self.status[subject][query] = 'TOP'
                else:

                    if query not in self.status[subject]:

                        orthologTable = orthologTable.append(row,ignore_index=True)

                        self.status[subject][query] = 'TOP'
                    else:

                        self.status[subject][query] = 'BDBH/TOP'


        return orthologTable

    '''
    The subset function returns a subset of the table which contains all the
    entries which could be potential orthologs in the blastp results.
    '''
    def subset(self,table,query):

        #Subset the table based on each of the parameters
        comparisons = table.loc[table['qseqid'] == query]
        comparisons = comparisons.loc[comparisons['evalue'] <= self.evalue]
        comparisons = comparisons.loc[comparisons['pident'] >= self.pident]
        comparisons = comparisons.loc[comparisons['scov'] >= self.scov]

        comparisons = comparisons.sort_values(['evalue','pident','scov'], ascending=[True,False,False])

        return comparisons


    def printTable(self,table,output):

        for index,row in table.iterrows():

            status = self.status[row['qseqid']][row['sseqid']]

            output.write('{}\t{}\n'.format(rowToString(row),status))

'''
Class Paralogs evaluates potential paralogs within any given genome.
'''
class Paralogs:

    '''
    Constructor for the Paralogs class.

    Parameters:
        -genome: the genome object for which paralogs are being evaluated
        -outdir: the path to the output directory
        -evalue: the maximum evalue to be considered for parology
        -ident: the minimum percent identity to be considered for parolgy
        -scov: the minimum subject coverage to be considered for parology
    '''
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

        print('Evaluating {} for Paralogs'.format(self.g1.name))

        #Run blast
        table = self.runBlastP(self.g1)

        # write paralog data to a file
        self.writeParalogs(table)

    '''
    Runs blastp using a genome and it's associated BLAST database.
    '''
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

    '''
    Loads the output of blastp into PANDAS tables
    '''
    def loadToTable(self,results):

        # load raw data into pandas table
        table = pd.read_table(results,header=None,sep='\t',names=self.columns)

        return self.processTable(table)

    '''
    Process the table by adding the column 'scov' (subject coverage) and calculating
    the subject coverage. Also process the subject ids to not include 'ref' in
    the name.
    '''
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

    '''
    The subset function returns a subset of the table which contains all the
    entries which could be potential orthologs in the blastp results.
    '''
    def subset(self, table):
        comparisons = table
        comparisons = comparisons.loc[comparisons['evalue'] <= self.evalue]
        comparisons = comparisons.loc[comparisons['ident'] >= self.ident]
        comparisons = comparisons.loc[comparisons['scov'] >= self.scov]

        # sort
        comparisons = comparisons.sort_values(['evalue','ident', 'scov'], ascending = [True, False, False])

        return comparisons

    '''
    Writes the paralogs to the output file.
    '''
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
'''
Class Intersection finds the protiens that have orthologs in all of the genomes.
All results are outputed as the ID of a protein in the reference genome.
'''
class Intersection:

    '''
    Constructor for Intersection class.

    Parameters:
        -genomes: a list of all the genome objects
        -orthologs: a list of all the ortholog objects
    '''
    def __init__(self,genomes,orthologs,outdir):

        self.orthologs = orthologs
        self.ref = orthologs[0].g1
        self.genomes = genomes
        self.outfile = outdir + '/results/intersection.txt'

        self.orthologDict = {}
        self.refComp = []
        self.common = set()

        print('Finding Intersection of Genomes')

        self.loadGenomes(orthologs)
        self.findCommon()
        self.intersect(self.orthologs)

    '''
    Creates a 2D dictionary which allows any ortholog object to be found easily
    when trying to find the intersection of all the genomes.
    '''
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

    '''
    Creates an initial set of sequence ids
    '''
    def findCommon(self):

        #get only the ortholog objects that contain the reference genome
        self.refComp = self.orthologDict[self.ref.name]

        #create a set of the reference proteins from the first ortholog object
        orthologDF = self.refComp.values()[0].orthologTable

        #print(orthologDF)

        refList = orthologDF["sseqid"].tolist()
        self.common = set(refList)


    '''
    Find the intersection of the proteins that are found as orthologs
    '''
    def intersect(self,orthologs):

        #flag to indicate whether gene is shared
        shared = 0


        for genome in self.genomes:

            if genome.name != self.ref.name:

                ortholog = self.orthologDict[self.ref.name][genome.name]

                possibleCommon = self.getRefOrthologs(self.ref,ortholog)

                self.common.intersection_update(possibleCommon)

                if not self.common:

                    print('No Common Orthologs')

                    return None

        #Write the results to a file

        intersectFile = open(self.outfile,'w')

        intersectFile.write('\n'.join(list(self.common)))

        intersectFile.close()

    '''
    Returns a set of protein ids that are of the genome provided
    '''
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


'''
A Function to turn a row in a PANDAS dataframe to a string
'''
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

    Intersection(genomes,orthologs,outdir)
