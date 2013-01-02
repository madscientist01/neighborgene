#!/usr/bin/python
# -*- coding: utf-8 -*-
# NeighborGene
# Scripts for the examination of neighborhood genes of bacterial genome
# Require locally installed Hmmer 3.0, CD-HIT, MySQL, pymysql
#
import re
import os
import pymysql
import argparse
import sys
import jinja2
import operator
from orf import ORF
from orfdrawer import ORFDrawer
from phmmerclust import PhmmerSearch
from colorcycler import ColorCycler
from hmmscan import Hmmer
from cdhitcluster import CDHitSearch
from collections import defaultdict
import subprocess
import SimpleHTTPServer
import SocketServer
import webbrowser

def list_has_duplicate_items(L):
    return len(L) > len(set(L))

def get_duplicate_items(L):
    D = defaultdict(int)
    for k in L:
        D[k] += 1
    return [k for (k, v) in D.items() if v > 1]

def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        raise argparse.ArgumentError("{0} does not exist".format(x))
    return x

def makeMultiFasta(fileList,newFileName):
    '''
    Concatenate Multiple FASTA file into one multiFasta File.
    fileList  should be List
    '''
    buffer = []
    for file in fileList:
        if os.path.exists(file):
            f = open(file,'r')
            buffer = buffer + f.readlines()
            f.close()
    if len(buffer)>0:
        f=open(newFileName,'w')
        f.writelines(buffer)
        f.close()
        return(True)
    else:
        return(False)

def esl_sfetch(multifasta, path, extractlist, extractedFile):
    '''
    Extract multifasta using easel-sfetch utility in hmmer with fasta headers
    '''
    if os.path.exists(multifasta):
        indexfilename = multifasta+".ssi"
        #
        #if index for esl_sfetch is not founded, generated it.
        #
        if not os.path.exists(indexfilename):
             p = subprocess.Popen([
                                    'esl-sfetch',
                                    '--index',
                                    multifasta,
                                    ], stdout=subprocess.PIPE)
             p_stdout = p.stdout.read()
             successString = "SSI index written to file "+ indexfilename
             try:
                 p_stdout.find(successString)
             except ValueError:
                 print "problem in indexing."
                 sys.exit()
             else:
                 pass
             finally:
                 pass
        #
        # Run esl-sfetch
        #
        p = subprocess.Popen([
                            'esl-sfetch',
                            '-f',
                            multifasta,
                            extractlist
                            ], stdout=subprocess.PIPE)
        p_stdout = p.stdout.read()
        #
        # Mutifasta will be extracted as STDOUT.
        # Save it
        #
        f = open(extractedFile,'w')
        f.write(p_stdout)
        f.close()
        #
        # Then split file as seperate FASTAs.
        # 
        lines = p_stdout.split('\n')
        capture = False
        fileNames = {}
        descriptions = {}
        organisms = {}
        refseqs = []
        buffer = []
        gi = re.compile('>gi\|(\S+)\|ref\|(\S+)\| (.*) \[(.*)\]')
        # gid = ''
        desc=''
        organism = ''
        id = ''
        refseq = ''
        capture = False
        for line in lines:
            match = gi.match(line)
            if match:
                if capture:
                    print "Extracting "+refseq + '.fasta'
                    fw = open(path+refseq + '.fasta', 'w')
                    fw.writelines(buffer)
                    fw.close
                    capture = False
                    refseq = ''
                    buffer = []
                    fileNames[id] = path+refseq + '.fasta'
                    descriptions[id]=desc
                    organisms[id]=organism
                    refseqs.append(id)

                id = match.group(2).strip()
                capture = True
                refseq = id
                desc = match.group(3)
                
            if capture:
                buffer.append(line+'\n')

        if capture:
            print "Extracting "+refseq + '.fasta'
            fw = open(path+refseq + '.fasta', 'w')
            fw.writelines(buffer)
            fw.close
            capture = False
            refseq = ''
            buffer = []
            fileNames[id] = path+refseq + '.fasta'
            descriptions[id]=desc
            refseqs.append(id)
            return (refseqs,descriptions,fileNames)
    else:
        print "{0} is not found. stop.".format(multifasta)
        sys.exit()


def extractRefSeq(data):
    '''
    Extract Refseq id from FASTA header
    '''
    gi = re.compile('gi\|(\S+)\|ref\|(\S+)\|')
    refseq = ''
    mat = gi.match(data)
    if mat:
        refseq = mat.group(2)
    else:
        refseq = data
    return refseq

def extractRepresentive(multifasta,saveName,overwrite=False):
    '''
    Extract Longest sequence of FASTA record in MultiFASTA file.
    Used for the representive sequences in Clustered FASTA file.
    '''
    capture = False
    length = {}
    buffers = {}
    buffer = []
    if os.path.exists(multifasta) or overwrite:
        f = open(multifasta)
        gi = re.compile('>gi\|(\S+)\|ref\|(\S+)\| (.*) \[(.*)\]')
        refseq = ''
        while True:
            line = f.readline()
            if not line:
                break
            match = gi.match(line)
            if match:
                if capture:
                    buffers[refseq]=buffer 
                    length[refseq]=len(''.join(buffer))
                    capture = False
                    refseq = ''
                    buffer = []

                capture = True
                refseq = match.group(2)    
            if capture:
                buffer.append(line)
        if capture:
            buffers[refseq]=buffer 
            length[refseq]=len(''.join(buffer))
        longest = max(length, key=length.get)
        fw = open(saveName, 'w')
        fw.writelines(buffers[longest])
        fw.close
    return ()


class NeighborGene(object):
    '''
    Process NeighborHood Genes and Perform Clustering
    '''
    def __init__(self,**kwargs):
        self.score = kwargs.get('score',300)
        self.annotationdb = kwargs.get('annotationdb')
        self.queryFasta = kwargs.get('queryfasta')
        self.overwrite = kwargs.get('overwrite',False)
        self.fastaseq = kwargs.get('fastaseq')
        self.outputFile = kwargs.get('outputfile')
        self.linkFormat = kwargs.get('linkFormat')
        self.clusteringMode = kwargs.get('clusteringmode','cdhit')
        self.folder = kwargs.get('folder')
        if self.folder:
            basename = os.path.basename(self.queryFasta)
            self.path = os.path.splitext(basename)[0]
        else:
            self.path = ""
        if not self.clusteringMode in ['cdhit', 'phmmer']:
            self.clusteringMode = 'cdhit'
        self.clusters = []
        self.dataList = []
        self.clusterNames = []
        self.clusterColors = []
        self.sourceList = []
        self.clusterDic = {}
        self.startDic = []
        self.endDic= []
        self.lengthDic = []
        self.organismDic = []
        self.annotations = []
        
    def clusterInfoInjection(self,clusters,originalSearchAccession):
        '''
        Afer clustering using phmmer, cluster informations are injected in this method.
        '''
        clusterNo ={}

        n=1
        self.clusters=clusters
        color = ColorCycler(initColor=30) 
        if len(self.clusters)<10:
            formatText = "Cluster {num:01d}"
        elif len(self.clusters)<100:
            formatText = "Cluster {num:02d}"
        elif len(self.clusters)<1000:
            formatText = "Cluster {num:03d}"
        else:
            formatText = "Cluster {num:04d}"

        for cluster in self.clusters:
            clusterName = formatText.format(num=n)
            clusterColor = color.getColor(clusterName)
            self.clusterNames.append(clusterName)
            self.clusterColors.append(clusterColor)
            for member in cluster:
                self.clusterDic[member] = n-1
                if member in originalSearchAccession:
                    if clusterName in clusterNo:
                        clusterNo[clusterName]+=1
                    else:
                        clusterNo[clusterName]=1 
                for (organismName, organism, start, end, length, orfs) in self.dataList:
                    for orf in orfs:
                        if member == orf.name:
                            orf.accession = clusterName
                            orf.color = clusterColor
            n+=1

        return (max(clusterNo,key=clusterNo.get))

    def clusterHMMscan(self):
        '''
        using hmmscan, find out Pfam Hits of the FIRST member of clusters
        '''
        representiveFASTA = []
        clusterDir=self.path+"cluster"
        if not(os.path.isdir(clusterDir)):
            print "Generating Directory {0}. All of representive cluster files will be saved there.".format(clusterDir)
            os.mkdir(clusterDir)
        
        if len(self.clusters)<10:
            formatText = "{num:01d}"
        elif len(self.clusters)<100:
            formatText = "{num:02d}"
        elif len(self.clusters)<1000:
            formatText = "{num:03d}"
        else:
            formatText = "{num:04d}"

        for i in range(len(self.clusters)):
            fileList = []
            [fileList.append(self.path+member+".fasta") for member in self.clusters[i]]
            combinedFastaName = self.path+"Cluster"+formatText.format(num=i+1)+".fasta"                                                            
            makeMultiFasta(fileList, combinedFastaName)
            saveName = clusterDir+"/ClusterRep"+formatText.format(num=i+1)+".fasta"
            extractRepresentive(combinedFastaName,saveName,overwrite=False)
            representiveFASTA.append(saveName)
        
        representedFastaName = self.path+"representive.fasta"
        makeMultiFasta(representiveFASTA, representedFastaName)
        hmmer = Hmmer(file=representedFastaName,
                      db="Pfam-A.hmm",
                      evalue=0.01,
                      threshold="No",
                      overwrite=False,
                      excludeRedundancy=False)
        if hmmer.runLocal():
            for i in range(len(self.clusters)):
                annotation = {}
                self.annotations.append(annotation)
            for hit in hmmer.features['domain']:
                if not hit.acc in self.annotations[self.clusterDic[extractRefSeq(hit.query)]]:
                    self.annotations[self.clusterDic[extractRefSeq(hit.query)]][hit.acc]=hit.desc.strip()
        return()

    def run(self):
        '''
        Main Pipeline of NeighborGene.py
        1. Using given FASTA file (self.queryFASTA), run Phmmer against sequence db(self.fastaseq) to extract candidate list
        2. Generate HMM with initial hit from Phmmer (hmmbuild in HMMER3.0) and run hmmSearch
        3. Using range query of mySQL, extract ORF located between upperBound (bp) and lowerBound location
        4. Using FASTA record fetched from mySQL, extract portions of FASTA sequences using esl-sFetch.
        5. Extracted sequence will be clustered using Phmmer (sensitive to low sequence homology, but extremely slow) or CD-HIT 
        6. Representive sequence from the cluster is annotated with Phmmer search against Pfam Database 
        7. Table of cluster and cluster diagram and ORF distributions will be drawn as SVG and saved as HTML
        '''
        db = self.annotationdb
        upperBound =10000
        lowerBound = 10000
        tempFile = "selected.fasta"
        memberLink = "http://www.ncbi.nlm.nih.gov/protein/{0}?report=genpept"
        PORT = 8000
        # linkFormat = 'http://localhost:'+str(PORT)+'/{0}'
        linkFormat =  self.linkFormat 
        dataFile = self.outputFile      
        allFiles=[] 
        if len(self.path)>0 and not(os.path.isdir(self.path)):
            print "Generating Directory {0}. All of generated files will be saved there.".format(self.path)
            os.mkdir(self.path)

        if len(self.path)>0 and self.path[len(self.path)-1]!="/":
            self.path = self.path+"/"
        # generate query hmm using iterative phmmer - hmmsearch.
        phmmer = PhmmerSearch(file=self.queryFasta,
                              db=self.fastaseq,
                              align=self.path+self.queryFasta+'.sto',
                              outputFile=self.path+self.queryFasta+'.hmmer',
                              score=self.score,
                              overwrite=self.overwrite)
        phmmer.runLocal()
        if len(phmmer.hits)==0:
	    print "No hit found. current score threshold is {0}. Consider decrease threshold no.".format(self.score)
	    sys.exit()
        for hit in phmmer.hits:
            print hit.target, hit.query, hit.desc, hit.evalue
        # generated align.sto file is converted as hmm using hmmbuild 
        phmmer.hmmbuild()

        # hmmsearch using generated hmm
        hmmSearch = PhmmerSearch(file=phmmer.hmmfile,
                                 db=self.fastaseq,
                                 score=self.score,
                                 algorithm='hmmsearch',
                                 outputFile=self.path+self.queryFasta+'.hmmer',
                                 path = self.path,
                                 overwrite=self.overwrite)
        hmmSearch.runLocal()
        originalSearchAccession = []
        for hit in hmmSearch.hits:
            originalSearchAccession.append(hit.target)
        # using sqlite db, fetch sequence near initial hmmer hit.
        print "Search neighborhood genes.."
        
        conn = pymysql.connect(user='root', host='localhost', db=db)
        c=conn.cursor()
        conn.text_factory = str
        extract = []
        accessions = []
        sourceList = []

        for hit in hmmSearch.hits:
            t=[hit.target]
            #
            # First query exact phmmer match (seed) from gff table
            #
            c.execute('SELECT gff.organism, gff.start, gff.end, gff.direction, gff.protein,organism.organismName FROM gff, organism WHERE gff.organism=organism.organism and protein=%s',t)
            for row in c.fetchall():
                (organism, start, end, direction, id,organismName) = row
                # search ORF located between lowerbound (bp) < gene of intereste < upperbound(bp)
                searchStart = start - lowerBound
                searchEnd = end + upperBound
                f = [searchStart, searchEnd,organism]
                orfData = []
                #
                # SQL query between regions and specific organism
                #
                query = "SELECT gff.organism, gff.start, gff.end, gff.direction, gff.protein, accession.accession FROM gff,accession WHERE gff.start BETWEEN %s AND %s AND gff.organism=%s AND gff.protein=accession.protein" 
                
                c.execute(query,f)
                for newrow in c.fetchall():    
                    (source, start, end, direction, name, accession) = newrow
                    
                    extract.append(name)
                    accessions.append(accession+"\n")
                    link = memberLink.format(name)
                    orf = ORF(source=source,
                              start=start,
                              end=end,
                              direction=direction,
                              link=link,
                              name=name,
                              organism = organismName
                              )

                    orfData.append(orf) 
                    sourceList.append(source)
                
                start = 1e+100000
                end = 0
                for orf in orfData:
                    if start > orf.start:
                        start = orf.start
                    if end < orf.end:
                        end = orf.end

                length = end - start
                
                self.dataList.append((organismName, organism, start, end, length, orfData))

                self.sourceList.append(source)
                self.startDic.append(start)
                self.endDic.append(end)
                self.lengthDic.append(length)
                self.organismDic.append(organismName)
        #
        # extract FASTA file based on the refseq id from self.fastaseq (FASTA database) using esl-sfetch
        #

        c.close()
        conn.close()
        extractlist = self.path+'extractlist.txt'
        f=open(extractlist,'w')
        f.writelines(accessions)
        f.close()
        (refseqs, descriptions,fileNames)=esl_sfetch(self.fastaseq, self.path, extractlist, self.path+tempFile)
                    
        allFiles+=fileNames.values()
        for i in range(len(self.dataList)):
            orfs = self.dataList[i][5]
            for orf in orfs:
                orf.description = descriptions[orf.name]
                orf.file = fileNames[orf.name]
       
        #
        # Clustering with phmmer or CD-Hit
        #

        if self.clusteringMode =="phmmer":
            phmmer = PhmmerSearch(file=self.path+tempFile,
                                  db=self.path+tempFile,
                                  evalue=1e-5,
                                  outputFile=self.path+tempFile+'.hmmer',
                                  overwrite=self.overwrite)
            phmmer.clustering()
            clusters=phmmer.clusters
        else:
            cdhit =CDHitSearch(file=self.path+tempFile,
                               thershold=0.4, 
                               overwrite=self.overwrite)
            cdhit.runLocal()
            clusters = cdhit.clusters
        #
        # Inject clustering info into self.dataList
        #
        standardCluster = self.clusterInfoInjection(clusters,originalSearchAccession)
        self.clusterHMMscan()

        # draw results
        # Sort self.dataList based on the name of organism
        #
        self.dataList.sort(key=operator.itemgetter(0))
        orfdraw = ORFDrawer(outputSVG=self.path+self.outputFile+".svg", 
                            outputHTML=self.path+self.outputFile,
                            results=self.dataList,
                            startDic=self.startDic,
                            endDic=self.endDic,
                            lengthDic=self.lengthDic,
                            organismDic=self.organismDic,
                            linkAddress=linkFormat.format(self.path+dataFile),
                            standardCluster = standardCluster,
                            source=self.sourceList,
                            path = self.path,
                            standardDirection = '+'
                            )  
        # orfdraw.drawSVG()
        (svgFileNames, svgContent, svgLabels)=orfdraw.drawMultiSVG()      
        #
        # Table Generation using Jinja2
        #  
        svgList = []

        for i in range(len(self.dataList)):
            organism = self.dataList[i][1]
            svg = {"id":organism,
                   "content":svgContent[i],
                   "label":svgLabels[i]
                 }
            svgList.append(svg)

        clusterTableSVG=orfdraw.drawClusterTable(self.clusters,self.clusterNames,self.clusterColors)     
        clusterTableLabel=orfdraw.drawClusterLabel(self.clusters,self.clusterNames)
        link = "<a href='http://www.ncbi.nlm.nih.gov/protein/{0}?report=genpept'>{1}</a>"   
        clusterlink = "<a id='{0}' name='{0}'>{0}</a>"
        annotationLink = "<a href='http://pfam.sanger.ac.uk/family/{0}'>{1}</a>"
        jinja_environment = \
        jinja2.Environment(loader=jinja2.FileSystemLoader(os.path.dirname(__file__)))
        template = jinja_environment.get_template('neighborgenetable.html')

        # header for the hits information table
        headers = ['Cluster', '#','Pfam Hits', 'Descriptions']
        
        # html rendering and save into HTML file.
        clusterInfo = []
        clusterWithoutAnnotation = []
        for i in range(len(clusters)):

            clust = self.clusterNames[i].replace(" ","")
            clusterName = clusterlink.format(clust)
            clusterMemberNo = len(clusters[i])
            desc = []
            d = []
            others=[]
            no=1
            #
            # Clean up Pfam annotation and Descriptions
            #
            for member in clusters[i]:
                if descriptions[member] in d:
                    others.append(link.format(member,int(no)))
                    no+=1
                else:
                    d.append(descriptions[member])
                    desc.append(link.format(member,descriptions[member]))      
            pfam = []

            for (pfamid, pfamdesc) in self.annotations[i].items():
                pfam.append(annotationLink.format(pfamid,pfamdesc))

            tableDescriptions = ' , '.join(desc)
            if len(others)>0:
                tableDescriptions = tableDescriptions+", others("+",".join(others)+")"
            annotationDescriptions = " ,".join(pfam)
            if len(annotationDescriptions)==0:
                clusterWithoutAnnotation.append(clust)

            cluster =  {
                        "clusterName":clusterName,
                        "clust":clust,
                        "clusterMemberNo":clusterMemberNo,
                        "annotationDescriptions":annotationDescriptions,
                        "tableDescriptions":tableDescriptions
                        }
            clusterInfo.append(cluster)
        
        for line in clusterWithoutAnnotation:
            print line
        # time = datetime.datetime.fromtimestamp(os.path.getmtime(hhblits.hhrfile))
        # print os.path.basename(hhblits.outputAlignFile)

        params = {
                    "filename":self.fastaseq,
                    "searchdb":self.annotationdb,
                    "clustering":self.clusteringMode,
                    "standardCluster":standardCluster.replace(' ','')
        }
        
        #
        # Render HTML Report and save it.
        #

        t = template.render(params = params, headers = headers, clusters=clusterInfo, svgList = svgList, clusterTableSVG = clusterTableSVG, clusterTableLabel = clusterTableLabel)
        f = open(self.path + self.outputFile, 'w')
        f.write(t)
        f.close()
      
        # run simpleHTTPServer and load html report in default Browser
        Handler = SimpleHTTPServer.SimpleHTTPRequestHandler
        httpd = SocketServer.TCPServer(("", PORT), Handler)
        print "serving at port", PORT
        webbrowser.open(linkFormat.format(self.path+self.outputFile))
        httpd.serve_forever()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-q',
        '--query_sequence',
        dest='file',
        help='file to process',
        required=True
        )
    parser.add_argument(
        '-d',
        '--database_sequence',
        action='store',
        dest='dbseq',
        help='database sequence',
        required=True
        )
    parser.add_argument(
        '-f',
        action='store_false',
        dest='folder',
        default=True,
        help='save result in the current folder. as default, it will generate new folder and save data into it.',
        )
    parser.add_argument(
        '-a',
        '--annotation_database',
        action='store',
        dest='annodb',
        default='allbacteria',
        help='annotation db',
        )
    parser.add_argument(
        '-s',
        '--score_cutoff',
        action='store',
        dest='score',
        default=100,
        type=float,
        help='HMMER score cutoff (Higher score means high stringency)',
        )
    parser.add_argument(
        '-o',
        '--output_file',
        action='store',
        dest='outputfile',
        default='output.html',
        help='Output HTML file',
        )
    parser.add_argument(
        '-l',
        '--link',
        action='store',
        dest='link',
        default='http://localhost:8000/{0}',
        help='output file link address',
        )
    parser.add_argument(
        '-w',
        '--overwrite',
        action='store_true',
        dest='overwrite',
        default=False,
        help='ignore current running results'
    )
    parser.add_argument(
        '-p',
        '--phmmer_clustering',
        action='store_true',
        dest='phmmerclustering',
        default=False,
        help='Use phmmer clustering (sensitive but much slower than default CD-HIT clustering)',
    )
    results = parser.parse_args() 
    if results.phmmerclustering:
        method = 'phmmer'
    else:
        method = 'cdhit'   
    neighbor = NeighborGene(score=results.score,
                            overwrite=results.overwrite, 
                            queryfasta=results.file,
                            fastaseq=results.dbseq,
                            clusteringmode=method, 
                            annotationdb = results.annodb,
                            outputfile=results.outputfile,
                            linkFormat=results.link,
                            folder = results.folder)
    neighbor.run()
