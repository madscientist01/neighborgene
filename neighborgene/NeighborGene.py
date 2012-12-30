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
from orf import ORF
from orfdrawer import ORFDrawer
from phmmerclust import PhmmerSearch
from colorcycler import ColorCycler
from htmltable import HTMLTable
from svgList import SVGList
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
                organism = match.group(4)
                
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
            organisms[id]=organism
            refseqs.append(id)
            return (refseqs,descriptions, organisms, fileNames)
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

def extractMultiFasta(multifasta,path,searchId,overwrite=False):
    '''
    extract given RefSeq Ids(searchId) as individual FASTA file from given multifasta file.
    Parse out refseq Id, descriptions and organism name from the FASTA header and return
    as dictionary refseqs,descriptions and organisms (Key is Refseq Id)
    '''
    capture = False
    fileNames = {}
    descriptions = {}
    organisms = {}
    refseqs = []
    buffer = []
    gi = re.compile('>gi\|(\S+)\|ref\|(\S+)\| (.*) \[(.*)\]')
    for i in searchId:
        if os.path.exists(path+i+'.fasta'):
            run = False
        else:
            run = True
    if run or overwrite:
        f = open(multifasta)
        gi = re.compile('>gi\|(\S+)\|ref\|(\S+)\| (.*) \[(.*)\]')
        # gid = ''
        desc=''
        organism = ''
        id = ''
        refseq = ''
        while True:
            line = f.readline()
            if not line:
                break
            match = gi.match(line)
            if match:
                if capture:
                    print refseq + '.fasta'
                    fw = open(path+refseq + '.fasta', 'w')
                    fw.writelines(buffer)
                    fw.close
                    fileNames[id] = path+refseq + '.fasta'
                    descriptions[id]=desc
                    organisms[id]=organism
                    refseqs.append(id)
                    capture = False
                    refseq = ''
                    buffer = []

                id = match.group(2).strip()
                if id in searchId:
                    capture = True
                    refseq = id
                    desc = match.group(3)
                    organism = match.group(4)             
            if capture:
                buffer.append(line)

        if capture:            
            print "Extracting "+refseq + '.fasta'
            fw = open(path+refseq + '.fasta', 'w')
            fw.writelines(buffer)
            fw.close
            capture = False
            refseq = ''
            buffer = []
            fileNames[id] = path+refseq + '.fasta'
    else:
        print "Files are exist."
        for i in searchId:
            fileNames[i]=path+i+'.fasta'
            f = open(fileNames[i])
            line = f.readline()
            match = gi.match(line)
            if match:
                id = match.group(2)
                descriptions[id]=match.group(3)
                organisms[id]=match.group(4)
                refseqs.append(id)             
    return (refseqs,descriptions, organisms, fileNames)

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
        for cluster in self.clusters:
            clusterName = "Cluster {num:02d}".format(num=n)
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
                for specie in self.dataList:
                    for orf in specie:
                        if member == orf.name:
                            orf.accession = clusterName
                            orf.color = clusterColor
            n+=1

        return (max(clusterNo,key=clusterNo.get))

    def findLength(self):
        '''     (FIX)
        Find start and end point of ORFs
        '''
      
        for i in range(len(self.dataList)):
           
            start = 1e+100000
            end = 0
            for orf in self.dataList[i]:
                if start > orf.start:
                    start = orf.start
                if end < orf.end:
                    end = orf.end
                organism = orf.organism
            length = end - start
            self.startDic.append(start)
            self.endDic.append(end)
            self.lengthDic.append(length)
            self.organismDic.append(organism)
        return()

    def clusterHMMscan(self):
        '''
        using hmmscan, find out Pfam Hits of the FIRST member of clusters
        '''
        representiveFASTA = []
        clusterDir=self.path+"cluster"
        if not(os.path.isdir(clusterDir)):
            print "Generating Directory {0}. All of representive cluster files will be saved there.".format(clusterDir)
            os.mkdir(clusterDir)
        
        for i in range(len(self.clusters)):
            fileList = []
            [fileList.append(self.path+member+".fasta") for member in self.clusters[i]]
            combinedFastaName = self.path+"Cluster{num:03d}.fasta".format(num=i+1)
            makeMultiFasta(fileList, combinedFastaName)
            saveName = clusterDir+"/ClusterRep{num:03d}.fasta".format(num=i+1)
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
        Main Pipeline
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
            c.execute('SELECT gff.organism, gff.start, gff.end, gff.direction, gff.protein FROM gff WHERE protein=%s',t)
            for row in c.fetchall():
                (organism, start, end, direction, id) = row
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
                              name=name)

                    orfData.append(orf) 
                    sourceList.append(source)
                    
                self.dataList.append(orfData)
                self.sourceList.append(source)
                #
                    # extract FASTA file based on the refseq id from self.fastaseq (FASTA database)
                    #
        c.close()
        conn.close()
        extractlist = self.path+'extractlist.txt'
        f=open(extractlist,'w')
        f.writelines(accessions)
        f.close()

        # (refseqs, descriptions,organisms,fileNames)=extractMultiFasta(self.fastaseq, self.path, extract,self.overwrite)
        (refseqs, descriptions,organisms,fileNames)=esl_sfetch(self.fastaseq, self.path, extractlist, self.path+tempFile)
                    

        allFiles+=fileNames.values()
        for i in range(len(self.dataList)):
            for orf in self.dataList[i]:
                orf.description = descriptions[orf.name]
                orf.organism = organisms[orf.name]
                orf.file = fileNames[orf.name]
       #
       # Concatenate all of hit as single mutifasta
       #
        # makeMultiFasta(allFiles, self.path+tempFile)

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

        standardCluster = self.clusterInfoInjection(clusters,originalSearchAccession)
        self.findLength()
        self.clusterHMMscan()

        # draw results
        orfdraw = ORFDrawer(outputSVG=self.path+self.outputFile+".svg", 
                            results=self.dataList,
                            startDic=self.startDic,
                            endDic=self.endDic,
                            lengthDic=self.lengthDic,
                            organismDic=self.organismDic,
                            linkAddress=linkFormat.format(dataFile),
                            standardCluster = standardCluster,
                            source=self.sourceList,
                            path = self.path,
                            standardDirection = '+'
                            )  
        orfdraw.drawSVG()
        (svgFileNames, svgContent)=orfdraw.drawMultiSVG()      
        #
        # Table Generation
        #
        header = ['Cluster', '#','Pfam Hits', 'Descriptions']
        table = HTMLTable(header=header)
        # Setup for the DataTable
        # Table content generations
        table.scriptcontent="""           
                $('#listtable').dataTable({
                    "sDom": "<'row'<'span8'l>r>t<'row'<'span8'i><'span8'p>>",
                     "iDisplayLength": 50,
                     "aoColumnDefs": [
                     { "sWidth": "60px", "aTargets": [ 0 ] },
                     { "sWidth": "20px", "aTargets": [ 1 ] },                
                     ]                                    
                    });
        """
        table.style ="""                  
            <style>
            table{
                font-family: "Arial",Sans-Serif;
                font-size: 12px;
                margin: 40px;
                width:1000px;
                text-align: left;
                border-collapse: collapse;  
                }
            tr.conditionalRowColor
            {

            }
                
             td.conditionalRowColor
            {
                background-color:#FFEEEE;
            }

            .scrollable {
            height: 100%;
            overflow: auto;
            }
            div.head {
                width:800px;
                font-family: Sans-Serif;
                font-size: 14px;
                border:3px solid #EEEEEE;
                border-radius: 10px;
                padding: 10px;
                align :center;
                background-color: #FFFFFF;
                }
           div.dataTables_length label {
                width: 460px;
                float: left;
                text-align: left;
            }
             
            div.dataTables_length select {
                width: 75px;
            }
             
            div.dataTables_filter label {
                float: right;
                width: 460px;
            }
             
            div.dataTables_info {
                padding-top: 8px;
            }
             
            div.dataTables_paginate {
                float: right;
                margin: 0;
            }
             
            table {
                clear: both;
            } 
            </style>
        """
        
        svgList = SVGList(svgEmbeddingTemplate='<div class="domain" id="{0}">{1}</div>')
        for i in range(len(self.dataList)):
            specieid = self.sourceList[i]
            svgList.svgContentFill([specieid, svgContent[i]])      
        svg = svgList.svgEmbedContent
        
        link = "<a href='http://www.ncbi.nlm.nih.gov/protein/{0}?report=genpept'>{1}</a>"   
        clusterlink = "<a id='{0}' name='{0}'>{0}</a>"
        annotationLink = "<a href='http://pfam.sanger.ac.uk/family/{0}'>{1}</a>"
        
        for i in range(len(clusters)):
            clust = self.clusterNames[i].replace(" ","")
            clusterName = clusterlink.format(clust)
            clusterMemberNo = len(clusters[i])
            desc = []
            for member in clusters[i]:
                desc.append(link.format(member,descriptions[member]))      
            pfam = []
            for (pfamid, pfamdesc) in self.annotations[i].items():
                pfam.append(annotationLink.format(pfamid,pfamdesc))

            tableDescriptions = ' , '.join(desc)
            annotationDescriptions = " ,".join(pfam)
            if clust !=standardCluster.replace(" ",""):
                table.tableContentFill([clusterName, clusterMemberNo, annotationDescriptions, tableDescriptions])
            else:
                table.tableContentMarkFill([clusterName, clusterMemberNo, annotationDescriptions, tableDescriptions])       
        
        svg2=orfdraw.drawClusterTable(self.clusters,self.clusterNames,self.clusterColors)     
        table.extra = '<div class="span12 scrollable">'+svg2+'</div><div class="span12 scrollable">' + svg + '</div>'  
        table.tableGenerate(self.path+dataFile)
        # run simpleHTTPServer and load html report in default Browser
        Handler = SimpleHTTPServer.SimpleHTTPRequestHandler
        httpd = SocketServer.TCPServer(("", PORT), Handler)
        print "serving at port", PORT
        webbrowser.open(linkFormat.format(self.path+dataFile))
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
