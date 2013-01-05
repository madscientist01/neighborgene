#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import re
import subprocess

def extractRefSeqId(multifasta):
    '''
    Extract Longest sequence of FASTA record in MultiFASTA file.
    Used for the representive sequences in Clustered FASTA file.
    '''
    ids = []
    if os.path.exists(multifasta):
        f = open(multifasta)
        gi = re.compile('>gi\|(\S+)\|ref\|(\S+)\| (.*) \[(.*)\]')
        while True:
            line = f.readline()
            if not line:
                break
            if line[0:1]==">":
                match = gi.match(line)
                if match:
                    ids.append(match.group(2))
        f.close()       
    return (ids)


class usearch(object):
    '''
    USEARCH (http://www.drive5.com/usearch/) Wrapper
    Using cd-hit, cluster set of fasta file depend on sequence similarity.
    '''
    def __init__(self, **kwargs):
        self.file = kwargs.get('file')
        self.hits = []
        self.threshold = kwargs.get('threshold',0.4)
        self.clustersDic = {}
        self.clusters = []
        self.cluster = []
        self.representive = {}
        self.overwrite = kwargs.get('overwrite',False)
        self.clusteredFASTA = ''
        
    def runLocal(self):
        #
        # Using Hmmscan in locally installed Hmmer3 packages, find locations of domains in the Fasta sequence and store
        #
        filesincluster = []
        if os.path.exists(self.file):
            allfilelist = extractRefSeqId(self.file)
           
            self.clusteredFASTA = self.file+"clustered.fasta"
            usearchOutFile = self.file + ".uc"
            print 'Clustring using usearch..'.format(self.file)
            if not os.path.exists(usearchOutFile) or self.overwrite:
                    # Run local phmmer
                p = subprocess.Popen([
                    'usearch',
                    '-cluster_fast',
                    self.file,
                    '-uc',
                    usearchOutFile,
                    '-id',
                    str(self.threshold),
                    '--centroids',
                    self.clusteredFASTA
                    ], stdout=subprocess.PIPE)
                p_stdout = p.stdout.read()
            if os.path.exists(usearchOutFile):
                f = open(usearchOutFile, 'r')
                data = f.readlines()
                f.close()
                cluster=[]
                for line in data:
                    split = line.split()
                    representive = split[0]
                    clusterNo = split[1]
                    accession = split[8]
                    refseq = accession.split('|')[3]
                    if (representive == "S") and clusterNo not in self.clustersDic:
                        cluster = []
                        cluster.append(refseq)
                        self.clustersDic[clusterNo]=cluster
                        filesincluster.append(refseq)

                    elif representive =="H":
                        self.clustersDic[clusterNo].append(refseq)
                        filesincluster.append(refseq)
                self.clusters = [cluster[1] for cluster in sorted(self.clustersDic.iteritems(), key=lambda x: len(x[1]), reverse=True)]
                #
                # usearch ignore very small proteins. so they are grouped as a new cluster
                #
                deletedcluster = []
                if len(allfilelist) <> len(filesincluster):
                    for file in allfilelist:
                        if not file in filesincluster:
                            deletedcluster.append(file)
                self.clusters.append(deletedcluster)

                return True
            else:
                print 'Problem in usearch running. check the parameters'
                sys.exit()
        else:
            print '{0} file is not exist'.format(self.file)
            sys.exit()
            return False

if __name__ == '__main__':
    
    u = usearch(file='selected.fasta', threshold=0.4)
    u.runLocal()
    for i in range(len(u.clusters)):
        print i, u.clusters[i]
    print len(u.clusters)