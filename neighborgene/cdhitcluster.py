#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import re
import os
import subprocess

def giRefDictionary(multifasta):
    '''
    Make dictionary with gi id as key, refseq id as value. Used for the parsing of CD-Hit Results
    '''
    giRefDic = {}
    if os.path.exists(multifasta):
        f = open(multifasta)
        gi = re.compile('>gi\|(\S+)\|ref\|(\S+)\| (.*) \[(.*)\]')
        while True:
            line = f.readline()
            if not line:
                break
            match = gi.match(line)
            if match:
                if not match.group(1) in giRefDic:    
                    giRefDic[match.group(1)]=match.group(2)
    else:
        print "File is not exist."
        return(None)             
    return (giRefDic)

class CDHitSearch(object):
    '''
    CD-HIT (http://bioinformatics.org/cd-hit/) Wrapper
    Using cd-hit, cluster set of fasta file depend on sequence similarity.
    '''
    def __init__(self, **kwargs):
        self.file = kwargs.get('file')
        self.word = kwargs.get('word', 2)
        self.hits = []
        self.color = kwargs.get('color')
        self.threshold = kwargs.get('threshold',0.4)
        self.clusters = []
        self.giRefDic = {}
        self.overwrite = kwargs.get('overwrite',False)
        self.clusteredFASTA = ''
        giRefDic = giRefDictionary(self.file)
        if giRefDic :
            self.giRefDic = giRefDic 

    def runLocal(self):
        #
        # Using Hmmscan in locally installed Hmmer3 packages, find locations of domains in the Fasta sequence and store
        #
        if os.path.exists(self.file):
            self.clusteredFASTA = self.file+"clustered.fasta"
            cdhitoutFile = self.clusteredFASTA + ".clstr"
            print 'Clustring using cd-hitr..'.format(self.file)
            if not os.path.exists(cdhitoutFile) or self.overwrite:
                    # Run local phmmer
                p = subprocess.Popen([
                    'cd-hit',
                    '-i',
                    self.file,
                    '-o',
                    self.clusteredFASTA,
                    '-c',
                    str(self.threshold),
                    '-n',
                    str(self.word)
                    ], stdout=subprocess.PIPE)
                p_stdout = p.stdout.read()
            headerRegex = re.compile("^>(Cluster \d+)")
            contentRegex = re.compile("^(\d+)\s+(\d+)aa, >gi\|(\d+)\|")     
            if os.path.exists(cdhitoutFile):
                f = open(cdhitoutFile, 'r')
                data = f.readlines()
                f.close()
                cluster=[]
                clusterName = ""
                for line in data:
                    match = headerRegex.match(line)
                    if match :
                        if len(clusterName)>0:
                            self.clusters.append(cluster)
                            cluster=[]
                        clusterName= match.group(1)
                    else:
                        match2 = contentRegex.match(line)
                        if match2:
                            # print match2.group(1), match2.group(2), match2.group(3)
                            cluster.append(self.giRefDic[match2.group(3)])
                self.clusters.append(cluster)
                return True
            else:
                print 'Problem in CD-Hit running. check the parameters'
                sys.exit()
        else:
            print '{0} file is not exist'.format(self.file)
            sys.exit()
            return False

if __name__ == '__main__':
    
    cdhit = CDHitSearch(file='selected.fasta', threshold=0.4)
    cdhit.runLocal()
    for i in range(len(cdhit.clusters)):
        print i, cdhit.clusters[i]