#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import re
import os
import subprocess

def extractRefSeq(data):
    gi = re.compile('gi\|(\S+)\|ref\|(\S+)\|')
    refseq = ''
    mat = gi.match(data)
    if mat:
        refseq = mat.group(2)
    else:
        refseq = data
    return refseq

class PhmmerSearchHit(object):

    def __init__(self, **kwargs):
        self.target = kwargs.get('target')
        self.query = kwargs.get('query')
        self.desc = kwargs.get('desc')
        self.evalue = kwargs.get('evalue')

class PhmmerSearch(object):
    '''
    Using phmmer algorithm of hmmer3, cluster set of fasta file depend on sequence similarity.
    '''
    def __init__(self, **kwargs):
        self.algorithm = kwargs.get('algorithm','phmmer')
        if not self.algorithm in ['phmmer','hmmsearch']:
            print "Algorithm is set as phmmer"
            self.algorithm = 'phmmer' 
        self.file = kwargs.get('file')
        self.db = kwargs.get('db')
        self.align = kwargs.get('align','align.sto')
        self.cutoff = kwargs.get('evalue', 1e-5)
        self.score = kwargs.get('score',0)
        self.hits = []
        self.color = kwargs.get('color')
        self.threshold = kwargs.get('cut_ga')
        self.overwrite = kwargs.get('overwrite',False)
        self.clusters = []
        self.hmmfile = self.align+".hmm"
        self.outputFile = kwargs.get('outputFile')


    def hmmbuild(self):
        if os.path.exists(self.align) or self.overwrite:
            print "generate HMM.."
            p = subprocess.Popen([
                        'hmmbuild',
                        self.hmmfile,
                        self.align,
                        ], stdout=subprocess.PIPE)
            p_stdout = p.stdout.read()
            if os.path.exists(self.hmmfile):
                return True
            else:
                print 'Problem in Hmmbuild running. check the parameters'
                sys.exit()

    def clustering(self):   
        self.runLocal()
        query = []
        subject = []
        included = []
        for hit in self.hits:
            query.append(hit.query)
            subject.append(hit.target)
            
        for i in range(len(query)):

            for cluster in self.clusters:
                if query[i] in cluster and not subject[i] in cluster and not subject[i] in included:
                    cluster.append(subject[i])
                    included.append(subject[i])
                    
            if not query[i] in included and not subject[i] in included:
                if query[i]!=subject[i]:
                    newCluster = []
                    newCluster.append(query[i])
                    newCluster.append(subject[i])
                    self.clusters.append(newCluster)
                    included.append(query[i])
                    included.append(subject[i])
                else:
                    newCluster = []
                    newCluster.append(query[i])
                    self.clusters.append(newCluster)
                    included.append(query[i])
                   
        for i in range(len(query)):
            if not query[i] in included:
                newCluster = []
                newCluster.append(query[i])
                included.append(query[i])
                self.clusters.append(newCluster)
        # n=1
        # no=0
        # for cluster in self.clusters:
        #     print "Cluster {0}".format(n)
        #     for member in cluster:
        #         print member
        #         no+=1
        #     n+=1
        # print no

    def runLocal(self):

        #
        # Using Hmmscan in locally installed Hmmer3 packages, find locations of domains in the Fasta sequence and store
        #

        if os.path.exists(self.file):
            hmmeroutFile = self.outputFile
            print 'Search {0} using {1}..'.format(self.file, self.algorithm)
            #
            # hmmscan --domtbout outputFile --cut_ga, DBfile, queryFile
            #
            # phmmer search against himself
            #
            if not os.path.exists(hmmeroutFile) or self.overwrite:
                if self.threshold == 'cut_ga':
                    # Run local phmmer
                    p = subprocess.Popen([
                        self.algorithm,
                        '-A',
                        str(self.align),
                        '-T',
                        str(self.score),
                        '--tblout',
                        hmmeroutFile,
                        '--cut_ga',
                        self.file,
                        self.db,
                        ], stdout=subprocess.PIPE)
                else:
                    p = subprocess.Popen([
                        self.algorithm,
                        '-A',
                        str(self.align),
                        '-T',
                        str(self.score),
                        '--tblout',
                        hmmeroutFile,
                        self.file,
                        self.db,
                        ], stdout=subprocess.PIPE)
                p_stdout = p.stdout.read()
            if os.path.exists(hmmeroutFile):
                f = open(hmmeroutFile, 'r')
                data = f.readlines()
                f.close()
                #
                # Parse table formatted domain datafile generated by '--domtblout'.
                #
                for line in data:
                    if line[0:8] == '# target':
                        desclocation = line.find('description')
                    if line[0:1] != '#':
                        splited = line.split()
                        if float(splited[4]) < self.cutoff:
                            target = extractRefSeq(splited[0])
                            query = extractRefSeq(splited[2])
                            evalue = splited[4]
                            desc = line[desclocation]
                            hit = PhmmerSearchHit(
                                target=target,
                                query=query,
                                desc=desc,
                                evalue=evalue,
                                )
                            self.hits.append(hit)
                return True
            else:
                print 'Problem in Hmmer3 running. check the parameters'
                sys.exit()
        else:
            print '{0} file is not exist'.format(self.file)
            sys.exit()
            return False

    
  