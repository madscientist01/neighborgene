#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# HHblits wrapper
# MadScientist http://madscientist.wordpress.com
#
# The MIT License
#
# Copyright (c) 2012 Suk Namgoong (suk.namgoong@gmail.com)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

from hhblitshit import HHblitsHit
from hhblitdrawer import HHBlitsDrawer
from abstractsequenceobject import AbstractSequenceObject
import subprocess
import os
import datetime
import sys
import argparse
import jinja2
import SimpleHTTPServer
import SocketServer
import webbrowser
import re
import glob

class HHblits(AbstractSequenceObject):

    '''
    HHblits wrapper class
    '''

    def __init__(self, **kwargs):

        super(HHblits, self).__init__(**kwargs)
        self.path = kwargs.get('path', '')
        self.cutoff = kwargs.get('evalue', 0.0000000001)
        self.threshold = kwargs.get('threshold')
        self.overwrite = kwargs.get('overwrite')
        self.dbDictionary = kwargs.get('dbDictionary')
        self.db = kwargs.get('db')
        self.psipredLocation = kwargs.get('psipredLocation')
        self.psipredDataLocation = kwargs.get('psipredDataLocation')
        self.cpu = kwargs.get('cpu', 4)
        self.folder = kwargs.get('folder')
        if self.folder:
            basename = os.path.basename(self.file)
            self.path = os.path.splitext(basename)[0] + '_hh'
        else:
            self.path = ''
        self.psipred = ''
        self.psipredConf = ''
        self.iteration = kwargs.get('iteration`')
        hits = []
        self.features['hhblits'] = hits
        self.tier = {}
        self.tier[0] = 'hhblits'
        secondary = []
        self.features['psipred'] = secondary
        self.tier[1] = 'psipred'
        self.svg = ''
        self.hhrfile=''
        self.outputAlignFile=''
        self.outputAlignFileFASTA=''
        
    def addHit(
        self,
        start,
        state,
        end,
        currentstate,
        ):

        if currentstate != 'C':
            hit = HHblitsHit(name=currentstate, queryStart=start,
                             queryEnd=end)
            hit.label = False
            hit.border = False
            hit.startshow = False
            hit.endshow = False
            hit.gradient = False
            if currentstate == 'H':
                color = 'lawngreen'
            else:
                color = 'red'
            hit.color = color
            self.features['psipred'].append(hit)

    def hitName(self, name):
        db = self.db
        if db == 'uniprot':
            matched = re.match('\S+\s([A-Za-z0-9\-\+:, ]+)OX=', name)
            if matched:
                outname = matched.group(1)
                extra = {}
            else:
                outname = name[:10]  # and labelXPos > labelXEnd:
                extra = {}
        if db == 'pdb':
            matched = \
                re.match('(\S{4})_(\S) ([A-Za-z0-9_,\'/\+:\.\-\(\) ]+); ([A-Za-z0-9,_/\'\.\(\)\+:\- ]+);(.*)\{(.*)\}'
                         , name)
            if matched:
                outname = matched.group(1).upper()
                chain = matched.group(2)
                description = matched.group(3)
                description2 = matched.group(4)
                matched2 = re.search('([0-9.]+)A', matched.group(5))
                if matched2:
                    resolution = matched2.group(1)
                else:
                    resolution = ''
                specie = matched.group(6)

                extra = {
                    'accession': outname,
                    'description': description,
                    'description2': description2,
                    'chain': chain,
                    'specie': specie,
                    'resolution': resolution,
                    }
            else:
                outname = name[:10]
                extra = {}
        return (outname, extra)

    def runLocal(self):
        '''
        Using Hmmscan in locally installed Hmmer3 packages, find locations of domains in the Fasta sequence and store
        '''

        if os.path.exists(self.dbDictionary[self.db] + '.cs219'):

            # Generate folder for data save

            if len(self.path) > 0 and not os.path.isdir(self.path):
                print 'Generating Directory {0}. All of generated files will be saved there.'.format(self.path)
                os.mkdir(self.path)

            if len(self.path) > 0 and self.path[len(self.path) - 1] \
                != '/':
                self.path = self.path + '/'

            self.outputAlignFile = self.path + self.file + '.a3m'
            self.outputAlignFileFASTA = self.path + self.file \
                + '.aligned.fasta'
            basename = os.path.basename(self.path + self.file)
            self.hhrfile = self.path + os.path.splitext(basename)[0] + '.hhr'

            if not os.path.exists(self.hhrfile) or self.overwrite:
                print 'Running HHblits of {0} using local HHblits'.format(self.file)
                print 'First pass of HHblits'

                #
                # In the first run of hhblits, search Uniprot20 database and generate 'seed' alignment
                # supplemented with psipred secondary structure predictions. (-addss options)
                #

                p = subprocess.Popen([
                    'hhblits',
                    '-d',
                    self.dbDictionary['uniprot'],
                    '-addss',
                    '-psipred',
                    self.psipredLocation,
                    '-psipred_data',
                    self.psipredDataLocation,
                    '-i',
                    self.file,
                    '-oa3m',
                    self.path + 'temp.a3m',
                    '-n',
                    '1'
                    ], stdout=subprocess.PIPE)

                p_stdout = p.stdout.read()
                print 'Second pass of HHblits'  
                p = subprocess.Popen([
                    'hhblits',
                    '-d',
                    self.dbDictionary['uniprot'],
                    '-addss',
                    '-psipred',
                    self.psipredLocation,
                    '-psipred_data',
                    self.psipredDataLocation,
                    '-i',
                    self.path+'temp.a3m',
                    '-oa3m',
                    self.path +'temp2.a3m',
                    '-n',
                    str(self.iteration),
                    '-norealign',
                    '-maxfilt',
                    '30000'
                    ], stdout=subprocess.PIPE)
                p_stdout = p.stdout.read()
                print 'Third pass of HHblits'
                #
                # choice of db (PDB, pfam, uniprot20, nr20, scop) will be searched with generated a3m alignment
                # from the previous HHblit step.
                #

                p = subprocess.Popen([
                    'hhblits',
                    '-o',
                    self.hhrfile,
                    '-d',
                    self.dbDictionary[self.db],
                    '-addss',
                    '-psipred',
                    self.psipredLocation,
                    '-psipred_data',
                    self.psipredDataLocation,
                    '-i',
                    self.path + 'temp2.a3m',
                    '-oa3m',
                    self.outputAlignFile,
                    '-n',
                    str(self.iteration),
                    '-e',
                    str(self.cutoff),
                    '-E',
                    str(self.cutoff)

                    ], stdout=subprocess.PIPE)
                p_stdout = p.stdout.read()
              
            # a3m alignment was converted as fasta using reformat.pl scripts in
            # HHSUITE packages
            #

            if not os.path.exists(self.outputAlignFile) or self.overwrite:
                p = subprocess.Popen(['reformat.pl', 'a3m', 'fas',
                        self.outputAlignFile, self.outputAlignFileFASTA],
                        stdout=subprocess.PIPE)
                p_stdout = p.stdout.read()

            #
            # Parsing hhr file and add informations into self.feature['hhblits']
            #

            if os.path.exists(self.hhrfile):
                f = open(self.hhrfile, 'r')

                filecontents = f.readlines()
                f.close()
                startRegex = re.compile("^No (\d+)")
                capture = False
                buffers = {}
                buffer = []

                for line in filecontents:
                    matchRegex = startRegex.match(line)
                    if matchRegex:
                        id = matchRegex.group(1)
                        if capture:
                            buffers[int(id)] = buffer
                            buffer = []
                        else:
                            capture = True
                    if capture:
                        buffer.append(line)

                #
                # Regex for the Parsing
                #

                nameRegex = re.compile('>(.*)$')
                paramRegex = \
                    re.compile('^Probab=([0-9\.]+)\s+E-value=([0-9e\.\+\-]+)\s+Score=([0-9\.]+)\s+Aligned_cols=([0-9\.]+)\s+Identities=([0-9\.]+)%\s+Similarity=([0-9\.]+)\s+Sum_probs=([0-9\.]+)'
                               )
                queryRegex = \
                    re.compile("^Q \S+\s+(\d+)\s+([A-Z\-]+)\s+(\d+)")
                targetRegex = \
                    re.compile("^T \S+\s+(\d+)\s+([A-Z\-]+)\s+(\d+)")
                confidenceRegex = re.compile('^Confidence([0-9 ]+)$')
                ssRegex = re.compile("Q ss_pred\s+([cCHhEe-]+)")
                dsspRegex = re.compile("T ss_dssp\s+([HBEGITSC-]+)")

                for id in sorted(buffers.iterkeys()):
                    buffer = buffers[id]
                    querySequence = ''
                    queryStart = 999999999
                    queryEnd = 0
                    targetSequence = ''
                    targetStart = 999999999
                    targetEnd = 0
                    confidence = ''
                    ssSequence = ''
                    dsspSequence = ''

                    for line in buffer:
                        namematch = nameRegex.match(line)
                        if namematch:
                            name = namematch.group(1)
                        match = paramRegex.match(line)
                        if match:
                            probability = match.group(1)
                            evalue = match.group(2)
                            score = match.group(3)
                            alignedCols = match.group(4)
                            identity = match.group(5)
                            similarity = match.group(6)
                            sumProbs = match.group(7)
                        match = queryRegex.match(line)
                        if match and line[2:11] != 'Consensus':
                            if queryStart > int(match.group(1)):
                                queryStart = int(match.group(1))
                            querySequence = querySequence \
                                + match.group(2)
                            queryEnd = match.group(3)

                        match = targetRegex.match(line)
                        if match and line[2:11] != 'Consensus':
                            if targetStart > int(match.group(1)):
                                targetStart = int(match.group(1))
                            targetSequence = targetSequence \
                                + match.group(2)
                            targetEnd = match.group(3)

                        match = ssRegex.match(line)
                        if match:
                            ssSequence = ssSequence + match.group(1)

                        match = dsspRegex.match(line)
                        if match:
                            dsspSequence = dsspSequence + match.group(1)

                        match = confidenceRegex.match(line)
                        if match:
                            confidence = confidence + match.group(1)[12:92]

                    if len(dsspSequence) > 0:
                        dsspSequence = re.sub('[ST]', 'C', dsspSequence)
                        dsspSequence = re.sub('[GI]', 'H', dsspSequence)
                        dsspSequence = re.sub('[B]', 'E', dsspSequence)

                    #
                    # Based on the e-value, color of bar in the digram will be changed.
                    #

                    probabilityfloat = float(probability)
                    if probabilityfloat == 100:
                        color = 'crimson'
                    elif probabilityfloat > 95:
                        color = 'red'
                    elif probabilityfloat > 90:
                        color = 'orange'
                    elif probabilityfloat > 80:
                        color = 'yellow'
                    elif probabilityfloat > 70:
                        color = 'lawngreen'
                    elif probabilityfloat > 50:
                        color = 'green'
                    elif probabilityfloat > 30:
                        color = 'blue'
                    elif probabilityfloat > 20:
                        color = 'violet'
                    else:
                        color = 'grey'
                    (processedName, extra) = self.hitName(name)
                    hhblits = HHblitsHit(
                        name=processedName,
                        probability=float(probability),
                        evalue=float(evalue),
                        score=float(score),
                        alignedCols=int(alignedCols),
                        identity=float(identity),
                        similarity=float(similarity),
                        sumProbs=float(sumProbs),
                        querySequence=querySequence,
                        queryStart=int(queryStart),
                        queryEnd=int(queryEnd),
                        querySecondary=ssSequence,
                        targetSequence=targetSequence,
                        targetStart=int(targetStart),
                        targetEnd=int(targetEnd),
                        confidence=confidence,
                        color=color,
                        gradient=True,
                        targetDSSPSecondary=dsspSequence,
                        )
                    if 'description' in extra:
                        hhblits.description = extra['description']
                    if 'description2' in extra:
                        hhblits.description2 = extra['description2']
                    if 'specie' in extra:
                        hhblits.specie = extra['specie']
                    if 'resolution' in extra:
                        hhblits.resolution = extra['resolution']
                    if self.db == 'pdb':
                        hhblits.labelLink = \
                            'http://www.rcsb.org/pdb/explore/explore.do?structureId={0}'.format(processedName)
                    self.features['hhblits'].append(hhblits)

            if os.path.exists(self.outputAlignFile):

                #
                # Read Secondary Structure Prediction of Query FASTA. a3m files contains
                # all of sequence informations, but we just need Secondary Structure Predictions
                #

                state = {'C': 'Coil', 'H': 'Helix', 'E': 'Sheet'}
                f = open(self.outputAlignFile)
                alignFile = f.readlines()
                f.close()
                psipred = \
                    '>ss_pred PSIPRED predicted secondary structure'
                psipredconf = '>ss_conf PSIPRED confidence values'
                if re.match(psipred, alignFile[0]):
                    self.psipred = alignFile[1]
                if re.match(psipredconf, alignFile[2]):
                    self.psipredConf = alignFile[3]
                currentstate = self.psipred[0:1]
                start = 1
                end = 1
                for i in range(1, len(self.psipred)):
                    if currentstate != self.psipred[i:i + 1]:
                        self.addHit(start, state, end, currentstate)
                        start = i + 1
                    currentstate = self.psipred[i:i + 1]
                    end = i + 1
                self.addHit(start, state, end, currentstate)
        else:
            print 'HHblits db is not found'.format(self.db)
            sys.exit()
            return False


def main(results):
    '''
    main worker method for the pipeline.
    '''
    files = results.files
    if len(files)==0:
        files = glob.glob("*.fasta")

    for singleFile in files:
        if os.path.exists(singleFile):
            hhblits = HHblits(  # define location of HHSUITE DB. Based on the location of actual DB, these should be changed.
                file=singleFile,
                db=results.db,
                dbDictionary={'pdb': '/Users/suknamgoongold/hhsuite/db/pdb70_06Oct12',
                              'uniprot': '/Users/suknamgoongold/hhsuite/db/uniprot20_2012_03',
                              'pfamA': '/Users/suknamgoongold/hhsuite/db/pfamA_v26.0_06Dec11'
                              },
                psipredLocation=results.psipredLocation,
                psipredDataLocation=results.psipredDataLocation,
                evalue=results.evalue,
                iteration=results.iteration,
                overwrite=results.overwrite,
                folder=results.folder,
                )
            hhblits.runLocal()

            draw = HHBlitsDrawer(hhblitResult=hhblits, outputSVG='hhblits.svg',
                                 fixWidth=True, canvasWidth=1200)
            hitmap = draw.drawSVG()
            for hit in hhblits.features['hhblits']:
                (svgfile, hit.svg) = draw.singleAlignDraw(hit)

            #
            # Table Generation using jinja2
            #

            jinja_environment = \
                jinja2.Environment(loader=jinja2.FileSystemLoader(os.path.dirname(__file__)))
            template = jinja_environment.get_template('table.html')

            # header for the hits information table

            headers = [
                'PDB ID',
                'Resolution(A)',
                'Description',
                'e-value',
                'Score',
                'Probability',
                'Aligned Residues',
                'Identity(%)',
                ]

            # html rendering and save into HTML file.
            (name,ext) = os.path.splitext(singleFile)
            outFileName = name+".html"
            time = datetime.datetime.fromtimestamp(os.path.getmtime(hhblits.hhrfile))
            print os.path.basename(hhblits.outputAlignFile)
            t = template.render(results= results, time = time, headers=headers, 
                                hits=hhblits.features['hhblits'], hitmap=hitmap, 
                                hhrfile = os.path.basename(hhblits.hhrfile), 
                                alignment=os.path.basename(hhblits.outputAlignFileFASTA))
            f = open(hhblits.path + outFileName, 'w')
            f.write(t)
            f.close()
        else:
            "{0} is not exist. Skipped.".format(singleFile)

    PORT=8000
    Handler = SimpleHTTPServer.SimpleHTTPRequestHandler
    httpd = SocketServer.TCPServer(("", PORT), Handler)
    print "serving at port", PORT
    webbrowser.open("http://localhost:8000/"+hhblits.path + outFileName)
    httpd.serve_forever()



if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-q', '--query', dest='files', default=[], nargs='+',
                        help='file to search')
    parser.add_argument('-d', '--db', dest='db', default='uniprot',
                        help='HHblits DB')
    parser.add_argument('-f', action='store_false', dest='folder',
                        default=True,
                        help='save result in the current folder. as default, it will generate new folder and save data into it.'
                        )
    parser.add_argument('-psipred', dest='psipredLocation',
                        default='/Users/Shared/blast-2.25/ncbi-blast-2.2.25+/bin'
                        , help='Location of psipred')
    parser.add_argument('-psipred_data', dest='psipredDataLocation',
                        default='/Users/Shared/blast-2.25/ncbi-blast-2.2.25+/bin/data'
                        , help='Location of psipred data')
    parser.add_argument('-e', '--evalue_cutoff', dest='evalue',
                        default=0.0000000001, help='E-value cutoff')
    parser.add_argument('-i', '--iterations', dest='iteration',
                        default=2, help='Number of Iteration')
    parser.add_argument(
        '-w',
        '--overwrite',
        action='store_true',
        dest='overwrite',
        default=False,
        help='ignore current running results',
        )
    parser.add_argument('-p', '--path', dest='path', default='',
                        help='path')
    results = parser.parse_args()
    main(results)
    