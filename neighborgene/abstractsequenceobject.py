from fetchutil import readFasta, readAccession
import sys

class AbstractSequenceObject(object):
     def __init__(self, **kwargs):

        self.file = kwargs.get('file')
        self.db = kwargs.get('db')
        fasta = readFasta(self.file)
        if fasta:
            (self.name, self.sequence) = fasta
            accession = readAccession(self.file)
            if accession:
                (self.source, self.accession) = accession
            else:
                self.source = None
                self.accession = None
            self.length = len(self.sequence)
            self.features = {}
        else:
            print "problem in {0} file.".format(self.file)
            sys.exit()
           
        