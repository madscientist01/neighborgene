#!/usr/bin/python
# -*- coding: utf-8 -*-

class HHblitsHit(object):

    '''
    HHblits hitdata container
    '''
    def __init__(self, **kwargs):
        self.name = kwargs.get('name')
        self.evalue = kwargs.get('evalue')
        self.probability = kwargs.get('probability')
        self.score = kwargs.get('score')
        self.alignedCols = kwargs.get('alignedCols')
        self.identity = kwargs.get('identity')
        self.similarity = kwargs.get('similarity')
        self.sumProbs = kwargs.get('sumProbs')
        self.querySequence = kwargs.get('querySequence')
        self.queryStart = kwargs.get('queryStart')
        self.queryEnd = kwargs.get('queryEnd')
        self.targetSequence=kwargs.get('targetSequence')
        self.targetStart = kwargs.get('targetStart')
        self.targetEnd = kwargs.get('targetEnd')
        self.confidence = kwargs.get('confidence')
        self.querySecondary = kwargs.get('querySecondary')
        self.targetDSSPSecondary = kwargs.get('targetDSSPSecondary')         
        self.tier = kwargs.get('tier',0)
        self.exclude = kwargs.get('exlude',False)
        self.color =kwargs.get('color')
        self.label = kwargs.get('label',True)
        self.labelLink = kwargs.get('labellink')
        self.border = kwargs.get('border',True)
        self.startshow = kwargs.get('startshow',True)
        self.endshow = kwargs.get('endshow',True)
        self.gradient = kwargs.get('gradient',False)
        self.description = ""
        self.description2 =""
        self.specie=""
        self.resolution = None


