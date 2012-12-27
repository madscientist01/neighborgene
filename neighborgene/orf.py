#!/usr/bin/python
# -*- coding: utf-8 -*-


class ORF(object):

    """
....Hmmer domain hits data container

....name : Name of hit. It is used for the label display in svgdrawer.py
....acc : Accession #
....desc : Description
....evalue : evalue of hit
...."""

    def __init__(self, **kwargs):
        self.name = kwargs.get('name')
        self.source = kwargs.get('source')
        self.file = kwargs.get('desc')
        self.start = int(kwargs.get('start'))
        self.end = int(kwargs.get('end'))
        self.direction = kwargs.get('direction')
        self.organism = kwargs.get('organism')
        self.link = kwargs.get('link')
        self.accession = kwargs.get('accession')
        self.description = ""
        self.exclude = False
        self.color = None
        self.label = True
        self.border = True
        self.startshow = True
        self.endshow = True
        self.gradient = False


