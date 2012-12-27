#!/usr/bin/python
# -*- coding: utf-8 -*-

#
# Simple HTML Table Generator
#
# Written by Madscientist (http://madscientist.wordpress.com, https://github.com/madscientist01/)
#
class SVGList(object):

    def __init__(self,**kwargs):

        self.header = '<br>'
        self.footer = '<br>'
        self.svgEmbeddingTemplate = kwargs.get('svgEmbeddingTemplate')
        self.svgEmbedContent = ''

    def svgContentFill(self, svgContent):
        fill = self.svgEmbeddingTemplate.format(*svgContent)
        self.svgEmbedContent = self.svgEmbedContent + fill

