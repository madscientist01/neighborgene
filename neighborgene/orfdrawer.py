#!/usr/bin/python
# -*- coding: utf-8 -*-
import xml.etree.cElementTree as ET
import re
from operator import attrgetter
from colorcycler import blackOrWhite

def splitsentence(word):
    wordlist = re.split("(\W+)", word)
    wordcount = len(wordlist)
    splited = []
    divpoint = wordcount // 2
    if wordcount > 2:
        splited.append(' '.join(wordlist[:divpoint]))
        splited.append(' '.join(wordlist[divpoint + 1:]))
    else:
        splited.append(word)
    return splited

class ORFDrawer(object):

    '''
    Class for the Hmmer domain orfs
    '''

    def __init__(self, **kwargs):
        self.outputHTML = kwargs.get('outputHTML')
        self.outputSVG = kwargs.get('outputSVG')
        self.results = kwargs.get('results')
        self.titlemode = kwargs.get('titlemode', False)
        self.scaleFactor = kwargs.get('scaleFactor', 0.2)
        self.startDic = kwargs.get('startDic')
        self.endDic = kwargs.get('endDic')
        self.lengthDic =kwargs.get('lengthDic')
        self.organismDic = kwargs.get('organismDic')
        self.linkAddress = kwargs.get('linkAddress')
        self.standardCluster = kwargs.get('standardCluster')
        self.standardid = kwargs.get('standardid')
        self.source = kwargs.get('source')
        self.path = kwargs.get('path')
        self.standardDirection = ""
        self.xDeltaDic = []
        self.clusterTable =  'clusterTable.svg'
        return

    def saveSVG(self, filename, doc):

        '''
         Save SVG based on the ElementTreeDoc
        '''
        f = open(filename, 'w')
        f.write('<?xml version="1.0" standalone="no"?>\n')
        f.write('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"\n')
        f.write('"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n')
        f.write(ET.tostring(doc))
        f.close()
        return ()

    def setStandard(self):
        standardStart = 0
        orfs = self.results[standardStart]
        for orf in orfs:
            if orf.accession == self.standardCluster:
                standardStart = orf.start - self.startDic[standardStart]
                self.standardDirection = orf.direction
                              
        for i in range(len(self.results)):
            self.xDeltaDic.append(0)
            for orf in self.results[i]:
                if orf.accession == self.standardCluster:
                    if self.standardDirection == orf.direction:
                        orfstart = orf.start - self.startDic[i]
                        self.xDeltaDic[i]=standardStart - orfstart
                        # print orfstart, standardStart
                    else:
                        orfend = self.endDic[i] - orf.end 
                        self.xDeltaDic[i]=standardStart - orfend 
        return (max([abs(x) for x in self.xDeltaDic]))

    def drawSVG(self):
        '''
        Draw SVG based on the orf annotation
        '''
        x = 50
        y = 70
        rightMargin = 100
        fontSize = 16
        boxHeight = 20
        maxLength = max(self.lengthDic)
        maxTitleLength = len(max(self.organismDic, key=lambda x:len(x)))
        self.setStandard()
        extraLeftMargin =0
        if min(self.xDeltaDic)<0:
            extraLeftMargin =abs(min(self.xDeltaDic))
        else:
            extraLeftMargin = -abs(min(self.xDeltaDic))
        extraRightMargin = max(self.xDeltaDic)
        canvasWidth = int((maxLength+extraLeftMargin+extraRightMargin) * self.scaleFactor)
        leftMargin = int(maxTitleLength * fontSize * 0.65)+int(extraLeftMargin*self.scaleFactor)+50
        effectiveWidth = canvasWidth - leftMargin - rightMargin
        conversion = float(effectiveWidth) / float(maxLength)
        yDelta = 60
        canvasHeight = len(self.results) * yDelta + 100

        
        # doc is elementTree container for svg

        doc = ET.Element('svg', width=str(canvasWidth),
                         height=str(canvasHeight), version='1.2',
                         xmlns='http://www.w3.org/2000/svg')
        doc.attrib['xmlns:xlink'] = 'http://www.w3.org/1999/xlink'

        # Color and Gradient Assignment

        defs = self.colorDef()
        doc.append(defs)

        # Draw several domain arcorfecture in single SVG
        
        for i in range(len(self.results)):
            specie = self.results[i]
            doc = self.singleLabel(
                i,
                specie,
                doc,
                x,
                y,
                leftMargin,
                fontSize,
                conversion,
                boxHeight,
                )
            doc = self.singleSVG(
                i,
                specie,
                doc,
                x,
                y,
                leftMargin,
                fontSize,
                conversion,
                boxHeight,
                )
            y += yDelta

        self.saveSVG(self.outputSVG, doc)
        return '\n'.join(ET.tostringlist(doc))

    def drawMultiSVG(self):
        '''
        Draw multiple SVG
        '''
        
        x = 50
        y = 25
        rightMargin = 100
        fontSize = 16
        boxHeight = 20
        maxLength = max(self.lengthDic)
        maxTitleLength = len(max(self.organismDic, key=lambda x:len(x)))
        # Prepare for alignments
        self.setStandard()
        extraLeftMargin =0
        if min(self.xDeltaDic)<0:
            extraLeftMargin =abs(min(self.xDeltaDic))
        else:
            extraLeftMargin = -abs(min(self.xDeltaDic))
        extraRightMargin = max(self.xDeltaDic)
        canvasWidth = int((maxLength+extraLeftMargin+extraRightMargin) * self.scaleFactor)
        leftMargin = int(extraLeftMargin*self.scaleFactor)+50
        effectiveWidth = canvasWidth - leftMargin - rightMargin
        conversion = float(effectiveWidth) / float(maxLength)
        labelWidth = int(maxTitleLength * fontSize * 0.65)
        canvasHeight = 50
        svgFileNames = {}
        svgContent = {}
        svgContentLabel = {}
        #
        # find alignment
        #
        for i in range(len(self.results)):
            specie = self.results[i]
            doc = ET.Element('svg', width=str(canvasWidth),
                             height=str(canvasHeight), version='1.2',
                             xmlns='http://www.w3.org/2000/svg')
            doc.attrib['xmlns:xlink'] = 'http://www.w3.org/1999/xlink'
            # Color and Gradient Assignment
            defs = self.colorDef()
            doc.append(defs)
            doc = self.singleSVG(
                i,
                specie,
                doc,
                x,
                y,
                leftMargin,
                fontSize,
                conversion,
                boxHeight,
                )
            organismName = re.sub("\W", "_", self.organismDic[i])
            svgFileName = self.path+organismName + '.svg'
            self.saveSVG(svgFileName, doc)
            label = ET.Element('svg', width=str(labelWidth),
                             height=str(canvasHeight), version='1.2',
                             xmlns='http://www.w3.org/2000/svg')
            label.attrib['xmlns:xlink'] = 'http://www.w3.org/1999/xlink'
            # Color and Gradient Assignment
            label = self.singleLabel(
                i,
                specie,
                label,
                x,
                y,
                leftMargin,
                fontSize,
                conversion,
                boxHeight,
                )
            svgFileNames[i] = svgFileName
            svgContent[i] = ET.tostring(doc)
            svgContentLabel[i] = ET.tostring(label)
        return (svgFileNames, svgContent, svgContentLabel)

    def colorDef(self):
         
        defs = ET.Element('defs')
        gradientList = []
        for orfs in self.results:
            for orf in orfs:
                if orf.gradient:
                    gradientid = 'gradient_' + orf.color
                    if not gradientid in gradientList:
                        gradientList.append(gradientid)
                        gradient = ET.Element(
                            'linearGradient',
                            id=gradientid,
                            x1='0%',
                            y1='-20%',
                            x2='0%',
                            y2='120%',
                            )
                        stop1 = ET.Element('stop', offset='0%',
                                style='stop-color:worfe;stop-opacity:1'
                                )
                        stop2 = ET.Element('stop', offset='40%',
                                style='stop-color:' + orf.color
                                + ';stop-opacity:1')
                        stop3 = ET.Element('stop', offset='60%',
                                style='stop-color:' + orf.color
                                + ';stop-opacity:1')
                        stop4 = ET.Element('stop', offset='100%',
                                style='stop-color:worfe;stop-opacity:1'
                                )
                        gradient.append(stop1)
                        gradient.append(stop2)
                        gradient.append(stop3)
                        gradient.append(stop4)
                        defs.append(gradient)
        return defs

    def singleSVG(self, i, orfs, doc, x, y, leftMargin, fontSize, conversion, boxHeight):
        '''
         Draw Protein Text
         Return ElementTree Doc object as SVG container 
        '''
         # Draw Line
        # print specieid, self.xDeltaDic[specieid], int(self.xDeltaDic[specieid]*conversion)
        line = ET.Element(
            'line',
            x1=str(leftMargin+int(self.xDeltaDic[i]*conversion)),
            y1=str(y),
            x2=str(leftMargin + int((self.xDeltaDic[i]+self.lengthDic[i]) * conversion)),
            y2=str(y),
            style='stroke:rgb(200,200,200);stroke-width:4',
            )

        doc.append(line)
        invertedMode = False 
        for orf in orfs:
            if orf.accession == self.standardCluster and orf.direction != self.standardDirection:
                invertedMode = True
         # Start and End Amino Acid Number
        if invertedMode:
            startText = str(self.endDic[i])
            endText = str(self.startDic[i])
        else:
            startText = str(self.startDic[i])
            endText = str(self.endDic[i])
        
        start = ET.Element('text', x=str(int (leftMargin - len(startText)*fontSize*0.5+int(self.xDeltaDic[i]*conversion))), y=str(y), 
                            fill='black', 
                            style='font-family:Sans-Serif;font-size:13px;text-anchor:right;dominant-baseline:middle'
                           )
        start.text = startText
        doc.append(start)
        end = ET.Element('text', x=str(leftMargin + int((self.xDeltaDic[i]+self.lengthDic[i]) * conversion)), y=str(y), 
                        fill='black',
                        style='font-family:Sans-Serif;font-size:13px;text-anchor:left;dominant-baseline:middle'
                        )
        end.text = endText
        doc.append(end)

        orfs.sort(key = attrgetter('start'))
        for orf in orfs:
            if not orf.exclude:
                color = orf.color
                if orf.border:
                    if orf.accession ==self.standardCluster:
                        border = ';stroke-width:3;stroke:red'
                    else:   
                        border = ';stroke-width:1;stroke:black'
                else:
                    border = ' '
                if orf.gradient:
                    style = 'fill:url(#' + 'gradient_' + orf.color + ')' + border
                else:
                    style = 'fill:' + color + border

                 #
                 # Draw rectanglar domains
                 #

                labelfont = 13
                rectHeight = boxHeight
                if len(str(orf.accession))*labelfont*0.5 > (orf.end-orf.start)*conversion:
                    labelYPos = y - int(boxHeight*0.8)
                    labelfont = 11
                    labelColor = 'black'
                    
                else:
                    labelYPos = int(y - 0.8 * boxHeight + boxHeight)
                    labelColor = blackOrWhite(orf.color)
                    


                numberYPos = y + boxHeight
                if invertedMode:
                    if orf.direction =="+":
                        direction = "-"
                    else:
                        direction = "+"
                    orfstart = self.lengthDic[i]+self.startDic[i]-orf.end+self.xDeltaDic[i]
                    orfend = self.lengthDic[i]+self.startDic[i]-orf.start+self.xDeltaDic[i]
                else:
                    orfstart = orf.start-self.startDic[i]+self.xDeltaDic[i]
                    orfend = orf.end-self.startDic[i]+self.xDeltaDic[i]
                    direction = orf.direction
                arrowhead = 10
                arrowheight = 7
                
                if direction == "+":
                    xPos = leftMargin + orfstart * conversion
                    yPos = y 
                
                    width = int((orfend - orfstart) * conversion)
                    height = rectHeight
                    arrow = ET.Element('polygon',
                                       points='{0},{1} {2},{3} {4},{5} {6},{7} {8},{9} {10},{11} {12},{13} {14},{15} {16},{17}'.format(
                        str(xPos), str(yPos),
                        str(xPos), str(yPos-arrowheight),                        
                        str(xPos + width-arrowhead), str(yPos-arrowheight),
                        str(xPos + width-arrowhead), str(yPos- height * 0.5),
                        str(xPos + width), str(yPos),
                        str(xPos + width-arrowhead), str(yPos + height * 0.5),
                        str(xPos + width-arrowhead), str(yPos+arrowheight),
                        str(xPos), str(yPos+arrowheight),
                        str(xPos), str(yPos),
                        ), style=style)
                else:
                    xPos = leftMargin + int(orfstart * conversion)
                    yPos = y
                    width = int((orfend - orfstart) * conversion)
                    height = rectHeight
                    arrow = ET.Element('polygon',
                                       points='{0},{1} {2},{3} {4},{5} {6},{7} {8},{9} {10},{11} {12},{13} {14},{15} {16},{17}'.format(
                        str(xPos), str(yPos),
                        str(xPos + arrowhead), str(yPos-height*0.5),                        
                        str(xPos + arrowhead), str(yPos-arrowheight),
                        str(xPos + width), str(yPos- arrowheight),
                        str(xPos + width), str(yPos),
                        str(xPos + width), str(yPos + arrowheight),
                        str(xPos + arrowhead), str(yPos+arrowheight),
                        str(xPos+arrowhead), str(yPos+height * 0.5),
                        str(xPos), str(yPos),
                        ), style=style)

                if orf.link:
                    link = ET.Element('a')
                    link.attrib['xlink:href'] = orf.link
                    link.append(arrow)
                    doc.append(link)
                else:
                    doc.append(arrow)

                 #
                 # Draw Domain Label
                 #

                labelXPos = leftMargin + int((orfstart + (orf.end - orf.start)*0.5 )* conversion)
                if orf.label:  # and labelXPos > labelXEnd:

                    textLabel = ET.Element('text', x=str(labelXPos), y=str(labelYPos),
                            fill=labelColor, style='font-family:Sans-Serif;font-size:'
                            + str(labelfont) + 'px;text-anchor:middle')
                    textLabel.text = orf.accession

                    if orf.link:
                        link = ET.Element('a')
                        link.attrib['xlink:href'] = orf.link
                        link.append(textLabel)
                        doc.append(link)
                    else:
                        doc.append(textLabel)

                 #
                 # Draw start and end aa numbers of the ORF
                 #
                 # Adjust location of number, based on the ORF length
                 #

                numberFont = 9
               
                if leftMargin + orfstart*conversion + numberFont*0.6*len(str(orfstart))>\
                    leftMargin + orfend*conversion - len(str(orfend))*numberFont*0.5:
                    orf.endshow = False
                    orf.startshow=False
               
                if orf.link:
                    startEndLink = ET.Element('a')
                    startEndLink.attrib['xlink:href'] = orf.link

                if orf.startshow:
                    orfStartText = ET.Element('text', x=str(leftMargin
                            + int(orfstart* conversion)),
                            y=str(numberYPos),
                            fill='black',
                            style='font-family:Arial;font-size:'
                             + str(numberFont)
                            + 'px;text-anchor:left;dominant-baseline:top'
                            )
                    if invertedMode:
                        orfStartText.text = str(orf.end)
                    else:
                        orfStartText.text = str(orf.start)
                    if orf.link:
                        startEndLink.append(orfStartText)
                        doc.append(startEndLink)
                    else:
                        doc.append(orfStartText)

                if orf.endshow:
                    orfEndText = ET.Element('text', x=str(leftMargin
                            + int(orfend  * conversion)
                            - (len(str(orfend))-1) * numberFont),
                            y=str(numberYPos),
                            fill='black',
                            style='font-family:Arial;font-size:'
                             + str(numberFont)
                            + 'px;text-anchor:right;dominant-baseline:top'
                            )
                    if invertedMode:
                        orfEndText.text = str(orf.start)
                    else:
                        orfEndText.text = str(orf.end)
                    if orf.link:
                        startEndLink.append(orfEndText)
                        doc.append(startEndLink)
                    else:
                        doc.append(orfEndText)
        return doc


    def singleLabel(self, i, orfs, doc, x, y, leftMargin, fontSize, conversion, boxHeight):
        '''
         Draw Protein Text
         Return ElementTree Doc object as SVG container 
        '''
        id = self.source[i]
        if len(id) > 0:

            linkAddress = 'http://www.ncbi.nlm.nih.gov/nuccore/{0}'
            
            link = ET.Element('a')
            link.attrib['xlink:href'] = \
                linkAddress.format(id)

        text = ET.Element('text', x=str(fontSize),
                          y=str(y),
                          fill='black',
                          style='font-family:Sans-Serif;font-weight:bold;font-size:16px;text-anchor:left;dominant-baseline:middle'
                          )
        text.text = self.organismDic[i]
        text.attrib['id'] = id

        if len(id) > 0:
            link.append(text)
            doc.append(link)
        else:
            doc.append(text)
        return doc

    def drawClusterTable(self,clusters,clusterNames,clusterColors):
        '''
        Draw SVG based on the orf annotation
        '''
        # Find Maximum Length of ORF (as nucleotide. Not amino acid!) in each Clusters
        clusterLengthMax = {}
        for specie in self.results:
            for orf in specie:
                length = orf.end - orf.start
                if not orf.accession in clusterLengthMax:
                    clusterLengthMax[orf.accession]=length
                elif length>clusterLengthMax[orf.accession]:
                    clusterLengthMax[orf.accession]=length

        y = 70
        rightMargin = 100
        fontSize = 16
        boxHeight = 20
        boxWidth = 20
        # maxTitleLength = len(max(self.organismDic, key=lambda x:len(x)))
        maxTitleLength = 0
        leftMargin = int(maxTitleLength * fontSize * 0.65)
        canvasWidth = int(leftMargin+len(clusters)*boxWidth+rightMargin)
        # effectiveWidth = canvasWidth - leftMargin - rightMargin
        yDelta = 25
        canvasHeight = len(self.results) * yDelta + 100
        # doc is elementTree container for svg

        doc = ET.Element('svg', width=str(canvasWidth),
                         height=str(canvasHeight), version='1.2',
                         xmlns='http://www.w3.org/2000/svg')
        doc.attrib['xmlns:xlink'] = 'http://www.w3.org/1999/xlink'
        # Color and Gradient Assignment
        defs = self.colorDef()
        doc.append(defs)
        clusterlink = self.linkAddress+'#{0}'
        # Draw Top Cluster Labels
        for i in range(len(clusters)):
            x = leftMargin + i* boxWidth
            clusterColor = clusterColors[i]
            clusterName = clusterNames[i]
            link = ET.Element('a')
            link.attrib['xlink:href'] = clusterlink.format(clusterName.replace(" ", ""))
            label = ET.Element('text', x=str(x),
                          y=str(y-20),
                          fill=clusterColor,
                          style='font-family:Sans-Serif;font-weight:bold;font-size:16px;text-anchor:left;dominant-baseline:middle'
                          )
            label.text = str(i+1)
            link.append(label)
            doc.append(link)
        # Draw several domain architure in single SVG
        for j in range(len(self.results)):
            
            specie = self.results[j]
            
            for i in range(len(clusters)):
                x = leftMargin + i* boxWidth
                clusterColor = clusterColors[i]
                clusterName = clusterNames[i]
                for orf in specie:
                    if orf.accession == clusterName:
                        ratio = float(orf.end-orf.start)/float(clusterLengthMax[orf.accession])
                        if orf.border:
                            if orf.accession ==self.standardCluster:
                                border = ';stroke-width:3;stroke:red'
                            else:   
                                border = ';stroke-width:1;stroke:black'
                        else:
                            border = ''
                       
                        style = border
                        link=ET.Element('a')
                        link.attrib['xlink:href']=orf.link
                        box = ET.Element('rect', x=str(x),
                                      y=str(y),
                                      width = str(int(boxWidth*ratio)),
                                      height = str(boxHeight),
                                      fill=clusterColor,
                                      style=style
                                      )
                        link.append(box)
                        doc.append(link)
            y += yDelta

        self.saveSVG(self.path+self.clusterTable, doc)
        return '\n'.join(ET.tostringlist(doc))

    def drawClusterLabel(self,clusters,clusterNames):
        '''
        Draw SVG based on the orf annotation
        '''
        # Find Maximum Length of ORF (as nucleotide. Not amino acid!) in each Clusters
        y = 70
        fontSize = 16
        boxHeight = 20
        maxTitleLength = len(max(self.organismDic, key=lambda x:len(x)))
        leftMargin = int(maxTitleLength * fontSize * 0.65)
        canvasWidth = int(leftMargin)
        # effectiveWidth = canvasWidth - leftMargin - rightMargin
        yDelta = 25
        canvasHeight = len(self.results) * yDelta + 100
        # doc is elementTree container for svg

        doc = ET.Element('svg', width=str(canvasWidth),
                         height=str(canvasHeight), version='1.2',
                         xmlns='http://www.w3.org/2000/svg')
        doc.attrib['xmlns:xlink'] = 'http://www.w3.org/1999/xlink'
        # Color and Gradient Assignment
        # Draw several domain architure in single SVG
        
        for j in range(len(self.results)):
            
            specieid = self.source[j]
            linkAddress = 'http://www.ncbi.nlm.nih.gov/nuccore/{0}'
            link = ET.Element('a')
            link.attrib['xlink:href'] = linkAddress.format(specieid)
            text = ET.Element('text', x=str(fontSize),
                              y=str(int(y+boxHeight/2)),
                              fill='black',
                              style='font-family:Sans-Serif;font-weight:bold;font-size:16px;text-anchor:left;dominant-baseline:bottom'
                              )
            text.text = self.organismDic[j]
            text.attrib['id'] = specieid
            link.append(text)
            doc.append(link)
            y += yDelta
        
        self.saveSVG(self.path+"label"+self.clusterTable, doc)     
        return (ET.tostring(doc))


