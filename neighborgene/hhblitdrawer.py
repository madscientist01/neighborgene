#!/usr/bin/python
# -*- coding: utf-8 -*-
import xml.etree.cElementTree as ET
import re

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

class HHBlitsDrawer(object):

    '''
    Class for the Hmmer domain hits
    '''

    def __init__(self, **kwargs):
        self.outputHTML = kwargs.get('outputHTML')
        self.outputSVG = kwargs.get('outputSVG')
        self.hhblitResult = kwargs.get('hhblitResult')
        self.scaleFactor = kwargs.get('scaleFactor', 1.5)
        self.fixWidth = kwargs.get('fixWidth',False)
        self.canvasWidth = kwargs.get('canvasWidth',1200)

    def colorDef(self):

        colors = [
            'aliceblue',
            'antiquewhite',
            'aqua',
            'aquamarine',
            'azure',
            'beige',
            'bisque',
            'black',
            'blanchedalmond',
            'blue',
            'blueviolet',
            'brown',
            'burlywood',
            'cadetblue',
            'chartreuse',
            'chocolate',
            'coral',
            'cornflowerblue',
            'cornsilk',
            'crimson',
            'cyan',
            'darkblue',
            'darkcyan',
            'darkgoldenrod',
            'darkgray',
            'darkgreen',
            'darkgrey',
            'darkkhaki',
            'darkmagenta',
            'darkolivegreen',
            'darkorange',
            'darkorchid',
            'darkred',
            'darksalmon',
            'darkseagreen',
            'darkslateblue',
            'darkslategray',
            'darkslategrey',
            'darkturquoise',
            'darkviolet',
            'deeppink',
            'deepskyblue',
            'dimgray',
            'dimgrey',
            'dodgerblue',
            'firebrick',
            'floralwhite',
            'forestgreen',
            'fuchsia',
            'gainsboro',
            'ghostwhite',
            'gold',
            'goldenrod',
            'gray',
            'green',
            'greenyellow',
            'grey',
            'honeydew',
            'hotpink',
            'indianred',
            'indigo',
            'ivory',
            'khaki',
            'lavender',
            'lavenderblush',
            'lawngreen',
            'lemonchiffon',
            'lightblue',
            'lightcoral',
            'lightcyan',
            'lightgoldenrodyellow',
            'lightgray',
            'lightgreen',
            'lightgrey',
            'lightpink',
            'lightsalmon',
            'lightseagreen',
            'lightskyblue',
            'lightslategray',
            'lightslategrey',
            'lightsteelblue',
            'lightyellow',
            'lime',
            'limegreen',
            'linen',
            'magenta',
            'maroon',
            'mediumaquamarine',
            'mediumblue',
            'mediumorchid',
            'mediumpurple',
            'mediumseagreen',
            'mediumslateblue',
            'mediumspringgreen',
            'mediumturquoise',
            'mediumvioletred',
            'midnightblue',
            'mintcream',
            'mistyrose',
            'moccasin',
            'navajowhite',
            'navy',
            'oldlace',
            'olive',
            'olivedrab',
            'orange',
            'orangered',
            'orchid',
            'palegoldenrod',
            'palegreen',
            'paleturquoise',
            'palevioletred',
            'papayawhip',
            'peachpuff',
            'peru',
            'pink',
            'plum',
            'powderblue',
            'purple',
            'red',
            'rosybrown',
            'royalblue',
            'saddlebrown',
            'salmon',
            'sandybrown',
            'seagreen',
            'seashell',
            'sienna',
            'silver',
            'skyblue',
            'slateblue',
            'slategray',
            'slategrey',
            'snow',
            'springgreen',
            'steelblue',
            'tan',
            'teal',
            'thistle',
            'tomato',
            'turquoise',
            'violet',
            'wheat',
            'white',
            'whitesmoke',
            'yellow',
            'yellowgreen',
            ]
        hmmcolors = {}
        colorIndex = 10
        defs = ET.Element('defs')
        gradientList = []
       
        for hit in self.hhblitResult.features['hhblits']:
            if not hit.color:
                if hit.name in hmmcolors:
                    hit.color = hmmcolors[hit.name]
                else:
                    hmmcolors[hit.name] = colors[colorIndex]
                    hit.color = colors[colorIndex]
                    if colorIndex < len(colors) - 1:
                        colorIndex += 1
                    else:
                        colorIndex = 0
            if hit.gradient:
                gradientid = 'gradient_' + hit.color
                if not gradientid in gradientList:
                    gradientList.append(gradientid)
                    gradient = ET.Element(
                        'linearGradient',
                        id=gradientid,
                        x1='0%',
                        y1='-10%',
                        x2='0%',
                        y2='110%',
                        )
                    stop1 = ET.Element('stop', offset='0%',
                            style='stop-color:white;stop-opacity:1'
                            )
                    stop2 = ET.Element('stop', offset='40%',
                            style='stop-color:' + hit.color
                            + ';stop-opacity:1')
                    stop3 = ET.Element('stop', offset='60%',
                            style='stop-color:' + hit.color
                            + ';stop-opacity:1')
                    stop4 = ET.Element('stop', offset='100%',
                            style='stop-color:white;stop-opacity:1'
                            )
                    gradient.append(stop1)
                    gradient.append(stop2)
                    gradient.append(stop3)
                    gradient.append(stop4)
                    defs.append(gradient)
        return defs

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

    def drawSVG(self):
        '''
        Draw SVG based on the hmmer domaim
        '''
        x = 50
        y = 70
        leftMargin = 100
        rightMargin = 100
        fontSize = 16
        boxHeight = 20
        length  = self.hhblitResult.length
        if self.fixWidth:
            canvasWidth = self.canvasWidth
        else:
            canvasWidth = int(length * self.scaleFactor)
        effectiveWidth = canvasWidth - leftMargin - rightMargin
        conversion = float(effectiveWidth) / float(length)
        # If protein name width is bigger than leftMargin, wrote Label at top of proteins (self.titlemode=True)
        yDelta = 50

        canvasHeight = len(self.hhblitResult.features['hhblits']) * yDelta + 100

        # doc is elementTree container for svg

        doc = ET.Element('svg', width=str(canvasWidth),
                         height=str(canvasHeight), version='1.2',
                         xmlns='http://www.w3.org/2000/svg')
        doc.attrib['xmlns:xlink'] = 'http://www.w3.org/1999/xlink'
        defs = self.colorDef()
        doc.append(defs)
        # Draw Ruler
        line = ET.Element(
            'line',
            x1=str(leftMargin),
            y1=str(y),
            x2=str(leftMargin + int(length * conversion)),
            y2=str(y),
            style='stroke:rgb(0,0,0);stroke-width:4',
            )
        doc.append(line)
        step = 100
        stepNo = length // step
        for i in range(1,stepNo+1):
            ruler = ET.Element(
                'line',
                x1=str(leftMargin+int(i*step*conversion)),
                y1=str(y),
                x2=str(leftMargin+int(i*step*conversion)),
                y2=str(y-5),
                style='stroke:rgb(0,0,0);stroke-width:1',
                )
            label = ET.Element(
                'text',
                x=str(leftMargin+int(i*step*conversion)-int(len(str(i*step))*fontSize*0.3)),
                y=str(y-10),
                style='font-family:Sans-Serif;font-weight:bold;font-size:16px;text-anchor:center;dominant-baseline:top',
                fill='black'
                )
            label.text = str(i*step)         
            doc.append(ruler)
            doc.append(label)
        startLabel = ET.Element(
            'text',
            x=str(leftMargin-10),
            y=str(y+5),
            style='font-family:Sans-Serif;font-weight:bold;font-size:16px;text-anchor:right;dominant-baseline:top',
            fill='black'
            )
        startLabel.text = "1"         
        doc.append(startLabel)
        endLabel = ET.Element(
            'text',
            x=str(leftMargin+int(length*conversion)+5),
            y=str(y+5),
            style='font-family:Sans-Serif;font-weight:bold;font-size:16px;text-anchor:center;dominant-baseline:top',
            fill='black'
            )
        endLabel.text = str(length)         
        doc.append(endLabel)
        # Draw Secondary Structure       
        doc = self.secondaryFeaturesSVG(self.hhblitResult.features['psipred'], doc, 'psipred', x, y+15, leftMargin, fontSize, conversion, int(boxHeight/3),0,False)
        doc = self.secondaryFeaturesSVG(self.hhblitResult.features['hhblits'], doc, '', x, y+30, leftMargin, fontSize, conversion, boxHeight,30,True)       
        self.saveSVG(self.outputSVG, doc)

        return '\n'.join(ET.tostringlist(doc))

    def secondaryFeaturesSVG(self, hits, doc, label, x, y, leftMargin, fontSize, conversion, boxHeight, yDelta,rounded):
        '''
         Draw Secondary Structure
         Return ElementTree Doc object as SVG container 
        '''
        featureLabel = ET.Element(
            'text',
            x=str(leftMargin-int(len(label)*fontSize*0.5)),
            y=str(y+5),
            style='font-family:Sans-Serif;font-weight:bold;font-size:'+str(fontSize)+'px;text-anchor:right;dominant-baseline:top',
            fill='black'
            )
        featureLabel.text = label        
        doc.append(featureLabel)

        for hit in hits:
            if not hit.exclude:
                if hit.color:
                    color = hit.color
                else:
                    color = 'blue'
                if hit.border:
                    border = ';stroke-width:1;stroke:black'
                else:
                    border = ''
                if hit.gradient:
                    style = 'fill:url(#' + 'gradient_' + hit.color \
                        + ')' + border
                else:
                    style = 'fill:' + color + border
                xPos = leftMargin + hit.queryStart * conversion
                yPos = y  
                width = int((hit.queryEnd - hit.queryStart) * conversion)
                if rounded:
                    rx = int(boxHeight/2)
                    ry = int(boxHeight/2)
                else:
                    rx = 0
                    ry = 0

                rect = ET.Element(
                        'rect',
                        x=str(xPos),
                        y=str(yPos),
                        rx=str(rx),
                        ry=str(ry),
                        width=str(width),
                        height=str(boxHeight),
                        style=style,
                        )
                doc.append(rect)

                if hit.startshow:
                    xPos = int(leftMargin+hit.queryStart * conversion - int(len(str(hit.queryStart))*fontSize*0.6))
                    yPos = int(y + boxHeight*0.7)
                    startLabel = ET.Element(
                        'text',
                        x=str(xPos),
                        y=str(yPos),
                        style='font-family:Sans-Serif;font-weight:bold;font-size:'+str(fontSize)+'px;text-anchor:left;dominant-baseline:bottom',
                        fill='black'
                        )
                    startLabel.text = str(hit.queryStart)
                    doc.append(startLabel)

                if hit.endshow:
                    xPos = int(leftMargin+hit.queryEnd * conversion+5)
                    yPos = int(y + boxHeight*0.7)
                    endLabel = ET.Element(
                        'text',
                        x=str(xPos),
                        y=str(yPos),
                        style='font-family:Sans-Serif;font-weight:bold;font-size:'+str(fontSize)+'px;text-anchor:left;dominant-baseline:bottom',
                        fill='black'
                        )
                    endLabel.text = str(hit.queryEnd)
                    doc.append(endLabel)
              
                if hit.label: 
                    name = hit.name
                    if len(name) * fontSize * 0.6  < (hit.queryEnd- hit.queryStart) * conversion:
                        xPos = leftMargin + int((hit.queryStart + (hit.queryEnd - hit.queryStart) * 0.5) * conversion)
                        align = 'middle'
                    else:
                        xPos = int(leftMargin + hit.queryEnd*conversion+10+len(str(hit.queryEnd))*fontSize*0.6)
                        align = 'left'
                    yPos = int(y + boxHeight*0.7)
                    textLabel = ET.Element('text',
                            x=str(xPos), y=str(yPos),
                            fill='black',
                            style='font-family:Sans-Serif;font-size:'
                             + str(fontSize)
                            + 'px;text-anchor:'+align)
                    
                    textLabel.text = name

                    if hit.labellink:
                        link = ET.Element('a')
                        link.attrib['xlink:href'] = hit.labellink
                        link.append(textLabel)
                        doc.append(link)
                    else:
                        doc.append(textLabel)


                y+=yDelta
        return doc

    def sequencesDraw(self):

        sequenceDetailFileNames = {}
        sequenceDetailSVGContent = {}
        for hmmer in self.hhblitResult:
            defs = self.colorDef()
            (svgFileName, svgContent) = self.singleSequenceDraw(hmmer,
                    defs)
            sequenceDetailFileNames[hmmer.name] = svgFileName
            sequenceDetailSVGContent[hmmer.name] = \
                ET.tostring(svgContent)

        return (sequenceDetailFileNames, sequenceDetailSVGContent)
   

    def singleAlignDraw(self, hit):
        '''
            Draw Single Sequence and Domain, Secondary Structure.
        '''
        yStart = 30
        y = yStart
        leftMargin = 150
        fontSize = 14
        columnWidth = 60
        canvasWidth = 800
        yDelta = 50
        yDiff = 15
        conversionFactor = 0.7
         # and labelXPos > labelXEnd:
        name = hit.name
        svgfilename = name+".svg"
        length = len(hit.querySequence)
        canvasHeight = int((length / float(columnWidth)
                           + 1.5) * yDelta)

        # doc is elementTree container for svg
        doc = ET.Element('svg', width=str(canvasWidth),
                         height=str(canvasHeight), version='1.2',
                         xmlns='http://www.w3.org/2000/svg')
        doc.attrib['xmlns:xlink'] = 'http://www.w3.org/1999/xlink'
        cursor = 0
        queryStart = hit.queryStart
        targetStart = hit.targetStart

        while cursor < length:
            (doc,queryStart)=self.sequenceDraw(queryStart, hit.querySequence, hit.targetSequence, columnWidth, cursor, length, fontSize, doc, y, conversionFactor, "Query", leftMargin)
            (doc,targetStart)=self.sequenceDraw(targetStart, hit.targetSequence, hit.querySequence, columnWidth, cursor, length, fontSize, doc, y+yDiff, conversionFactor, "Target", leftMargin)
            doc=self.secondaryDraw(hit.querySecondary, columnWidth, cursor, length, fontSize, doc, y-yDiff, conversionFactor, "SS psipred", leftMargin)
            cursor += columnWidth
            y += yDelta
            queryStart+=1
            targetStart+=1

        name = hit.querySequence[cursor:length]
        aminoacids = ET.Element('text', x=str(leftMargin), y=str(y),
                                fill='black',
                                style='font-family:Courier;font-size:'
                                + str(fontSize)
                                + 'px;text-anchor:left;dominant-baseline:middle'
                                )
        aminoacids.text = name
        doc.append(aminoacids)

        self.saveSVG(svgfilename, doc)
        return (svgfilename, doc)

    def sequenceDraw(self, sequenceStart, sequence, comparisonSequence, columnWidth, cursor, length, fontSize, doc, y, conversionFactor, name, leftMargin):
        
        
        sequenceName= ET.Element('text', x=str(int(leftMargin*0.7)),
                             y=str(y), fill='black',
                             style='font-family:Arial;font-size:'
                             + str(fontSize)
                             + 'px;text-anchor:end;dominant-baseline:middle'
                             )
        sequenceName.text = str(name)
        doc.append(sequenceName)

        numbers = ET.Element('text', x=str(leftMargin - fontSize),
                             y=str(y), fill='black',
                             style='font-family:Courier;font-size:'
                             + str(fontSize)
                             + 'px;text-anchor:end;dominant-baseline:middle'
                             )
        numbers.text = str(cursor + sequenceStart)
        doc.append(numbers)

        text = sequence[cursor : cursor + columnWidth]
        compareText = comparisonSequence[cursor : cursor + columnWidth]
        width = fontSize * 0.7
        height = fontSize *1.1
        realWidth = len(text.replace("-",""))
        for i in range(len(text)):
            aa = text[i]
            caa = compareText[i]
            color = self.aminoToColor(aa)
            compcolor = self.aminoToColor(caa)
            xPos = str(int(leftMargin + fontSize * i * conversionFactor))
                
            if (color == compcolor):
                style = 'fill:'+color+';fill-opacity:0.5'
                rect = ET.Element(
                    'rect',
                    x=str(xPos),
                    y=str(int(y-fontSize/2)),
                    width=str(width),
                    height=str(height),
                    style=style,
                    )
                doc.append(rect)

            aminoacids = ET.Element('text', x=xPos, y=str(y),
                fill='black',
                style='font-family:Courier;font-size:'
                + str(fontSize)
                + 'px;text-anchor:left;dominant-baseline:middle'
                )
            aminoacids.text = aa
            doc.append(aminoacids)

        if cursor + columnWidth < length:

            endnumbers = ET.Element('text', x=str(leftMargin
                    + (columnWidth + 5) * fontSize
                    * conversionFactor), y=str(y), fill='black',
                    style='font-family:Courier;font-size:'
                    + str(fontSize)
                    + 'px;text-anchor:end;dominant-baseline:middle')
            endnumbers.text = str(cursor + realWidth+sequenceStart)
            doc.append(endnumbers)

        return(doc,sequenceStart+realWidth-columnWidth)

    def secondaryDraw(self, sequence, columnWidth, cursor, length, fontSize, doc, y, conversionFactor, name, leftMargin):
               
        sequenceName= ET.Element('text', x=str(int(leftMargin*0.7)),
                             y=str(y), fill='black',
                             style='font-family:Arial;font-size:'
                             + str(fontSize)
                             + 'px;text-anchor:end;dominant-baseline:middle'
                             )
        sequenceName.text = str(name)
        doc.append(sequenceName)
        text = sequence[cursor : cursor + columnWidth]
        width = fontSize * 0.7
        height = fontSize *0.5
        for i in range(len(text)):
            aa = text[i]
            if aa in ['H', 'h']:
                color = 'lawngreen'
            elif aa in ['E','e']:
                color = 'red'
            if aa in ['H','h', 'E','e']:
                xPos = str(int(leftMargin + fontSize * i * conversionFactor))      
                style = 'fill:'+color
                rect = ET.Element(
                    'rect',
                    x=str(xPos),
                    y=str(int(y-fontSize/2)),
                    width=str(width),
                    height=str(height),
                    style=style,
                    )
                doc.append(rect)
        return(doc)

    def aminoToColor(self, aa):
        if aa in ['L','V','I','F','W','M','Y']:
            color = 'lawngreen'
        if aa in ['D','E','N']:    
            color = 'red'
        if aa in ['K','R','H','Q']:
            color = 'blue'
        if aa in ['P']:
            color = 'pink'
        if aa in ['G','-','X']:
            color = 'white'    
        if aa in ['S','T','C','A']:
            color = 'yellow'

        return color

   


