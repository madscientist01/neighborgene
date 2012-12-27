#!/usr/bin/python
# -*- coding: utf-8 -*-
import jinja2
import os 

jinja_environment = jinja2.Environment(loader=jinja2.FileSystemLoader(os.path.dirname(__file__)))

template = jinja_environment.get_template('table.html')

headers=['PDB ID', 'Description', 'e-value','score','identity']

t = template.render(headers=headers)

f = open("out.html",'w')
f.write(t)
f.close()