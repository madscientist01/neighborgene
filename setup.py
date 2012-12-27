#!/usr/bin/python

from setuptools import setup

setup(name='neighborgene',
      version='0.1',
      description='NeighborGene : Find Neighborhood gene',
      author='Suk Namgoong',
      author_email='suk.namgoong@gmail.com',
      url='https://github.com/madscientist01/neighborgene',
      platforms='any',
      install_requires=['jinja2','pymysql'],
      packages=['neighborgene'],
      scripts=['neighborgene/NeighborGene.py','neighborgene/hhblits.py']
      )
