#!/usr/bin/python
# -*- coding: utf-8 -*-
import sqlite3
import os
from phmmerclust import PhmmerSearch


if __name__ == '__main__':

    db = "vibrio.db"
    upperbound =10000
    lowerbound = 10000
    phmmer = PhmmerSearch(file='VopC.fasta',db='vibrio.fasta', align='VopC.sto', evalue=1e-50,overwrite=False)
    phmmer.runLocal()
    for hit in phmmer.hits:
        print hit.target, hit.query, hit.desc, hit.evalue
    phmmer.hmmbuild()
    hmmsearch = PhmmerSearch(file=phmmer.hmmfile, db='vibrio.fasta', score=300, algorithm='hmmsearch', overwrite=False)
    hmmsearch.runLocal()
    if os.path.exists(db):
        conn = sqlite3.connect(db)
        c=conn.cursor()
    for hit in hmmsearch.hits:
        # print hit.target, hit.query, hit.desc, hit.evalue
        t=[hit.target]
        for row in c.execute('SELECT * FROM gff WHERE protein=?',t):
            (organism, start, end, direction, id) = row
            searchStart = start - lowerbound
            searchEnd = end + upperbound
            f = [searchStart, searchEnd,organism]
            for newrow in c.execute('SELECT * FROM gff WHERE start BETWEEN ? AND ? AND organism=?',f):
                (organism, start, end, direction, id) = row
                print organism, start, end, direction, id