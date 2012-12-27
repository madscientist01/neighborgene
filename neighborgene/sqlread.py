#!/usr/bin/python

import pymysql
cnx = pymysql.connect(user='root', host='localhost', db='allbacteria')
cursor = cnx.cursor()
adddata = ( "SELECT protein, other FROM allbacteria"
            "WHERE protein IN (%s)")
lists = ["YP_002539552.1", "YP_001542603.1", "NP_396389.2", "NP_355941.2","NP_800868.1", "NP_800869.1", "NP_800870.1", "NP_800871.1", "NP_800872.1", "NP_800873.1", "NP_800874.1", "NP_800875.1", "NP_800876.1", "NP_800877.1", "NP_800878.1", "NP_800879.1", "NP_800880.1", "NP_800881.1", "NP_800881.1.list", "NP_800882.1", "NP_800883.1", "NP_800884.1", "NP_800885.1", "NP_800886.1", "NP_800887.1"]
in_ids = ", ".join(map(lambda x: '%s', lists))
adddata = adddata % ('%s', in_ids)
rows = cursor.execute(adddata)

for row in rows:
	(protein, other) = row
	print protein, other

cursor.close()
cnx.close()