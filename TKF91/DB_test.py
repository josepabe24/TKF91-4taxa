# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 09:23:41 2021

@author: josh
"""

from pyArango.connection import *

# creating connection
conn=Connection(username='root', password='josh123')

# creating database
db=conn.createDatabase(name="tree")
db

# creating a collection
anc_collection=db.createCollection(name="ancestors")

# creating document

doc1=anc_collection.createDocument()
doc1["A"]=left

doc2=anc_collection.createDocument()
doc2["node1"]=anc_left

doc1._key="A"
doc1.save()
doc2._key="node1"
doc2.save()

db.createCollection(name="ChildOf", className = 'Edges')

AQL="""INSERT { _from: "ancestors/A", _to: "ancestors/node1" } INTO ChildOf Return NEW"""
qry_result=db.AQLQuery(AQL)


anc_query = """
// First find the start node, i.e., A
FOR c IN ancestors
    FILTER c._key== "A"
    // Then start a Graph traversal from that start node
    FOR v IN 1..1 OUTBOUND c ChildOf
    RETURN v[v._key]
"""

query_result = db.AQLQuery(anc_query, rawResults=True)
for doc in  query_result:
    print(doc)
    print()