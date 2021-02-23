# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 09:23:41 2021

@author: josh
"""

from neo4j import GraphDatabase


# making connection
uri = "bolt://localhost:7687"
driver = GraphDatabase.driver(uri, auth=("neo4j", "josh123"), encrypted = False)

# session=driver.session()

# result=session.run("MATCH (n) RETURN n")

for record in result:
    print(record)

block=m_1.tolist()
iblock=left.tolist()
def create_child(tx, block):
    tx.run("CREATE (a:child {A: $A, B:$B,C:$C, D:$D} )",
           A=block[0],B=block[1],C=block[2],D=block[3])

with driver.session() as session:
    session.write_transaction(create_child, block)   
    
    
    
def create_inode(tx, iblock):
    tx.run("MATCH (c:child)"
        "CREATE (c)-[:child_of]->(:node {L: $A, R:$B})",
           A=iblock[0],B=iblock[1])

with driver.session() as session:
    session.write_transaction(create_inode, iblock)   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
driver.close()
