# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 10:41:44 2021

@author: josh
"""

from Bio import SeqIO
import numpy as np



# reading and storing the MSA
msa=[]
with open(r"C:\Daten\TKF91\script\TKF91\Input\msa.fasta", "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        msa.append(str(record.seq))
        
# converting to 1's and 0's
orig="ACGT-"
new="11110"

change=str.maketrans(orig,new)
msa_indel=[list(item.translate(change)) for item in msa]

# converting to a numpy array
msa_indel=np.array(msa_indel,dtype=int)
        
# blocking the msa
def main_blocking(msa_indel):
    """
    Splits the MSA into blocks based on the presence of ancestral link
    """
    col_sums=np.sum(msa_indel,0)
    blk_idx=list(np.where(col_sums>=3))
    
    im_block=[]
    m_block=[]
    
    for i in range(len(blk_idx[0])):
        pos=blk_idx[0][i]
        
        if i==0:
            if pos==i:
                 im_block.append(np.zeros((4,1),dtype=int))
            else:
                 im_block.append(msa_indel[:,:pos])
        elif i==len(blk_idx[0])-1:
            m_block.append(msa_indel[:,pos:])
        else:
            m_block.append(msa_indel[:,blk_idx[0][i-1]:pos])
            
    return [im_block,m_block]
        
        
# test main_blocking(msa_indel)
[im_block,m_block]=main_blocking(msa_indel)
    
    

