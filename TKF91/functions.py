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
            m_block.append(msa_indel[:,blk_idx[0][i-1]:pos])
            m_block.append(msa_indel[:,pos:])
        else:
            m_block.append(msa_indel[:,blk_idx[0][i-1]:pos])
            
    return [im_block,m_block]
        
        
# test main_blocking(msa_indel)
[im_block,m_block]=main_blocking(msa_indel)


# creating internal ancestors

# for mortal blocks
# get the block
m_1=m_block[0]
n=m_1.shape[1]
# split into left and right subtree
left=m_1[0:2,:]
right=m_1[2:,:]

# finding position of ancestral characters
def m_ancestral_generator(block):
    """
    return ancestoral possibilities
    """
    anc_lst=[]
    n=block.shape[1]
    if n==1:
        anc_lst.append(np.array([1],dtype=int))
    
    else:
        c_sum=np.sum(block[:,1:],0)
        anc_idx=np.where(c_sum==2)[0]
        anc_idx
        
        
        if anc_idx.size==0:
            for i in range(2**(n-1),2**(n)):
                anc=np.array(list(format(i,'b')),dtype=int)
                anc.shape=(1,n)
                anc_lst.append(anc)
        else:
            
            for i in range(2**(n-1),2**(n)):
                anc=np.array(list(format(i,'b')),dtype=int)
                
                anc_pos=np.where(anc[1:]==1)[0]
                
                if np.array_equal(anc_pos, anc_idx):
                    anc.shape=(1,n)
                    anc_lst.append(anc)
    return anc_lst

anc_left=m_ancestral_generator(left)
anc_right=m_ancestral_generator(right)



    

    

