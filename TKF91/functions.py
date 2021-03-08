# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 10:41:44 2021

@author: josh
"""

from Bio import SeqIO
import numpy as np
from itertools import product



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
m_1=m_block[1]
# n=m_1.shape[1]

im_1=im_block[0]
# split into left and right subtree
left=m_1[0:2,:]
right=m_1[2:,:]

left_im=im_1[0:2,:]
right_im=im_1[2:,:]



# ancestral generator for mortal block

def m_ancestral_generator(block):
    """
    return ancestoral possibilities for mortal block
    """
    anc_lst=[]
    n=block.shape[1]
    if n==1:
        anc_lst.append(np.array([1],dtype=int))
    
    else:
        c_sum=np.sum(block[:,1:],0)
        anc_idx=np.where(c_sum==2)[0]
        anc_check=np.where(c_sum==1)[0]
        # anc_idx
        
        
        if anc_idx.size==0:
            for i in range(2**(n-1),2**(n)):
                anc=np.array(list(format(i,'b')),dtype=int)
                anc.shape=(1,n)
                if len(anc_check)>0:
                    count=0
                    for pos in anc_check:
                        if anc[0][pos+1]==1 or anc[0][pos]==1:
                            count+=1
                            if count==len(anc_check):
                                anc_lst.append(anc)
                        else:
                            break
                else:
                    anc_lst.append(anc)
                    
                
        else:
            
            for i in range(2**(n-1),2**(n)):
                anc=np.array(list(format(i,'0{}b'.format(n))),dtype=int)
                
                anc_pos=np.where(anc[1:]==1)[0]
                
                if len(np.setdiff1d(anc_idx,anc_pos))==0:
                    anc.shape=(1,n)
                    if len(anc_check)>0:
                        count=0
                        for pos in anc_check:
                            if anc[0][pos+1]==1 or anc[0][pos]==1:
                                count+=1
                                # print(count,count==len(anc_check))
                                if count==len(anc_check):
                                    anc_lst.append(anc)
                            else:
                                break
                    else:
                        anc_lst.append(anc)
                    
    return anc_lst

anc_left=m_ancestral_generator(left)
anc_right=m_ancestral_generator(right)

# ancestral generator for immortal block
def im_ancestral_generator(block):
    """
    return ancestoral possibilities for immortal block
    """
    anc_lst=[]
    n=block.shape[1]
    if n==1:
        anc_lst.append(np.array([1],dtype=int))
        anc_lst.append(np.array([0],dtype=int))
    else:
        c_sum=np.sum(block[:,1:],0)
        anc_idx=np.where(c_sum==2)[0]
        
        if anc_idx.size==0:
            for i in range(2**(n)):
                anc=np.array(list(format(i,'0{}b'.format(n))),dtype=int)
                anc.shape=(1,n)
                anc_lst.append(anc)
        else:
            
            for i in range(0,2**(n)):
                anc=np.array(list(format(i,'b')),dtype=int)
                
                anc_pos=np.where(anc[1:]==1)[0]
                
                if np.array_equal(anc_pos, anc_idx):
                    anc.shape=(1,n)
                    anc_lst.append(anc)
    return anc_lst

anc_left_im=im_ancestral_generator(left_im)
anc_right_im=im_ancestral_generator(right_im)
############################################################################
# probabilities of ancestral fates
def beta(lamda, mu, t):
    beta_val=(1-np.exp((lamda-mu)*t))/(mu-(lamda*(np.exp((lamda-mu)*t))))
    return beta_val

def p_k(lamda, mu, t, k):
    b=beta(lamda, mu, t)
    prob=np.exp(-1*mu*t)*(1-(lamda*b))*((lamda*b)**(k-1))
    return prob
                        
def p_prime_k(lamda, mu, t, k):
    b=beta(lamda, mu, t)
    prob=(1-(np.exp(-mu*t))-(mu*b))*(1-(lamda*b))*((lamda*b)**(k-1))
    return prob

def p_prime_0(lamda, mu, t):
    b=beta(lamda, mu, t)
    return (mu*b)

def p_dprime_k(lamda, mu, t, k):
    b=beta(lamda, mu, t)
    prob=(1-(lamda*b))*((lamda*b)**(k-1))
    return prob
############################################################################
# subblocking
def sub_blocking(block,anc):
    """
    Splits the block into subblocks based on the presence of ancestral link
    """
    # col_sums=np.sum(msa_indel,0)
    blk_idx=list(np.where(anc==1))[0]
    
    im_block=[]
    blocks=[]
    
    if len(blk_idx)==0:
        im_block.append(block[:,:])
    
    elif len(blk_idx)==1:
        pos=blk_idx[0]
        if pos==0:
            blocks.append(block[:,:])
        else:
            imblock.append(block[:,:pos])
            blocks.append(block[:,pos:])
            
        
    else:
        for i in range(len(blk_idx)):
            pos=blk_idx[i]
            
            if i==0:
                if pos==i:
                     # im_block.append(np.zeros((4,1),dtype=int))
                     continue
                else:
                     im_block.append(block[:,:pos])
            elif i==len(blk_idx)-1:
                blocks.append(block[:,blk_idx[i-1]:pos])
                blocks.append(block[:,pos:])
            else:
                blocks.append(block[:,blk_idx[i-1]:pos])
            
    return [im_block,blocks]

# likelihood mortal block
def mblock_lK(mblock_lst,lamda, mu, t):
    
    prob=1
    
    for mblock in mblock_lst:
        
        children=np.sum(mblock,1)
            
        for i in range(2):
            k=children[i]
            print(k)
            if k==0:
                prob=prob*p_prime_0(lamda, mu, t)
                print(prob, "del")
            else:
                if mblock[i,0]==0:
                    prob=prob*p_prime_k(lamda, mu, t, k)
                    print(prob, "del and birth")
                else:
                    prob=prob*p_k(lamda, mu, t, k)
                    print(prob, "sur and birth")
                
        prob=(lamda/mu)*prob
        print(prob)
    
    return prob
        
# likelihood immortal block
def imblock_lK(imblock,lamda, mu, t):
    
    prob=1
    
    children=np.sum(imblock,1)
            
    for i in range(2):
        k=children[i]+1
        print(k)
        prob=prob*p_dprime_k(lamda, mu, t,k)
        print(prob, "immortal")
        
    prob=(1-(lamda/mu))*prob
    print(prob)
    
    return prob                
                
    
######################################################################   

# creating dictionary to store LK
internal_dict=dict()
root_dict=dict()    

# for mortal blocks
# get the block
sel_block=m_block[1]

# split into left and right subtree
left=sel_block[0:2,:]
right=sel_block[2:,:]

# ancestral generation
anc_left=m_ancestral_generator(left)
anc_right=m_ancestral_generator(right)

# finding combinations
anc_comb=list(product(anc_left, anc_right))

block_lk=0

for item in anc_comb:
    sel_lft=item[0][0]
    sel_rgt=item[1][0]
    # print(sel_lft)
    lst_anc_l=sel_lft.tolist()
    keyl=str(lst_anc_l)
    
    if keyl in internal_dict.keys():
        lblk_lk=internal_dict[keyl]
        print("dictionary")
    else:
        sub_left=sub_blocking(left,sel_lft)
        lblk_lk=mblock_lK(sub_left[1],0.4, 0.6, 0.5)
        internal_dict[keyl]=lblk_lk
    lst_anc_r=sel_rgt.tolist()  
    keyr=str(lst_anc_r)
    
    if keyr in internal_dict.keys():
        rblk_lk=internal_dict[keyr]
        print("dictionary")
    else:
        sub_right=sub_blocking(right,sel_rgt)
        rblk_lk=mblock_lK(sub_right[1],0.4, 0.6, 0.5)
        internal_dict[keyr]=rblk_lk
        
    prod=lblk_lk*rblk_lk
    
    intrnl_anc=np.array([sel_lft,sel_rgt])
    print(intrnl_anc)
    root_anc=m_ancestral_generator(intrnl_anc)
    print(root_anc)
    
    prod_1=0
    for anc in root_anc:
        print(anc)
        key2=anc.tolist()
        key2.extend(lst_anc_l)
        key2.extend(lst_anc_r)
        key2=str(key2)
        if key2 in root_dict.keys():
            lblk_lk=root_dict[key2]
            print("dictionary")
        else:
            cols=anc.shape[1]
            anc.shape=(cols,)
            sub_blocks=sub_blocking(intrnl_anc,anc)
            lblk_lk=mblock_lK(sub_blocks[1],0.4, 0.6, 0.5)
            root_dict[key2]=lblk_lk
            
        prod_1+=lblk_lk
        print(prod_1)
        
    block_lk+=(prod*prod_1)
    
    
    

