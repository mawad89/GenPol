#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 17:55:02 2020

@author: awad
"""
import os
def fasta_sep(fasta,out=os.getcwd()):
    a=list(open(fasta))
    a=''.join(a)
    a=a.split('>')
    del a[0]
    for base in a:
        b=open(out+'/'+base[2]+'.fa','w')
        b.writelines('>'+base)
        b.close()
    return('Done')
fasta='genome.fa'
out=os.getcwd()
fasta_sep(fasta,out)
