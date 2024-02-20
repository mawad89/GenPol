#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 00:00:31 2020

@author: awad
"""

import os
def polish(sdi,fa,output=os.getcwd(),name='genome'):
    a=list(open(sdi))
    f=(list(open(fa)))
    b=''.join(f[1:]).replace('\n','')
    contig=f[0]
    c=['del']
    for base in b:
        c.append(base)
    prob,cc=[],[]
    n=0
    for base in a:
        d=base.split('\t')
        s=int(d[1])
        m=int(d[2])
        r=d[3]
        alt=d[4]
        for bas in range(n,s):
            cc.append(c[bas])
        if m==0 and r==c[s]:
            cc.append(alt.replace('\t',''))
        elif m<0 and r==c[s:(m*-1)+1]:
            cc.append(alt.replace('\t',''))
        elif m>1:
            cc.append(alt.replace('\t',''))
        else:
            cc.append(c[s])
            prob.append('error\t'+str(s)+'\t'+c[s]+'\n')
        n=s+1
    for base in range(n,len(c)):
        cc.append(c[base])
    seq=''.join(cc[1:]).replace('.','')
    new_fasta=open(output+'/'+name+'_polished.fa','w')
    new_fasta.write(contig+seq+'\n')
    problematic=open(output+'/'+name+'_problematic.txt','w')
    problematic.writelines(''.join(prob))
    return()
    
import argparse
parser = argparse.ArgumentParser(prog="polishing",usage='%(prog)s -h  [options] -s <sdi_file>  -r <_genome> ',description='genome_polishing_software',)
parser.add_argument("-s",nargs=1,type=str,required=True,metavar='variant calling sdi file',dest="sdi")
parser.add_argument("-f",nargs=1,type=str,required=True,metavar='reference fastafile',dest="fa")
parser.add_argument("-o",nargs=1, type=str, default=[os.getcwd()], metavar='output path',dest="output")
parser.add_argument("-n",nargs=1, type=str, default=['genome'], metavar='output file name',dest="name")
args = parser.parse_args()
sdi=args.sdi[0]
fa=args.fa[0]
output=args.output[0]
name=args.name[0]
polish(sdi=sdi,fa=fa,output=output,name=name)