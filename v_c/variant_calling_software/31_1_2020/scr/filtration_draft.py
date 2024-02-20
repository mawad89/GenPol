#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 15:12:58 2020

@author: awad
"""

a=list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/illumina_polishing/hirsuta/new_10x/new_polishing/0q/secondround/0qm/1_filtered.sdi'))
n=0
b=[]
v=[]
while n<len(a):
    base=a[n]
    bas=a[n-1]
    if base.split('\t')[4]=='N':
        b.append(a[n])
        if int(bas.split('\t')[1])>=int(base.split('\t')[1])-5 and bas.split('\t')[2]=='-1':
            m=a.index(bas)
            while (int(a[m].split('\t')[1])==int(a[m-1].split('\t')[1]) or int(a[m-1].split('\t')[1])-1) and a[m].split('\t')[2]=='-1':
                v.append(a[m])
                m=m-1
    n=n+1
for base in b:
    q=a.index(base)
    del a[q]
for base in v:
    q=a.index(base)
    del a[q]
bb=b[:]
b=[]
n=0
while n<len(bb):
    base=bb[n]
    if base.split('\t')[3]=='N' and base.split('\t')[4]=='N':   
        b.append(bb[n])
    n=n+1
c=[]
n=0
while n<len(b)-1:
    base=b[n]
    bas=b[n+1]
    bst=base.split('\t')[1]
    while int(base.split('\t')[1])+1==int(bas.split('\t')[1]) and n<len(b)-2:
        n=n+1
        base=b[n]
        bas=b[n+1]
    bnd=base.split('\t')[1]
    calc=int(bnd)-int(bst)
    c.append(base.split('\t')[0]+'\t'+bst+'\t'+str(calc)+'\t'+'N'*calc+'\t'+'N'*calc+'\n')
    n=n+1
d=[]
for base in c:
    d.append(base.split('\t')[1])