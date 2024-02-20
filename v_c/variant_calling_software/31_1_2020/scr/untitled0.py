#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 17:13:21 2020

@author: awad
"""
'''
a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom1/hifi_new_20/genome.fa'))).split('>')
b=[a[4],a[3],a[2]]
n=0
c=[]
for base in b:
    c.append('>Chr1 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr1.fa','w')
d.writelines(c)
d.close()

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom2/hifi_merge_2/genome.fa'))).split('>')
b=[a[1],a[2],a[4],a[3]]
n=0
c=[]
for base in b:
    c.append('>Chr2 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr2.fa','w')
d.writelines(c)
d.close()

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom3/hifim/genome.fa'))).split('>')
b=[a[10]]
c=[]
n=0
for base in b:
    c.append('>Chr3 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr3.fa','w')
d.writelines(c)
d.close()

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom4/hifim/genome.fa'))).split('>')
b=[a[8],a[21]]
c=[]
n=0
for base in b:
    c.append('>Chr4 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr4.fa','w')
d.writelines(c)
d.close()


a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom5/hifi_merge/genome.fa'))).split('>')
b=[a[25],a[142]]
c=[]
n=0
for base in b:
    c.append('>Chr5 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr5.fa','w')
d.writelines(c)
d.close()

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom6/hifim/genome.fa'))).split('>')
aa=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom6/hifim/chr6_small_arm.fa'))).split('>')
b=[aa[1],a[14]]
c=[]
n=0
for base in b:
    c.append('>Chr6 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr6.fa','w')
d.writelines(c)
d.close()

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom7/hifi_merge/genome.fa'))).split('>')
b=[a[15]]
c=[]
n=0
for base in b:
    c.append('>Chr7 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr7.fa','w')
d.writelines(c)
d.close()

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom8/hifim/genome.fa'))).split('>')
b=[a[6],a[12]]
c=[]
n=0
for base in b:
    c.append('>Chr8 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr8.fa','w')
d.writelines(c)
d.close()

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom9/hifi_merge_2/genome.fa'))).split('>')
b=[a[1]]
c=[]
n=0
for base in b:
    c.append('>Chr9 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr9.fa','w')
d.writelines(c)
d.close()

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom10/hifi_merge/chr10.fa'))).split('>')
b=[a[1]]
c=[]
n=0
for base in b:
    c.append('>Chr10 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr10.fa','w')
d.writelines(c)
d.close()

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom11/hifim/genome.fa'))).split('>')
b=[a[9]]
c=[]
n=0
for base in b:
    c.append('>Chr11 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr11.fa','w')
d.writelines(c)
d.close()

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom12/hifim/genome.fa'))).split('>')
b=[a[11],a[65]]
c=[]
n=0
for base in b:
    c.append('>Chr12 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr12.fa','w')
d.writelines(c)
d.close()
a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom13/hifim/genome.fa'))).split('>')
b=[a[2]]
c=[]
n=0
for base in b:
    c.append('>Chr13 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr13.fa','w')
d.writelines(c)
d.close()

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom14/hifim/genome.fa'))).split('>')
b=[a[1]]
c=[]
n=0
for base in b:
    c.append('>Chr14 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr14.fa','w')
d.writelines(c)
d.close()

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom15/hifim/genome.fa'))).split('>')
b=[a[32]]
c=[]
n=0
for base in b:
    c.append('>Chr15 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr15.fa','w')
d.writelines(c)
d.close()

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom16/hifi_mix_3/genome.fa'))).split('>')
b=[a[3]]
c=[]
n=0
for base in b:
    c.append('>Chr16 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr16.fa','w')
d.writelines(c)
d.close()

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom17/hifi_merge/genome.fa'))).split('>')
b=[a[19]]
c=[]
n=0
for base in b:
    c.append('>Chr17 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr17.fa','w')
d.writelines(c)
d.close()

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom18/hifi_mix_2/genome.fa'))).split('>')
b=[a[46],a[5]]
c=[]
n=0
for base in b:
    c.append('>Chr18 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr18.fa','w')
d.writelines(c)
d.close()

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom19/hifi_merge_2/genome.fa'))).split('>')
b=[a[1],a[2]]
c=[]
n=0
for base in b:
    c.append('>Chr19 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr19.fa','w')
d.writelines(c)
d.close()

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom20/hifi_new_20/genome.fa'))).split('>')
b=[a[2]]
c=[]
n=0
for base in b:
    c.append('>Chr20 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr20.fa','w')
d.writelines(c)
d.close()

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom21/canu_mix_2/canu_mix_2.contigs.fasta'))).split('>')
b=[a[1]]
c=[]
n=0
for base in b:
    c.append('>Chr21 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr21.fa','w')
d.writelines(c)
d.close()

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chrom22/hifi_mix_3/genome.fa'))).split('>')
b=[a[52]]
c=[]
n=0
for base in b:
    c.append('>Chr22 Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/Chr22.fa','w')
d.writelines(c)
d.close()
'''

a=''.join(list(open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/chromx/hifim/genome.fa'))).split('>')
b=[a[1],a[27],a[2]]
c=[]
n=0
for base in b:
    c.append('>ChrX Scaffold_'+str(n+1)+'\tLen='+str(len(''.join(base.split('\n')[1:])))+'\n'+'\n'.join(base.split('\n')[1:]))
    n=n+1
c=''.join(c)
d=open('/netscratch/dep_tsiantis/grp_gan/awad/todo/Human/hifi/human_genome/ChrX.fa','w')
d.writelines(c)
d.close()




