#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 16:11:43 2020

@author: awad
"""
import os
def polish(sdi,fa,output=os.getcwd(),name='genome'):
    a=list(open(sdi))
    n=0
    while n<len(a):
        base=a[n]
        if base.split('\t')[3]=='N' or base.split('\t')[4]=='N':
            del a[n]
        else:
            n=n+1
    b={}
    b['-1']=[]
    b['1']=[]
    b['0']=[]
    for base in a:
        b[base.split('\t')[2]]+=[base]
    n=0 #next while loop will delated in filtration update
    while n<len(b['0'])-1:
        if int(b['0'][n].split('\t')[1])==int(b['0'][n+1].split('\t')[1]):
            del b['0'][n]
        else:
            n=n+1
    c={}
    c['-1']=[]
    c['1']=[]
    c['0']=[]
    for base in b.keys():
        d=b[base]
        d.append('.\t-1\t.\t.\t.\t.\t.\n')
        n=0
        while n<len(d)-1:
            e=int(d[n].split('\t')[1])
            f=int(d[n+1].split('\t')[1])
            if f!=e+1:
                k=d[n].split('\t')
                if base=='0':
                    k[2]='0'+str(len(k[3]))
                    k[6]='0\n'
                elif base=='1':
                    k[2]='+1'
                    k[6]='+\n'
                elif base=='-1':
                    k[2]='-1'
                    k[6]='-\n'
                k[5]=k[1]
                c[base]+=['\t'.join(k)]
                n=n+1
            else:
                g=True
                h=e
                l=d[n].split('\t')[3]
                m=d[n].split('\t')[4]
                while g==True:
                    i=d[n].split('\t')
                    j=d[n+1].split('\t')
                    if int(j[1])!=int(i[1])+1:
                        g=False
                        if base=='0':
                            i[2]='0'+str(len(l))
                            i[5]=';'.join(map(str,list(range(h,h+len(l)))))
                            i[6]='0\n'
                            #i[5]=map(str,i[5])
                            #i[5]=';'.join(i[5])
                        elif base=='-1':
                            i[2]='-'+str(len(l))
                            i[5]=';'.join(map(str,list(range(h,h+len(l)))))
                            i[6]='-\n'
                            #i[5]=map(str,i[5])
                            #i[5]=';'.join(i[5])
                        elif base=='1':
                            i[2]='+'+str(len(m))
                            i[5]=';'.join(map(str,list(range(h,h+len(m)))))
                            i[6]='+\n'
                            #i[5]=map(str,i[5])
                            #i[5]=';'.join(i[5])
                        c[base]+=[i[0]+'\t'+str(h)+'\t'+i[2]+'\t'+l+'\t'+m+'\t'+i[5]+'\t'+i[6]]
                        n=n+1
                    else:
                        if l!='.':
                            l=l+j[3]
                        else:
                            l=l
                        if m!='.':
                            m=m+j[4]
                        else:
                            m=m
                        n=n+1
    d=c['0']+c['-1']+c['1']
    d=sorted(d, key=lambda i:(int(i.split('\t')[1]),int(i.split('\t')[6].replace('0\n','1').replace('+\n','0').replace('-\n','2')),int(i.split('\t')[2].replace('0','').replace('+','').replace('-',''))))
    try:
        del a,b,base,c,e,f,g,h,i,j,k,l,m,n,sdi
    except:
        a,b,base,c,e,f,g,h,i,j,k,l,m,n,sdi='','','','','','','','','','','','','','',''
        del a,b,base,c,e,f,g,h,i,j,k,l,m,n,sdi
    a=d
    del d
    f=(list(open(fa)))
    b=''.join(f[1:]).replace('\n','')
    contig=f[0].replace('>','').split(' ')[0].split('\t')[0]
    c=['del']
    for base in b:
        c.append(base)
    prob,cc,dd=[],[],[]
    o=0
    n=0
    for base in a:
        d=base.split('\t')
        s=int(d[1])
        ro=d[2][0]
        m=int(d[2][1:])
        r=d[3].replace('\t','')
        alt=d[4].replace('\t','')
        for bas in range(n,s):
            cc.append(c[bas])
        if ro=='0' and r==''.join(c[s:s+m]):
            cc.extend(map(str,alt))
            #dd.append('origenal\t'+str(s)+'\tSNP\t'+str(m)+'\tsupposed\t'+str(o+(s))+'\trealy\t'+str(len(cc)-m)+'\n')
            n=s+m
        elif ro=='-' and r==''.join(c[s:s+m]):
            cc.extend(list(alt*m))
            n=s+m
            #dd.append('origenal\t'+str(s)+'\tdeletion\t'+str(m)+'\tsupposed\t'+str(o+(s))+'\trealy\t'+str(len(cc)-m)+'\n')
        elif ro=='+':
            cc.extend(map(str,alt))
            #dd.append('origenal\t'+str(s)+'\tinsertion\t'+str(m)+'\tsupposed\t'+str(o+(s))+'\trealy\t'+str(len(cc)-m)+'\n')
            o=o+m
            n=s
        else:
            cc.extend(c[s:s+m])
            prob.append('error\t'+str(s)+'\t'+c[s]+'\n')
    
    for base in range(n,len(c)):
        cc.append(c[base])
    seq=''.join(cc[1:]).replace('.','')
    new_fasta=open(output+'/'+name+'_polished.fa','w')
    new_fasta.write('>'+contig.replace('\n','')+'\t'+str(len(seq))+'\n'+seq+'\n')
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