#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 16:59:16 2019

@author: awad
"""
#back to not that to understand
import os
from rc import rc
def variant_caller(bed_file,ref_genome,draft,output_path=os.getcwd(),output_name='var_call_awad',chr_spliter=' ',contig_spliter=' ',bse=1):
    a=list(open(bed_file))
    d=list(open(ref_genome))
    d=''.join(d)
    d=d.split('>')
    del d[0]
    e=list(open(draft))
    e=''.join(e)
    e=e.split('>')
    del e[0]
    s=0
    z={}
    y={}
    x={}
    for base in a:
        b=base.split('\t')
        chr_name=b[0]
        chr_pos=int(b[1])
        contig_name=b[3]
        quality=b[4]
        strand=b[5]
        cigar=b[6]
        contig_pos=0
        c,n,m=[],0,0
        while n<len(cigar):
            base=cigar[n]
            if base == 'D' or base == 'H' or base == 'I' or base == 'M' or base == 'N' or base == 'P' or base == 'S' or base == 'x' or base == '=':
                c.append(cigar[m:n+1])
                m=n+1
                n=n+1
            else:
                n=n+1
        for base in d:
            if chr_name == base.split('\n')[0].split(chr_spliter)[0]:
                ch1=base
        ch1=''.join(ch1.split('\n')[1:])
        for base in e:
            if contig_name == base.split('\n')[0].split(contig_spliter)[0]:
                contig=base
        contig=''.join(contig.split('\n')[1:])
        if strand=='-':
            contig=rc(contig)
        n=0
        f=['awad\tawad\t\tawad\t\tawad\t\tawad\t\tawad\t\tawad\t\tawad\t\tawad\t\tawad']
        ff=[]
        while n<len(c):
            if c[n][-1]=='S' or c[n][-1]=='H':
                contig_pos=contig_pos+int(c[n][:-1])
            elif c[n][-1]=='D' or c[n][-1]=='N':
                m=int(c[n][:-1])
                delation=''
                insertion='.'
                r=0
                while r<m:
                    delation+=ch1[chr_pos+r]
                    r=r+1
                f.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+delation+'\t'+insertion+'\t'+quality+'\tbass\t'+str(contig_pos+bse)+'\t'+'-'+str(m)+'\t'+strand+'\n')
                ff.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+delation+'\t'+insertion+'\t'+quality+'\tbass\t'+str(contig_pos+bse)+'\t'+'-'+str(m)+'\t'+strand+'\n')
                print(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+delation+'\t'+insertion+'\t'+quality+'\tbass\t'+str(contig_pos+bse)+'\t'+'-'+str(m)+'\t'+strand+'\n')
                chr_pos=chr_pos+m
            elif c[n][-1]=='I' or c[n][-1]=='P':
                m=int(c[n][:-1])
                delation='.'
                insertion=''
                r=0
                while r<m:
                    insertion+=contig[contig_pos+r]
                    r=r+1
                f.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+delation+'\t'+insertion+'\t'+quality+'\tbass\t'+str(contig_pos+bse)+'\t'+'+'+str(m)+'\t'+strand+'\n')
                ff.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+delation+'\t'+insertion+'\t'+quality+'\tbass\t'+str(contig_pos+bse)+'\t'+'+'+str(m)+'\t'+strand+'\n')
                print(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+delation+'\t'+insertion+'\t'+quality+'\tbass\t'+str(contig_pos+bse)+'\t'+'+'+str(m)+'\t'+strand+'\n')
                contig_pos=contig_pos+m
            elif c[n][-1]=='M' or c[n][-1]=='=':
                chr_seg=ch1[chr_pos:chr_pos+int(c[n][:-1])]
                contig_seg=contig[contig_pos:contig_pos+int(c[n][:-1])]
                if chr_seg!=contig_seg:
                    q=0
                    while q<len(chr_seg):
                        if chr_seg[q]==contig_seg[q]:
                            qqq=0
                            while q<len(chr_seg) and chr_seg[q]==contig_seg[q]:
                                q=q+1
                                qqq=qqq+1
                            chr_pos=chr_pos+qqq
                            contig_pos=contig_pos+qqq
                        else:
                            qq=0
                            while q<len(chr_seg) and chr_seg[q]!=contig_seg[q]:
                                qq=qq+1
                                q=q+1
                            m=qq
                            ref=''
                            con=''
                            r=0
                            while r<m:
                                con+=contig[contig_pos+r]
                                ref+=ch1[chr_pos+r]
                                r=r+1
                            if str(chr_pos+bse) ==f[-1].split('\t')[1]:
                                new=f[-1].split('\t')
                                new_n=new[:3]+[ref]+[new[4]+con]+new[5:]
                                f.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+ref+'\t'+con+'\t'+quality+'\tbass\t'+str(contig_pos+bse)+'\t'+'*'+str(m)+'\t'+strand+'\n')
                                ff[-1]='\t'.join(new_n)
                                print(ff[-1])
                            elif str(contig_pos+bse) ==f[-1].split('\t')[7]:
                                new=f[-1].split('\t')
                                new_n=new[:3]+[new[3]+ref]+[con]+new[5:]
                                f.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+ref+'\t'+con+'\t'+quality+'\tbass\t'+str(contig_pos+bse)+'\t'+'*'+str(m)+'\t'+strand+'\n')
                                ff[-1]='\t'.join(new_n)
                                print(ff[-1])
                            else:
                                f.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+ref+'\t'+con+'\t'+quality+'\tbass\t'+str(contig_pos+bse)+'\t'+'*'+str(m)+'\t'+strand+'\n')
                                ff.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+ref+'\t'+con+'\t'+quality+'\tbass\t'+str(contig_pos+bse)+'\t'+'*'+str(m)+'\t'+strand+'\n')
                            print(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+ref+'\t'+con+'\t'+quality+'\tbass\t'+str(contig_pos+bse)+'\t'+'*'+str(m)+'\t'+strand+'\n')
                            contig_pos=contig_pos+qq
                            chr_pos=chr_pos+qq
                            qq=0
                else:
                    chr_pos=chr_pos+int(c[n][:-1])
                    contig_pos=contig_pos+int(c[n][:-1])
            elif c[n][-1]=='X':
                m=int(c[n][:-1])
                ref=''
                con=''
                r=0
                while r<m:
                    con+=contig[contig_pos+r]
                    ref+=ch1[chr_pos+r]
                    r=r+1
                f.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+ref+'\t'+con+'\t'+quality+'\tbass\t'+str(contig_pos+bse)+'\t'+'*'+str(m)+'\t'+strand+'\n')
                ff.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+ref+'\t'+con+'\t'+quality+'\tbass\t'+str(contig_pos+bse)+'\t'+'*'+str(m)+'\t'+strand+'\n')
                print(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+ref+'\t'+con+'\t'+quality+'\tbass\t'+str(contig_pos+bse)+'\t'+'*'+str(m)+'\t'+strand+'\n')
                contig_pos=contig_pos+m
                chr_pos=chr_pos+m
            n=n+1
        if f[0]=='awad\tawad\t\tawad\t\tawad\t\tawad\t\tawad\t\tawad\t\tawad\t\tawad\t\tawad':
            del f[0]
        fff=ff[:]
        fff.append('0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0')
        nn=0
        while nn<len(fff)-1:
            dd=fff[nn].split('\t')
            ee=fff[nn+1].split('\t')
            if (int(dd[1])+len(dd[3])==int(ee[1])) and (int(dd[7])+len(dd[4])==int(ee[7])):
                while (int(dd[1])+len(dd[3])==int(ee[1])) and (int(dd[7])+len(dd[4])==int(ee[7])):
                    new_s=dd[:3]+[dd[3]+ee[3]]+[dd[4]+ee[4]]+dd[5:8]+[str(len(dd[4]+ee[4])-len(dd[3]+ee[3]))]+[dd[9]]
                    if new_s[8][0]!='+' and new_s[8][0]!='-' and new_s[8][0]!='*':
                        new_s[8]='+'+new_s[8]
                    fff[nn+1]='\t'.join(new_s)
                    del fff[nn]
                    dd=fff[nn].split('\t')
                    ee=fff[nn+1].split('\t')
            else:
                nn=nn+1
        if fff[-1]=='0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0':
            del fff[-1]
        z[s]=f
        y[s]=ff
        x[s]=fff
        s=s+1
    
    h=open(output_path+'/'+output_name+'.cvf','w')
    for key in z.keys():
        q=z[key]
        for base in q:
            h.writelines(base)
    h.close()
    h=open(output_path+'/'+output_name+'_gather.cvf','w')
    for key in y.keys():
        q=y[key]
        for base in q:
            h.writelines(base)
    h.close()
    h=open(output_path+'/'+output_name+'_final.cvf','w')
    for key in x.keys():
        q=x[key]
        for base in q:
            h.writelines(base)
    h.close()    
    return(s)