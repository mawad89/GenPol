#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 11:38:17 2019

@author: awad
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 23:24:31 2019

@author: awad
"""

import os
from rc import rc
'''
bed_file='/netscratch/dep_tsiantis/grp_gan/awad/todo/18_2_2019_vc/canu_nano/canu_sort_1.bed'
ref_genome='/biodata/dep_tsiantis/common/chigenome/chi_v1.fa'
draft='/netscratch/dep_tsiantis/grp_gan/awad/todo/9_1_2019_new_nano/canu_all/canu_assemble/canu_assemble.contigs.fasta'
output_path='/netscratch/dep_tsiantis/grp_gan/awad/todo/18_2_2019_vc/canu_nano'
chr_spliter=' '
contig_spliter=' '
bse=1
'''


def variant_caller(bed_file,ref_genome,draft,output_path=os.getcwd(),output_name='var_call_awad',chr_spliter=' ',contig_spliter=' ',bse=1,qty=30):
    a=list(open(bed_file))
    d=list(open(ref_genome))
    d=''.join(d)
    d=d.split('>')
    del d[0]
    e=list(open(draft))
    e=''.join(e)
    if e[0]=='>':
        e=e.split('>')
        typ='fa'
    else:
        e=e.split('@')
        typ='fq'
    del e[0]
    s=0
    z={}
    y={}
    x={}
    aa=a[:]
    a=[]
    for base in aa:
        if int(base.split('\t')[4])>qty-1:
            a.append(base)
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
                ch1=base.upper()
        ch1=''.join(ch1.split('\n')[1:])
        for base in e:
            if contig_name == base.split('\n')[0].split(contig_spliter)[0]:
                contig=base.upper()
        if typ=='fa':
            contig=''.join(contig.split('\n')[1:])
        else:
            ch1=''.join(ch1.split('\n')[1:])
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
                f.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+delation+'\t'+insertion+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse)+'\t'+'-'+str(m)+'\t'+strand+'\n')
                ff.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+delation+'\t'+insertion+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse)+'\t'+'-'+str(m)+'\t'+strand+'\n')
                print(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+delation+'\t'+insertion+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse)+'\t'+'-'+str(m)+'\t'+strand+'\n')
                chr_pos=chr_pos+m
            elif c[n][-1]=='I' or c[n][-1]=='P':
                m=int(c[n][:-1])
                delation='.'
                insertion=''
                r=0
                while r<m:
                    insertion+=contig[contig_pos+r]
                    r=r+1
                f.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+delation+'\t'+insertion+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse)+'\t'+'+'+str(m)+'\t'+strand+'\n')
                ff.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+delation+'\t'+insertion+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse)+'\t'+'+'+str(m)+'\t'+strand+'\n')
                print(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+delation+'\t'+insertion+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse)+'\t'+'+'+str(m)+'\t'+strand+'\n')
                contig_pos=contig_pos+m
            elif c[n][-1]=='M' or c[n][-1]=='=':
                chr_seg=ch1[chr_pos:chr_pos+int(c[n][:-1])]
                contig_seg=contig[contig_pos:contig_pos+int(c[n][:-1])]
                if chr_seg!=contig_seg:
                    q=0
                    while q<len(chr_seg):
                        if q<len(contig_seg) and chr_seg[q]==contig_seg[q]:#add the first one need to review
                            qqq=0
                            while q<len(chr_seg) and chr_seg[q]==contig_seg[q]:
                                q=q+1
                                qqq=qqq+1
                            chr_pos=chr_pos+qqq
                            contig_pos=contig_pos+qqq
                        else:
                            qq=0
                            while q<len(chr_seg) and q<len(contig_seg) and chr_seg[q]!=contig_seg[q]:#add the second one need to review
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
                                f.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+ref+'\t'+con+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse)+'\t'+'*'+str(m)+'\t'+strand+'\n')
                                ff[-1]='\t'.join(new_n)
                                print(ff[-1])
                            elif str(contig_pos+bse) ==f[-1].split('\t')[7]:
                                new=f[-1].split('\t')
                                new_n=new[:3]+[new[3]+ref]+[con]+new[5:]
                                f.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+ref+'\t'+con+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse)+'\t'+'*'+str(m)+'\t'+strand+'\n')
                                ff[-1]='\t'.join(new_n)
                                print(ff[-1])
                            else:
                                f.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+ref+'\t'+con+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse)+'\t'+'*'+str(m)+'\t'+strand+'\n')
                                ff.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+ref+'\t'+con+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse)+'\t'+'*'+str(m)+'\t'+strand+'\n')
                            print(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+ref+'\t'+con+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse)+'\t'+'*'+str(m)+'\t'+strand+'\n')
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
                f.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+ref+'\t'+con+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse)+'\t'+'*'+str(m)+'\t'+strand+'\n')
                ff.append(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+ref+'\t'+con+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse)+'\t'+'*'+str(m)+'\t'+strand+'\n')
                print(chr_name+'\t'+str(chr_pos+bse)+'\t'+contig_name+'\t'+ref+'\t'+con+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse)+'\t'+'*'+str(m)+'\t'+strand+'\n')
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
                    new_s=dd[:3]+[dd[3]+ee[3]]+[dd[4]+ee[4]]+dd[5:8]+[str(len(dd[4]+ee[4])-len(dd[3]+ee[3]))]+dd[9:]
                    if len(new_s[3])>1 and new_s[3][-1]=='.' and len(new_s[4])>1 and new_s[4][-1]=='.':
                        new_s[3]=new_s[3].replace('.','')
                        new_s[4]=new_s[4].replace('.','')
                    elif len(new_s[3])>1 and new_s[3][-1]=='.':
                        new_s[3]=new_s[3].replace('.','')
                        new_s[8]=str(len(dd[4]+ee[4])-len(dd[3]+ee[3]))
                    elif len(new_s[4])>1 and new_s[4][-1]=='.':
                        new_s[4]=new_s[4].replace('.','')
                        new_s[8]=str(len(dd[4]+ee[4])-len(dd[3]+ee[3]))
                    if new_s[8][0]!='+' and new_s[8][0]!='-' and new_s[8][0]!='*':
                        new_s[8]='+'+str(new_s[8])
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
    h=open(output_path+'/'+output_name+'_first.vc','w')
    for key in z.keys():
        q=z[key]
        for base in q:
            h.writelines(base)
    h.close()
    h=open(output_path+'/'+output_name+'_gather.vc','w')
    for key in y.keys():
        q=y[key]
        for base in q:
            h.writelines(base)
    h.close()
    h=open(output_path+'/'+output_name+'_with_r.vc','w')
    for key in x.keys():
        q=x[key]
        for base in q:
            h.writelines(base)
    h.close()    
    return(s)

import argparse
parser = argparse.ArgumentParser(prog="variant_caller",usage='%(prog)s -h  [options] -b <bed_file>  -r <ref_genome> -g <draft>',description='Variant caller',)
parser.add_argument("-b",nargs=1,type=str,required=True,metavar='mapping bed file',dest="bed")
parser.add_argument("-r",nargs=1,type=str,required=True,metavar='referance genome fasta',dest="ref")
parser.add_argument("-g",nargs=1,type=str,required=True,metavar='draft genome fasta',dest="draft")
parser.add_argument("-o",nargs=1, type=str, default=[os.getcwd()], metavar='output path',dest="output")
parser.add_argument("-f",nargs=1, type=str, default=['var_call_awad'], metavar='output file name',dest="name")
parser.add_argument("-m",nargs=1, type=str, default=[' '], metavar='chromosome spilter',dest="ch")
parser.add_argument("-c",nargs=1, type=str, default=[' '], metavar='contig spilter',dest="cont")
parser.add_argument("-s",nargs=1, type=int, default=[0], metavar='base 0 or 1',dest="base")
parser.add_argument("-q",nargs=1, type=int, default=[30], metavar='mapping quality',dest="qty")
args = parser.parse_args()
bed=args.bed[0]
ref=args.ref[0]
draft=args.draft[0]
output=args.output[0]
name=args.name[0]
ch=args.ch[0]
cont=args.cont[0]
base=args.base[0]
qty=args.qty[0]
variant_caller(bed_file=bed,ref_genome=ref,draft=draft,output_path=output,output_name=name,chr_spliter=ch,contig_spliter=cont,bse=base,qty=qty)
