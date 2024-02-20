#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 11:53:44 2020

@author: awad
"""

import os
from rc import rc
'''
bed_file='/netscratch/dep_tsiantis/grp_gan/awad/todo/elegans/comp/new_try/scaffolding/final/final_version/quiver_4/m.bed'
ref_genome='/netscratch/dep_tsiantis/grp_gan/awad/todo/elegans/comp/new_try/scaffolding/final/final_version/quiver_4/m.fa'
draft='/netscratch/dep_tsiantis/grp_gan/awad/todo/elegans/comp/new_try/mapping_corrected_first_pacbio/corrected_first_pacbio_pacbio_fly/contig_44_awad.bam.read.fq'
output_path='/netscratch/dep_tsiantis/grp_gan/awad/todo/elegans/comp/new_try/scaffolding/final/final_version/quiver_4'
output_name='m_awad_new_try'
chr_spliter=' '
contig_spliter=' '
bse=1
qty=30
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
    x={}
    y={}
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
            contig=contig.split('\n')[1]
        if strand=='-':
            contig=rc(contig)
        n=0
        f=['awad\tawad\t\tawad\t\tawad\t\tawad\t\tawad\t\tawad\t\tawad\t\tawad\t\tawad\n']
        while n<len(c):
            if c[n][-1]=='S' or c[n][-1]=='H':
                contig_pos=contig_pos+int(c[n][:-1])
            elif c[n][-1]=='D' or c[n][-1]=='N':
                m=int(c[n][:-1])
                delation=''
                insertion='.'
                r=0
                while r<m:
                    delation=ch1[chr_pos+r]
                    f.append(chr_name+'\t'+str(chr_pos+bse+r)+'\t'+contig_name+'\t'+delation+'\t'+insertion+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse+r)+'\t'+'-'+'1'+'\t'+strand+'\n')
                    print(chr_name+'\t'+str(chr_pos+bse+r)+'\t'+contig_name+'\t'+delation+'\t'+insertion+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse+r)+'\t'+'-'+'1'+'\t'+strand+'\n')
                    r=r+1
                chr_pos=chr_pos+r
            elif c[n][-1]=='I' or c[n][-1]=='P':
                m=int(c[n][:-1])
                delation='.'
                insertion=''
                r=0
                while r<m:
                    insertion=contig[contig_pos+r]
                    f.append(chr_name+'\t'+str(chr_pos+bse+r)+'\t'+contig_name+'\t'+delation+'\t'+insertion+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse+r)+'\t'+'+'+'1'+'\t'+strand+'\n')
                    print(chr_name+'\t'+str(chr_pos+bse+r)+'\t'+contig_name+'\t'+delation+'\t'+insertion+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse+r)+'\t'+'+'+'1'+'\t'+strand+'\n')
                    r=r+1
                contig_pos=contig_pos+r
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
                                con=contig[contig_pos+r]
                                ref=ch1[chr_pos+r]
                                f.append(chr_name+'\t'+str(chr_pos+bse+r)+'\t'+contig_name+'\t'+ref+'\t'+con+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse+r)+'\t'+'*1'+'\t'+strand+'\n')
                                print(chr_name+'\t'+str(chr_pos+bse+r)+'\t'+contig_name+'\t'+ref+'\t'+con+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse+r)+'\t'+'*1'+'\t'+strand+'\n')
                                r=r+1              
                            chr_pos=chr_pos+qq
                            contig_pos=contig_pos+qq
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
                    con=contig[contig_pos+r]
                    ref=ch1[chr_pos+r]
                    f.append(chr_name+'\t'+str(chr_pos+bse+r)+'\t'+contig_name+'\t'+ref+'\t'+con+'\t'+quality+'\tSl_1\t'+str(contig_pos+bse+r)+'\t'+'*1'+'\t'+strand+'\n')
                    r=r+1
                contig_pos=contig_pos+r
                chr_pos=chr_pos+r
            n=n+1
        z[s]=f
        s=s+1

    h=open(output_path+'/'+output_name+'_first.vc','w')
    for key in z.keys():
        q=z[key]
        for base in q:
            h.writelines(base)
    h.close()
    #h=open(output_path+'/'+output_name+'_covarage.vc','w')
    #for key in cccc.keys():
    #    h.writelines(str(key)+'\t'+str(cccc[key])+'\n')
    #h.close()
                    
    
    return(s)
#a=variant_caller(bed_file=bed_file,ref_genome=ref_genome,draft=draft,output_path=output_path,output_name='new_try_awad',chr_spliter=' ',contig_spliter=' ',bse=bse,qty=qty)

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

