#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 17:02:47 2019

@author: awad
"""
import os

def variant_run(bed_file,ref,draft,bed_out_name='map',variant_out_name='calling_variant',output_path=os.getcwd(),software=os.path.dirname(os.path.realpath(__file__))+'/scr/variant_caller_28_1_2020.py',filterate=os.path.dirname(os.path.realpath(__file__))+'/scr/filter.py',chr_spliter="' '",contig_spliter="' '",bse=1,qual=30,q='multicore20',core=4,job_name='Awad_v_c',submit=False):
    a=list(open(bed_file))
    n,m,r,h,b=0,0,1,1,[]
    if output_path[-1]=='/':
       output_path=output_path[:-1] 
    os.makedirs(output_path+'/'+variant_out_name)
    output_path=output_path+'/'+variant_out_name
    d=open(output_path+'/'+'run.sh','w')
    d.writelines('#!/bin/sh -login\n')
    z=open(output_path+'/'+'file_assembly.sh','w')
    z.writelines('#!/bin/sh -login\ncat *first.vc > '+variant_out_name+'_first_all.vc\ncat *gather.vc > '+variant_out_name+'_gather_all.vc\ncat *with_r.vc > '+variant_out_name+'_with_r_all.vc')
    z.close()
    z=open(output_path+'/'+'filter.sh','w')
    z.writelines('#!/bin/sh -login\nname='+variant_out_name+'\nv_c_file='+variant_out_name+'_with_r_all.vc'+'\nsoftware='+filterate+'\noutput='+output_path+'\n')
    z.writelines('python $software -v $v_c_file -f $name -o $output\n')
    z.close()
    while n<len(a):
        while m<100000 and n<len(a):#review
            b.append(a[n])
            n=n+1
            m=m+1
        m=0
        c=open(output_path+'/'+bed_out_name+'_'+str(r)+'.bed','w')
        e=open(output_path+'/'+'c_v'+'_'+str(r)+'.sh','w')
        e.writelines('#!/bin/sh -login\nref='+ref+'\ndraft='+draft+'\nname='+variant_out_name+'_'+str(h)+'\nsoftware='+software+'\nqt='+str(qual)+'\nbed='+bed_out_name+'_'+str(r)+'.bed'+'\noutput='+output_path+'\nchr_spliter='+chr_spliter+'\ncontig_spliter='+contig_spliter+'\nbse='+str(bse)+'\n')
        e.writelines('python $software -b $bed -r $ref -g $draft -o $output -f $name -s $bse -q $qt\n')
        e.close()
        d.writelines('bsub -q '+q+' -n '+str(core)+' -J '+job_name+' sh '+'c_v'+'_'+str(r)+'.sh\n')
        h=h+1
        for base in b:
            c.writelines(base)
        c.close()
        r=r+1
        b=[]
    d.writelines('bsub -q normal -J gathering -w "done('+"'"+job_name+"'"+')" sh file_assembly.sh')
    d.writelines('\nbsub -q normal -w "done('+"'"+'gathering'+"'"+')" sh filter.sh')
    d.close()
    if submit==True:
        os.system('cd '+output_path)
        return(os.system("sh run.sh"))
    return ('done')


import argparse
parser = argparse.ArgumentParser(prog="variant_caller",usage='%(prog)s -h  [options] -b <bed_file>  -r <ref_genome> -g <draft>',description='Variant caller',)
parser.add_argument("-b",nargs=1,type=str,required=True,metavar='mapping bed file',dest="bed")
parser.add_argument("-r",nargs=1,type=str,required=True,metavar='referance genome fasta',dest="ref")
parser.add_argument("-g",nargs=1,type=str,required=True,metavar='draft genome fasta',dest="draft")
parser.add_argument("-o",nargs=1, type=str, default=[os.getcwd()], metavar='output path',dest="output")
parser.add_argument("-f",nargs=1, type=str, default=['var_call_awad'], metavar='output file name',dest="name")
parser.add_argument("-s",nargs=1, type=int, default=[1], metavar='base 0 or 1',dest="base")
parser.add_argument("-q",nargs=1, type=int, default=[30], metavar='mapping quality',dest="qty")
parser.add_argument("-i",nargs=1, type=str, default=['normal'], metavar='Grid_queue',dest="ch")
parser.add_argument("-n",nargs=1, type=str, default=['variant_calling'], metavar='Grid_job_name',dest="gn")
parser.add_argument("-t",nargs=1, type=bool, default=[False], metavar='Grid_job_submitiom',dest="submit")
args = parser.parse_args()
bed=args.bed[0]
ref=args.ref[0]
draft=args.draft[0]
output=args.output[0]
name=args.name[0]
ch=args.ch[0]
base=args.base[0]
qty=args.qty[0]
gn=args.gn[0]
submit=args.submit[0]
variant_run(bed_file=bed,ref=ref,draft=draft,variant_out_name=name,output_path=output,bse=base,qual=qty,q=ch,job_name=gn,submit=submit)
