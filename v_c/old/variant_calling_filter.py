#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 18:01:10 2019

@author: awad
"""
import os
def dect_collect(dect,key,mother_dect,region_dect,region):
    dect[key]=[]
    for base in mother_dect[key]:
        cc='\t'.join(base.split('\t')[3:5])
        if cc==region_dect[key][region]:
            dect[key]=base
            break
    return(dect)
    
def filterate(variant_list,files=False,output_path=os.getcwd(),output_name='var_call_awad'):
    #if type(variant_list)==file: 
    a=list(open(variant_list))
    while '\n' in a:
        del a[a.index('\n')]
    #elif type(variant_list)==dict:
    #    a=[]
     #   for key in variant_list.keys():
      #      for base in variant_list:
       #         a.append(base)
   # elif type(variant_list)==list:
    #    a=variant_list
    b={}
    for base in a:
        cc='\t'.join(base.split('\t')[:2])
        b[cc]=[]
    for base in a:
        cc='\t'.join(base.split('\t')[:2])
        b[cc]+=[base]
        
    c,d={},{}
    for base in b.keys():
        if len(b[base])>1:
            c[base]=b[base]
        else:
            d[base]=b[base]
    e,f={},{}
    for base in c.keys():
        e[base]=[]
        f[base]=[]
        for ba in c[base]:
            e[base]+=['\t'.join(ba.split('\t')[3:5])]
        for ba in set(e[base]):
            f[base]+=[ba+'\t'+str(e[base].count(ba))+'\t'+str(float(e[base].count(ba))/len(e[base]))]
        f[base]=sorted(f[base],key= lambda i:float(i.split('\t')[-1]),reverse= True)
    g,h,i,j,k,l,m,n,o={},{},{},{},{},{},{},{},{}
    gg,hh,ii,jj,ll,mm,oo={},{},{},{},{},{},{}
    for key in f.keys(): 
        if float(f[key][0].split('\t')[-1])==1:
            g[key]=e[key][0]
            gg=dect_collect(gg,key,c,e,0)
        elif float(f[key][0].split('\t')[-1])>0.5:
            h[key]=e[key][0]
            hh=dect_collect(hh,key,c,e,0)
        elif float(f[key][0].split('\t')[-1])==0.5 and float(f[key][1].split('\t')[-1])==0.5:
            i[key]=e[key]
            ii=dect_collect(ii,key,c,e,0)
            oo=dect_collect(oo,key,c,e,1)
        elif float(f[key][0].split('\t')[-1])==0.5:
            n[key]=e[key][0]
            ii=dect_collect(ii,key,c,e,0)
        elif float(f[key][0].split('\t')[-1])==float(f[key][1].split('\t')[-1]) and float(f[key][0].split('\t')[-1])>0.349:
            j[key]=e[key][:2]
            jj=dect_collect(jj,key,c,e,0)
            oo=dect_collect(oo,key,c,e,1)
        elif float(f[key][0].split('\t')[-1])>0.349:
            k[key]=e[key][0]
            jj=dect_collect(jj,key,c,e,0)
        elif len(f[key])==3 and float(f[key][0].split('\t')[-1])>0.329:
            l[key]=e[key]
            ll=dect_collect(ll,key,c,e,0)
            o=dect_collect(o,key,c,e,1)
            o=dect_collect(o,key,c,e,2)
        else:
            m[key]=e[key]
            mm=dect_collect(mm,key,c,e,0)
    p=[]
    for base in d.keys():
        p.append(d[base][0])
    for base in gg.keys():
        p.append(gg[base].replace('Sl_1','r_1'))
    for base in hh.keys():
        p.append(hh[base].replace('Sl_1','r_>50'))
    for base in ii.keys():
        p.append(ii[base].replace('Sl_1','r_50'))
    for base in jj.keys():
        p.append(jj[base].replace('Sl_1','r_<50'))
    for base in ll.keys():
        p.append(ll[base].replace('Sl_1','r_33'))
    for base in mm.keys():
        p.append(mm[base].replace('Sl_1','r_<33'))
    q=[]
    for base in oo.keys():
        q.append(oo[base].replace('Sl_1','r_50'))
    for base in o.keys():
        q.append(o[base].replace('Sl_1','r_33'))
    p=sorted(p,key= lambda i:(i.split('\t')[0],int(i.split('\t')[1])))
    q=sorted(q,key= lambda i:(i.split('\t')[0],int(i.split('\t')[1])))
    if files ==True:
        h=open(output_path+'/'+output_name+'_final.vc','w')
        for base in p:
            h.writelines(base)
        h.close()
        h=open(output_path+'/'+output_name+'_problematic.vc','w')
        for base in q:
            h.writelines(base)
        h.close()
    return(p,q)
    
import argparse
parser = argparse.ArgumentParser(prog="variant_caller_filter",usage='%(prog)s -h  [options] -v <variant file> ',description='Variant caller filter',)
parser.add_argument("-v",nargs=1,type=str,required=True,metavar='variant file',dest="var")
parser.add_argument("-o",nargs=1, type=str, default=[os.getcwd()], metavar='output path',dest="output")
parser.add_argument("-f",nargs=1, type=str, default=['var_call_awad'], metavar='output file name',dest="name")
parser.add_argument("-m",nargs=1, type=str, default=[' '], metavar='chromosome spilter',dest="ch")
parser.add_argument("-c",nargs=1, type=str, default=[' '], metavar='contig spilter',dest="cont")
parser.add_argument("-s",nargs=1, type=int, default=[0], metavar='base 0 or 1',dest="base")
args = parser.parse_args()
var=args.var[0]
output=args.output[0]
name=args.name[0]
ch=args.ch[0]
cont=args.cont[0]
base=args.base[0]
filterate(var,files=True,output_path=output,output_name=name)