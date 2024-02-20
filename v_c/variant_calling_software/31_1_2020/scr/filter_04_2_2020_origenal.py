#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 13:06:20 2020

@author: awad
"""
import os
def filteration(bed,vc_file,output_path=os.getcwd(),output_name='var_call_awad'):
    a=list(open(bed))
    b=list(open(vc_file))
    z={}
    x,y={},{}
    for base in b:
        z[base.split('\t')[0]]=[]
    for base in b:
        z[base.split('\t')[0]]+=[base]
    aaa=[]
    bbb=[]
    for base in a:
        aaa.append(int(base.split('\t')[1]))
        bbb.append(int(base.split('\t')[2]))
    
    cccc={}
    for nnn in range(0,max(bbb)+1):
        cccc[nnn]=0
    ddd=[]
    for nnn in range(0,len(bbb)):
        ddd.append((aaa[nnn],bbb[nnn]))
    
    for base in ddd:
        for nnn in range(base[0],base[1]+1):
            cccc[nnn]+=1
    ccc={}
    for nnn in range(0,max(bbb)+1):
        ccc[nnn]=[]
    del z['awad']
    for key in z.keys():
        q=z[key]
        for base in q:
            try:
                ccc[int(base.split('\t')[1])]+=[base]
            except:
                awad=0
    del q
    eee={}
    for nnn in range(1,max(bbb)+1):
        eee[nnn]=[['.','A',0],['.','T',0],['.','C',0],['.','G',0],['.','A',0],['.','T',0],['.','C',0],['.','G',0],['.','.',0],['N','N',0]]
    chr_name=z.keys()[0]#del it in the update
    
    for base in ccc.keys():
        q=ccc[base]
        #chr_name=base#del it in the update
        for ba in q:
            if ba.split('\t')[3]=='.':
                if ba.split('\t')[4]=='A':
                    eee[base][0][-1]+=1
                elif ba.split('\t')[4]=='T':
                    eee[base][1][-1]+=1
                elif ba.split('\t')[4]=='C':
                    eee[base][2][-1]+=1
                elif ba.split('\t')[4]=='G':
                    eee[base][3][-1]+=1
                else:
                    eee[base][9][-1]+=1
            else:
                if ba.split('\t')[4]=='A':
                    eee[base][4][-1]+=1
                    eee[base][4][0]=ba.split('\t')[3]
                    eee[base][5][0]=ba.split('\t')[3]
                    eee[base][6][0]=ba.split('\t')[3]
                    eee[base][7][0]=ba.split('\t')[3]
                    eee[base][8][0]=ba.split('\t')[3]
                elif ba.split('\t')[4]=='T':
                    eee[base][5][-1]+=1
                    eee[base][4][0]=ba.split('\t')[3]
                    eee[base][5][0]=ba.split('\t')[3]
                    eee[base][6][0]=ba.split('\t')[3]
                    eee[base][7][0]=ba.split('\t')[3]
                    eee[base][8][0]=ba.split('\t')[3]
                elif ba.split('\t')[4]=='C':
                    eee[base][6][-1]+=1
                    eee[base][4][0]=ba.split('\t')[3]
                    eee[base][5][0]=ba.split('\t')[3]
                    eee[base][6][0]=ba.split('\t')[3]
                    eee[base][7][0]=ba.split('\t')[3]
                    eee[base][8][0]=ba.split('\t')[3]
                elif ba.split('\t')[4]=='G':
                    eee[base][7][-1]+=1
                    eee[base][4][0]=ba.split('\t')[3]
                    eee[base][5][0]=ba.split('\t')[3]
                    eee[base][6][0]=ba.split('\t')[3]
                    eee[base][7][0]=ba.split('\t')[3]
                    eee[base][8][0]=ba.split('\t')[3]
                elif ba.split('\t')[4]=='.':
                    eee[base][8][-1]+=1
                    eee[base][4][0]=ba.split('\t')[3]
                    eee[base][5][0]=ba.split('\t')[3]
                    eee[base][6][0]=ba.split('\t')[3]
                    eee[base][7][0]=ba.split('\t')[3]
                    eee[base][8][0]=ba.split('\t')[3]
                else:
                    eee[base][9][-1]+=1
    ff=[]
    for base in eee.keys():
        qq=eee[base]
        qqq=cccc[base-1]
        if qqq==0:
            qqq=1
        qqqq=qqq-(qq[0][-1]+qq[1][-1]+qq[2][-1]+qq[3][-1]+qq[4][-1]+qq[5][-1]+qq[6][-1]+qq[7][-1]+qq[8][-1]+qq[9][-1])
        if qqqq==0:
            qqqq=-1
        for ba in qq:
            if ba[0]=='.' and ba[2]!=0:
                ff.append(chr_name+'\t'+str(base)+'\t'+str(ba[0])+'\t'+str(ba[1])+'\t+1\t'+str(qqq)+'\t'+str(qqqq)+'\t'+str(ba[2])+'\t'+str(float(float(qqqq)/float(qqq)))+'\t'+str(float(float(ba[2])/float(qqq)))+'\t'+str(float(float(ba[2])/float(qqqq)))+'\n')
            else:
                if ba[1]=='.' and ba[2]!=0:
                    ff.append(chr_name+'\t'+str(base)+'\t'+str(ba[0])+'\t'+str(ba[1])+'\t-1\t'+str(qqq)+'\t'+str(qqqq)+'\t'+str(ba[2])+'\t'+str(float(float(qqqq)/float(qqq)))+'\t'+str(float(float(ba[2])/float(qqq)))+'\t'+str(float(float(ba[2])/float(qqqq)))+'\n')
                elif ba[1]!='.' and ba[2]!=0:
                    ff.append(chr_name+'\t'+str(base)+'\t'+str(ba[0])+'\t'+str(ba[1])+'\t*1\t'+str(qqq)+'\t'+str(qqqq)+'\t'+str(ba[2])+'\t'+str(float(float(qqqq)/float(qqq)))+'\t'+str(float(float(ba[2])/float(qqq)))+'\t'+str(float(float(ba[2])/float(qqqq)))+'\n')
    
    x[chr_name]=ff
    fff=[]
    
    for key in x.keys():
        q=x[key]
        for base in q:
            if abs(float(base.split('\t')[-1][:-1]))>1:#7.2.2020
                fff.append(base)
    ffff={}
    for bas in fff:
        ffff[bas.split('\t')[1]]=[]
    for bas in fff:
        ffff[bas.split('\t')[1]]+=[bas]
    fff=[]
    for bas in ffff.keys():
        if len(ffff[bas])==1:
            fff.append(ffff[bas][0])
        else:
            fffff=sorted(ffff[bas], key=lambda i:abs(float(i.split('\t')[-1])), reverse=True)
            fff.append(fffff[0])
    y[chr_name]=fff
    h=open(output_path+'/'+output_name+'_final.vc','w')
    for key in x.keys():
        q=x[key]
        for base in q:
            h.writelines(base)
    h.close()
    h=open(output_path+'/'+output_name+'_filtered.vc','w')
    hh=open(output_path+'/'+output_name+'_filtered.sdi','w')
    for key in y.keys():
        q=y[key]
        q=sorted(q, key=lambda i:int(i.split('\t')[1]))
        for base in q:
            h.writelines(base)
            hh.writelines(base.split('\t')[0]+'\t'+base.split('\t')[1]+'\t'+base.split('\t')[4].replace('+','').replace('*1','0')+'\t'+base.split('\t')[2]+'\t'+base.split('\t')[3]+'\t*\t*\n')
    h.close()
    hh.close()
    return(x,y)
import argparse
parser = argparse.ArgumentParser(prog="variant_caller",usage='%(prog)s -h  [options] -b <bed_file>  -r <ref_genome> -g <draft>',description='Variant caller',)
parser.add_argument("-b",nargs=1,type=str,required=True,metavar='mapping bed file',dest="bed")
parser.add_argument("-v",nargs=1,type=str,required=True,metavar='vc_file',dest="vc")
parser.add_argument("-o",nargs=1, type=str, default=[os.getcwd()], metavar='output path',dest="output")
parser.add_argument("-f",nargs=1, type=str, default=['var_call_awad'], metavar='output file name',dest="name")
args = parser.parse_args()
bed=args.bed[0]
vc_file=args.vc[0]
output=args.output[0]
name=args.name[0]
filteration(bed=bed,vc_file=vc_file,output_path=output,output_name=name)