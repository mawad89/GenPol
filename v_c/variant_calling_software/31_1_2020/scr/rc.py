#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 15:01:29 2018

@author: awad
"""

def rc(seq):
    seq=seq.upper()
    a=[]
    for base in seq:
        if base=='A':
            a.append('T')
        elif base=="T":
            a.append("A")
        elif base=='C':
            a.append('G')
        elif base=='G':
            a.append("C")
        elif base=="W":
            a.append("W")
        elif base=="S":
            a.append("S")
        elif base=="M":
            a.append("K")
        elif base=="K":
            a.append("M")
        elif base=="R":
            a.append("Y")
        elif base=="Y":
            a.append('R')
        elif base=="B":
            a.append("V")
        elif base=="V":
            a.append("B")
        elif base=="D":
            a.append("H")
        elif base=="H":
            a.append("D")
        else:
            a.append("N")
    a.reverse()
    return (''.join(a))