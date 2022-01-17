#!/usr/bin/env python
'''
Created on 15/01/2018

@author: sium
'''

import yaml
import sys
from re import search


#hpv_name=search('(HPV[0-9]+)',sys.argv[1]).group(1)

"""
with open(sys.argv[1]) as txt:
    content=yaml.load(txt.read().replace("\t"," "))


    for stat,value in content.items():
        print("{stat}\t{value}\t{strain}\t{file}".format(stat=stat,value=value,strain=hpv_name,file=sys.argv[1].strip(".statistics.txt")))
"""

content=yaml.load(sys.stdin.read().replace("\t"," "))

for stat in sorted(content.keys()):
    print("{stat}\t{value}".format(stat=stat, value=content[stat]))












