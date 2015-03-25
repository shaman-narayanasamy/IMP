#!/usr/bin/env python
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-i", "--infile", type=str,
        help="Supply a BRITE htext functional hierarchy flatfile. If none is supplied the script attempts to download one from http://www.kegg.jp/kegg-bin/download_htext?htext=br08901.keg&format=htext&filedir=")
parser.add_argument("-o", "--outfile", type=str,
        help="Write tab-delimited pathway hierarchy to outfile. Defaults to stdout")
parser.add_argument("-I", "--include", type=str,
        help="Include only these pathways in the lookup")
parser.add_argument("-E", "--exclude", type=str,
        help="Exclude these pathways in the lookup")
args = parser.parse_args()

import sys

def ReadFile(f):
    d = {}
    hin = open(f, 'r')
    for line in hin: d[line.rstrip()] = ""
    hin.close()
    return d.keys()

import csv

def ParseBRITE(data, include, exclude):
    pathways = {}
    for line in data:
        try: 
            if not line[0] in ["A","B","C"]: continue
        except IndexError: continue
        if line[0] == "A": A = line.split(">")[1].split("<")[0]
        if line[0] == "B": B = (line.lstrip("B")).lstrip(" ")
        if line[0] == "C":
            line = line.lstrip("C")
            line = line.lstrip(" ")
            pwy = line[0:5]
            if len(include) > 0 and not pwy in include: continue
            if pwy in exclude: continue
            name = (line.replace(pwy,"")).lstrip()
            pathways[pwy] = {"name": name, "hier": "|".join([A,B])}
    return pathways

exclude = []
if args.exclude: exclude = ReadFile(args.exclude)
include = []
if args.include: include = ReadFile(args.include)

## Read the hierarchy flatfile
data = []
if args.infile: 
    hin = open(args.infile, 'r')
    data = hin.read()
else:
    import urllib2
    f = urllib2.urlopen("http://www.kegg.jp/kegg-bin/download_htext?htext=br08901.keg&format=htext&filedir=")
    data = f.read()

data = data.split("\n")
pathways = ParseBRITE(data, include, exclude)

if args.outfile: hout = open(args.outfile, 'w')
else: hout = sys.stdout
houtcsv = csv.writer(hout, delimiter = '\t')

for p, d in pathways.iteritems(): houtcsv.writerow([p,d["name"],d["hier"]])
hout.close()
    

