#!/usr/bin/env python

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-o", "--outfile", type=str,
        help="Write tab-delimited enzyme to pathway table to outfile. Defaults to stdout")
parser.add_argument("-I", "--include", type=str,
        help="Include only these enzymes in the lookup")
parser.add_argument("-E", "--exclude", type=str,
        help="Exclude these enzymes in the lookup")
args = parser.parse_args()

import sys,csv, urllib2

def ReadFile(f):
    d = {}
    hin = open(f, 'r')
    for line in hin: d[line.rstrip()] = ""
    hin.close()
    return d.keys()

from bioservices.kegg import KEGG, KEGGParser
k = KEGG()
p = KEGGParser()

exclude = []
if args.exclude: exclude = ReadFile(args.exclude)
if args.include: enzymes = ReadFile(args.include)
else: enzymes = p.enzymeIds

ecs = {}

if args.outfile: hout = open(args.outfile, 'w')
else: hout = sys.stdout
houtcsv = csv.writer(hout, delimiter = '\t')

for ec in enzymes:
    ec = ec.replace("ec:","")
    if ec in exclude: continue
    l = []
    try: result = k.get(ec)
    except urllib2.HTTPError: continue
    parsed = p.parse(result)
    ## Check if the enzyme is obsolete
    if "Obsolete" in parsed["entry"]: continue
    #sys.stderr.write(ec+"\n")
    ## Check that the enzyme has a pathway key
    try: parsed["pathway"]
    except KeyError: 
        ## Try to find pathways by orthology
        try: parsed["orthology"]
        except KeyError: continue
        try: ko = parsed["orthology"].rsplit()[0]
        except AttributeError: 
            kos = parsed["orthology"].keys()
            ## KOs here are only potential KOs because of occasional bad parsing
            foundKO = False
            for ko in kos:
                try: int(ko[1:])
                except ValueError: continue
                if ko[0] == "K": 
                    foundKO = True
                    break
            if foundKO: result = k.get(ko)
            else: continue  
        parsed = p.parse(result)
        try: parsed["pathway"]
        except KeyError: 
            ## Try brite if no pathway
            try:brite = parsed["brite"]
            except KeyError: continue
            paths = []
            for item in brite:
                try: part = item.split(":")[-1]
                except IndexError: continue
                part = part.replace("]","")
                path = part[2:]
                try: int(path)
                except ValueError: continue
                paths.append(path)
            if len(paths)==0: continue
            parsed["pathway"] = {}
            for path in paths: parsed["pathway"][path] = ""


    ## If only one pathway the pathway ID needs to be handled differently
    try: 
        parsed["pathway"].keys()
        pathways = parsed["pathway"].keys()
    except AttributeError: pathways = [parsed["pathway"].rsplit()[0]]
    for pathway in pathways: 
        pathway = pathway[2:]
        try: int(pathway)
        except ValueError: continue
        l.append(pathway)
    l = list(set(l))
    for pathway in l: houtcsv.writerow([pathway,ec])
hout.close()

