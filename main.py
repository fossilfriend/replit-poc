
import json
import pprint
import sys
import csv
from collections import OrderedDict
from utils import qw, VepJsonParser

def extract_consequences(annotation, alleles, ctype):
    ''' extract consequences '''
    result = {}
    consequences = annotation[ctype + '_consequences']
    for a in alleles:
        result[a] = [conseq for conseq in consequences if conseq['variant_allele'] == a]

    return result
    

pp = pprint.PrettyPrinter(indent=4)
            
# alleles = ['C', 'G']
alleles = ['A']
with open('multiple_conseq.json') as fh:
    for line in fh:
        annotation = json.loads(line.rstrip())
        tconseq = extract_consequences(annotation, alleles, 'transcript')
        pp.pprint(tconseq)


parser = VepJsonParser("ranks.txt", verbose=False)
pp.pprint(parser.get_consequence_map())
# parser.set_annotation()