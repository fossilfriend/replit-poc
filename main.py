
import json
import pprint
import sys


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


a = ['a', 'b', 'c']
b = ['c', 'a', 'b']

