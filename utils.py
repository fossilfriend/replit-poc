

from collections import Counter, OrderedDict
import csv

def qw(s, returnTuple=False):
    '''
    mimics perl's qw function
    usage: qw('a b c') will yield ['a','b','c']
    returnTuple: return a tuple if true, otherwise return list
    '''
    if returnTuple:
        return tuple(s.split())
    else:
        return s.split()

'''
utils for parsing VEP JSON output
'''
# pylint: disable=line-too-long,invalid-name


class VepJsonParser(object):
    ''' class to organize utils for parsing VEP JSON output '''

    def __init__(self, rankingFileName, verbose=True):
        self.__codingConsequences = qw('synonymous_variant missense_variant inframe_insertion inframe_deletion stop_gained stop_lost stop_retained_variant start_lost frameshift_variant coding_sequence_variant')
        self._verbose = verbose
        self._consequenceRankings = self.parse_ranking_file(rankingFileName)
        self._annotation = None


    def is_coding_consequence(self, conseqStr):
        ''' check term against list of coding consequences and return 
        True if found '''
        matches = [value for value in conseqStr.split(',') 
                       if value in self.__codingConsequences]
        if len(matches > 0):
            return True
        else:
            return False


    def parse_ranking_file(self, fileName):
        ''' parse ranking file and save as dictionary lookup '''
        if self._verbose:
            1
            # warning("Parsing ranking file:,", fileName)

        result = OrderedDict()
        with open(fileName, 'r') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                conseq = row['consequence']
                result[conseq]  = {}
                result[conseq]['coding'] = self.is_coding_consequence(conseq)
                result[conseq]['rank'] = row['adsp_ranking']
                result[conseq]['impact'] = row['adsp_impact']
                
        return result

    
    def find_matching_consequence(self, terms):
        ''' match list of consequences against those in the
        ranking file and return ranking info associated with match;
        throw error if not found '''
        if len(terms) == 0:
            return self.get_conseq_rank(terms[0])

        for conseqStr in self._consequenceRankings:
            conseqList = conseqStr.split(',')
            if is_equivalent_list(terms, conseqList):
                return self.get_conseq_rank(conseqStr) 
            else:
                raise IndexError('Consequence combination ' + ','.join(terms)
                                     + ' not found in ranking file.')


    def assign_consequence_rank(self, conseq):
        ''' assemble consequence terms into ordered comma delim string
        and lookup in ranking table; add rank to consequence annotation;
        backup old impact information
        '''
        terms = conseq['consequence_terms']
        matchingConseq = self.find_matching_consequence(terms) 

        conseq['vep_impact'] =  conseq['impact']
        conseq['impact'] = matchingConseq['impact']
        conseq['rank'] = matchingConseq['rank']
        conseq['is_coding'] = matchingConseq['is_coding']

        return conseq
        

    def __verify_annotation(self):
        ''' check that annotation is set '''
        assert self._annotation is not None, \
          "DEBUG - must set value of _annotation in the VEP parser to access it"



    def rank_and_sort_consequences(self):
        ''' takes an array of consequences (too which adsp ranking
        has been applied, and reorders) '''
        pass

    # =========== modifiers ==================
    def set_annotation(self, annotation):
        ''' set the annotation json '''
        self._annotation = annotation

    
    def set(self, key, value):
        ''' set a value to annotation json '''
        self.__verify_annotation()
        self._annotation[key] = value


    # =========== accessors ==================
    def coding_consequences(self):
        ''' return coding variants'''
        return self.__codingConsequences

    def consequence_rank_map(self):
        ''' return consequence rankings '''
        return self._consequenceRankings

    def get_conseq_rank(self, conseq):
        ''' return value from consequence rank map for the specified
        consequence '''
        if conseq in self._consequenceRankings:
            return self._consequenceRankings[conseq]
        else:
            raise IndexError('Consequence ' + conseq + ' not found in ranking file.')


    def get_ranked_consequences(self, alleles, conseqType):
        ''' extract consequences and apply ranking;
        convert from list to list of dicts, keyed on allele'''
        result = {}
        consequences = self._annotation[conseqType + '_consequences']
        for a in alleles:
            result[a] = [self.assign_consequence_rank(conseq) for conseq in consequences
                         if conseq['variant_allele'] == a]

        return result


    def get_frequencies(self):
        ''' extract frequencies from colocated_variants section '''
        self.__verify_annotation()

        cv = self._annotation['colocated_variants']
        if len(cv) > 1:
            raise NotImplementedError("More than one co-located variant; can't extract frequencies")
        return cv[0]['frequencies'] if 'frequencies' in cv[0] else None


    def get(self, key): 
        ''' get the annotation value associated with the key '''
        self.__verify_annotation()

        if key == 'frequencies':
            return self.get_frequencies()
        else:
            return self._annotation[key]


    def get_annotation(self):
        ''' return updated annotation '''
        return self._annotation
            

# ======================
# Functions
# ======================


def get_most_severe_consequence(conseq):
    ''' returns first element of the conseq array '''
    return conseq[0]



# --------------
# helpers

def is_equivalent_list(list1, list2):
    ''' test if two lists contain the same elements;
    order does not matter'''
    if Counter(list1) == Counter(list2):
        return True
    return False
