'''
File: NullSeq_Functions.py
Author: Sophia Liu
Date: May 16, 2016
Description: This script will generate random nucleotide sequences with a
             given GC content either based on a nucleotide or amino acid
             sequence or amino acid usage probabitliy
'''

import collections
from Bio.Seq import Seq
from Bio.Data import CodonTable
import random
import math
import bisect
import numpy as np
import pandas as pd
import os.path

def Codons_for_AA(n):
    ''' Determines the codons used to code an amino acid given a codon table

    Parameters
    ----------
    n : int
        the codon table index according to NCBI

    Returns
    -------
    CodonforAA : dictionary
        Key : str
            Amino acids in Single letter notation (caps)
        Values : str
            Three letter codon (in caps)
    '''


    GivenCodonTable = CodonTable.unambiguous_dna_by_id[n]
    nucleotides = ['A', 'T', 'C', 'G']
    CodonforAA = {}
    for first in nucleotides:
        for second in nucleotides:
            for third in nucleotides:
                Codon = first + second + third
                if Codon not in CodonTable.unambiguous_dna_by_id[n].stop_codons:
                    if GivenCodonTable.forward_table[Codon] in CodonforAA.keys():
                        CodonforAA[GivenCodonTable.forward_table[Codon]].append(Codon)
                    else:
                        CodonforAA[GivenCodonTable.forward_table[Codon]] = [Codon]
                else:
                    pass
    return CodonforAA

def AAUsage_from_csv(filename):
    ''' Gets AA usage dictionary from csv filename

    Parameters:
    -----------
    filename : str
        path to csv file

    Returns:
    --------
    AAusage : dictionary
        Keys : str
            AA abrv
        Values : float
            P(AA)

    AAUsagedf : pandas dataframe
    '''
    AAUsageDf = pd.read_csv(filename)
    return AAUsageDf

def get_Equal_AA_Prob(n):
    '''Determines the Equal AA frequency according to a codon table

    Parameters:
    -----------
    n :
        codon table

    Returns:
    --------
    AAdict : dictionary
        Keys : str
            amino acid
        Values : float
            Probability of using the amino acid
    '''

    AAdict = Codons_for_AA(n)
    EqualAAdf = pd.DataFrame()
    AAlist = AAdict.keys()
    freqlist = [1/len(AAlist)]*len(AAlist)
    EqualAAdf['AA'] = AAlist
    EqualAAdf['f'] = freqlist
    return EqualAAdf

def df_to_dict(df, n):
    ''' Translates pandas dataframe to dictionary

    Parameters:
    ----------
    df : pandas DataFrame

    Returns : dictionary
        Keys : str
            AA
        Values : float
            P(AA)
    '''
    
    AAUsageDict = {}
    valuelist = df.values
    for item in valuelist:
        AAUsageDict[item[0]] = item[1]
    
    AAUsageDict = clean_AA_dict(AAUsageDict, n)

    return AAUsageDict

def cdf(weights):
    ''' Translates weights to cdf

    Parameters:
    -----------
    weights : list
        list of weights

    Returns:
    -------
    result : list
        cdf of the weights
    '''

    total = sum(weights)
    result = []
    cumsum = 0
    for w in weights:
        cumsum += w
        result.append(cumsum / total)
    return result

def choice(population, weights):
    '''Chooses from the population with probablities
        determined by the weights

    Parameters:
    -----------
    population : list
        list of item to choose from

    weights : list
        weights of each item in population
    '''

    assert len(population) == len(weights)
    cdf_vals = cdf(weights)
    x = random.random()
    idx = bisect.bisect(cdf_vals, x)
    return population[idx]

def get_Random_AA_Seq(AAprob, length):
    '''Creates a random amino acid sequence based on a given AA frequency

    Parameters:
    -----------
    AAprob : dictionary
        Keys : str
            Amino acid
        Values : float
            Probability of using the corresponding amino acid

    length : int
        Number of Amino acids in the sequence

    Returns:
    --------
    AAseq : str
        random AA sequence drawn from a specified amino acid Probability and length

    Notes:
    ------
    Start/Stop codon not included
    '''

    AAlist = []
    Problist = []

    for AA, Prob in AAprob.items():
        AAlist.append(AA)
        Problist.append(Prob)

    AAseq = ''

    for i in range(length):
        AAseq += choice(AAlist, Problist)

    return AAseq

def get_AA_Count(seq, n, Nuc=True):
    '''Determines the number of each Amino acids used in the sequence

    Parameters:
    -----------
    seq : str
        The sequence to be evaluated
        if Nuc is True:
            seq is the nucleotide sequence
        if Nuc is False:
            seq is the amino acid sequence

    n : int
        Codon table to be used

    Nuc : bool (optional)
        True:
            evaulates seq as a nucleotide sequence
        False:
            evaulates seq as a Amino acid sequence

    Returns:
    --------
    AAUsageDict : dictionary
        Keys : str
            Amino acid symbol
        Values : float
            the number of the corresponding AA in the sequence

    Notes:
    -----
    Does not include first AA
    '''

    if Nuc:
        # print(seq)
        AA = str(Seq(seq).translate(table=n)[1:-1])
    else:
        AA = seq[1:]
    AAUsageDict = collections.Counter(AA)

    return AAUsageDict

def get_AA_Freq(seq, n, nucleotide=True):
    '''Determines the frequency of each Amino acids used in the sequence

    Parameters:
    -----------
    seq : str
        The sequence to be evaluated
        if nucleotide is True:
            seq is the nucleotide sequence
        if nucleotide is False:
            seq is the amino acid sequence

    n : int
        Codon table to be used

    nucleotide : bool (optional)
        True:
            evaulates seq as a nucleotide sequence
        False:
            evaulates seq as a Amino acid sequence

    Returns:
    --------
    AAUsageDict : dictionary
        Keys : str
            Amino acid symbolget_AA_Count
        Values : float
            the frequency of the corresponding amino acid
    '''

    AAUsageDict = get_AA_Count(seq, n, Nuc=nucleotide)
    for AA in AAUsageDict.keys():
        if nucleotide:
            AAUsageDict[AA] = AAUsageDict[AA]/((len(seq)-6)/3)
        else:
            AAUsageDict[AA] = AAUsageDict[AA]/(len(seq)-1)

    AAUsageDict = clean_AA_dict(AAUsageDict, n)
    return AAUsageDict



def clean_AA_dict(AAUsageDict, n):
    '''
    Recently added function to add in amino acids with zero frequency and ensure only usage of the 
    standard alphabet.

    --AJH
    '''
    CforAA = Codons_for_AA(n) 
    for aa in CforAA.keys():
        if aa not in AAUsageDict.keys():
            AAUsageDict[aa] = 0.0
    
    assert np.isclose(1.0, np.sum(list(AAUsageDict.values()))), 'The amino acid frequencies from your file do not sum to 1. Exiting'
    
    assert set(list(CforAA.keys())) == set(list(AAUsageDict.keys())), 'Problematic amino acid usage frequency dictionary, '\
            'I found {} amino acids instead of the {} included in the indicated translation table (which, for the record is {})'.format(len(AAUsageDict.keys()), len(CforAA.keys()), n)
    
    return AAUsageDict


def get_ATCG_Count(seq):
    '''Determines the number of each nucleotide in the sequence

    Parameters:
    -----------
    seq : str
        nucleotide sequence to be evaluated

    Returns:
    --------
    ATCGcount : dictionary
        Keys : str
            A, T, C, G
        Values : float
            number of corresponding nucleotide in the sequence

    Note:
    -----
    input sequence should include start and stop Codons
    calculated values do not take into accout the start and stop codons
    '''

    ATCGcount = collections.Counter([seq[x] for x in range(3, len(seq)-3)])

    return ATCGcount

def get_GC_Count(seq):
    '''Determines the number GC nucleotides in the sequence

    Parameters:
    -----------
    seq : str
        nucleotide sequence to be evaluated

    Returns:
    --------
    GCcount : float
        the number of GC nucleotides in the sequence

    Note:
    -----
    input sequence should include start and stop Codons
    calculated values do not take into accout the start and stop codons
    '''

    ATCGcount = get_ATCG_Count(seq)
    GCcount = ATCGcount['G'] + ATCGcount['C']

    return GCcount

def get_GC_Freq(seq):
    '''Determines the frequency of GC nucleotide in the sequence

    Parameters:
    -----------
    seq : str
        nucleotide sequence to be evaluated

    Returns:
    --------
    GCfreq : float
        the frequency of GC nucleotides in the sequence

    Note:
    -----
    input sequence should include start and stop Codons
    calculated values do not take into accout the start and stop codons
    '''

    GCfreq = get_GC_Count(seq)/(len(seq)-6)

    return GCfreq

def GC_Content_for_Codons(n):
    ''' Determines the number of G/C nucletides for all codons in the
    codon table

    Parameters
    ----------

    n : int
        The codon table index according to NCBI

    Returns
    -------
    CodonGCContent : dictionary
        Key : str
            Three letter codon (in caps)
        Values : int
            Number of G/C nucleotides in the codon
    '''

    CodonforAA = Codons_for_AA(n)
    CodonGCContent = collections.defaultdict(int)
    for AA in CodonforAA:
        for codon in CodonforAA[AA]:
            codonusage = collections.Counter(codon)
            CodonGCContent[codon] = codonusage['G'] + codonusage['C']

    return CodonGCContent

def minmax_GC_for_AA(n):
    ''' Determines the maximum or minimum number
    of G/C nucleotide for a given AAsequence

    Parameters:
    -----------
    n : int
        NCBI translation table

    Returns:
    --------
    minGCDict : dictionary
        Keys : str
            AA
        Values : int
            smallest number of G/C nucleotides

    maxGCDict : dictionary
        Keys : str
            AA
        Values : int
            largest number of G/C nucleotides
    '''

    CodonforAA = Codons_for_AA(n)
    CodonGCContent = GC_Content_for_Codons(n)
    maxGCDict = {}
    minGCDict = {}
    for AA in CodonforAA.keys():
        templist = [CodonGCContent[codon] for codon in CodonforAA[AA]]
        maxGCDict[AA] = max(templist)
        minGCDict[AA] = min(templist)

    return (minGCDict, maxGCDict)

def get_maxmin_GC_count(AAfreq, n):
    '''Determines the maximum and minimum GC ratio for
    specified AA compostion

    Parameters:
    -----------
    AAfreq : dictionary
        Keys : str
            AA
        Values : float
            frequency of AA in sequence

    n : int
        NCBI translation table

    Returns:
    --------
    low : float
        lowest possible GC ratio

    high : float
        highest possible GC ratio
    '''

    (minGCDict, maxGCDict) = minmax_GC_for_AA(n)
    high = 0
    low = 0
    for AA in AAfreq:
        high += AAfreq[AA]*maxGCDict[AA]/3
        low += AAfreq[AA]*minGCDict[AA]/3

    return (low, high)

def evaluate_possibility(AAfreq, GC, l, n):
    '''Determines whether is it possible to obtain the GC content
    for the given AA compostion

    Parameters:
    -----------
    AAfreq : dictionary
        Keys : str
            AA
        Values : float
            frequency of AA in sequence

    GC : float
        desired GC content [0, 100]

    l : int
        length of sequence

    n : int
        NCBI translation table

    Returns:
    --------
    bool
        True : if the GC content is obtainable
        False: if the GC content is not obtainable
    '''

    (L, H) = get_maxmin_GC_count(AAfreq, n)
    if GC/100 > L and GC/100 < H:
        return True
    else:
        return False

def get_beta(given, length, CforAA, setAAProb, n):
    '''Determines the value of beta given the GC content of the sequence

    Parameters:
    -----------
    given : float
        GC content can be speficied (%)

    length : float
        number of AA in the sequnece

    CforAA : dictionary
        the codon table being used
        key : str
            AA
        Values : list
            list of codons that can be used to code for the AA

    setAAProb : dictionary
        the specified AA usage probability
        key : str
            AAprob
        Values : float
            P(AA)

    Returns:
    -------
    beta : float
        the constant in the equation P = exp(-beta*N)/Z
    '''

    m = given*length*3

    rl = -20
    rr = 20

    avl = compute_average_GC(rl, length, CforAA, setAAProb, n) - m
    avr = compute_average_GC(rr, length, CforAA, setAAProb, n) - m

    if avl*avr > 0:
        if avl < 0:
            return rl
        else:
            return rr
    else:

        rm = (rl+rr)/2
        avv = compute_average_GC(rm, length, CforAA, setAAProb, n) - m

        while math.fabs(avv) > 1e-7:
            if avv*avl > 0:
                rl = rm
            else:
                rr = rm
            rm = 0.5*(rl+rr)
            avv = compute_average_GC(rm, length, CforAA, setAAProb, n) - m

    return rm

def Probability_Given_Beta(beta, CforAA, n):
    ''' Gives the probability of using each codon given beta.
        The probabilites are determined using the bolzmann distribtion
        where P = exp(-beta*N)/Z
        N is the number of G/C nucleotides in the codon

    Parameters:
    ----------
    beta : float
        The constant beta used in the calulations of the probability

    n : int
        The codon table index according to NCBI

    Returns:
    -------
    probs : dictionary
        Key : str
            Codons  (caps)
        Values : float
            The probability of using the codon (P(codon|AA))
    '''

    probs = {}
    CodonforAA = CforAA
    CodonGC = GC_Content_for_Codons(n)

    for AA in CodonforAA:
        Z = sum([math.exp(-beta*CodonGC[codon]) for codon in CodonforAA[AA]])
        for codon in CodonforAA[AA]:
            probs[codon] = math.exp(-beta*CodonGC[codon])/Z
    return probs

def compute_average_GC(b, length, CforAA, setAAProb, n, given='beta'):
    '''Determines the average GC of the given sequence given codon
       probablities or a given beta

    Parameters:
    -----------
    b : float or dictionary
        When calulating for a given beta b is a float.
        When calulating for a codon distribution is a dictionary
            key : str
                codon
            Values : float
                P(codon)

    length : float
        number of AA in the sequnece

    CforAA : dictionary
        the codon table being used
        key : str
            AA
        Values : list
            list of codons that can be used to code for the AA

    setAAProb : dictionary
        the specified AA usage probability
        key : str
            AAprob
        Values : float
            P(AA)

    given : str (optional)
        When given is 'beta' the fuction calculates average GC for beta.
        When given is 'probability' the the fuction calculates average GC
        for a given codon usage proability

    Returns:
    --------
    sum(GCnumber) : float
        The expected number GC nucleotides for the given sequence given
        the conditions
    '''

    CodonGC = GC_Content_for_Codons(n)
    GCnumber = []

    if given == 'beta':
        Probability = Probability_Given_Beta(b, CforAA, n)
        CodonforAA = CforAA
        AAcount = {}
        for AA in setAAProb:
            AAcount[AA] = setAAProb[AA]*length

        for AA in CodonforAA:
            GCnumber.append(AAcount[AA]*sum([Probability[codon]*CodonGC[codon] \
                                             for codon in CodonforAA[AA]]))

    return sum(GCnumber)

def Normalize_Probabilites(rawprob, n):
    ''' Normalizes the Probabilities of synonymus codons so that
        it sums to 1

        Parameters:
        -----------
        rawprob : dictionary
            Keys : str
                three letter codon CodonFreq
            Values :
                raw Probabilities associated with using that codonProb
        n : int
            codon table

        Returns:
        --------
        NormProb : dictionary
            Keys : str
                AA
            Values : tuple
            [0] : list of codons
            [1] : associated Normalized Probabilities
    '''

    CodonforAA = Codons_for_AA(n)
    NormProb = {}

    for AA, Codonlist in CodonforAA.items():
        rawproblist = []
        for codon in Codonlist:
            rawproblist.append(rawprob[codon])
        s = sum(rawproblist)
        NormProb[AA] = (Codonlist, [i/s for i in rawproblist])

    return NormProb

def parse_fastafile(filename):
    ''' parses fasta files into list of sequence

    Parameters:
    -----------
    filename : str
        path to fasta file

    Returns:
    --------
    seqlist : list
        list of sequences in the fasta file
    '''

    with open(filename, 'r') as f:
        lines = f.readlines()

    seqlist = []
    for l in lines:
        if l[0] == '>':
            pass
        else:
            seqlist.append(l.split('\n')[0])

    return seqlist[0]

def get_Random_Nuc_Seq(AAfreq=None, AAseq=None, nseq=None, GC=None, length=None, n=None):
    '''Gets list of Random nucleotide sequence given specified Parameters

    Parameters:
    -----------
    AAfreq : dictionary (opt)
        The amino acid frequency of the random sequence will be drawn from
        If not specified, equal AA usage is assumed
        Keys : str
            Amino acid
        Values : float
            amino acid usage frequency

    AAseq : str (opt)
        when specified, the random sequences have exactally the same AA usage as the sequence
        with start and stop codons

    nseq : int
        the number of random sequences generated

    GC : float
        GC ratio of the sequence

    length : int
        length of the sequence (number of amino acids)

    n : int
        codon table

    Returns:
    --------
    SeqListBeta : list
        list of randomly generated nucleotide sequences
    '''

    CodonforAA = Codons_for_AA(n)
    def CC(tupl):
        return 'ATG' + ''.join(tupl) + 'TGG'
    AAsequence = None

    if AAfreq is None and AAseq is None:
        AAFrequency = get_Equal_AA_Prob(n)
    elif AAseq is None and AAfreq is not None:
        AAFrequency = AAfreq
    else:
        AAsequence = AAseq
        AAFrequency = get_AA_Freq(AAseq, n, nucleotide=False)
        length = len(AAseq)-1

    b = get_beta(GC, length, CodonforAA, AAFrequency, n)
    codonProb = Probability_Given_Beta(b, CodonforAA, n)
    betaProb = Normalize_Probabilites(codonProb, n)

    if AAsequence is None:
        pofC = []
        pofAA = []
        codonlist = []

        for AA, clistp in betaProb.items():
            pofAA = pofAA + [AAFrequency[AA]]*len(CodonforAA[AA])
            pofC = pofC + clistp[1]
            codonlist = codonlist + clistp[0]

        prob = np.array(pofC)*np.array(pofAA)
        prob = prob/sum(prob)
        SeqListBeta = np.random.choice(codonlist, size=(nseq, length), replace=True, p=prob)
        SeqListBeta = np.apply_along_axis(CC, 1, np.array(SeqListBeta))

        return SeqListBeta

    else:

        SeqListBeta = [['ATG']*nseq]
        for AA in AAsequence[1:]:
            RandCforAA = np.random.choice(betaProb[AA][0], size=(nseq),
                                          replace=True, p=betaProb[AA][1])
            SeqListBeta.append(list(RandCforAA))
        SeqListBeta = np.apply_along_axis(CC, 0, np.array(SeqListBeta[1:]))
        return SeqListBeta




