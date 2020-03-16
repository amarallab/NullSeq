'''
File: NullSeq.py
Author: Sophia Liu
Date: May 16, 2016
Description: This script will generate random nucleotide sequences with a
             given GC content either based on a nucleotide or amino acid
             sequence or amino acid usage probabitliy
'''

from Bio.Seq import Seq
import NullSeq_Functions as NS
import argparse
import os.path

def main(outfile, number=None, TT=None, AAfile=None,
         seqfile=None, l=None, gc=None, ES=False):

    #Intro
    print()
    print('*************************')
    print('*   Welcome to NullSeq  *')
    print('*************************')

    N = TT
    of = outfile
    if AAfile is not None:
        GC = gc
        if AAfile.split('.')[-1] == 'csv':
            AAUsage = NS.df_to_dict(NS.AAUsage_from_csv(AAfile))
            length = l
            AASequence = None
            operatingmode = 'AA Usage Frequency'
            #B, O, U, X and Z are often used for ambiguous amino acid calls.
            # This warns the user when these are included in AAfile as they are not accepted later.
            bad_AA = ('B', 'O', 'U', 'X', 'Z')
            for AA in bad_letters:
              if AA in AAUsage.keys():
                print(AA,' is an invalid Amino Acid code for this program, please remove.')
                raise ValueError('Invalid Amino Acid code in AA Usage Probabilities')
        else:
            AASequence = NS.parse_fastafile(AAfile)
            AAUsage = NS.get_AA_Freq(AASequence, N, nucleotide=False)
            operatingmode = 'Existing Sequence - AA'
            if l == None:
                length = len(AASequence)-1
            else:
                length = l
            if ES:
                pass
            else:
                AASequence = None

    elif seqfile is not None:
        NucleotideSeq = NS.parse_fastafile(seqfile)
        AASequence = str(Seq(NucleotideSeq).translate(table=N)[0:-1])
        AAUsage = NS.get_AA_Freq(AASequence, N, nucleotide=False)
        operatingmode = 'Exisitng Sequence - Nucleotide'
        if gc is None:
            gc = NS.get_GC_Freq(NucleotideSeq)*100
        else:
            pass

        if l is None:
            length = len(AASequence)-1
        else:
            length = l

        if ES:
            pass
        else:
            AASequence = None


    print('\n***********************************************************************')
    print('Translation Table:\t%s' %N)
    print('Operating Mode:\t\t%s' %operatingmode)
    print('AAFreq :')
    for AA in AAUsage:
        print('\t\t' + AA, '\t', round(AAUsage[AA], 4))

    if AASequence is None:
        pass
    else:
        splitAA = [AASequence[i:i+50] for i in range(0, len(AASequence), 50)]
        print('AAseq :\t\t', splitAA[0])
        for j in range(len(splitAA)):
            if j != 0:
                print('\t\t', splitAA[j])
            else:
                pass

    print('GC content :\t\t%s' % gc)
    print('Length :\t\t%s' % length)
    print('Number of sequence :\t%s' % number)
    print('***********************************************************************')

    if NS.evaluate_possibility(AAUsage, gc, length, N):
        seqlist = NS.get_Random_Nuc_Seq(AAfreq=AAUsage, AAseq=AASequence,
                                        GC=gc/100, length=length, n=N, nseq=number)
        with open(of, 'w') as f:
            for i in range(number):
                f.write('>%s\n' %str(i+1))
                f.write(seqlist[i])
                f.write('\n')
    else:
        print('The GC content selected is not possible for this AA composition.\n \
                Please select a new GC content and retry!')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make Random Sequences')
    parser.add_argument('-m', action='store_true',
                        help='use exact AA mode')
    parser.add_argument('-n', type=int, metavar='number', default=100,
                        help='number of random sequences generated')
    parser.add_argument('-l', type=int, metavar='length', default = 200,
                        help='Length of random sequence')
    parser.add_argument('--TT', type=int, metavar='TransTable', default=11,
                        choices=[1, 2, 3, 4, 5,
                                 6, 9, 10, 11, 12,
                                 13, 14, 16, 21, 22,
                                 23, 24, 25],
                        help='NCBI translation Table')
    parser.add_argument('--seq', help='location of sample sequence')
    parser.add_argument('--AA', help='location of csv file specifiying Amino acid usage frequency')
    parser.add_argument('--GC', type=float,
                        help='GC content of desired sequence [0, 100]. \
                        If none is given, takes GC content from seq')
    parser.add_argument('-o', default='NullSeq_Output.fasta',
                        help='Output file (fasta format)')
    args = parser.parse_args()


    if args.n < 1:
        raise ValueError('Number of generated random sequences must be greater than 0')
    else:
        pass

    if args.AA is not None:
        if os.path.isfile(args.AA):
            if args.AA.split('.')[-1] == 'csv':
                if args.GC is None or args.l is None:
                    raise ValueError('Missing Parameters')
                else:
                    main(args.o, AAfile=args.AA, number=args.n,
                         TT=args.TT, gc=args.GC, l=args.l)
            elif args.AA.split('.')[-1] == 'fasta':
                if args.GC is None:
                    raise ValueError('Missing Parameters')
                else:
                    if args.m:
                        main(args.o, AAfile=args.AA, number=args.n,
                             TT=args.TT, gc=args.GC, l=args.l, ES=args.m)
                    else:
                        main(args.o, AAfile=args.AA, number=args.n,
                             TT=args.TT, gc=args.GC, l=args.l, ES=args.m)
            else:
                raise ValueError('Improper File')
        else:
            raise ValueError('Improper File')
    else:
        if args.seq is None:
            raise ValueError('Input Proper AA usage or nucleotide sequence file')
        else:
            if os.path.isfile(args.seq) and args.seq.split('.')[-1] == 'fasta':
                if args.m:
                    main(args.o, number=args.n, seqfile=args.seq,
                         l=args.l, gc=args.GC, TT=args.TT, ES=args.m)
                else:
                    main(args.o, number=args.n, seqfile=args.seq,
                         l=args.l, gc=args.GC, TT=args.TT, ES=args.m)
            else:
                raise ValueError('Improper Input File')
