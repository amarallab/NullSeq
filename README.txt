#############################
#                           #
#          Nullseq          #
#                           #
#############################


Parameters to run nullseq:

$ python nullseq.py [-h] [-m] [-n number] [-l length] [--TT TransTable] [--seq SEQ] [--AA AA] [--GC GC] [-o O]

-m                  Use exact primary AA sequence from input file

-n number           Number of random sequences generated            < Default: 1 >

-l length           Length of random sequence (# of AAs)            < Optional >
                    * Must be specified if running script with
                      --AA csv input.
                    * Otherwise uses length of input seqeunce

--TT TransTable     NCBI translation table                          < Default: 11 >

--seq SEQ           Path to FASTA file with nulceotide sequences    < Optional >
                    * either SEQ or AA must be specified
                    * Sequences must include Stop and Start codons

--AA AA             Path to FASTA file with AA sequence or          < Optional >
                    csv with AA usage probabilities
                    * either SEQ or AA must be specified

--GC GC             GC content of random sequence                   < Optional >
                    * Must be specified if running script using
                      --AA input.
                    * Otherwise uses GC of input nucleotide seqeunce
                    * [0,100]

-o O                Output file Path                                < Default: NullSeq_Output.fasta >

Example:

1. Generating Random Sequences From a Known Nuclotide Sequence:

    $ python nullseq.py -n 10 --seq test_Nseq.fasta

        --GC and -l can be used to specify GC content and
        length of random sequence

        The GC content will be targeted at 50%

        The random sequences will be 500 amino acids in length

    $ python nullseq.py -n 10 --seq test_Nseq.fasta --GC 50 -l 500

        -m indicates the use of the exact primary amino acid sequence
        of the input sequence

    $ python nullseq.py -m -n 10 --seq test_Nseq.fasta --GC 50


2. Generating Random Sequences from AA Usage Probabilities:

    * --GC -l must be specified

    $ python nullseq.py -n 10 --AA AAUsage.csv --GC 50 -l 500

        Generates random sequence according to primary amino acid
        usage probabilities in AAUsage.csv


3. Generating Random Sequences from Primary AA Sequences:

    * --GC must be specified

    $ python nullseq.py -n 10 --AA test_AAseq.fasta --GC 50

        -l can be used to specify length of random sequence

    $ python nullseq.py -n 10 --AA test_AAseq.fasta --GC 50 -l 500

        -m indicates the use of the exact primary amino acid sequence
        of the input sequence

    $ python nullseq.py -m -n 10 --AA test_AAseq.fasta --GC 50

Please Cite:
Liu SS, Hockenberry AJ, Lancichinetti A, Jewett MC, Amaral LAN (2016) NullSeq: A Tool for Generating Random Coding Sequences with Desired Amino Acid and GC Contents. PLoS Comput Biol 12(11): e1005184. doi:10.1371/journal.pcbi.1005184

