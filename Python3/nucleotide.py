"""
Background
----------
DNA is made up of two strands that are sequences of 4 types of bases,
each of which is assigned a single character code.

    A = Adenine
    C = Cytosine
    G = Guanine
    T = Thymine

DNA sequences can be represented as sequences of these 4 letters (n.b. this is
an oversimplification). These represent the bases along just one strand as the
other can be easily inferred from the fact that the bases on the two strands
pair up with A<->T and C<->G.)

Our understanding of what is in a DNA sequence is built-up incrementally
by looking at repeated internal sequences. This is built up from the statistics
and indexes of short repeats.

The Problem
-----------
Our simple problem is, given a string that represents a DNA sequence, calculate
a table of how many times each sub-sequence of length-N nucleotides occurs in it.
For example, the table for the sub-sequences of length 2 in the tiny DNA
sequence 'aacact' would be:

    { "aa":1, "ac":2, "ca":1, "ct":1 }
"""

# From https://benchmarksgame-team.pages.debian.net/benchmarksgame/download/knucleotide-input.txt
# THREE Homo sapiens frequency, extracted for workshop to duck i/o issue.
DNA_Sequence='''
aacacttcaccaggtatcgtgaaggctcaagattacccagagaacctttgcaatataaga
atatgtatgcagcattaccctaagtaattatattctttttctgactcaaagtgacaagcc
ctagtgtatattaaatcggtatatttgggaaattcctcaaactatcctaatcaggtagcc
atgaaagtgatcaaaaaagttcgtacttataccatacatgaattctggccaagtaaaaaa
tagattgcgcaaaattcgtaccttaagtctctcgccaagatattaggatcctattactca
tatcgtgtttttctttattgccgccatccccggagtatctcacccatccttctcttaaag
gcctaatattacctatgcaaataaacatatattgttgaaaattgagaacctgatcgtgat
tcttatgtgtaccatatgtatagtaatcacgcgactatatagtgctttagtatcgcccgt
gggtgagtgaatattctgggctagcgtgagatagtttcttgtcctaatatttttcagatc
gaatagcttctatttttgtgtttattgacatatgtcgaaactccttactcagtgaaagtc
atgaccagatccacgaacaatcttcggaatcagtctcgttttacggcggaatcttgagtc
taacttatatcccgtcgcttactttctaacaccccttatgtatttttaaaattacgttta
ttcgaacgtacttggcggaagcgttattttttgaagtaagttacattgggcagactcttg
acattttcgatacgactttctttcatccatcacaggactcgttcgtattgatatcagaag
ctcgtgatgattagttgtcttctttaccaatactttgaggcctattctgcgaaatttttg
ttgccctgcgaacttcacataccaaggaacacctcgcaacatgccttcatatccatcgtt
cattgtaattcttacacaatgaatcctaagtaattacatccctgcgtaaaagatggtagg
ggcactgaggatatattaccaagcatttagttatgagtaatcagcaatgtttcttgtatt
aagttctctaaaatagttacatcgtaatgttatctcgggttccgcgaataaacgagatag
attcattatatatggccctaagcaaaaacctcctcgtattctgttggtaattagaatcac
acaatacgggttgagatattaattatttgtagtacgaagagatataaaaagatgaacaat
tactcaagtcaagatgtatacgggatttataataaaaatcgggtagagatctgctttgca
attcagacgtgccactaaatcgtaatatgtcgcgttacatcagaaagggtaactattatt
aattaataaagggcttaatcactacatattagatcttatccgatagtcttatctattcgt
tgtatttttaagcggttctaattcagtcattatatcagtgctccgagttctttattattg
ttttaaggatgacaaaatgcctcttgttataacgctgggagaagcagactaagagtcgga
gcagttggtagaatgaggctgcaaaagacggtctcgacgaatggacagactttactaaac
caatgaaagacagaagtagagcaaagtctgaagtggtatcagcttaattatgacaaccct
taatacttccctttcgccgaatactggcgtggaaaggttttaaaagtcgaagtagttaga
ggcatctctcgctcataaataggtagactactcgcaatccaatgtgactatgtaatactg
ggaacatcagtccgcgatgcagcgtgtttatcaaccgtccccactcgcctggggagacat
gagaccacccccgtggggattattagtccgcagtaatcgactcttgacaatccttttcga
ttatgtcatagcaatttacgacagttcagcgaagtgactactcggcgaaatggtattact
aaagcattcgaacccacatgaatgtgattcttggcaatttctaatccactaaagcttttc
cgttgaatctggttgtagatatttatataagttcactaattaagatcacggtagtatatt
gatagtgatgtctttgcaagaggttggccgaggaatttacggattctctattgatacaat
ttgtctggcttataactcttaaggctgaaccaggcgtttttagacgacttgatcagctgt
tagaatggtttggactccctctttcatgtcagtaacatttcagccgttattgttacgata
tgcttgaacaatattgatctaccacacacccatagtatattttataggtcatgctgttac
ctacgagcatggtattccacttcccattcaatgagtattcaacatcactagcctcagaga
tgatgacccacctctaataacgtcacgttgcggccatgtgaaacctgaacttgagtagac
gatatcaagcgctttaaattgcatataacatttgagggtaaagctaagcggatgctttat
ataatcaatactcaataataagatttgattgcattttagagttatgacacgacatagttc
actaacgagttactattcccagatctagactgaagtactgatcgagacgatccttacgtc
gatgatcgttagttatcgacttaggtcgggtctctagcggtattggtacttaaccggaca
ctatactaataacccatgatcaaagcataacagaatacagacgataatttcgccaacata
tatgtacagaccccaagcatgagaagctcattgaaagctatcattgaagtcccgctcaca
atgtgtcttttccagacggtttaactggttcccgggagtcctggagtttcgacttacata
aatggaaacaatgtattttgctaatttatctatagcgtcatttggaccaatacagaatat
tatgttgcctagtaatccactataacccgcaagtgctgatagaaaatttttagacgattt
ataaatgccccaagtatccctcccgtgaatcctccgttatactaattagtattcgttcat
acgtataccgcgcatatatgaacatttggcgataaggcgcgtgaattgttacgtgacaga
gatagcagtttcttgtgatatggttaacagacgtacatgaagggaaactttatatctata
gtgatgcttccgtagaaataccgccactggtctgccaatgatgaagtatgtagctttagg
tttgtactatgaggctttcgtttgtttgcagagtataacagttgcgagtgaaaaaccgac
gaatttatactaatacgctttcactattggctacaaaatagggaagagtttcaatcatga
gagggagtatatggatgctttgtagctaaaggtagaacgtatgtatatgctgccgttcat
tcttgaaagatacataagcgataagttacgacaattataagcaacatccctaccttcgta
acgatttcactgttactgcgcttgaaatacactatggggctattggcggagagaagcaga
tcgcgccgagcatatacgagacctataatgttgatgatagagaaggcgtctgaattgata
catcgaagtacactttctttcgtagtatctctcgtcctctttctatctccggacacaaga
attaagttatatatatagagtcttaccaatcatgttgaatcctgattctcagagttcttt
ggcgggccttgtgatgactgagaaacaatgcaatattgctccaaatttcctaagcaaatt
ctcggttatgttatgttatcagcaaagcgttacgttatgttatttaaatctggaatgacg
gagcgaagttcttatgtcggtgtgggaataattcttttgaagacagcactccttaaataa
tatcgctccgtgtttgtatttatcgaatgggtctgtaaccttgcacaagcaaatcggtgg
tgtatatatcggataacaattaatacgatgttcatagtgacagtatactgatcgagtcct
ctaaagtcaattacctcacttaacaatctcattgatgttgtgtcattcccggtatcgccc
gtagtatgtgctctgattgaccgagtgtgaaccaaggaacatctactaatgcctttgtta
ggtaagatctctctgaattccttcgtgccaacttaaaacattatcaaaatttcttctact
tggattaactacttttacgagcatggcaaattcccctgtggaagacggttcattattatc
ggaaaccttatagaaattgcgtgttgactgaaattagatttttattgtaagagttgcatc
tttgcgattcctctggtctagcttccaatgaacagtcctcccttctattcgacatcgggt
ccttcgtacatgtctttgcgatgtaataattaggttcggagtgtggccttaatgggtgca
actaggaatacaacgcaaatttgctgacatgatagcaaatcggtatgccggcaccaaaac
gtgctccttgcttagcttgtgaatgagactcagtagttaaataaatccatatctgcaatc
gattccacaggtattgtccactatctttgaactactctaagagatacaagcttagctgag
accgaggtgtatatgactacgctgatatctgtaaggtaccaatgcaggcaaagtatgcga
gaagctaataccggctgtttccagctttataagattaaaatttggctgtcctggcggcct
cagaattgttctatcgtaatcagttggttcattaattagctaagtacgaggtacaactta
tctgtcccagaacagctccacaagtttttttacagccgaaacccctgtgtgaatcttaat
atccaagcgcgttatctgattagagtttacaactcagtattttatcagtacgttttgttt
ccaacattacccggtatgacaaaatgacgccacgtgtcgaataatggtctgaccaatgta
ggaagtgaaaagataaatat
'''.replace( '\n', '' )

def r(dict, el):
    if type(el) == str:
        if el not in dict:
            dict[el] = 0
        dict[el] += 1
    return dict

def get_list_of_string_slices(dna, length):
    substrings = []
    for i in range(len(dna) - length + 1):
        s = dna[i:i+length]
        substrings.append(s)
    return substrings

def tableOfSubsequences( *, dna, length ):

    """
    Compute the frequency distribution of subsequences of a given length in a
    stretch of DNA.
    :param dna: The DNA sequence we are interested in.
    :param length: The length of the subsequences to generate.
    :return: A dictionary mapping subsequence to count (i.e. a bag)
    """
    from functools import reduce
    # pairs = list(filter(lambda x: len(x) == 2, splitter(dna)) )
    # print(pairs)
    ##should be able to do with reduce - come back to this
    substrings = get_list_of_string_slices(dna,length =length)
    return reduce(r,substrings ,{})
    # return reduce(r,dna, {}) #functools.reduce(func, iter, [initial_value])

if __name__ == '__main__':
    print (tableOfSubsequences(dna = DNA_Sequence,length= 4 ))
