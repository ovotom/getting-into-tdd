from unittest import TestCase, main
from nucleotide import *
# import tdd as tdd

class TestTableOfSubsequences(TestCase):
    def test_pairs_2_length(self):
        seq = 'aaa'
        ans = {'aa': 2}
        ret = tableOfSubsequences(dna =seq, length = 2)
        print(ret)
        self.assertEqual(ans, ret)

    def test_longer_seq(self):
        seq = 'aacact'
        ans = { "aa":1, "ac":2, "ca":1, "ct":1 }
        ret = tableOfSubsequences(dna =seq, length = 2)
        self.assertEqual(ans, ret)
    pass

class Testget_list_of_string_slices(TestCase):
    def test_slice(self):
        dna = 'aa'
        self.assertEqual(get_list_of_string_slices(dna,length =2),['aa'])
    def test_longer_slice(self):
        dna = 'aaa'
        self.assertEqual(get_list_of_string_slices(dna,length =2),['aa','aa'])


if __name__ == '__main__':
    main()