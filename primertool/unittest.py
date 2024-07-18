import unittest
from primertool import functions
from primertool import exceptions as ex


class PrimertoolTest(unittest.TestCase):
    def test_query_ucsc_database(self):
        # SETUP
        query = f'SELECT chrom, strand, name2, exonCount, cdsStart, cdsEnd, exonStarts, exonEnds ' \
                f'FROM refGene WHERE name="NM_000451"'
        expected_result = [('chrX', '+', 'SHOX', 6, 630897, 644636, b'624343,630465,634617,640820,640998,644390,', b'624602,631174,634826,640878,641087,646823,'),
                           ('chrY', '+', 'SHOX', 6, 630897, 644636, b'624343,630465,634617,640820,640998,644390,', b'624602,631174,634826,640878,641087,646823,')]
        # TEST
        result = functions.query_ucsc_database('hg38', query)
        self.assertEqual(result, expected_result, 'Yielded unexpected result')

    def test_get_gene_information(self):
        # SETUP
        nm_number = 'NM_000451'
        expected_result = {'nm_number': 'NM_000451',
                           'chromosome': 'chrY',
                           'strand': '+',
                           'name': 'SHOX',
                           'exoncount': 6,
                           'cds_start': 630897,
                           'cds_end': 644636,
                           'exon_starts': [624343, 630465, 634617, 640820, 640998, 644390],
                           'exon_ends': [624602, 631174, 634826, 640878, 641087, 646823]}
        # TEST
        result = functions.get_gene_information('hg38', nm_number)
        self.assertEqual(result, expected_result, 'Yielded unexpected result')


if __name__ == '__main__':
    unittest.main()
