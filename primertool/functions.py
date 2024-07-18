"""
This module contains static auxiliary functions that are used in the primertool package.
"""
from typing import Tuple
import genomepy
import re
import requests
import hgvs.parser
import hgvs.assemblymapper
import hgvs.variantmapper
import hgvs.sequencevariant
import hgvs.location
import hgvs.dataproviders.uta
import hgvs.exceptions
import logging
from . import exceptions
from .insilicopcr import InSilicoPCR
from .ucsc_database import query as ucsc_query

logger = logging.getLogger(__package__)


def remove_whitespaces(func: callable) -> callable:
    """Decorator to remove whitespaces from all string arguments before passing them to the decorated function.

    Args:
        func (callable): The function to decorate
    Returns:
        The decorated function
    """

    def wrapper(*args, **kwargs):
        args = [arg.replace(' ', '') if isinstance(arg, str) else arg for arg in args]
        kwargs = {key: value.replace(' ', '') if isinstance(value, str) else value for key, value in kwargs.items()}
        return func(*args, **kwargs)

    return wrapper


def get_gene_information(genome_assembly: str, nm_number: str) -> dict:
    """Retrieve gene information from RefSeq database for a given NM number.

    Args:
        genome_assembly (str): The genome assembly to use ('hg38' or 'hg19')
        nm_number (str): The NM number of the gene to retrieve information for
    Returns:
        A dictionary containing the gene information
    """
    logger.info(f'Collecting gene information from RefSeq database for {nm_number}.')
    query = f'SELECT chrom, strand, name2, exonCount, cdsStart, cdsEnd, exonStarts, exonEnds ' \
            f'FROM refGene WHERE name="{nm_number}"'

    query_results = ucsc_query(genome_assembly, query, local=True)

    if not query_results:
        raise exceptions.PrimertoolInputError(f'Could not find gene information for {nm_number} in RefSeq database')

    if len(query_results) > 1:
        for entry in query_results:
            if len(entry[0]) < 6:
                result = entry
    else:
        result = query_results[0]

    gene = dict(nm_number=nm_number,
                chromosome=result[0],
                strand=result[1],
                name=result[2],
                exoncount=result[3],
                cds_start=result[4],
                cds_end=result[5],
                exon_starts=[int(x) for x in result[6].decode("utf-8").split(',') if x],
                exon_ends=[int(x) for x in result[7].decode("utf-8").split(',') if x]
                )
    return gene


def calculate_targets(target_start: int, target_end: int, primer_bases: int) -> dict:
    """Defining the sequence start and end positions.

    Primer bases are the number of bases to each side of the target sequence in which primer3 looks for a possible
    primer. These are added to sequence start and end respectively to calculate the sequence start and end positions.
    The size range is a list with the length of the target sequence and target sequence with the primer bases added.
    Primer3 should design primers in this size range.

    Args:
        target_start (int): Start position of the target sequence
        target_end (int): End position of the target sequence
        primer_bases (int): Number of bases to each side of the target sequence in which primer3 looks for a possible primer
    Returns:
        A dictionary containing the sequence start and end positions, the number of bases to each side of the target
        sequence in which primer3 looks for a possible primer, the length of the target sequence, and the size range for
        primer3.
    """
    seq_start = target_start - primer_bases
    seq_end = target_end + primer_bases
    target_base = primer_bases
    target_length = target_end - target_start
    size_range = [target_length, target_length + int(primer_bases / 2)]

    target_info = dict(target_start=target_start,
                       target_end=target_end,
                       seq_start=seq_start,
                       seq_end=seq_end,
                       target_base=target_base,
                       target_length=target_length,
                       size_range=size_range)
    return target_info


def mask_snps(genome: genomepy.Genome, chromosome: str, seq_start: int, seq_end: int, genome_assembly: str) -> str:
    """Mask common SNP positions with an N in the sequence.

    Using the get_snps() function to retrieve common SNPs in the given sequence from UCSC. Extract the genomic sequence
    between the given positions via genomepy and replace the bases at common SNP positions with an N.

    Args:
        genome (genomepy.Genome): genomepy genome object
        chromosome (str): chromosome name as string (e.g. 'chr1')
        seq_start (int): start position of the sequence
        seq_end (int): end position of the sequence
        genome_assembly (str): genome assembly as string (e.g. 'hg38')
    Returns:
        The sequence with common SNPs masked as a string
    """
    sequence = str(genome[chromosome][seq_start:seq_end]).upper()

    snps = get_snps(chromosome, seq_start, seq_end, genome_assembly)
    if snps is None:
        seq_snps = sequence
    else:
        snp_seq = list(sequence)
        for j in snps:
            snp_seq[j - 2] = 'N'
        seq_snps = ''.join(snp_seq)

    return seq_snps


def get_snps(chromosome: str, seq_start: int, seq_end: int, genome_assembly: str) -> list or None:
    """Retrieve common SNPs from the UCSC database.

    Query the UCSC database with the chromosome and position data to find common SNPs in the sequence, which need to
    be masked. SNPs are returned as list.

    Args:
        chromosome (str): chromosome name as string (e.g. 'chr1')
        seq_start (int): start position of the sequence
        seq_end (int): end position of the sequence
        genome_assembly (str): genome assembly as string (e.g. 'hg38')
    Returns:
        List of common SNPs in the sequence
    """
    snpquery = (f"SELECT chromStart "
                f"FROM snp150Common "
                f"WHERE chrom='{chromosome}' AND class='single' AND chromEnd BETWEEN'{seq_start}'AND'{seq_end}'")
    snp = ucsc_query(genome_assembly, snpquery, local=True)
    snps = []
    if snp is not None:
        for idx, i in enumerate(snp):
            snps.append(seq_end - int(str(snp[idx])[1:-2]))

    return snps


def mutalyzer_error_handler(response: dict) -> None or exceptions.PrimertoolInputError:
    """Checks for errors in the mutalyzer response and raises an exception if there is an error.

    Args:
        response (dict): The response from the mutalyzer API
    Raises:
        PrimertoolInputError: If there is an error in the response
    """
    if 'message' in response and 'custom' in response:
        # If the entries at the top level of the response are message and custom, there is a problem with the input
        logger.info(response['message'])

        # Handle infos and errors
        if 'infos' in response['custom']:
            for info in response['custom']['infos']:
                logger.info(f'{info["code"]}: {info["details"]}')  # print all infos
        if 'errors' in response['custom']:
            for error in response['custom']['errors']:
                logger.error(f'{error["code"]}: {error["details"]}')  # print all errors

            error_code = response['custom']['errors'][0]['code']
            error_message = response['custom']['errors'][0]['details']

            if error_code == 'EPARSE':
                raise exceptions.PrimertoolInputError('There is an error in the given mutation', error_code,
                                                      error_message)
            elif error_code == 'ERETR':
                raise exceptions.PrimertoolInputError('The given NM number has an error and could not be found',
                                                      error_code,
                                                      error_message)
            elif error_code == 'ENOINTRON':
                raise exceptions.PrimertoolInputError('The given NM number has an error and could not be found',
                                                      error_code,
                                                      error_message)
            elif error_code == 'ESYNTAXUC':
                raise exceptions.PrimertoolInputError(error_code, error_message)
            else:
                raise exceptions.PrimertoolInputError('There was a problem with the input. ', error_code, error_message)


def parse_mutation(mutation: str) -> hgvs.sequencevariant.SequenceVariant or exceptions.PrimertoolInputError:
    """Parse mutation with hgvs parser.

    Gene names in brackets are removed from the variant (eg: eg: NM_003165.6(STXBP1):c.1702G>A).
    The mutation is then parsed using hgvs.parser and a hgvs tree object (see
    https://hgvs.readthedocs.io/en/stable/key_concepts.html#variant-object-representation) is returned.
    Parses coding and genomic variants.

    Args:
        mutation (str): mutation in HGVS nomenclature
    Raises:
        PrimertoolInputError: If there is an error in the mutation nomenclature
    Returns:
        hgvs.sequencevariant.SequenceVariant
    """
    mutation = re.sub(r'\([^)]*\)', '', mutation)

    try:
        hp = hgvs.parser.Parser()
        hgvs_mutation = hp.parse_hgvs_variant(mutation)
    except hgvs.exceptions.HGVSParseError as msg:
        raise exceptions.PrimertoolInputError(
            f'Could not parse the input. There is a problem with the hgvs nomenclature: {msg}')
    return hgvs_mutation


def split_nm(nm_number: str) -> Tuple[str, int]:
    """Split variant.ac from hgvs object into transcript and version number

    If the transcript number includes a version number, split at '.'. If not, set the version number to 1.

    Args:
        nm_number (str): NM number (e.g. NM_000451.3)
    Returns:
        Tuple of transcript and version number
    """
    nm_split = nm_number.split('.')
    transcript = nm_split[0]

    if len(nm_split) == 1:
        version = 1
    else:
        version = int(nm_split[1])

    return transcript, version


def find_sequence_positions(exon_starts: list, exon_ends: list, exon_count: int, strand: str,
                            mutation_position: hgvs.location.Interval) -> dict:
    """Establish if variant position is in an exon.

    The hgvs.parser object has the start and end position of the variant (in case of an indel rather than a SNP). Using
    these positions and the lists of exon starts and ends, determine if the variant is located in an exon.

    The information is returned as a dictionary. If the gene is based on the - strand, the exon number needs to be
    inverted.

    Args:
        exon_starts (list): List of exon start positions
        exon_ends (list): List of exon end positions
        exon_count (int): Number of exons
        strand (str): Strand of the gene ('+' or '-')
        mutation_position (hgvs.location.Interval): Start and end position of the variant
    Returns:
        Dictionary containing the exon number, start and end position of the variant, the length of the variant, and a
        boolean indicating if the variant is in an exon
    """
    mut_start = int(str(mutation_position.start))
    mut_end = int(str(mutation_position.end))
    exon_number = 0
    exon_len = 0
    is_in_exon = False

    for exon in range(exon_count):
        if exon_starts[exon] <= mut_start and mut_end <= exon_ends[exon]:
            is_in_exon = True
            exon_len = exon_ends[exon] - exon_starts[exon]
            if strand == '-':
                exon_number = exon_count - exon
            else:
                exon_number = exon + 1

    return dict(exon_number=exon_number, mut_start=mut_start, mut_end=mut_end, mut_length=mut_end - mut_start,
                is_in_exon=is_in_exon, exon_len=exon_len)


def filter_unique_primers(primer3_dict: dict) -> Tuple[dict, bool]:
    """Check primer uniqueness and remove all primers with multiple binding sites, so that only the uniquely binding
    ones remain.
    Also returns a flag if all primers were invalid (i.e. not uniquely binding), but only if there was at least one
    primer to begin with.

    Args:
        primer3_dict (dict): Primer3 output dictionary
    Returns:
        Tuple of the filtered primer3_dict and a flag indicating if all primers were invalid
    """
    pre_filter_primer_count = primer3_dict['PRIMER_PAIR_NUM_RETURNED']

    idx = 0
    while idx < primer3_dict['PRIMER_PAIR_NUM_RETURNED']:  # iterate over all generated primers
        forward_primer = primer3_dict[f'PRIMER_LEFT_{idx}_SEQUENCE']
        reverse_primer = primer3_dict[f'PRIMER_RIGHT_{idx}_SEQUENCE']
        if not InSilicoPCR(forward_primer, reverse_primer).is_uniquely_binding():
            # delete primerpair with index idx from primers
            logger.info("Purging primer pair")
            primer3_dict = purge_primer_pair(primer3_dict, idx)
            # don't increment idx because primer idx was removed and replaced by primer idx+1
        else:
            idx = idx + 1

    post_filter_primer_count = primer3_dict['PRIMER_PAIR_NUM_RETURNED']

    # If >0 primers were found, but all were invalid (i.e. not uniquely binding), set flag
    all_primers_invalid_flag = pre_filter_primer_count > 0 and post_filter_primer_count == 0

    return primer3_dict, all_primers_invalid_flag


def purge_primer_pair(primer3_dict: dict, index: int) -> dict:
    """Removes primer pair of given index from primer3_dict and updates the dictionary so that the remaining data stays
    consistent. (I.e. reset indices so enumeration does not have gaps and update)

    Args:
        primer3_dict (dict): Primer3 output dictionary
        index (int): Index of primer pair to remove
    Returns:
        Primer3 output dictionary with the primer pair of the given index removed
    """
    # Remove primerpair with given index (remove every item from the dict where the index appears in the key)
    for key, value in primer3_dict.copy().items():
        if re.match(f'.*{index}.*', key):
            del primer3_dict[key]

    # Reduce index of all primers with larger index than the one we removed
    for idx in range(index + 1, primer3_dict['PRIMER_PAIR_NUM_RETURNED']):
        for key, value in primer3_dict.copy().items():
            # if the corresponding index is found
            if re.findall(str(idx), key):
                # Rename keys
                new_key = reduce_numbers_in_string(key)
                primer3_dict[new_key] = primer3_dict.pop(key)

    # Count down number of returned primers
    primer3_dict['PRIMER_LEFT_NUM_RETURNED'] -= 1
    primer3_dict['PRIMER_RIGHT_NUM_RETURNED'] -= 1
    primer3_dict['PRIMER_PAIR_NUM_RETURNED'] -= 1

    return primer3_dict


def reduce_numbers_in_string(input_string: str) -> str:
    """Takes an input string and reduces any (positive) integer in it by one.

    Args:
        input_string (str): String to reduce numbers in
    Returns:
        String with all numbers reduced by one
    """
    # Define a regex pattern to match numbers
    pattern = r"\d+"

    # Find all occurrences of the pattern in the input string
    matches = re.findall(pattern, input_string)

    # Decrement each matched number by 1
    decremented_values = [str(int(match) - 1) for match in matches]

    # Replace the original numbers with the decremented values
    for i, match in enumerate(matches):
        input_string = input_string.replace(match, decremented_values[i], 1)

    return input_string


def correct_intronic_variant(response: dict, variant: str) -> str or exceptions.PrimertoolIntronicPositionError:
    """
    Handle the case where the given variant is intronic. If the offset is <= 5, drop the offset. Otherwise,
    raise an exception.

    Args:
        response (dict): The response from the mutalyzer API
        variant (str): The variant to correct
    Raises:
        PrimertoolIntronicPositionError: If the offset is too large to be corrected automatically (>)
    Returns:
        The corrected variant
    """
    offset = response['custom']['corrected_model']['variants'][0]['location']['offset']['value']

    if offset <= 5:
        # remove offset from variant
        variant = re.sub(r'\+\d+', '', variant)

        # correct possible sequence mismatch
        response = requests.get(f'https://mutalyzer.nl/api/normalize/{variant}?only_variants=false')

        if response.json()['custom']['errors'][0]['code'] == 'ESEQUENCEMISMATCH':
            # find the correct base (A, G, C, T)
            replacement_letter = re.search(r'found ([AGCT]) instead',
                                           response.json()['custom']['errors'][0]['details']).group(1)

            # The replacement of the substituent base is only done to avoid the A>A case, which leads to errors later on
            base_partner = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

            # replace the char right before '>' in variant string with the replacement letter
            variant = re.sub(r'([A-Z])>([A-Z])', f'{replacement_letter}>{base_partner[replacement_letter]}', variant)

        return variant

    else:
        raise exceptions.PrimertoolIntronicPositionError(
            f"The given variant is intronic and the offset is too large to be corrected automatically. Please use the "
            f"genomic position instead. "
        )
