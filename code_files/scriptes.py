from pilot_cancer import *


def analyze_single_chromosome(chromosome_data_file, chrom_num,
                              reference):
    """
    This function will analyze a single chromosome,
    create interval tables for:
    Window sizes 20, 30, 50
    Error sizes 5%, 10%, 15%
    for all the permutations of the values above, we will create
    a line in the final plot
    """
    for i in range(2):


