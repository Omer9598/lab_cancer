import os

from interval_analyze import *
from dict_analyzer import *
from file_analyzer import *
import sys


def process_child_file(file_path, reference_type, window_size, error_size):
    """
    Process a child file and return the processed dictionary.
    """
    child_dict = create_and_filter_dictionary(file_path, reference_type)
    windowed_dict = process_dict(child_dict, reference_type, window_size, error_size)
    interval_list = create_intervals(windowed_dict)
    return interval_list


def calc_coverage(interval_list, chrom_num):
    """
    This function will calculate the coverage of an interval list given
    The coverage is the number of base pairs in all intervals, divided
    by the whole chromosome
    """
    chromosome_sizes = {
        1: 249250621, 2: 243199373, 3: 198022430, 4: 191154276,
        5: 180915260, 6: 171115067, 7: 159138663, 8: 146364022, 9: 141213431,
        10: 135534747, 11: 135006516, 12: 133851895, 13: 115169878,
        14: 107349540, 15: 102531392, 16: 90354753, 17: 81195210,
        18: 78077248, 19: 59128983, 20: 63025520, 21: 48129895,
        22: 51304566
    }
    interval_coverage_sum = 0
    for interval in interval_list:
        interval_coverage_sum += interval["end"] - interval["start"]

    return interval_coverage_sum / chromosome_sizes[chrom_num]


def create_tables_and_plots(input_file, reference_type, save_directory, invert,
                            window_size, error_size):
    """
    This function will create interval table from the given family.txt file
    """
    # Creating the common cancer variants dict
    common_cancer_variants_dict = (
        create_common_cancer_genes_dict(r"data_files/BROCA.genes.tsv"))

    if invert:
        # Inverting the file and saving the new path
        file_to_split = invert_reference_genome_haplotype(input_file, "data_files")
        path_to_save_interval_table = save_directory + "/inverted_interval_tables"
        path_to_save_interval_plots = save_directory + "/inverted_interval_plots"
    else:
        file_to_split = input_file
        path_to_save_interval_table = save_directory + "/interval_tables"
        path_to_save_interval_plots = save_directory + "/interval_plots"

    # Splitting the file to separate chromosome files
    split_file_to_chromosomes(file_to_split,
                              save_directory + "/chromosomes")

    chromosome_coverage_dict = {}
    # creating interval table for each chromosome
    for chrom_num in range(1, 23):
        interval_list = single_chromosome_process(
            save_directory + f"/chromosomes/chromosome_{chrom_num}.txt",
            reference_type, path_to_save_interval_table,
            path_to_save_interval_plots, False, chrom_num,
            window_size, error_size)
        update_cancer_variant_dict(interval_list,
                                   common_cancer_variants_dict)
        # Adding the interval coverage of the current chromosome
        chromosome_coverage_dict[chrom_num] = calc_coverage(interval_list,
                                                            chrom_num)

    merge_haplotype_tables(path_to_save_interval_table, chromosome_coverage_dict)
    write_common_genes_to_file(path_to_save_interval_table,
                               common_cancer_variants_dict)


def single_chromosome_process(input_path, reference_type,
                              output_directory_tables,
                              output_directory_plots,
                              inverted,
                              chromosome_number,
                              window_size, error_size):
    """
    This function will process a single chromosome given, creating an
    interval table, and a plot.
    """
    if inverted:
        file_to_process = (
            invert_reference_genome_haplotype(input_path, output_directory_tables))
    else:
        file_to_process = input_path
    num_of_children = open_and_split_children_files(file_to_process, chromosome_number)

    interval_children_list = []
    for i in range(1, num_of_children + 1):
        file_path = f"child_{i}_chrom_{chromosome_number}.txt"
        interval_list = process_child_file(file_path, reference_type,
                                           window_size, error_size)
        interval_children_list.append(interval_list)
        # Delete temp child file
        os.remove(file_path)

    shared_interval_list = shared_interval(interval_children_list)
    create_table(shared_interval_list, output_directory_tables)

    plot_title = f'chromosome {chromosome_number} interval plot'
    # plot_interval(shared_interval_list, plot_title,
    #               save_dir=output_directory_plots)

    return shared_interval_list


def main():
    args = sys.argv
    if len(args) not in [7, 9]:
        print("Invalid number of arguments.\n "
              "for all chromosomes: \n"
              "input_file reference inverted(0 or 1) window_size error_size"
              " output_directory \n"
              "for a single chromosome: \n"
              "input_file reference inverted(0 or 1) window_size error_size"
              " output_directory_tables output_directory_plots"
              " chromosome_number ")
        sys.exit(1)

    input_file = args[1]
    reference = args[2]
    inverted = bool(int(args[3]))
    window_size = int(args[4])
    error_size = int(args[5])

    # Whole genome process
    if len(args) == 7:
        output_directory = args[6]
        # Running the code on the given arguments
        create_tables_and_plots(input_file, reference, output_directory,
                                inverted, window_size, error_size)

    # One chromosome process
    if len(args) == 9:
        output_directory_tables = args[6]
        output_directory_plots = args[7]
        chromosome_number = int(args[8])

        single_chromosome_process(input_file,
                                  reference,
                                  output_directory_tables,
                                  output_directory_plots,
                                  inverted,
                                  chromosome_number,
                                  window_size, error_size)


if __name__ == '__main__':
    main()

    # single_chromosome_process("family3/chromosomes/chromosome_17.txt",
    #                           "parent",
    #                           "temp", "temp",
    #                           True, 17)

    # # analyzing chromosome 13 only
    # num_of_children = open_and_split_children_files("data_files/HR7.chr13.genotypes.tab")
    #
    # interval_children_list = []
    # for i in range(1, num_of_children + 1):
    #     file_path = f'{"child_"}{i}{".txt"}'
    #     interval_list = process_child_file(file_path, SIBLING_REFERENCE)
    #     interval_children_list.append(interval_list)
    #
    # shared_interval_list = shared_interval(interval_children_list)
    # # plot_title = f'Chromosome {13} interval plot'
    # # plot_interval(shared_interval_list, plot_title, save_dir="family3")
    #
    # create_table(shared_interval_list, "family3")

    # # Analyzing family1 - call the function after preprocess
    # create_tables_and_plots(r"data_files/preprocess.genotypes.generation1.txt",
    #                         PARENT_REFERENCE, r"family1", True)

    # # Analyzing family2 - after preprocess
    # create_tables_and_plots(r"data_files/HR3.genotypes.tab", SIBLING_REFERENCE,
    #                         "family2", True)

# num_of_children = open_and_split_children_files(save_directory +
#                                                 f"/chromosomes/chromosome_{chrom_num}.txt")
#
# interval_children_list = []
# for i in range(1, num_of_children + 1):
#     file_path = f'{"child_"}{i}{".txt"}'
#     interval_list = process_child_file(file_path, reference_type)
#     interval_children_list.append(interval_list)
#
# shared_interval_list = shared_interval(interval_children_list)
# # plot_title = f'Chromosome {chrom_num} interval plot'
# # plot_interval(shared_interval_list, plot_title, save_dir=path_to_save_interval_plots)
#
# create_table(shared_interval_list, path_to_save_interval_table)
# update_cancer_variant_dict(shared_interval_list,
#                            common_cancer_variants_dict)
