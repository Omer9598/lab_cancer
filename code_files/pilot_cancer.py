import os

from code_files.interval_analyze import *
from code_files.dict_analyzer import *
from code_files.file_analyzer import *
from code_files.test_scripts import *
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
        create_common_cancer_genes_dict("/Users/dahansarah/PycharmProjects/lab_cancer_new/data_files/BROCA.genes.tsv"))

    if invert:
        # Inverting the file and saving the new path
        file_to_split = invert_reference_genome_haplotype(input_file, save_directory)
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
            path_to_save_interval_plots, invert, chrom_num,
            window_size, error_size)
        update_cancer_variant_dict(interval_list,
                                   common_cancer_variants_dict)
        # Adding the interval coverage of the current chromosome
        chromosome_coverage_dict[chrom_num] = calc_coverage(interval_list,
                                                            chrom_num)

    merge_haplotype_tables(path_to_save_interval_table, chromosome_coverage_dict,
                           window_size, error_size, invert)
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
    file_to_process = input_path
    if inverted:
        # Avoid inverting the file when inverted is True
        pass
    num_of_children, children_filenames = open_and_split_children_files(file_to_process)
    interval_children_list = []
    for i in range(1, num_of_children + 1):
        child_filename = children_filenames[i - 1]
        interval_list = process_child_file(child_filename, reference_type,
                                           window_size, error_size)
        interval_children_list.append(interval_list)
        # Delete the last child file after processing
        os.remove(child_filename)

    shared_interval_list = shared_interval(interval_children_list)
    create_table(shared_interval_list, output_directory_tables, window_size,
                 error_size, inverted)

    # todo - uncomment when finish
    # plot_title = f'chromosome {chromosome_number} interval'
    # plot_interval(shared_interval_list, plot_title,
    #               save_dir=output_directory_plots)

    return shared_interval_list


def analyze_single_chromosome(chromosome_data_file, chrom_num, reference, output_directory):
    """
    This function will analyze a single chromosome,
    create interval tables for:
    Window sizes 20, 30, 50
    Error sizes 5%, 10%, 15%
    for all the permutations of the values above, we will create
    two plots - one with error percent vs coverage and another with
    window size vs coverage
    """
    window_size_dict = {}
    error_coverage_dict = {}
    window_coverage_dict = {}

    # Creating all the interval tables
    for i in range(2):
        for window_size in [20, 30, 50]:
            for error in [0.95, 0.9, 0.85]:
                interval_list = single_chromosome_process(
                    chromosome_data_file, reference, output_directory, output_directory,
                    i, chrom_num, window_size, window_size * error)
                # Updating the window size and error dicts
                key = f'chrom_{chrom_num}_window_{window_size}_error_{error}'
                window_size_dict[key] = [window_size, error]
                error_coverage_dict[key] = calc_coverage(interval_list, chrom_num)

                # Add window size and coverage to window_coverage_dict
                if window_size not in window_coverage_dict:
                    window_coverage_dict[window_size] = []
                window_coverage_dict[window_size].append(calc_coverage(interval_list, chrom_num))

        plot_chromosome_analyze(chrom_num, error_coverage_dict, window_size_dict, i)


def plot_chromosome_analyze(chrom_num, error_coverage_dict, window_size_dict, inverted):
    # Plot for Error Percent vs Coverage
    plt.figure(figsize=(12, 6))
    # Iterate through window sizes and plot lines for each
    for window_size in [20, 30, 50]:
        window_data = [(values[1], error_coverage_dict[key])
                       for key, values in window_size_dict.items()
                       if values[0] == window_size]
        error_percents, coverages = zip(*window_data)
        plt.plot(error_percents, coverages, marker='o', label=f'Window Size {window_size}')
    plt.title(f'Error Percent vs Coverage for Chromosome {chrom_num}')
    plt.xlabel('Error Percent')
    plt.ylabel('Coverage')
    plt.legend()
    plt.grid(True)
    plot_path = "temp_script/error_coverage_chr{}inverted{}.png".format(chrom_num, inverted)
    plt.savefig(plot_path)


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


def single_run(file, output_dir1, output_dir2, chrom, inverted, window_size, error_size):
    new_file = preprocess_file(file, output_dir1)
    if inverted:
        new_file = invert_reference_genome_haplotype(new_file, output_dir1)
    single_chromosome_process(new_file, 'parent', output_dir2, output_dir2, inverted, chrom, window_size, error_size)


if __name__ == '__main__':

    # single_run('/Users]/dahansarah/PycharmProjects/lab_cancer_new/tests/family1/chromosomes/chromosome_8.txt',
    #            '/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1',
    #            '/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1/chrom_8_analyse',
    #            8, 1, 20, 16)
    #
    # print('good')
    #
    # single_run('/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1/chromosomes/chromosome_8.txt',
    #            '/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1',
    #            '/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1/chrom_8_analyse',
    #            8, 0, 20, 16)
    #
    # print('good')
    #
    # check_right_coverage("/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1/real.shared.tsv",
    #                          f"/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1/chrom_8_analyse/table_8_window_20_error_16_inverted_False.txt",
    #                          f"/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1/chrom_8_analyse/table_8_window_20_error_16_inverted_True.txt",
    #                          f"/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1/all_coverage_results/chrom_8_coverage_results/chrom_8_window_20_error_16")


    # pre_file = preprocess_file("/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1/chromosomes/chromosome_1.txt", "/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1")
    # print(pre_file)
    # new_file = invert_reference_genome_haplotype('/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1/processed_chromosome_1.txt',
    #                                              '/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1')
    # print(new_file)
    # single_chromosome_process('/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1/inverted_processed_chromosome_1.txt', "parent",
    #                           "/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1/chrom_1_analyse",
    #                           "/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1/chrom_1_analyse",
    #                           1, 1, 20, 17)

    # check_right_coverage("/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1/real.shared.tsv",
    #                          f"/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1/chrom_1_analyse/table_1_window_20_error_17_inverted_False.txt",
    #                          f"/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1/chrom_1_analyse/table_1_window_20_error_17_inverted_True.txt",
    #                          f"/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1/all_coverage_results/chrom_1_coverage_results/chrom_1_window_20_error_17")

    # single_chromosome_process("tests/family1/inverted_chromosomes/chromosome_22.txt", "parent",
    #                           "tests/family1/chrom_22_analyze", "tests/family1/chrom_22_analyze",
    #                           1, 22, 20, 16)
    # main()
    # create_tables_and_plots("test_data_files/simulated.family.genotypes.tsv", "parent",
    #                         "tests/family1", 1, 50, 48)

    # errors = {16: 20, 18: 20, 19: 20, 40: 50, 45: 50, 48: 50, 90: 100, 95: 100,
    #           98: 100, 145: 150, 190: 200}

    # for error, window in errors.items():
    #     check_right_coverage("tests/family1/real.shared.tsv",
    #                          f"tests/family1/chrom_1_analyze/table_1_window_{window}_error_{error}_inverted_False.txt",
    #                          f"tests/family1/chrom_1_analyze/table_1_window_{window}_error_{error}_inverted_True.txt",
    #                          f"tests/family1/all_coverage_results/chrom_1_coverage_results/chrom_1_window_{window}_error_{error}")
    #
    # merge_coverage_files("tests/family1/all_coverage_results/chrom_22_coverage_results",
    #                      "tests/family1/all_coverage_results/chrom_22_coverage_results/chrom_22_merged")

    # for error, window in errors.items():
    #     single_chromosome_process("tests/family1/chromosomes/chromosome_1.txt", "parent",
    #                               "tests/family1/chrom_1_analyze", "tests/family1/chrom_1_analyze",
    #                               0, 1, window, error)
    # plot_f1_score('/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1/all_coverage_results/chrom_22_coverage_results/chrom_22_merged')
    # plot_coverage('/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1/all_coverage_results/chrom_22_coverage_results/chrom_22_merged')

    # single_chromosome_process("tests/family1/inverted_chromosomes/chromosome_22.txt", "parent",
    #                           "tests/family1/chrom_22_analyze", "tests/family1/chrom_22_analyze",
    #                           1, 22, window, error)

    # print(calculate_coverage({'1': [(0, 200)]},
    #                    {'1': [(10, 20), (45, 300)]}))

    # preprocess_file("test_data_files/GP_3siblings.HET.tab", "test_data_files")
    # create_tables_and_plots("/Users/dahansarah/PycharmProjects/lab_cancer_new/test_data_files/simulated.family.genotypes.tsv",
    #                         "parent",
    #                         "/Users/dahansarah/PycharmProjects/lab_cancer_new/tests/family1",
    #                         1, 10, 8)

    # single_chromosome_process("tests/family1/chromosomes/chromosome_22.txt", "parent",
    #                           "tests/family1/inverted_tables", "tests/family1/inverted_plots",
    #                           1, 22, 20, 18)
    pass