from interval_analyze import *
from dict_analyzer import *
from file_analyzer import *
import sys


def plot_data(child_1_dict, child_2_dict, plot_title):
    """
    This function will plot data from two dictionaries (child_1 and child_2).
    The x-axis will be the chromosome position, and the y-axis will be the
    haplotype value of the current variant (position).
    Red dots will be added at positions 1.05 and 1.95 for child_1.
    """
    # Extract relevant information for child_1
    positions_child_1 = list(child_1_dict.keys())
    haplotypes_child_1 = [entry[2] for entry in child_1_dict.values()]

    # Extract relevant information for child_2
    positions_child_2 = list(child_2_dict.keys())
    haplotypes_child_2 = [(entry[2] - 0.015) for entry in child_2_dict.values()]

    # Create DataFrames for Plotly Express
    data_child_1 = {"Chromosome Position": positions_child_1,
                    "Haplotype": haplotypes_child_1}
    data_child_2 = {"Chromosome Position": positions_child_2,
                    "Haplotype": haplotypes_child_2}

    df_child_1 = pd.DataFrame(data_child_1)
    df_child_2 = pd.DataFrame(data_child_2)

    # Create an interactive scatter plot for child_1
    fig = px.scatter(df_child_1, x="Chromosome Position", y="Haplotype",
                     labels={'Chromosome Position': 'Chromosome Position',
                             'Haplotype': 'Haplotype'},
                     color_discrete_sequence=['blue'],  # Set color for child_1
                     title=plot_title)

    # Add data for child_2 to the same plot with red color
    fig.add_trace(px.scatter(df_child_2, x="Chromosome Position", y="Haplotype",
                             color_discrete_sequence=['red'],
                             # Set color for child_2
                             labels={'Haplotype': 'Haplotype (Child 2)'})
                  .data[0])

    # Adjust legend position
    fig.update_layout(showlegend=True)

    # Show the plot in an HTML window
    fig.show()


def process_child_file(file_path, reference_type):
    """
    Process a child file and return the processed dictionary.
    """
    child_dict = create_and_filter_dictionary(file_path, reference_type)
    windowed_dict = process_dict(child_dict, reference_type)
    interval_list = create_intervals(windowed_dict)
    return interval_list


def create_tables_and_plots(input_file, reference_type, save_directory, invert):
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

    # creating interval table for each chromosome
    for chrom_num in range(1, 23):
        interval_list = single_chromosome_process(
            save_directory + f"/chromosomes/chromosome_{chrom_num}.txt",
            reference_type, path_to_save_interval_table,
            path_to_save_interval_plots, False, chrom_num)
        update_cancer_variant_dict(interval_list, common_cancer_variants_dict)

    merge_haplotype_tables(path_to_save_interval_table, common_cancer_variants_dict)


def single_chromosome_process(input_path, reference_type,
                              output_directory_tables,
                              output_directory_plots,
                              inverted,
                              chromosome_number):
    """
    This function will process a single chromosome given.
    """
    if inverted:
        file_to_process = (
            invert_reference_genome_haplotype(input_path, output_directory_tables))
    else:
        file_to_process = input_path
    num_of_children = open_and_split_children_files(file_to_process)

    interval_children_list = []
    for i in range(1, num_of_children + 1):
        file_path = f'{"child_"}{i}{".txt"}'
        interval_list = process_child_file(file_path, reference_type)
        interval_children_list.append(interval_list)

    shared_interval_list = shared_interval(interval_children_list)
    create_table(shared_interval_list, output_directory_tables)

    plot_title = f'chromosome {chromosome_number} interval plot'
    # plot_interval(shared_interval_list, plot_title,
    #               save_dir=output_directory_plots)

    return shared_interval_list


def main():
    if len(sys.argv) not in [5, 7]:
        print("Invalid number of arguments")
        sys.exit(1)

    input_file = sys.argv[1]
    reference = sys.argv[2]
    inverted = bool(int(sys.argv[3]))

    # Whole genome process
    if len(sys.argv) == 5:
        output_directory = sys.argv[4]
        # Running the code on the given arguments
        create_tables_and_plots(input_file, reference, output_directory, inverted)

    # One chromosome process
    if len(sys.argv) == 7:
        output_directory_tables = sys.argv[4]
        output_directory_plots = sys.argv[5]
        chromosome_number = int(sys.argv[6])

        single_chromosome_process(input_file,
                                  reference,
                                  output_directory_tables,
                                  output_directory_plots,
                                  inverted,
                                  chromosome_number)


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
