from code_files.interval_analyze import *
from code_files.dict_analyzer import *
from code_files.file_analyzer import *


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
        file_to_split = invert_reference_genome_haplotype(input_file,
                                                          save_directory)
        path_to_save_interval_table = save_directory + "/inverted_interval_tables"
        path_to_save_interval_plots = save_directory + "/inverted_interval_plots"
    else:
        file_to_split = input_file
        path_to_save_interval_table = save_directory + "/interval_tables"
        path_to_save_interval_plots = save_directory + "/interval_plots"

    # splitting the file to separate chromosome files
    split_file_to_chromosomes(file_to_split,
                              save_directory + "/chromosomes")

    # creating interval table for each chromosome
    for chrom_num in range(1, 23):
        num_of_children = open_and_split_children_files(save_directory +
                                                        f"/chromosomes/chromosome_{chrom_num}.txt")

        interval_children_list = []
        for i in range(1, num_of_children + 1):
            file_path = f'{"child_"}{i}{".txt"}'
            interval_list = process_child_file(file_path, reference_type)
            interval_children_list.append(interval_list)

        shared_interval_list = shared_interval(interval_children_list)
        plot_title = f'Chromosome {chrom_num} interval plot'
        plot_interval(shared_interval_list, plot_title, save_dir=path_to_save_interval_plots)

        create_table(shared_interval_list, path_to_save_interval_table)

        update_cancer_variant_dict(shared_interval_list,
                                   common_cancer_variants_dict)

    merge_haplotype_tables(path_to_save_interval_table, common_cancer_variants_dict)


def main():
    # Analyzing family1 - call the function after preprocess
    create_tables_and_plots(r"data_files/preprocess.genotypes.generation1.txt",
                            PARENT_REFERENCE, r"family1", False)

    # Analyzing family2 - after preprocess
    create_tables_and_plots(r"data_files/HR3.genotypes.tab", SIBLING_REFERENCE,
                            "family2", True)


if __name__ == '__main__':
    main()

# shared_interval_list = [
#     {'chromosome': '2', 'start': 215590369, 'end': 3000000000000},  # True
#     {'chromosome': '17', 'start': 59756546, 'end': 8000000000000},  # False
#     {'chromosome': '11', 'start': 108093558, 'end': 108239826},  # True
#     {'chromosome': '13', 'start': 32889616, 'end': 32973808}  # False
# ]

    # # running on a single chromosome - 22
    # num_of_children = open_and_split_children_files(r"family2/chromosomes/chromosome_22.txt")
    #
    # interval_children_list = []
    # for i in range(1, num_of_children + 1):
    #     file_path = f'{"child_"}{i}{".txt"}'
    #     interval_list = process_child_file(file_path, SIBLING_REFERENCE)
    #     interval_children_list.append(interval_list)
    #
    # shared_interval_list = shared_interval(interval_children_list)
    # create_table(shared_interval_list, r"family2/inverted_interval_tables")
    #
    # plot_title = f'Chromosome 22 interval plot'
    # plot_interval(shared_interval_list, plot_title,
    #               save_dir="family2/interval_plots")
