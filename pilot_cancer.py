from interval_analyze import *
from dict_analyzer import *
from file_analyzer import *


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


def process_child_file(file_path):
    """
    Process a child file and return the processed dictionary.
    """
    child_dict = create_and_filter_dictionary(file_path)
    windowed_dict = process_dict(child_dict)
    interval_list = create_intervals(windowed_dict)
    return interval_list


def create_tables_and_plots(input_file):
    """
    This function will create interval table from the given family.txt file
    """
    invert_reference_genome_haplotype(input_file, r"generation1_inverted.txt")
    # preprocessing the all chromosome file
    preprocess_file(r"inverted.txt",
                    r"preprocess.genotypes.generation1.txt")
    # splitting the file to separate chromosome files
    split_file_to_chromosomes(r"preprocess.genotypes.generation1.txt",
                              r"genotypes_generation1_chromosomes")

    # creating interval table for each chromosome
    for chrom_num in range(1, 23):
        num_of_children = open_and_split_children_files \
            (f"genotypes_generation1_chromosomes/chromosome_{chrom_num}.txt")

        interval_children_list = []
        for i in range(1, num_of_children + 1):
            file_path = f'{"child_"}{i}{".txt"}'
            interval_list = process_child_file(file_path)
            interval_children_list.append(interval_list)

        shared_interval_list = shared_interval(interval_children_list)
        plot_title = f'Chromosome {chrom_num} interval plot'
        plot_interval(shared_interval_list, plot_title, save_dir='interval_plots')

        create_table(shared_interval_list, r"haplotype_interval_tables")

    # merging the tables into a single long table
    merge_haplotype_tables(r"haplotype_interval_tables")


def main():
    create_tables_and_plots(r"all_chromosomes_HR1.txt")


if __name__ == '__main__':
    main()
