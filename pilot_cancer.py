from itertools import islice
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
from tabulate import tabulate
from interval_analyze import *
from dict_analyzer import *


def preprocess_file(input_file_path, output_file_path):
    """
    This function will preprocess the file given:
    Changing "./." to 0/0, or "./1" to 0/1 etc, and save the result in a new
    file.
    """
    # Read the input file and process each line
    with open(input_file_path, 'r') as file:
        lines = file.readlines()

    processed_lines = []

    for line in lines:
        # Split the line into columns
        columns = line.strip().split('\t')

        # Process each column starting from the 5th column
        for i in range(4, len(columns)):
            # Replace "./." with 0/0, "./1" with 0/1, etc.
            columns[i] = columns[i].replace("./.", "0/0")\
                .replace("./1", "0/1").replace("1/.", "1/0")

        # Join the columns back into a line
        processed_line = '\t'.join(columns)
        processed_lines.append(processed_line)

    # Join the processed lines into a string
    result = '\n'.join(processed_lines)

    # Write the processed data to the output file
    with open(output_file_path, 'w') as output_file:
        output_file.write(result)


def open_and_split_file(file_path):
    """
    This function will open the file and split it into n files.
    Each file will contain column[0] column[1], column[4], and subsequent
    columns from 5 to the last column.
    The number of child files will be determined based on the available
    columns.
    """
    # Open the input file
    with open(file_path, 'r') as infile:
        # Get header columns
        header_columns = infile.readline().strip().split('\t')

        # Determine the number of child files based on available columns
        num_children = len(header_columns) - 5

        # Iterate through each child file
        for child_num in range(1, num_children + 1):
            child_filename = f'child_{child_num}.txt'

            # Open the child file for writing
            with open(child_filename, 'w') as child_file:
                # Write the header columns to the child file
                header_line = '\t'.join([header_columns[0], header_columns[1],
                                         header_columns[4],
                                         header_columns[child_num + 4]])
                child_file.write(header_line + '\n')

                # Iterate through each line in the input file
                infile.seek(0)  # Reset file pointer to the beginning
                next(infile)  # Skip the header line
                for line in infile:
                    # Split the line into columns
                    columns = line.strip().split('\t')

                    # Extract the desired columns for the current child file
                    child_line = f"{columns[0]}\t{columns[1]}\t{columns[4]}\t"\
                                 f"{columns[child_num + 4]}\n"

                    # Write the line to the child file
                    child_file.write(child_line)

    return num_children


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
    haplotypes_child_2 = [(entry[2] - 0.015) for entry in
                          child_2_dict.values()]

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
    fig.add_trace(
        px.scatter(df_child_2, x="Chromosome Position", y="Haplotype",
                   color_discrete_sequence=['red'],  # Set color for child_2
                   labels={'Haplotype': 'Haplotype (Child 2)'}).data[0])

    # Adjust legend position
    fig.update_layout(showlegend=True)

    # Show the plot in an HTML window
    fig.show()


def create_table(data_list):
    """
    This function will create the shared haplotype intervals table
    The table will be in a new .txt file, ordered in the following format - each
    row represents an interval, and will be as follows:
    chromosome - the chromosome number
    start - the starting position of the interval in the row's chromosome
    end - the ending position
    haplotype - the haplotype of the current interval
    """
    # Reorder the keys to make "chromosome" the first column
    data_list_reordered = [
        {k: entry[k] for k in ['chromosome', 'start', 'end', 'haplotype']} for
        entry in data_list]

    # Specify the file path
    file_path = 'haplotype.interval.table.txt'

    # Write the data to a text file without any table formatting
    with open(file_path, 'w') as file:
        for entry in data_list_reordered:
            file.write(f"{entry['chromosome']}\t{entry['start']}\t"
                       f"{entry['end']}\t{entry['haplotype']}\n")


def process_child_file(file_path):
    """
    Process a child file and return the processed dictionary.
    """
    child_dict = create_and_filter_dictionary(file_path)
    windowed_dict = process_dict(child_dict)
    interval_list = create_intervals(windowed_dict)
    return interval_list


def main():
    # preprocessing the all chromosome file
    preprocess_file(r"genotypes.generation1.txt",
                    r"preprocess.genotypes.generation1.txt")
    # Opening the file
    num_of_children = open_and_split_file(
        r"preprocess.genotypes.generation1.txt")
    # num_of_children = open_and_split_file(r"HR1.ch13.phased.tsv")
    interval_children_list = []
    for i in range(1, num_of_children + 1):
        file_path = f'{"child_"}{i}{".txt"}'
        interval_list = process_child_file(file_path)
        interval_children_list.append(interval_list)

    shared_interval_list = shared_interval(interval_children_list)

    # plot_interval(shared_interval_list,
    #               "shared haplotypes of child 1 and child 2")

    create_table(shared_interval_list)


if __name__ == '__main__':
    main()
