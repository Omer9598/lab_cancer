from itertools import islice
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
from tabulate import tabulate
from interval_analyze import *
from dict_analyzer import *


def open_and_split_file(file):
    """
    This function will open the file and split it into n files.
    Each file will contain column[1], column[4], and subsequent columns from
    5 to the last column.
    The number of child files will be determined based on the available
    columns.
    """
    # Open the input file
    with open(file, 'r') as infile:
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
                header_line = '\t'.join([header_columns[1], header_columns[4],
                                         header_columns[child_num + 4]])
                child_file.write(header_line + '\n')

                # Iterate through each line in the input file
                infile.seek(0)  # Reset file pointer to the beginning
                next(infile)  # Skip the header line
                for line in infile:
                    # Split the line into columns
                    columns = line.strip().split('\t')

                    # Extract the desired columns for the current child file
                    child_line = f"{columns[1]}\t{columns[4]}\t" \
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


def create_table(data_list, chromosome_num):
    # Add a new column for chromosome numbers to each dictionary in the list
    for entry in data_list:
        entry['chromosome'] = chromosome_num

    # Reorder the keys to make "chromosome" the first column
    data_list_reordered = [
        {k: entry[k] for k in ['chromosome', 'start', 'end', 'haplotype']} for
        entry in data_list]

    # Convert the reordered list of dictionaries to a table
    table = tabulate(data_list_reordered, headers="keys", tablefmt="pretty")

    # Specify the file path
    file_path = 'table.txt'

    # Write the table to a text file
    with open(file_path, 'w') as file:
        file.write(table)


def process_child_file(file_path):
    """
    Process a child file and return the processed dictionary.
    """
    child_dict = create_and_filter_dictionary(file_path)
    windowed_dict = process_dict(child_dict)
    interval_list = create_intervals(windowed_dict)
    return interval_list


def main():
    # Opening the file
    num_of_children = open_and_split_file(r"HR1.ch13.phased.tsv")
    interval_children_list = []
    for i in range(1, num_of_children + 1):
        file_path = f'{"child_"}{i}{".txt"}'
        interval_list = process_child_file(file_path)
        interval_children_list.append(interval_list)

    shared_interval_list = shared_interval(interval_children_list)

    plot_interval(shared_interval_list,
                  "shared haplotypes of child 1 and child 2")

    create_table(shared_interval_list, 13)


if __name__ == '__main__':
    main()
