from itertools import islice
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
from tabulate import tabulate

WINDOW_NUM = 20
WINDOW_ERR = 18
CHROMOSOME_NUM = 13


def check_heterozygous_parent(parent):
    if parent == '0|0' or parent == '1|1':
        return False
    return True


def check_homozygous_child(child):
    if child == '0|1' or child == '1|0':
        return False
    return True


def filter_dict(dictionary):
    """
    filtering the dict:
    the mother has to be heterozygous, for us to know which side she inherited
    and the child has to be homozygous, for us to know that he got one of his
    alleles from the mother
    :return: the filtered dictionary
    """
    filtered_dict = {}
    for key, value in dictionary.items():
        # Apply check_heterozygous_parent to the value associated with the key
        result_parent = check_heterozygous_parent(value[0])

        # Apply check_homozygous_child to the value associated with the key
        result_child = check_homozygous_child(value[1])

        # If both functions return True, add the key to the filtered_dict
        if result_parent and result_child:
            filtered_dict[key] = value

    return filtered_dict


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


def create_and_filter_dictionary(file_path):
    """
    This file will create dictionary from the given file, in the following
    format:
    {position: [mother, child], ...}
    """
    result_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            position = int(columns[0])
            values = [columns[1], columns[2]]
            result_dict[position] = values
    return filter_dict(result_dict)


def add_haplotype(my_dict):
    for position, values in my_dict.items():
        left_side_parent, right_side_parent = values[0].split('|')
        left_side_child, right_side_child = values[1].split('|')

        if left_side_parent == left_side_child:
            # Condition 1: If left sides are equal, add 1 to the list
            values.append(1)
        elif right_side_parent == left_side_child:
            # Condition 2: If the right side of the parent is equal to the
            # left side of the child, add 2 to the list
            values.append(2)


def add_confidence(my_dict):
    for position, values in my_dict.items():
        haplotype = values[-1]  # Get the haplotype for the current position
        count = 1
        count_10 = 1

        keys_iterator = iter(my_dict.keys())
        next_position = next(keys_iterator)

        # Skip positions until the current position in the outer loop
        while next_position != position:
            next_position = next(keys_iterator)

        for next_position in islice(keys_iterator, WINDOW_NUM):
            if count_10 == WINDOW_NUM:
                break
            else:
                count_10 += 1
                next_haplotype = my_dict[next_position][
                    -1]  # Get the haplotype for the next position
                # Check if haplotypes match
                if haplotype == next_haplotype:
                    # If they match, increment the count
                    count += 1
        values.append(count)


def process_dict(data_dict):
    """
    This function will process the dicts to be plotted
    """
    add_haplotype(data_dict)
    add_confidence(data_dict)
    filtered_dict = filter_low_score(data_dict)
    return filtered_dict


def filter_low_score(data_dict):
    """
    This function will filter the variants (keys in the dict) which their
    score is under 8
    :returns a new filtered dict
    """
    filtered_dict = dict()
    for key, value in data_dict.items():
        if value[3] > WINDOW_ERR:
            filtered_dict[key] = value
    return filtered_dict


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


def create_intervals(haplotype_dict):
    """
    This function will create a list for each child in the following format:
    [(interval number is the index) {start position: , end position: ,
    haplotype: (1 or 2)]
    each interval starts with the position of a variant from haplotype 1 or 2,
    and ends when the next variant is from the opposite haplotype, where a new
    interval will start
    """
    intervals = []
    current_interval = None

    for position, haplotype in haplotype_dict.items():
        if current_interval is None:
            # Start a new interval
            current_interval = {"start": position, "end": position,
                                "haplotype": haplotype[2]}
        elif haplotype[2] == current_interval["haplotype"]:
            # Continue the current interval
            current_interval["end"] = position
        else:
            # Start a new interval as haplotype changed
            intervals.append({"start": current_interval["start"],
                              "end": current_interval["end"],
                              "haplotype": current_interval["haplotype"]})
            current_interval = {"start": position, "end": position,
                                "haplotype": haplotype[2]}

    # Add the last interval
    if current_interval is not None:
        intervals.append({"start": current_interval["start"],
                          "end": current_interval["end"],
                          "haplotype": current_interval["haplotype"]})

    return intervals


def shared_interval(interval_lists):
    """
    This function will create a new list containing intervals that are shared
    in all the lists given, according to the haplotype
    """
    # Initialize shared_intervals with the intervals from the first list
    shared_intervals = interval_lists[0]

    # Iterate through the remaining lists
    for interval_list in interval_lists[1:]:
        # A temporary list to store shared intervals for the current list
        temp_shared_intervals = []

        # Iterate through each interval in the current list
        for interval_1 in shared_intervals:
            for interval_2 in interval_list:
                if (
                        interval_1["haplotype"] == interval_2["haplotype"]
                        and interval_1["start"] <= interval_2["end"]
                        and interval_1["end"] >= interval_2["start"]
                ):
                    # Calculate the intersection of intervals
                    start = max(interval_1["start"], interval_2["start"])
                    end = min(interval_1["end"], interval_2["end"])
                    temp_shared_intervals.append({"start": start, "end": end,
                                                  "haplotype": interval_1[
                                                      "haplotype"]})

        # Update shared_intervals with the current shared intervals
        shared_intervals = temp_shared_intervals

    return shared_intervals


def plot_interval(interval_list, plot_title):
    """
    This function plots intervals as straight lines, where each interval is
    represented by a line starting from the "start" to "end" keys on the
    x-axis,
    and the height determined by the "haplotype" key on the y-axis.
    """
    # Create a DataFrame for Plotly Express
    data = {"Start": [], "End": [], "Haplotype": [], "Interval": []}

    for i, interval in enumerate(interval_list):
        start_position = interval["start"]
        end_position = interval["end"]
        haplotype = interval["haplotype"]

        # Append data to the DataFrame
        data["Start"].extend([start_position, end_position])
        data["End"].extend([start_position, end_position])
        data["Haplotype"].extend([haplotype, haplotype])
        data["Interval"].extend([f'Interval {i + 1}', f'Interval {i + 1}'])

    # Create a DataFrame
    df = pd.DataFrame(data)

    # Create an interactive line plot using Plotly Express
    fig = px.line(df, x="Start", y="Haplotype", color="Interval",
                  labels={'Start': 'Chromosome Position',
                          'Haplotype': 'Haplotype'},
                  title=plot_title)

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

    create_table(shared_interval_list, CHROMOSOME_NUM)


if __name__ == '__main__':
    main()
