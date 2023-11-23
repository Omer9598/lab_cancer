from itertools import islice
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
# import matplotlib.pyplot as plt


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
    This function will open the file and split it to the 2 files - mother
    compared to child 1 and mother compared to child 2
    """
    # Open the input and output files
    with open(file, 'r') as infile, open('child_1.txt', 'w') as child_1, open(
            'child_2.txt', 'w') as child_2:
        # Iterate through each line in the input file
        for line in infile:
            # Split the line into columns
            columns = line.strip().split('\t')

            # Extract the desired columns for child_1
            child_1_line = f"{columns[1]}\t{columns[4]}\t{columns[5]}\n"

            # Write the line to the child_1 file
            child_1.write(child_1_line)

            # Extract the desired columns for child_2
            child_2_line = f"{columns[1]}\t{columns[4]}\t{columns[6]}\n"

            # Write the line to the child_2 file
            child_2.write(child_2_line)


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

        for next_position in islice(keys_iterator, 20):
            if count_10 == 20:
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
        if value[3] > 18:
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


def shared_interval(interval_list_1, interval_list_2):
    """
    This function will create a new list containing intervals that are shared
    in both lists given, according to the haplotype
    """
    shared_intervals = []

    for interval_1 in interval_list_1:
        for interval_2 in interval_list_2:
            if (
                    interval_1["haplotype"] == interval_2["haplotype"]
                    and interval_1["start"] <= interval_2["end"]
                    and interval_1["end"] >= interval_2["start"]
            ):
                # Calculate the intersection of intervals
                start = max(interval_1["start"], interval_2["start"])
                end = min(interval_1["end"], interval_2["end"])
                shared_intervals.append({"start": start, "end": end,
                                         "haplotype": interval_1["haplotype"]})

    return shared_intervals


def plot_interval(interval_list, plot_title):
    """
    This function plots intervals as straight lines, where each interval is
    represented by a line starting from the "start" to "end" keys on the x-axis,
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
        data["Interval"].extend([f'Interval {i+1}', f'Interval {i+1}'])

    # Create a DataFrame
    df = pd.DataFrame(data)

    # Create an interactive line plot using Plotly Express
    fig = px.line(df, x="Start", y="Haplotype", color="Interval",
                  labels={'Start': 'Chromosome Position',
                          'Haplotype': 'Haplotype'},
                  title=plot_title)

    # Show the plot in an HTML window
    fig.show()


def main():
    open_and_split_file(r"HR1.ch13.phased.tsv")

    # creating and filtering the child dicts
    child_1_dict = create_and_filter_dictionary('child_1.txt')
    child_2_dict = create_and_filter_dictionary('child_2.txt')

    # adding the haplotype and score values to the keys
    child_1_windowed_dict = process_dict(child_1_dict)
    child_2_windowed_dict = process_dict(child_2_dict)

    child_1_interval_list = create_intervals(child_1_windowed_dict)
    child_2_interval_list = create_intervals(child_2_windowed_dict)

    shared_interval_list = shared_interval(child_1_interval_list,
                                           child_2_interval_list)

    plot_interval(shared_interval_list,
                  "shared haplotypes of child 1 and child 2")

    # plotting the dicts
    # plot_data(child_1_windowed_dict, child_2_windowed_dict, 'children 1
    # and 2,'
    #                                                         ' chromosome
    #                                                         13 with filter')

    # plot_data(child_1_dict, child_2_dict, 'children 1 and 2,'
    #                                       ' chromosome 13 without filter')


if __name__ == '__main__':
    main()
