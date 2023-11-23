from itertools import islice
import plotly.express as px
import pandas as pd


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
            # Condition 2: If the right side of the parent is equal to the left side of the child, add 2 to the list
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

        for next_position in islice(keys_iterator, 50):
            if count_10 == 50:
                break
            else:
                count_10 += 1
                next_haplotype = my_dict[next_position][-1]  # Get the haplotype for the next position
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
        if value[3] > 45:
            filtered_dict[key] = value
    return filtered_dict


def plot_data(data_dict_to_plot, plot_title):
    """
    This function will plot the given dict:
    x axis will be the chromosome position
    y axis will be the haplotype value of the current variant (position)
    """
    # Extract relevant information
    positions = list(data_dict_to_plot.keys())
    haplotypes = [entry[2] for entry in data_dict_to_plot.values()]

    # Create a DataFrame for Plotly Express
    data = {"Chromosome Position": positions, "Haplotype": haplotypes}
    df = pd.DataFrame(data)

    # Create an interactive scatter plot
    fig = px.scatter(df, x="Chromosome Position", y="Haplotype",
                     labels={'Chromosome Position': 'Chromosome Position',
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
    child_1_to_plot_dict = process_dict(child_1_dict)
    child_2_to_plot_dict = process_dict(child_2_dict)

    # plotting the dicts
    plot_data(child_1_to_plot_dict, 'child 1, chromosome 13')
    plot_data(child_2_to_plot_dict, 'child 2, chromosome 13')
    print('x')


if __name__ == '__main__':
    main()
