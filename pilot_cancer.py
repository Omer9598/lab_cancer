from collections import Counter
import matplotlib.pyplot as plt


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


def process_rows_in_batches(data_dict, batch_size=10):
    """
    Creating a dict in the following format:
    key = position in the chromosome
    value = [alleles: (mother, child), haplotype: (1 or 2), Score: (0 to 10)]
    the haplotype is determined by the
    """
    result_dict = {}

    keys = list(data_dict.keys())

    for start in range(len(keys) - batch_size + 1):
        end = start + batch_size
        batch_keys = keys[start:end]

        if len(batch_keys) < 2:
            continue

        result_key = f'{batch_keys[0]}_{batch_keys[-1]}'
        left_parent_counter = Counter({"Score": 0})
        right_parent_counter = Counter({"Score": 0})

        for key in batch_keys:
            values = data_dict[key]

            if len(values) == 2:
                left_side_parent, right_side_parent = values[0].split('|')
                left_side_child, right_side_child = values[1].split('|')

                if left_side_parent == left_side_child:
                    left_parent_counter.update({"Score": 1})

                if right_side_parent == left_side_child:
                    right_parent_counter.update({"Score": 1})

        result_dict[result_key + 'left_parent'] = dict(left_parent_counter)
        result_dict[result_key + 'right_parent'] = dict(right_parent_counter)
    return result_dict


def filter_low_score(data_dict):
    """
    This function will filter the variants (keys in the dict) which their
    score is under 8
    :returns a new filtered dict
    """
    filtered_dict = dict()
    for key, value in data_dict.items():
        if value[3] > 8:
            filtered_dict[key] = value
    return filtered_dict


def plot_data(data_dict_to_plot):
    """
    This function will plot the given dict:
    x axis will be the chromosome position
    y axis will be the haplotype value of the current variant (position)
    """
    # Extract relevant information
    positions = list(data_dict_to_plot.keys())
    haplotypes = [entry[2] for entry in data_dict_to_plot.values()]

    # Create a range for the x-axis
    x_range = range(min(positions), max(positions) + 1)

    # Create the dot plot
    plt.scatter(positions, haplotypes, s=100, c='blue', marker='o')

    # Set plot labels and title
    plt.xlabel('Chromosome Position')
    plt.ylabel('Haplotype')
    plt.title('Dot Plot of Haplotypes')

    # Set x-axis ticks to cover the entire range
    plt.xticks(x_range)

    plt.show()


def main():
    open_and_split_file(r"HR1.ch13.phased.tsv")

    # creating and filtering the child dicts
    child_1_dict = create_and_filter_dictionary('child_1.txt')
    child_2_dict = create_and_filter_dictionary('child_2.txt')

    # adding the haplotype and score values to the keys
    child_1_batch_dict = process_rows_in_batches(child_1_dict)
    child_2_batch_dict = process_rows_in_batches(child_2_dict)


if __name__ == '__main__':
    main()