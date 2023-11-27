from itertools import islice
from collections import deque


WINDOW_NUM = 20
WINDOW_ERR = 18


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


def create_and_filter_dictionary(file_path):
    """
    This file will create dictionary from the given file, in the following
    format:
    {position: [mother, child], ...}
    """
    result_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            row = line.strip().split('\t')
            # skip the first line
            if row[0] == "CHROM":
                continue
            position = int(row[1])
            values = []
            # adding the inheritance values of each child
            for i in range(2, len(row)):
                values.append(row[i])
            # adding the chromosome number
            values.append(row[0])
            result_dict[position] = values
    return filter_dict(result_dict)


def add_haplotype(my_dict):
    """
    This function will add the haplotype of each variant to the child_dict
    """
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
    """
    This function will add the confidence value to the dict
    the confidence is the number of variants that similar to the current
    variant, in a given window - the size of the window is determined by
    WINDOW_NUM global parameter
    """
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
                next_haplotype = my_dict[next_position][-1]  # Get the haplotype for the next position
                # Check if haplotypes match
                if haplotype == next_haplotype:
                    # If they match, increment the count
                    count += 1
        values.append(count)


def add_confidence_efficient(data_dict):
    """
    This function will add the confidence value to the dict
    the confidence is the number of variants that are similar to the current
    variant, in a given window - the size of the window is determined by
    WINDOW_NUM global parameter
    """
    positions = list(data_dict.keys())  # positions already sorted

    # Starting the window - loop through the first WINDOW_NUM positions
    window_1_counter = 0
    window_2_counter = 0
    position_index = WINDOW_NUM - 1
    for i in range(WINDOW_NUM):
        cur_position = positions[i]
        haplotype = data_dict[cur_position][-1]
        if haplotype == 1:
            window_1_counter += 1
        if haplotype == 2:
            window_2_counter += 1

    for first_position_in_window, value_list in data_dict.items():
        cur_haplotype = value_list[-1]
        # add the confidence value to the dict and updating the window
        if cur_haplotype == 1:
            data_dict[first_position_in_window].append(window_1_counter)
            window_1_counter -= 1
        if cur_haplotype == 2:
            data_dict[first_position_in_window].append(window_2_counter)
            window_2_counter -= 1

        if position_index < len(positions):
            last_position_in_window = positions[position_index]
            last_position_haplotype = data_dict[last_position_in_window][-1]
            if last_position_haplotype == 1:
                window_1_counter += 1
            if last_position_haplotype == 2:
                window_2_counter += 1
            position_index += 1


def process_dict(data_dict):
    """
    This function will process the dicts to be plotted
    """
    add_haplotype(data_dict)
    add_confidence_efficient(data_dict)
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
        if value[-1] > WINDOW_ERR:
            filtered_dict[key] = value
    return filtered_dict
