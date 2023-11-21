from collections import Counter
# import matplotlib.pyplot as plt
import cancer_path



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
    :param dictionary:
    :return:
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


def create_dictionary(file_path):
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
    return result_dict


# def process_rows_in_batches(data_dict, batch_size=10):
#     """
#     We check the 4 different cases and create 4 dictionaries for each case.
#     the format of each dict is as follows:
#     {start_position - end_position, the score of key's window (10 mutations)}
#     the 4 dictionaries will be: left-left, left-right, right-left, right-right
#     :returns a dict that contains 4 dicts above for each 10-sized "window" of
#     mutations
#     """
#     result_dict = {}
#
#     keys = list(data_dict.keys())
#
#     for start in range(len(keys) - batch_size + 1):
#         end = start + batch_size
#         batch_keys = keys[start:end]
#
#         if len(batch_keys) < 2:
#             continue
#
#         result_key = f'{batch_keys[0]}'
#         left_parent_counter = Counter({"Score": 0})
#         right_parent_counter = Counter({"Score": 0})
#
#         for key in batch_keys:
#             values = data_dict[key]
#
#             if len(values) == 2:
#                 left_side_parent, right_side_parent = values[0].split('|')
#                 left_side_child, right_side_child = values[1].split('|')
#
#                 if left_side_parent == left_side_child:
#
#                     left_parent_counter.update({"Score": 1})
#
#                 if right_side_parent == left_side_child:
#                     right_parent_counter.update({"Score": 1})
#
#         result_dict[result_key] = dict(left_parent_counter)
#         result_dict[result_key] = dict(right_parent_counter)
#     return result_dict


def add_haplotype(data_dict):
    """
    We check the 4 different cases and create 4 dictionaries for each case.
    the format of each dict is as follows:
    {start_position - end_position, the score of key's window (10 mutations)}
    the 4 dictionaries will be: left-left, left-right, right-left, right-right
    :returns a dict that contains 4 dicts above for each 10-sized "window" of
    mutations
    """
    for position, values in data_dict.items():
        left_side_parent, right_side_parent = values[0].split('|')
        left_side_child, right_side_child = values[1].split('|')

        if left_side_parent == left_side_child:
            # Condition 1: If left sides are equal, add 1 to the list
            values.append(1)
        elif right_side_parent == left_side_child:
            # Condition 2: If the right side of the parent is equal to the left side of the child, add 2 to the list
            values.append(2)
    return data_dict


def add_confidence(data_dict):
    for position, values in data_dict.items():
        haplotype = values[-1]  # Get the haplotype for the current position
        count = 0  # Initialize the count

        for next_position in data_dict.keys():
            for i in range(10):

                # Check if the next position exists in the dictionary
                next_haplotype = data_dict[next_position][-1]  # Get the haplotype for the next position
                # Check if haplotypes match
                if haplotype == next_haplotype:
                    # If they match, increment the count
                    count += 1

            # Append the count to the list in the value for the current key
            values.append(count)
    return data_dict


def main():
    open_and_split_file(cancer_path.file_path)

    # creating the child dicts
    child_1_dict = create_dictionary('child_1.txt')
    child_2_dict = create_dictionary('child_2.txt')

    # filtering the dicts
    child_1_dict = filter_dict(child_1_dict)
    child_2_dict = filter_dict(child_2_dict)

    child_1_batch_dict = add_haplotype(child_1_dict)
    child_1_batch_dict = add_confidence(child_1_dict)
    child_2_batch_dict = add_haplotype(child_2_dict)
    child_2_batch_dict = add_confidence(child_2_dict)

if __name__ == '__main__':
    main()
