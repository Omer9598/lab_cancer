from itertools import islice


WINDOW_NUM = 20
WINDOW_ERR = 18
PARENT_REFERENCE = "parent"
SIBLING_REFERENCE = "sibling"


def check_heterozygous(parent):
    if parent == '0|0' or parent == '1|1':
        return False
    return True


def check_homozygous(child):
    if child == '0|1' or child == '1|0':
        return False
    return True


def filter_dict_parent_reference(dictionary):
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
        result_parent = check_heterozygous(value[0])

        # Apply check_homozygous_child to the value associated with the key
        result_child = check_homozygous(value[1])

        # If both functions return True, add the key to the filtered_dict
        if result_parent and result_child:
            filtered_dict[key] = value

    return filtered_dict


def filter_dict_sibling_reference(dictionary):
    """
    This function will filter the dict given, the same way as
    filter_dict_parent_reference function, but not filtering in case that both
    reference and non-reference are homozygous
    """
    filtered_dict = {}
    for key, value in dictionary.items():
        parent = value[0]
        child = value[1]

        result_parent = check_heterozygous(parent)
        result_child = check_homozygous(child)
        result_sibling = (parent == '0|0' and child == '1|1') or\
                         (child == '0|0' and parent == '1|1')

        if result_parent and result_child or result_sibling:
            filtered_dict[key] = (parent, child)

    return filtered_dict


def create_and_filter_dictionary(file_path, reference_type):
    """
    This function will create dictionary from the given file, in the following
    format:
    {position: [mother, child], ...}
    the filter depends on the reference type given - sibling reference will
    not filter homozygous siblings
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
    # Filter the dict according to reference type
    if reference_type == PARENT_REFERENCE:
        return filter_dict_sibling_reference(result_dict)
    if reference_type == SIBLING_REFERENCE:
        return filter_dict_sibling_reference(result_dict)


def add_haplotype_parent_reference(my_dict):
    """
    This function will add the haplotype of each variant to the child_dict
    """
    for position, values in my_dict.items():
        left_side_parent, right_side_parent = values[0].split('|')
        left_side_child, right_side_child = values[1].split('|')

        if left_side_parent == left_side_child:
            # If left sides are equal, add 1 to the list
            values.append(1)
        elif right_side_parent == left_side_child:
            # If the right side of the parent is equal to the left side of the
            # child, add 2 to the list
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


def process_dict(data_dict):
    """
    This function will process the dicts to be plotted
    """
    add_haplotype_parent_reference(data_dict)
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
        if value[-1] > WINDOW_ERR:
            filtered_dict[key] = value
    return filtered_dict
