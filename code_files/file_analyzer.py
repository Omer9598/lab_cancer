import os
import tempfile

import pandas as pd

replacements = {"./.": "0|0", "./1": "0|1", "1/.": "1|0", "1/1": "1|1",
                "1/0": "1|0", "0/1": "0|1", "0/0": "0|0"}


def preprocess_file(input_file_path, output_directory):
    """
    This function will preprocess the file given:
    Changing "./." to 0/0, or "./1" to 0/1 etc. and save the result in a new
    file.
    """
    with open(input_file_path, 'r') as file:
        lines = file.readlines()
    processed_lines = []
    header_index = None
    for i, line in enumerate(lines):
        if line.startswith("#CHROM"):
            header_index = i
            break
    if header_index is not None:
        header_line = lines[header_index].strip().replace("#", "")
        processed_lines.append(header_line)
        # Process lines after the header
        for line in lines[header_index + 1:]:
            columns = line.strip().split('\t')
            # Process each column starting from the 5th column
            for i in range(4, len(columns)):
                # Replace values based on the dictionary
                columns[i] = replacements.get(columns[i], columns[i])
            processed_line = '\t'.join(columns)
            processed_lines.append(processed_line)
        result = '\n'.join(processed_lines)
        input_file_name = os.path.basename(input_file_path)
        output_file_name = "processed_" + input_file_name
        output_file_path = os.path.join(output_directory, output_file_name)
        with open(output_file_path, 'w') as output_file:
            output_file.write(result)
        return output_file_path


def merge_haplotype_tables(input_directory, chromosome_coverage_dict,
                           window_size, error_size, invert):
    """
    This function will merge all the tables given in the input_directory
    The format will be as follows:
    CHROM START END HAPLOTYPE
    Specifying each interval's chromosome, location, and haplotype
    """
    merged_intervals = []
    for chrom_num in range(1, 23):
        file_path = (f"{input_directory}/table_{chrom_num}_"
                     f"window_{window_size}_error_{error_size}_inverted_"
                     f"{bool(invert)}.txt")
        with open(file_path, 'r') as file:
            for line in file:
                columns = line.strip().split()
                merged_intervals.append(columns)
        chrom_coverage = round(chromosome_coverage_dict[chrom_num] * 100, 1)
        columns.append(str(chrom_coverage) + "%")
    output_path_merged = f"{input_directory}/merged_haplotype_intervals.txt"
    with open(output_path_merged, 'w') as output_file:
        output_file.write("CHROM\tSTART\tEND\tHAPLOTYPE\tCERTAINTY\tCOVERAGE\n")
        # Write each merged interval to the output file
        for interval in merged_intervals:
            output_file.write('\t'.join(interval) + '\n')
    output_path_excel = f"{input_directory}/merged_haplotype_Excel.xlsx"
    convert_txt_to_excel(output_path_merged, output_path_excel)


def write_common_genes_to_file(output_directory, cancer_genes_dict):
    """
    Writes common cancer genes information to a .txt file and
    to Excel afterward in the specified output directory.
    """
    output_path_txt = f"{output_directory}/common_cancer_genes.txt"

    with open(output_path_txt, 'w') as output_file:
        output_file.write("Common cancer genes in shared intervals: \n")
        for variant_name, variant_info in cancer_genes_dict.items():
            output_file.write(f"{variant_name} - {variant_info[3]}\n")

    # Convert the common genes txt file to Excel
    df = pd.read_csv(output_path_txt, delimiter='\t')
    excel_output_path = f"{output_directory}/common_cancer_genes_Excel.xlsx"
    df.to_excel(excel_output_path, index=False)


def create_table(data_list, output_directory, window_size, error_size, inverted):
    """
    This function will create the shared haplotype intervals table
    The table will be in a new .txt file, ordered in the following format -
    each row represents an interval, and will be as follows:
    chromosome - the chromosome number
    start - the starting position of the interval in the row's chromosome
    end - the ending position
    haplotype - the haplotype of the current interval
    certainty_level - 1 for areas with the same haplotype,
    -1 for areas with opposing haplotypes
    the .txt file will be saved in the interval_tables directory
    """
    os.makedirs(output_directory, exist_ok=True)
    # Group data by chromosome
    grouped_data = {}
    for entry in data_list:
        chromosome = entry['chromosome']
        if chromosome not in grouped_data:
            grouped_data[chromosome] = []
        grouped_data[chromosome].append(entry)
    # Write data for each chromosome
    for chromosome, chromosome_data in grouped_data.items():
        data_list_reordered = [
            {k: entry[k] for k in ['chromosome', 'start', 'end', 'haplotype']}
            for entry in chromosome_data]
        # Calculate certainty level
        for entry in data_list_reordered:
            entry['certainty_level'] = -1 if (entry['haplotype'] ==
                                              chromosome_data[0]['haplotype']
                                              or entry['haplotype'] == 0) else 1
        file_path = os.path.join(output_directory,
                                 f'table_{chromosome}_window_{window_size}'
                                 f'_error_{error_size}_inverted_'
                                 f'{bool(inverted)}.txt')
        with open(file_path, 'w') as file:
            for entry in data_list_reordered:
                file.write(f"{entry['chromosome']}\t{entry['start']}\t"
                           f"{entry['end']}\t{entry['haplotype']}\t"
                           f"{entry['certainty_level']}\n")


def invert_reference_genome_haplotype(input_file, output_directory):
    """
    This function will create a new file from the file given, that will contain
    the inverted variants (if needed) of the reference genotype
    for example:
    reference - 0|1 child - 1|1, the reference will be inverted to 1|0
    reference - 0|1 child - 0|0, the reference won't be inverted
    If cases where the reference is homozygous or the child is heterozygous
    will be ignored
    In other words, the reference will always inherit the left side (haplotype =
    1)
    """
    inverted_lines = []
    with open(input_file, 'r') as file:
        header = file.readline().strip()
        inverted_lines.append(header)
        for line in file:
            columns = line.strip().split('\t')
            btn1, btn2 = columns[4], columns[5]
            ref_genotype = btn1.split('|')
            child_genotype = btn2.split('|')
            # Check if the conditions for inversion are met
            if ref_genotype[0] == '0' and ref_genotype[1] == '1' and\
                    child_genotype[0] == '1' and child_genotype[1] == '1':
                # Invert the reference genotype to 1|0
                columns[4] = '1|0'
            elif ref_genotype[0] == '1' and ref_genotype[1] == '0' and\
                    child_genotype[0] == '0' and child_genotype[1] == '0':
                # Invert the reference genotype to 0|1
                columns[4] = '0|1'
            inverted_lines.append('\t'.join(columns))
    os.makedirs(output_directory, exist_ok=True)
    output_file_name = "inverted_" + os.path.basename(input_file)
    output_file_path = os.path.join(output_directory, output_file_name)
    with open(output_file_path, 'w') as output_file:
        output_file.write('\n'.join(inverted_lines))
    return output_file_path


def open_and_split_children_files(file_path):
    """
    This function will open the file and split it into n files.
    Each file will contain column[0] column[1], column[4], and subsequent
    columns from 5 to the last column.
    The number of child files will be determined based on the available
    columns.
    """
    child_filenames = []
    with open(file_path, 'r') as infile:
        header_columns = infile.readline().strip().split('\t')
        # Determine the number of child files based on available columns
        num_children = len(header_columns) - 5
        for child_num in range(1, num_children + 1):
            with tempfile.NamedTemporaryFile(
                    mode='w+', delete=False) as child_file:
                child_filenames.append(child_file.name)
                # Write the header columns to the child file
                header_line = '\t'.join([header_columns[0], header_columns[1],
                                         header_columns[4],
                                         header_columns[child_num + 4]])
                child_file.write(header_line + '\n')
                infile.seek(0)  # Reset file pointer to the beginning
                next(infile)  # Skip the header line
                for line in infile:
                    columns = line.strip().split('\t')
                    # Extract the desired columns for the current child file
                    child_line = f"{columns[0]}\t{columns[1]}\t{columns[4]}\t" \
                                 f"{columns[child_num + 4]}\n"
                    child_file.write(child_line)
    return num_children, child_filenames


def split_file_to_chromosomes(input_file, output_directory):
    """
    This function will split the input_file to different files, according to the
    number of chromosomes (e.g. 23 chromosomes in the input file)
    The files the function creates will be placed in the output_directory
    """
    # Read the input file into a df, specifying datatype for 'CHROM' as str
    df = pd.read_table(input_file, sep='\t', dtype={'CHROM': str})
    grouped = df.groupby('CHROM')
    os.makedirs(output_directory, exist_ok=True)
    for chrom, group in grouped:
        output_file = os.path.join(output_directory, f'chromosome_{chrom}.txt')
        group.to_csv(output_file, sep='\t', index=False)


def convert_txt_to_excel(input_file, output_excel):
    df = pd.read_csv(input_file, sep='\t')
    df.to_excel(output_excel, index=False)
