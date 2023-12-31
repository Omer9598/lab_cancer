import os
import pandas as pd


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

    # Process the header line separately
    header_columns = lines[0].strip().split()
    processed_header = '\t'.join(header_columns)
    processed_lines.append(processed_header)

    for line in lines[1:]:
        # Split the line into columns
        columns = line.strip().split('\t')

        # Process each column starting from the 5th column
        for i in range(4, len(columns)):
            # Replace "./." with 0/0, "./1" with 0/1, etc.
            columns[i] = columns[i].replace("./.", "0|0")\
                .replace("./1", "0|1").replace("1/.", "1|0")\
                .replace("1/1", "1|1").replace("1/0", "1|0")\
                .replace("0/1", "0|1").replace("0/0", "0|0")

        # Join the columns back into a line
        processed_line = '\t'.join(columns)
        processed_lines.append(processed_line)

    # Join the processed lines into a string
    result = '\n'.join(processed_lines)

    # Write the processed data to the output file
    with open(output_file_path, 'w') as output_file:
        output_file.write(result)


def invert_reference_genome_haplotype(input_file, output_file):
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
            btn1, btn2 = columns[2], columns[3]

            # Process the reference and child genotypes
            ref_genotype = btn1.split('|')
            child_genotype = btn2.split('|')

            # Check if the conditions for inversion are met
            if ref_genotype[0] == '0' and ref_genotype[1] == '1' and child_genotype[0] == '1' and child_genotype[1] == '1':
                # Invert the reference genotype to 1|0
                columns[2] = '1|0'
                inverted_lines.append('\t'.join(columns))
            elif ref_genotype[0] == '1' and ref_genotype[1] == '0' and child_genotype[0] == '0' and child_genotype[1] == '0':
                # Invert the reference genotype to 0|1
                columns[2] = '0|1'
                inverted_lines.append('\t'.join(columns))

    # Write the inverted lines to the output file
    with open(output_file, 'w') as output_file:
        output_file.write('\n'.join(inverted_lines))


def merge_haplotype_tables(input_directory):
    """
    This function will merge all the tables given in the input_directory
    The format will be as follows:
    CHROM START END HAPLOTYPE
    Specifying each interval's chromosome, location, and haplotype
    """
    # Initialize a list to store the merged intervals
    merged_intervals = []

    # Loop through chromosomes 1 to 22
    for chrom_num in range(1, 23):
        # Read each haplotype interval table file
        file_path = f"{input_directory}/haplotype_interval_table_{chrom_num}.txt"

        try:
            with open(file_path, 'r') as file:
                # Skip the header line
                next(file)

                # Process each line in the file
                for line in file:
                    # Split the line into columns
                    columns = line.strip().split()

                    # Append the columns to the merged intervals list
                    merged_intervals.append(columns)
        except FileNotFoundError:
            print(f"File not found: {file_path}")

    # Write the merged intervals to a new file in the input directory
    output_path = f"{input_directory}/merged_haplotype_intervals.txt"
    with open(output_path, 'w') as output_file:
        # Write the header line
        output_file.write("CHROM\tSTART\tEND\tHAPLOTYPE\n")

        # Write each merged interval to the output file
        for interval in merged_intervals:
            output_file.write('\t'.join(interval) + '\n')


def create_table(data_list, output_directory):
    """
    This function will create the shared haplotype intervals table
    The table will be in a new .txt file, ordered in the following format - each
    row represents an interval, and will be as follows:
    chromosome - the chromosome number
    start - the starting position of the interval in the row's chromosome
    end - the ending position
    haplotype - the haplotype of the current interval
    the .txt file will be saved in the haplotype_interval_tables
    directory
    """
    # Group data by chromosome
    grouped_data = {}
    for entry in data_list:
        chromosome = entry['chromosome']
        if chromosome not in grouped_data:
            grouped_data[chromosome] = []
        grouped_data[chromosome].append(entry)

    # Write data for each chromosome
    for chromosome, chromosome_data in grouped_data.items():
        # Reorder the keys to make "chromosome" the first column
        data_list_reordered = [
            {k: entry[k] for k in ['chromosome', 'start', 'end', 'haplotype']} for
            entry in chromosome_data]

        file_path = os.path.join(output_directory,
                                 f'haplotype_interval_table_{chromosome}.txt')

        # Write the data to a text file without any table formatting
        with open(file_path, 'w') as file:
            for entry in data_list_reordered:
                file.write(f"{entry['chromosome']}\t{entry['start']}\t"
                           f"{entry['end']}\t{entry['haplotype']}\n")


def open_and_split_children_files(file_path):
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


def split_file_to_chromosomes(input_file, output_directory):
    """
    This function will split the input_file to different files, according to the
    number of chromosomes (e.g 23 chromosomes in the input file)
    The files the function creates will be placed in the output_directory
    """
    # Read the input file into a df, specifying datatype for 'CHROM' as str
    df = pd.read_table(input_file, sep='\t', dtype={'CHROM': str})

    # Group the DataFrame by the 'CHROM' column
    grouped = df.groupby('CHROM')

    # Iterate through each group and write to separate files
    for chrom, group in grouped:
        # Define the output file path for each chromosome
        output_file = os.path.join(output_directory, f'chromosome_{chrom}.txt')

        # Write the group to the output file
        group.to_csv(output_file, sep='\t', index=False)
