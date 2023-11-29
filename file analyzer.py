

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

