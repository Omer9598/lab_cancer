import os
import re
import numpy as np

from matplotlib import pyplot as plt


def parse_real_data(real_data_file):
    real_data = {}
    with open(real_data_file, 'r') as f:
        for line in f:
            chromosome, start, end = map(int, line.strip().split('\t'))
            real_data.setdefault(chromosome, []).append((start, end))
    return real_data


def parse_expected_data(expected_data_file):
    expected_data = {}
    with open(expected_data_file, 'r') as f:
        next(f)  # skip header
        for line in f:
            if len(line.strip().split('\t')) == 6:
                chromosome, start, end, haplotype, certainty, coverage = line.strip().split('\t')
            if len(line.strip().split('\t')) == 5:
                chromosome, start, end, haplotype, certainty = line.strip().split('\t')
            if certainty == '1':
                expected_data.setdefault(chromosome, []).append((int(start), int(end)))
    return expected_data


def calculate_coverage(real_data, expected_data):
    coverage = {}
    all_chrom_coverage = 0
    all_chrom_true = 0

    for chromosome in real_data.keys():
        real_intervals = real_data[chromosome]
        if str(chromosome) not in expected_data:
            continue
        expected_intervals = expected_data[str(chromosome)]
        total_coverage = 0
        for expected_interval in expected_intervals:
            interval_coverage = 0
            for real_interval in real_intervals:
                overlap_start = max(expected_interval[0], real_interval[0])
                overlap_end = min(expected_interval[1], real_interval[1])
                overlap_length = max(0, overlap_end - overlap_start)
                interval_coverage += overlap_length
            total_coverage += interval_coverage
        total_real_len = sum(end - start for start, end in real_intervals)
        total_expected_len = sum(end - start for start, end in expected_intervals)
        coverage_percentage = total_coverage / total_real_len * 100 if total_real_len > 0 else 0
        false_negative = (total_real_len - total_coverage) / total_real_len * 100
        false_positive = (total_expected_len - total_coverage) / total_expected_len * 100

        precision = coverage_percentage / (coverage_percentage + false_positive)
        recall = coverage_percentage / (coverage_percentage + false_negative)
        if recall + precision == 0:
            f1_score = 0
        else:
            f1_score = 200 * (precision * recall / (precision + recall))

        coverage[chromosome] = (coverage_percentage, false_negative, false_positive, f1_score)
        all_chrom_coverage += total_coverage
        all_chrom_true += total_real_len

    return coverage, (all_chrom_coverage / all_chrom_true)


def write_results(output_file, expected_coverage, inverted_coverage, total_coverage, total_coverage_inverted):
    # Define the width for each column
    col_widths = {"Chromosome": 12, "Right_Coverage": 15, "False_Negative": 15, "False_Positive": 15, "F1_score": 10}

    # Create a format string with appropriate column widths
    header_fmt = "{:<{Chromosome}}{:<{Right_Coverage}}{:<{False_Negative}}{:<{False_Positive}}{:<{F1_score}}\n"
    row_fmt = "{:<{Chromosome}}{:<{Right_Coverage}.2f}{:<{False_Negative}.2f}{:<{False_Positive}.2f}{:<{F1_score}.2f}\n"

    # Open the file and write the formatted content
    with open(output_file, 'w') as f:
        # Write the header
        f.write(header_fmt.format("Chromosome", "Right_Coverage", "False_Negative", "False_Positive", "F1_score",
            **col_widths))

        # Write the data rows
        for chromosome in sorted(expected_coverage.keys()):
            chromosome_str = str(chromosome)  # Ensure chromosome is a string
            f.write(row_fmt.format(chromosome_str, expected_coverage[chromosome][0], expected_coverage[chromosome][1],
                expected_coverage[chromosome][2], expected_coverage[chromosome][3], **col_widths))
            f.write(
                row_fmt.format(chromosome_str + "i", inverted_coverage[chromosome][0], inverted_coverage[chromosome][1],
                    inverted_coverage[chromosome][2], inverted_coverage[chromosome][3], **col_widths))
            f.write("\n")

        # Write the total coverage
        f.write(f"Total Coverage Regular: {100 * total_coverage:.2f}%\n")
        f.write(f"Total Coverage Inverted: {100 * total_coverage_inverted:.2f}%\n")


def check_right_coverage(real_data_file, expected_data_file, expected_inverted_data_file, output_file):
    real_data = parse_real_data(real_data_file)
    expected_data = parse_expected_data(expected_data_file)
    expected_inverted_data = parse_expected_data(expected_inverted_data_file)

    expected_coverage, total_coverage = calculate_coverage(real_data, expected_data)
    inverted_coverage, total_coverage_inverted = calculate_coverage(real_data, expected_inverted_data)

    write_results(output_file, expected_coverage, inverted_coverage,
                  total_coverage, total_coverage_inverted)


def merge_coverage_files(directory_path, output_file):
    merged_data = []

    # Iterate over all files in the directory
    for file_name in os.listdir(directory_path):
        file_path = os.path.join(directory_path, file_name)
        if not os.path.isfile(file_path):
            continue  # Skip directories

        # Extract chromosome, window, and error from the file name
        match = re.match(r'.*chrom_(\d+)_window_(\d+)_error_(\d+)', file_name)
        if match:
            chromosome, window, error = match.groups()
        else:
            print(f"Skipping file {file_name} as it doesn't match the pattern.")
            continue

        # Read data from the file
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('Chromosome'):
                    continue  # Skip header
                fields = line.strip().split('\t')
                if len(fields) == 5:  # Check if the line has expected format
                    _, right_coverage, false_negative, false_positive, f1_score = fields
                    merged_data.append((chromosome, window, error, f1_score, right_coverage))

    # Write merged data to output file
    with open(output_file, 'w') as f:
        f.write("chromosome\twindow\terror\tf1_score\tcoverage\n")
        for entry in merged_data:
            formatted_entry = '\t'.join('{:<10}'.format(str(item)) for item in entry)
            f.write(formatted_entry + '\n')


def read_data_from_file(file_path):
    data = []
    with open(file_path, 'r') as file:
        # Skip the header line
        next(file)
        for line in file:
            values = line.strip().split('\t')
            entry = {
                'chromosome': int(values[0]),
                'window': int(values[1]),
                'error': int(values[2]),
                'f1_score': float(values[3]),
                'coverage': values[4]
            }
            data.append(entry)
    return data


def plot_f1_score(file_path):
    data = read_data_from_file(file_path)
    x_values = []
    y_values = []
    colors = []

    for entry in data:
        entry['coverage'] = float(entry['coverage'].strip('%'))

    sorted_data = sorted(data, key=lambda x: (x['window'], x['error']))

    for entry in sorted_data:
        x_values.append(f"{entry['window']}-{entry['error']}")
        y_values.append(entry['f1_score'])
        colors.append('green' if colors.count('green') < colors.count('red') else 'red')

    plt.figure(figsize=(10, 6))
    plt.scatter(x_values, y_values, c=colors)
    plt.title('Plot window-error f1_score')
    plt.xlabel('Window-Error')
    plt.ylabel('F1 Score')
    plt.yticks(np.arange(0, 101, 5))
    plt.tight_layout()
    plt.show()


def plot_coverage(data):
    x_values = []
    y_values = []
    colors = []

    for entry in data:
        entry['coverage'] = float(entry['coverage'].strip('%'))

    sorted_data = sorted(data, key=lambda x: (x['window'], x['error']))

    for entry in sorted_data:
        x_values.append(f"{entry['window']}-{entry['error']}")
        y_values.append(entry['coverage'])
        colors.append('green' if colors.count('green') < colors.count('red') else 'red')

    plt.figure(figsize=(10, 6))
    plt.scatter(x_values, y_values, c=colors)
    plt.title('Plot window-error coverage')
    plt.xlabel('Window-Error')
    plt.ylabel('Coverage (%)')
    plt.yticks(np.arange(0, 101, 5))
    plt.tight_layout()
    plt.show()