import os
import re

import pandas as pd
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


def write_results(output_file, expected_coverage, inverted_coverage,
                  total_coverage, total_coverage_inverted):
    with open(output_file, 'w') as f:
        f.write("Chromosome\tRight_Coverage\tFalse_Negative\tFalse_Positive\tF1_score\n")
        for chromosome in sorted(expected_coverage.keys()):
            f.write(f"{chromosome}\t{expected_coverage[chromosome][0]:.2f}%\t{expected_coverage[chromosome][1]:.2f}%"
                    f"\t{expected_coverage[chromosome][2]:.2f}%\t{expected_coverage[chromosome][3]:.2f}\n")
            f.write(f"{chromosome}i\t{inverted_coverage[chromosome][0]:.2f}%\t{inverted_coverage[chromosome][1]:.2f}%"
                    f"\t{inverted_coverage[chromosome][2]:.2f}%\t{inverted_coverage[chromosome][3]:.2f}\n\n")

        f.write(f"Total Coverage Regular: {100 * total_coverage:.2f}%\n"
                f"Total Coverage Inverted: {100 * total_coverage_inverted:.2f}%")


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


def plot_f1_score(table_data):
    # Assuming table_data is a list of strings where each string represents a row in the table
    # Split each row into columns using tab as a separator
    data_rows = [row.split('\t') for row in table_data]

    # Extract column names from the first row
    column_names = data_rows[0]

    # Find the index of the "f1_score" column
    f1_score_index = column_names.index("f1_score")

    # Extract data for plotting
    points = [(int(row[1]), int(row[2]), float(row[f1_score_index])) for row in data_rows[1:]]
    points_sorted = sorted(points, key=lambda x: (x[0], x[1]))  # Sort by window and then error

    # Unpack sorted points
    window_error_labels = [f"{window}-{error}" for window, error, _ in points_sorted]
    f1_scores = [score for _, _, score in points_sorted]

    # Plotting
    plt.scatter(window_error_labels, f1_scores, color='green', s=10, edgecolors='black')
    plt.xlabel('Window-Error')
    plt.ylabel('F1 Score')
    plt.title('F1 Score for Different Window-Error Combinations')

    # Set y-axis limits to range from 0 to 100%
    plt.ylim(0, 100)

    # Add grid lines at every 5% interval on the y-axis
    plt.yticks(range(0, 101, 5))
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for better visibility
    plt.show()


def plot_coverage(file_path):
    # Lecture des données depuis le fichier
    df = pd.read_csv(file_path, delimiter='\t')

    # Création d'une liste de tuples pour chaque point
    points = list(zip(df['window'], df['error'], df['coverage']))

    # Tri des points en fonction de la fenêtre et de l'erreur
    points.sort(key=lambda x: (x[0], -x[1]))

    # Séparation des points triés en listes distinctes
    sorted_window, sorted_error, sorted_coverage = zip(*points)

    # Création du graphique
    plt.figure(figsize=(10, 6))
    plt.scatter(range(len(points)), sorted_coverage, marker='o')

    # Étiquetage des points avec les valeurs de la fenêtre et de l'erreur
    for i, txt in enumerate(zip(sorted_window, sorted_error)):
        plt.annotate(f'{txt[0]}-{txt[1]}', (i, sorted_coverage[i]), textcoords="offset points", xytext=(0, 5), ha='center')

    # Ajout des étiquettes d'axe et du titre
    plt.xlabel('Window-Error')
    plt.ylabel('Coverage (%)')
    plt.title('Graphique de Coverage en fonction de Window-Error')
    plt.xticks(range(len(points)))

    # Affichage du graphique
    plt.show()

