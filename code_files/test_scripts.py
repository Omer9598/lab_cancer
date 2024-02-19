

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
    for chromosome in real_data.keys():
        real_intervals = real_data[chromosome]
        expected_intervals = expected_data.get(chromosome, [])
        true_positive = sum(1 for interval in expected_intervals if any(start <= interval[0] and end >= interval[1] for start, end in real_intervals))
        false_positive = sum(1 for interval in expected_intervals if all(start > interval[1] or end < interval[0] for start, end in real_intervals))
        total_intervals = len(real_intervals)
        coverage[chromosome] = (true_positive / total_intervals * 100, false_positive)
    return coverage


def write_results(output_file, expected_coverage, inverted_coverage):
    with open(output_file, 'w') as f:
        f.write("Chromosome\tRight Coverage\tWrong Coverage\n")
        for chromosome in sorted(expected_coverage.keys()):
            f.write(f"{chromosome}\t{expected_coverage[chromosome][0]:.2f}%\t{expected_coverage[chromosome][1]}\n")
            f.write(f"{chromosome}_inv\t{inverted_coverage[chromosome][0]:.2f}%\t{inverted_coverage[chromosome][1]}\n")


def check_right_coverage(real_data_file, expected_data_file, expected_inverted_data_file, output_file):
    real_data = parse_real_data(real_data_file)
    expected_data = parse_expected_data(expected_data_file)
    expected_inverted_data = parse_expected_data(expected_inverted_data_file)

    expected_coverage = calculate_coverage(real_data, expected_data)
    inverted_coverage = calculate_coverage(real_data, expected_inverted_data)

    write_results(output_file, expected_coverage, inverted_coverage)

