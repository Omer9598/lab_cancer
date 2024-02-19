

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
        coverage[chromosome] = (coverage_percentage, false_negative, false_positive)
        all_chrom_coverage += total_coverage
        all_chrom_true += total_real_len
    return coverage, all_chrom_coverage / all_chrom_true


def write_results(output_file, expected_coverage, inverted_coverage,
                  total_coverage, total_coverage_inverted):
    with open(output_file, 'w') as f:
        f.write("Chromosome\tRight Coverage\tFalse Negative\tFalse Positive\n")
        for chromosome in sorted(expected_coverage.keys()):
            f.write(f"{chromosome}\t{expected_coverage[chromosome][0]:.2f}%\t{expected_coverage[chromosome][1]:.2f}%"
                    f"\t{expected_coverage[chromosome][2]:.2f}%\n")
            f.write(f"{chromosome}i\t{inverted_coverage[chromosome][0]:.2f}%\t{inverted_coverage[chromosome][1]:.2f}%"
                    f"\t{inverted_coverage[chromosome][2]:.2f}%\n\n")

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

