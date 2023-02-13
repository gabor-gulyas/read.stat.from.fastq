import argparse
import os
import gzip
import numpy as np
import pandas as pd
import statistics


def read_fastq(file):
    with open(file, "r") as fastq:
        while True:
            lines = [fastq.readline().strip() for i in range(4)]
            if not lines[0]:
                break
            yield lines


def process_fastq(file, out_dir):
    sample_name = os.path.splitext(os.path.basename(file))[0]
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, sample_name + ".tsv")

    with open(out_file, "w") as out:
        out.write("read_id\tread_length\tread_quality\taverage_base_quality\n")

        for i, lines in enumerate(read_fastq(file)):
            read_id, sequence, _, quality = lines
            read_length = len(sequence)
            read_quality = sum([ord(q) - 33 for q in quality]) / read_length
            average_base_quality = read_quality / read_length
            out.write(f"{read_id}\t{read_length}\t{read_quality:.2f}\t{average_base_quality:.2f}\n")

            if (i + 1) % 1000 == 0:
                print(f"Processed {i + 1} reads for sample {sample_name}")


def readstat_calculator(sample, read_lengths, read_qualities, avg_base_qualities):
    read_count = len(read_lengths)
    length_min = min(read_lengths)
    length_median = statistics.median(read_lengths)
    length_mean = statistics.mean(read_lengths)
    length_max = max(read_lengths)

    quality_min = min(read_qualities)
    quality_median = statistics.median(read_qualities)
    quality_mean = statistics.mean(read_qualities)
    quality_max = max(read_qualities)

    avg_base_quality_min = min(avg_base_qualities)
    avg_base_quality_median = statistics.median(avg_base_qualities)
    avg_base_quality_mean = statistics.mean(avg_base_qualities)
    avg_base_quality_max = max(avg_base_qualities)

    return [sample, read_count, length_min, length_median, length_mean, length_max,
            quality_min, quality_median, quality_mean, quality_max,
            avg_base_quality_min, avg_base_quality_median, avg_base_quality_mean, avg_base_quality_max]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_folder", required=True, help="Path to the folder containing fastq files")
    parser.add_argument("-o", "--output_folder", required=True,
                        help="Path to the folder where the results will be saved")
    args = parser.parse_args()

    for file in os.listdir(args.input_folder):
        if file.endswith((".fastq", ".fastq.gz")):
            file_path = os.path.join(args.input_folder, file)
            if file.endswith(".gz"):
                with gzip.open(file_path, "rt") as fastq:
                    process_fastq(fastq, args.output_folder)
            else:
                process_fastq(file_path, args.output_folder)

    print("Starting second for loop")

    datatable = pd.DataFrame(
        columns=["sample", "read_count", "length_min", "length_median", "length_mean", "length_max",
                 "quality_min", "quality_median", "quality_mean", "quality_max",
                 "avg_base_quality_min", "avg_base_quality_median", "avg_base_quality_mean", "avg_base_quality_max"])

    for file in os.listdir(args.output_folder):
        if file.endswith(".tsv"):
            print(file, "calc")
            file_path = os.path.join(args.output_folder, file)
            sample_name = os.path.splitext(os.path.basename(file))[0]
            try:
                data = pd.read_csv(file_path, sep='\t')
            except Exception as e:
                print("Error reading file:", file_path)
                print("Error:", e)
                continue
            print("Columns in data:", data.columns)
            try:
                out_data = readstat_calculator(sample_name, data["read_length"], data["read_quality"],
                                               data["average_base_quality"])
            except Exception as e:
                print("Error running readstat_calculator for file:", file_path)
                print("Error:", e)
            try:
                datatable.loc[file] = out_data
            except Exception as e:
                print("Error while printing data table: ", e)

    print(datatable)
    try:
        datatable.to_csv(os.path.join(args.output_folder, "read.stat.tsv"), index=False, sep='\t')
    except Exception as e:
        print("Error while writing output dataframe: ", e)


if __name__ == "__main__":
    main()
