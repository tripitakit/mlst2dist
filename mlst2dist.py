#!/usr/bin/python

"""
mlst2dist.py
This program calculate a distance (dissimilarity) matrix from a chewBBACA results_alleles.tsv table (https://github.com/B-UMMI/chewBBACA).

It converts all occurrences of LNF, NIPH, NIPHEM, PLOT3, PLOT5, ALM, ASM, LOTSC in 0 and transform all INF-X in the corresponding X inferred allele-calls.
Any call value that isn't an integer is converted to 0.
The program then calculates a matrix of pairwise dissimilarities using a Hamming Distance modified with correction for missing data (the allele-calls converted to 0s).

The distance algorithm implemented in mlst2dist.py is described in Galpern P, Manseau, M, Hettinga P, Smith K, and Wilson P. (2012) allelematch: an R package for identifying unique multilocus genotypes where genotype error and missing data may be present. Molecular Ecology Resources 12:771-778

The dissimilarity matrix output is saved to disk with default PHYLIP square matrix format, or optional TSV format.
"""


__author__ = "Patrick De Marta"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser(
        description="Converts a chewBBACA alleles table into a dissimilarity matrix.")

    parser.add_argument("input", help="chewBBACA alleles table .tsv")
    parser.add_argument("output", help="distance matrix output filename")
    parser.add_argument("-f", "--outfmt", choices=["CSV", "MEG", "PHY"], default="MEG",
                        help="format options for the output matrix [default: MEG (MEGAX format)]")
    args = parser.parse_args()
    return args


def sanitize(line):
    sanitized_line = line.replace("\n", "").split("\t")
    for index, a in list(enumerate(sanitized_line))[1:]:

        if "INF-" in a:
            sanitized_allele = a.replace("INF-", "")
        elif a in ["LNF", "PLOT3", "PLOT5", "NIPH", "NIPHEM", "ASM", "ALM", "LOTSC"]:
            sanitized_allele = "0"
        else:
            try:
                sanitized_allele = str(int(a))
            except:
                print(f"Unexpected allele call: {a}")
                sanitized_allele = "0"

        sanitized_line[index] = sanitized_allele.replace(" ", "")
    return sanitized_line


def transform_alleles(input):
    allele_matrix = []
    loci_vector = []
    samples_vector = []

    with open(input, "r") as alleles:
        for index, line in enumerate(alleles):
            if (index != 0):
                sanitized_line = sanitize(line)
                samples_vector.append(sanitized_line[0])
                allele_matrix.append(sanitized_line[1:])
            else:
                loci_vector = line.replace("\n", "").split("\t")[1:]
    return (samples_vector, loci_vector, allele_matrix)


def write_phylip(distance_matrix, samples, output_fname):
    nofsamp = len(samples)
    header = f"{nofsamp}\n"
    body = render_square_matrix(distance_matrix, samples, nofsamp)
    write_outfile(output_fname, header, body)


def write_tsv(distance_matrix, samples, output_fname):
    header = "\t".join(samples) + "\n"
    body = render_square_matrix(distance_matrix, samples, len(samples))
    write_outfile(output_fname, header, body)


def render_square_matrix(distance_matrix, samples, nofsamp):
    body = ""
    for i in range(nofsamp):
        samp_line = f"{samples[i].ljust(10)}"
        for j in range(nofsamp):
            samp_line += f"\t{distance_matrix[i][j]}"
        body += samp_line + "\n"
    return body


def write_mega(distance_matrix, samples, output_fname):
    header = build_mega_header(samples)
    body = build_mega_body(distance_matrix, samples)
    write_outfile(output_fname, header, body)


def build_mega_header(samples):
    header = f"#mega;\n!Title: <your-title-here>;\n!Format DataType=Distance DataFormat=LowerLeft NTaxa={len(samples)};\n\n"
    for i, sample in enumerate(samples):
        header += f"[{pad_mega_index(i+1)}] #{sample}\n"
    header += "\n"
    return header


def pad_mega_index(i):
    return str(i).rjust(2, " ")


def build_mega_body(distance_matrix, samples):
    body = ""

    ruler = "[     "
    for idx in list(range(len(samples))):
        ruler += f"       {idx + 1}"

    ruler += " ]\n"

    body = ruler + "[ 1]        \n"
    for i in range(1, len(distance_matrix)):
        row = ""
        for j in range(i):
            row += " " + str(distance_matrix[i][j])

        body += f"[{pad_mega_index(i+1)}]  " + row + "        \n"
    return body


def write_outfile(output_fname, header, body):
    outh = open(output_fname, "w")
    outh.write(header)
    outh.write(body)
    outh.close()


def make_matrix(input):
    samples, loci, alleles_matrix = transform_alleles(input)

    nloci = len(loci)
    nsamp = len(samples)

    distance_matrix = []

    for samp_i in range(nsamp):
        alleles_vector_a = alleles_matrix[samp_i]

        samp_line = []
        for samp_j in range(nsamp):
            alleles_vector_b = alleles_matrix[samp_j]
            n_matches = 0
            n_missing = 0

            for x in range(nloci):
                if (alleles_vector_a[x] == alleles_vector_b[x]):
                    n_matches += 1
                else:
                    if (alleles_vector_a[x] == 0 or alleles_vector_b[x] == 0):
                        n_missing += 1

            corrected_hamming_dist_a_b = round(
                1 - (n_matches/nloci) + (n_missing/(2*nloci)), 5)

            samp_line.append(corrected_hamming_dist_a_b)

        distance_matrix.append(samp_line)

    return(distance_matrix, samples)


def dispatch_output(distance_matrix, samples, outfmt, output):
    if outfmt == "PHY":
        write_phylip(distance_matrix, samples, output)
    if outfmt == "TSV":
        write_tsv(distance_matrix, samples, output)
    if outfmt == "MEG":
        write_mega(distance_matrix, samples, output)


def validate_infile(input):
    return os.path.exists(input)


def main():
    args = parse_args()
    input_fname = args.input
    output_fname = args.output
    output_format = args.outfmt

    if validate_infile(input_fname):
        distance_matrix, samples = make_matrix(input_fname)
        dispatch_output(distance_matrix, samples, output_format, output_fname)

        print(
            f"Done. The dissimilarity matrix has been saved in {output_fname}")

    else:
        raise RuntimeError('Failed to open the input file.')


if __name__ == "__main__":
    main()
