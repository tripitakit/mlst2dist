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
    parser.add_argument("-f", "--outfmt", choices=["PHYLIP", "TSV", "CSV"], default="PHYLIP",
                        help="format options for the output matrix, default to PHYLIP")
    args = parser.parse_args()
    return args


def sanitize(line):
    sanitized_line = line.replace("\n", "").split("\t")
    for index, a in list(enumerate(sanitized_line))[1:]:
        sanitized_allele = ""
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
    return(sanitized_line)


def transform_alleles_table(input):
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
    return(samples_vector, loci_vector, allele_matrix)


def make_matrix(input, output, outfmt):
    samples, loci, alleles_matrix = transform_alleles_table(
        input)

    nloci = len(loci)
    nsamp = len(samples)

    if outfmt == "PHYLIP":
        header = f"{nsamp}\n"
    if outfmt == "TSV":
        header = "\t" + "\t".join(samples) + "\n"

    matrix_output_body_text = ""
    for samp_i in range(nsamp):
        samp_a = samples[samp_i]
        alleles_vector_a = alleles_matrix[samp_i]

        if outfmt == "PHYLIP":
            samp_line = f"{samp_a.ljust(10)}"

        if outfmt == "TSV":
            samp_line = f"{samp_a}"

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

            corrected_hamming_dist_a_b = round(1 - (
                n_matches/nloci) + (n_missing/(2*nloci)), 5)

            samp_line += f"\t{corrected_hamming_dist_a_b}"
        matrix_output_body_text += f"{samp_line}\n"

    # write the output matrix file
    outh = open(output, "w")
    outh.write(header)
    outh.write(matrix_output_body_text)
    outh.close()


def validate_infile(input):
    return os.path.exists(input)


def main():
    inputs = parse_args()
    input = inputs.input
    output = inputs.output
    outfmt = inputs.outfmt

    if validate_infile(input):
        make_matrix(input, output, outfmt)
        print(f"Done. The dissimilarity matrix has been saved in {output}")

    else:
        raise RuntimeError('Failed to open the input file.')


if __name__ == "__main__":
    main()
