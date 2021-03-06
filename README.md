# mlst2dist: computes similarity matrices from chewBBACA MLST alleles tables

(https://github.com/B-UMMI/chewBBACA).

mlst2dist converts all occurrences of LNF, NIPH, NIPHEM, PLOT3, PLOT5, ALM, ASM, LOTSC in 0 and transforms all INF-X in the corresponding X inferred alleles.
Any other allele-call value left that isn't an integer is converted to 0.
The program then calculates a matrix of pairwise similarities using a Hamming Distance modified with correction for missing data (the allele-calls converted to 0s).

The dissimilarity matrix output is saved to disk with default MEGA lower-left matrix format (tested with MEGA v.11.0.10) or with optionals PHYLIP symmetric square and a raw TSV formats.

## Usage

    $ mlst2dist --help

    usage: mlst2dist.py [-h] [-f {MEG,PHY,TSV}] input output

    Converts a chewBBACA alleles table into a dissimilarity matrix.

    positional arguments:
      input                 chewBBACA alleles table .tsv
      output                distance matrix output filename

    optional arguments:
      -h, --help            show this help message and exit
      -f {MEG,PHY,TSV}, --outfmt {MEG,PHY,TSV}
                            format options for the output matrix [default: MEG (MEGA format)]

### Credits and References

* https://github.com/tseemann/cgmlst-dists as the inception of mlst2dist.py

* https://github.com/B-UMMI/chewBBACA/blob/master/CHEWBBACA/utils/Extract_cgAlleles.py as sanity check for missing-data identification 

* The distance algorithm implemented in mlst2dist.py is described in Galpern P, Manseau, M, Hettinga P, Smith K, and Wilson P. (2012) allelematch: an R package for identifying unique multilocus genotypes where genotype error and missing data may be present. Molecular Ecology Resources 12:771-778



Feedbacks and contributions are welcome!
