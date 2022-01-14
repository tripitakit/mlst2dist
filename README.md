# mlst2dist

This program calculate a distance (dissimilarity) matrix from a chewBBACA results_alleles.tsv table.

It converts all occurrences of LNF, NIPH, PLOT3, PLOT5  with 0, all INF-X in the corresponding X inferred allele-calls and anything else that isn't an integer value into 0; then calculates a matrix of pairwise dissimilarities using the Hamming distances corrected for missing data. The dissimilarity matrix output is saved to disk with default PHYLIP square matrix format, or optional TSV format.


    usage: mlst2dist.py [-h] [-f {PHYLIP,TSV}] input output

    Converts a chewBBACA alleles table into a dissimilarity matrix.

    positional arguments:
      input                 chewBBACA alleles table .tsv
      output                distance matrix output filename

    optional arguments:
      -h, --help            show this help message and exit
      -f {PHYLIP,TSV}, --outfmt {PHYLIP,TSV}
                            format options for the output matrix, default to PHYLIP


## Example

### Input alleles table (test_allele.tsv)

    FILE	G1	G2	G3	G4	G5	G6	G7	G8	G9
    S1	1	INF-2	3	2	1	5	NIPH	ASM	ALM
    S2	1	1	1	1	NIPH	5	1	NIPHEM	INF-1
    S3	INF-1	2	3	4	1	3	INF-1	LOTSC	1
    S4	1	LNF	2	4	1	3	1	LNF	2
    S5	1	2	ASM	2	1	3	2	1	LNF
    S6	2	INF-8	3	PLOT3	PLOT5	3	INF-2	INF-1	INF-1


### Dissimilarity matrix (PHYLIP squared matrix format)

    $ mlst2dist.py test_alleles.tsv test_dist.phy
    $ cat test_dist.phy
    6
    S1        	0.0	0.66667	0.44444	0.66667	0.44444	0.88889
    S2        	0.66667	0.0	0.55556	0.66667	0.88889	0.77778
    S3        	0.44444	0.55556	0.0	0.33333	0.55556	0.66667
    S4        	0.66667	0.66667	0.33333	0.0	0.66667	0.88889
    S5        	0.44444	0.88889	0.55556	0.66667	0.0	0.66667
    S6        	0.88889	0.77778	0.66667	0.88889	0.66667	0.0

### Dissimilarity matrix (TSV format)

    $ mlst2dist.py test_alleles.tsv test_dist.tsv --outfmt TSV
    $ cat test_dist.tsv

      S1	S2	S3	S4	S5	S6
    S1	0.0	0.66667	0.44444	0.66667	0.44444	0.88889
    S2	0.66667	0.0	0.55556	0.66667	0.88889	0.77778
    S3	0.44444	0.55556	0.0	0.33333	0.55556	0.66667
    S4	0.66667	0.66667	0.33333	0.0	0.66667	0.88889
    S5	0.44444	0.88889	0.55556	0.66667	0.0	0.66667
    S6	0.88889	0.77778	0.66667	0.88889	0.66667	0.0


### PHYLIP NEIGHBOR - UPGMA outfile

    $ cat outfile

      6 Populations

    Neighbor-Joining/UPGMA method version 3.697


    UPGMA method

    Negative branch lengths allowed


                +------------S1        
            +---2 
            !   +------------S5        
        +---3 
        !   !      +---------S3        
      +-4   +------1 
      ! !          +---------S4        
    --5 ! 
      ! +--------------------S2        
      ! 
      +----------------------S6        


    From     To            Length          Height
    ----     --            ------          ------
      5        4          0.04167         0.04167
      4        3          0.05556         0.09722
      3        2          0.06945         0.16667
      2     S1            0.22222         0.38889
      2     S5            0.22222         0.38889
      3        1          0.12500         0.22222
      1     S3            0.16667         0.38889
      1     S4            0.16667         0.38889
      4     S2            0.34722         0.38889
      5     S6            0.38889         0.38889


### PHYLIP NEIGHBOR - UPGMA outtree

    $ cat outtree
    ((((S1:0.22222,S5:0.22222):0.06945,(S3:0.16667,S4:0.16667):0.12500):0.05556,
    S2:0.34722):0.04167,S6:0.38889);
