## Output
All output files will be placed in a directory named `results` that will be created the first time you run bam2bakR. The output of greatest interest, the gzipped cB.csv file, will be in `results/cB/`. Columns that can be kept in the final cB (see the keepcols option in the config to choose among these options) are:

  * qname (read ID)
  * nA (number of As in read)
  * nC (number of Cs)
  * nT (number of Ts)
  * nG (number of Gs)
  * rname (chromosome name)
  * GF (ENSEMBL ID of gene to which read maps)
  * EF (ENSEMBL ID of gene to which read maps if read overlaps any exonic region)
  * XF (ENSEMBL ID of gene to which read maps if read only overlaps with exonic regions)
  * FR (Strandedness of read; F = forward, R = reverse. Only will have F if single-end sequencing)
  * sj (TRUE if read overlaps a splice junction)
  * ai (TRUE if read overlaps any intronic region)
  * io (TRUE if read exclusively overlaps an intronic region)
  * ei (TRUE if read maps to intronic and exonic regions)
  * TA (number of T-to-A mutations)
  * CA (number of C-to-A mutations)
  * GA (number of G-to-A mutations)
  * AT (number of A-to-T mutations)
  * CT (number of C-to-T mutations)
  * GT (number of G-to-T mutations)
  * NT (number of N-to-T mutations, where N is any nucleotide)
  * AC (number of A-to-C mutations)
  * TC (number of T-to-C mutations)
  * GC (number of G-to-C mutations)
  * NC (number of N-to-C mutations, where N is any nucleotide)
  * AG (number of A-to-G mutations)
  * TG (number of T-to-G mutations)
  * CG (number of C-to-G mutations)
  * NG (number of N-to-G mutations, where N is any nucleotide)
  * AN (number of A-to-N mutations, where N is any nucleotide)
  * TN (number of T-to-N mutations, where N is any nucleotide)
  * CN (number of C-to-N mutations, where N is any nucleotide)
  * GN (number of G-to-N mutations, where N is any nucleotide)

The tdf files to make color-coded tracks are in: `results/tracks/`.

Other output includes:

* Sorted and filtered bam files in `results/sf_reads/`
* HTseq output text and bam files in `results/htseq/`
* SNP calls in `results/snps/`
* .csv files with counts of all mutation types in `results/counts/`
* Scale factors calculated with edgeR in `results/normalization/`