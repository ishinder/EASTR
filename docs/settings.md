# Settings

## Required Arguments


* `--gtf` - Input GTF file containing transcript annotations.
* `--bed` - Input BED file with intron coordinates.
* `--bam` - Input BAM file or a TXT file containing a list of BAM files with
  read alignments.
!!! Warning
    Only one of the below input options (`GTF`, `BED`, or `BAM`) should be provided.

* `-r`, `--reference` - Reference FASTA genome used in alignment.
* `-i`, `--bowtie2_index` - Path to Bowtie2 index for the reference genome.

## Optional Arguments

- `--bt2_k` - Minimum number of distinct alignments found by bowtie2 for a
  junction to be considered spurious. **Default**: *10*.
- `-o` - Length of the overhang on either side of the splice junction. **Default**:
  *50*.
- `-a` - Minimum required anchor length in each of the two exons. **Default**: *7*.
- `--min_duplicate_exon_length`- Minimum length that a one-anchor alignment
  shift must meet or exceed to be considered as representing duplicated exons.
  It is used to differentiate between exon duplications and spurious splice
  alignments. **Default**: *27*.
- `--min_junc_score` - Minimum number of supporting spliced reads required per
  junction. **Default**: *1*.
- `--trusted_bed` - Path to a BED file path with trusted junctions, which will
  not be removed by EASTR.
- `--verbose` - Display additional information during BAM filtering, including
  the count of total spliced alignments and removed alignments.
- `--removed_alignments_bam` - Write removed alignments to a BAM file.
- `-p` - Number of parallel processes. **Default**: *1*.

## Minimap2 Parameters

- `-A` - Matching score. **Default**: *3*.
- `-B` - Mismatching penalty. **Default**: *4*.
- `-O` - Gap open penalty. **Default**: *[12, 32]*.
- `-E` - Gap extension penalty. **Default**: *[2, 1]*.
- `-k` - K-mer length for alignment. **Default**: *3*.
- `--scoreN` - Score of a mismatch involving ambiguous bases. **Default**: *1*.
- `-w` - Minimizer window size. **Default**: *2*.
- `-m` - Discard chains with chaining score. **Default**: *25*.

## Output Options

- `--out_original_junctions` - Write original junctions to the output file or
  directory.
- `--out_removed_junctions` - Write removed junctions to the output file or
  directory; the default output is to the terminal.
- `--out_filtered_bam` - Write filtered bams to the output file or directory.
- `--filtered_bam_suffix` - Suffix added to the name of the output BAM files.
  **Default**: *'_EASTR_filtered'*.

## Other arguments

- `-p` - Number of parallel processes. **Default**: *1*.
