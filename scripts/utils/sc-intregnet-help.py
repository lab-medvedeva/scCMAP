#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse

example_text = """Example usage:

START_INPUT_DIR="$INTREGNET_DIR/data/Dermal_fibroblast/input_files"
TARGET_INPUT_DIR="$INTREGNET_DIR/data/HSC/input_files"

sc-intregnet --tfs-in-combo 6 \\
             --start-sample-name "Dermal_fibroblast" \\
             --start-exp "$START_INPUT_DIR/Dermal_fibroblast_GSM1865616_exp.tsv" \\
             --start-h3k27ac "$START_INPUT_DIR/Dermal_fibroblast_H3K27ac_hg38.bed" \\
             --start-h3k4me3 "$START_INPUT_DIR/Dermal_fibroblast_H3K4me3_hg38.bed" \\
             --target-sample-name "HSC" \\
             --target-exp "$TARGET_INPUT_DIR/HSC_bulk_exp.tsv" \\
             --target-access-regs "$TARGET_INPUT_DIR/HSC_scATAC.bed" \\
             --target-h3k27ac "$TARGET_INPUT_DIR/Hematopoietic_Stem_Cells_H3K27ac_hg38.bed" \\
             --target-h3k4me3 "$TARGET_INPUT_DIR/Hematopoietic_Stem_Cells_H3K4me3_hg38.bed" \\
             --target-specific-chip-seq "$TARGET_INPUT_DIR/Hematopoietic_Stem_Cells_ChiPseq_hg38.bed" \\
             --custom-output-dir "$INTREGNET_DIR/data/HSC/run_on_cell_spec_chipseq" > ./logs.txt
"""

parser = argparse.ArgumentParser(
    prog="scINTREGNET",
    description="scINTREGNET - method for inferring transcription factor combinations for desired cell conversions",
    epilog=example_text,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser_required = parser.add_argument_group("Required named arguments")

parser.add_argument(
    "--sc-intregnet-path",
    type=str,
    required=True,
    help="Path to scINTREGNET core directory",
    default="~/sc-intregnet",
)
parser_required.add_argument(
    "--tfs-in-combo",
    type=int,
    required=True,
    help="Expected number of transcription factors in a combination",
)
parser_required.add_argument(
    "--start-sample-name",
    type=str,
    required=True,
    help="Name for initial/starting cell identity",
)
parser_required.add_argument(
    "--start-exp",
    type=str,
    required=True,
    help="Processed expression file for initial/starting cell identity",
)
parser.add_argument(
    "--start-h3k27ac",
    type=str,
    required=True,
    help="Active enhancer (H3K27ac) peak BED file for initial/starting cell identity",
)
parser.add_argument(
    "--start-h3k4me3",
    type=str,
    required=True,
    help="Active promoter (H3K4me3) peak BED file initial/starting cell identity",
)
parser_required.add_argument(
    "--target-sample-name",
    type=str,
    required=True,
    help="Name for target cell identity",
)
parser_required.add_argument(
    "--target-exp",
    type=str,
    required=True,
    help="Processed expression file for target cell identity. Supress '--target-exp-ranking' argument",
)
parser_required.add_argument(
    "--target-exp-ranking",
    type=str,
    required=False,
    help="List of DEGs in HGNC gene format ranked by Fold-Change. Supress '--target-exp' argument",
)
parser_required.add_argument(
    "--target-access-regs",
    type=str,
    required=True,
    help="Peak BED file from bulk/single-cell ATAC-seq or DNase-seq for target cell identity",
)
parser.add_argument(
    "--target-h3k27ac",
    type=str,
    required=True,
    help="Active enhancer (H3K27ac) peak BED file for target cell identity",
)
parser.add_argument(
    "--target-h3k4me3",
    type=str,
    required=True,
    help="Active promoter (H3K4me3) peak BED file for target cell identity",
)
parser.add_argument(
    "--target-specific-chip-seq",
    type=str,
    required=False,
    help="BED file for a cell-type specific ChiP-seq peaks from Chip-seq TFs experiments",
)
parser.add_argument(
    "--custom-output-dir",
    type=str,
    required=False,
    help="Path to a custom output directory for temporary files and final TFs combo predictions",
)

parser.parse_args()
