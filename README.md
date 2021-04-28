# scCMAP
Computational search for target genes for chemical reprograming of the cells using gene regulatory networks

scINTREGNET can be used to improve results of the [Connectivity Mapping (CMAP)](https://www.broadinstitute.org/connectivity-map-cmap) pipeline for drug repurposing by accounting of Transription-Regulation-Network for target cell population and providing core transcription factors as drug targets.

To start INTREGNET with the Connectivity mapping function, run cm_intregnet script, specifying: Name of work directory, Name of initial cell type, Name of a target cell type, Expression of initial cell, Expression of target cell, List of TFs associated with initial cell type, List of TFs associated with target cell type.

Example: 
./cm_intregnet.sh /home/iuliia/scINTREGNET Fibroblasts Neurons Fibroblast_aggregated_exp.tsv Neurons_agg.tsv Neurons_to_Fibroblasts_TFs_in_DEGs_promoters.txt Fibroblasts_to_Neurons_TFs_in_DEGs_promoters.txt

The final table includes chemicals, the signature of which mimics the cellular conversion from initial to target cell type. The chemicals are ranked by a score which represents an overlap between the obtained core transcription factors and signature genes.
