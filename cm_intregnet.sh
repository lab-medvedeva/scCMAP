./sc-intregnet.sh --sc-intregnet-path "$INTREGNET_DIR" --tfs-in-combo 1 --start-sample-name "${start-type}"
--start-exp "$CMAP_INPUT_DIR/${start-type}_aggregated_exp.tsv" --target-sample-name "${target_type}"
--target-exp-ranking "$CMAP_INPUT_DIR/${start-type}_to_${target_type}_TFs_in_DEGs_promoters.txt"



./sc-intregnet.sh --sc-intregnet-path "$INTREGNET_DIR" --tfs-in-combo 1             --start-sample-name "${target_type}"             --start-exp "$CMAP_INPUT_DIR/${target_type}_agg_exp.tsv"             --target-sample-name "${start-type}"             --target-exp-ranking "$CMAP_INPUT_DIR/Panc_to_${start-type}_TFs_in_DEGs_promoters.txt"


export TF_path="${TFS_PRED_DIR}/${initial_celltype}_to_${SAMPLE_NAME}_FinalVals.tsv"
python3 $WORKDIR/scripts/main/cmap.py --TF_path "$TF_path" > $TFS_PRED_DIR/cmap.csv
