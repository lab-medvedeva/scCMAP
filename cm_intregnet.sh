INTREGNET_DIR="$1"
CMAP_INPUT_DIR="$1/data/CM_exp_data"

./sc-intregnet.sh --sc-intregnet-path "$1" --tfs-in-combo 1             --start-sample-name "$2"             --start-exp "$CMAP_INPUT_DIR/$4"             --target-sample-name "$3"             --target-exp-ranking "$CMAP_INPUT_DIR/$6"

export TF_path1="$1/data/$2/tfs_predictions/$3_to_$2_FinalVals.tsv"

./sc-intregnet.sh --sc-intregnet-path "$1" --tfs-in-combo 1             --start-sample-name "$3"             --start-exp "$CMAP_INPUT_DIR/$5"             --target-sample-name "$2"             --target-exp-ranking "$CMAP_INPUT_DIR/$7"

export TF_path2="$1/data/$3/tfs_predictions/$2_to_$3_FinalVals.tsv"
python3 $1/scripts/main/cmap2.py --TF_path1 "$TF_path1" --TF_path2 "$TF_path2"> cmap3.csv
