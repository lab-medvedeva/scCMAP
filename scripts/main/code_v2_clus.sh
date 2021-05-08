#!/bin/bash

set -e

check_file_content() {
	local FILE_PATH=$1
	if [[ ! -s $FILE_PATH ]]; then
		echo "[EXIT]: UNABLE TO GENERATE TFs REGULATORY NETWORK (TRN)"
		echo "[REASON]: $FILE_PATH is empty"
		exit 1
	fi
}

if [ -s diff.txt ]; then
	rm -f empty.txt
	touch full.txt
else
	rm -f full.txt
	touch empty.txt
fi

## NOTE:
# This version computes the interaction and binding coordinates
# for core Enh/Pro and first neighbors (e.g. core 2 -> neighbor1 -> core 1)
# where first neighbor is also regulated by core TF

## UTILS

# # Convert to hg38
# python CrossMap.py bed GRCh37_to_GRCh38.chain.gz \
# 	"${target_H3K27ac}" \
# 	"${target_H3K27ac%.bed}_hg38.bed"

# python CrossMap.py bed GRCh37_to_GRCh38.chain.gz \
# 	"${target_H3K4me3}" \
# 	"${target_H3K4me3%.bed}_hg38.bed"

# target_H3K27ac="${target_H3K27ac%.bed}_hg38.bed"
# target_H3K4me3="${target_H3K4me3%.bed}_hg38.bed"

# grep -v -P "\t0.000000" ${target_access_regions}.tmp > $target_access_regions # remove closed regions

## INPUT

# target_H3K27ac
# target_H3K4me3
# target_access_regions
# gene_ranking_file
# coreTFs_file

## DEFAULT FILES

# ${AnimalTFDB_TF}
# ${Processed_GeneHancer_GRCh38}
# ${ChIPseq_for_networks}
# ${GRCh38_promoter}

## OUTPUT

# "${SAMPLE_NAME}_Shell1_Enh_Int_Acc_OnlyTF.txt"
# "${SAMPLE_NAME}_Shell1_Pro_Int_Acc_OnlyTF.txt"

## RUN

cd $OUTPUT_DIR

printf "\n Remove enhancers-gene association with low significance \n"

awk 'FNR==NR {a[$1]; next} ($(NF-1) in a)' $coreTFs_file ${Processed_GeneHancer_GRCh38} >"${SAMPLE_NAME}_Enhancers.tmp"
wc -l "${SAMPLE_NAME}_Enhancers.tmp"
awk '(NR>1) && ($8 > 10 ) ' "${SAMPLE_NAME}_Enhancers.tmp" >filtered_Enh.tmp # Remove those enhancers-gene association where score is less than 10
sort -k1,1 -k2,2n filtered_Enh.tmp >"sorted_${SAMPLE_NAME}_Enhancers.tmp"
wc -l "sorted_${SAMPLE_NAME}_Enhancers.tmp"

printf "\n Checking if target_H3K27ac is provided \n"
# if $target_H3K27ac is provided
if [[ -n $target_H3K27ac ]]; then
	bedtools intersect -a "sorted_${SAMPLE_NAME}_Enhancers.tmp" -b $target_H3K27ac >"${SAMPLE_NAME}_Act_Enhancers.tmp"
else
	cat "sorted_${SAMPLE_NAME}_Enhancers.tmp" >"${SAMPLE_NAME}_Act_Enhancers.tmp"
fi

wc -l "${SAMPLE_NAME}_Act_Enhancers.tmp"
sort -k1,1 -k2,2n "${SAMPLE_NAME}_Act_Enhancers.tmp" >"sorted_${SAMPLE_NAME}_Act_Enhancers.tmp"

printf "\n Looking for "entire" TFBS peak falling completely within Active Enhancer region \n"
# -F 1.0 will look for "entire" TFBS peak falling completely within Active Enhancer region.
bedtools intersect -a "sorted_${SAMPLE_NAME}_Act_Enhancers.tmp" -b ${ChIPseq_for_networks} -F 1.0 -wo -sorted >"${SAMPLE_NAME}_All_TFs_In_Core_Enh.tmp"
wc -l "${SAMPLE_NAME}_All_TFs_In_Core_Enh.tmp"

awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' $coreTFs_file "${SAMPLE_NAME}_All_TFs_In_Core_Enh.tmp" >"${SAMPLE_NAME}_coreTFs_file_In_CoreEnh.tmp"
wc -l "${SAMPLE_NAME}_coreTFs_file_In_CoreEnh.tmp"

#### Getting first neighbours

printf "\n Getting first neighbours within TRN \n"

awk 'FNR==NR {a[$1]; next} !($(NF-2) in a)' $coreTFs_file "${SAMPLE_NAME}_All_TFs_In_Core_Enh.tmp" >NC_TF_In_CoreEnh.tmp
wc -l NC_TF_In_CoreEnh.tmp
cat NC_TF_In_CoreEnh.tmp | awk '{print $12}' | sort -n | uniq >NC_Uniq_TFs_Binding_CoreEnh.tmp
wc -l NC_Uniq_TFs_Binding_CoreEnh.tmp
awk 'FNR==NR {a[$1]; next} ($(NF-1) in a)' NC_Uniq_TFs_Binding_CoreEnh.tmp ${Processed_GeneHancer_GRCh38} >NC_TFs_Enh.tmp
awk '(NR>1) && ($8 > 10)' NC_TFs_Enh.tmp >Filtered_NC_TFs_Enh.tmp # Remove those enhancers-gene association where score is less than 10
wc -l NC_TFs_Enh.tmp
wc -l Filtered_NC_TFs_Enh.tmp

# get active non-core TFs enhancers which are first neigbors

printf "\n Get active non-core TFs enhancers which are first neigbors \n"

sort -k1,1 -k2,2n Filtered_NC_TFs_Enh.tmp >Filtered_NC_TFs_Enh_Sorted.tmp

printf "\n Checking if target_H3K27ac is provided \n"
# if $target_H3K27ac is provided
if [[ -n $target_H3K27ac ]]; then
	bedtools intersect -a Filtered_NC_TFs_Enh_Sorted.tmp -b $target_H3K27ac >Filtered_NC_TFs_Act_Enh.tmp
else
	cat Filtered_NC_TFs_Enh_Sorted.tmp >Filtered_NC_TFs_Act_Enh.tmp
fi

wc -l Filtered_NC_TFs_Act_Enh.tmp
sort -k1,1 -k2,2n Filtered_NC_TFs_Act_Enh.tmp >Filtered_NC_TFs_Act_Enh_Sorted.tmp

# get TF bindings in non-core active enhancers

printf "\n Get TF bindings in non-core active enhancers \n"

bedtools intersect -a Filtered_NC_TFs_Act_Enh_Sorted.tmp -b ${ChIPseq_for_networks} -F 1.0 -wo -sorted >"${SAMPLE_NAME}_All_TFs_In_NonCore_Enh.tmp"
wc -l "${SAMPLE_NAME}_All_TFs_In_NonCore_Enh.tmp"
awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' $coreTFs_file "${SAMPLE_NAME}_All_TFs_In_NonCore_Enh.tmp" >"${SAMPLE_NAME}_Final_coreTFs_file_In_NonCore_Enh.tmp"
wc -l "${SAMPLE_NAME}_Final_coreTFs_file_In_NonCore_Enh.tmp"

# get non-core TF names whose enhancers are bound by core enhancers

printf "\n Get non-core TF names whose enhancers are bound by core enhancers \n"

cat "${SAMPLE_NAME}_Final_coreTFs_file_In_NonCore_Enh.tmp" | awk '{print $7}' | sort -n | uniq >NC_Uniq_TF_Names_Bound_By_CoreTF.tmp
wc -l NC_Uniq_TF_Names_Bound_By_CoreTF.tmp
awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' NC_Uniq_TF_Names_Bound_By_CoreTF.tmp "${SAMPLE_NAME}_All_TFs_In_Core_Enh.tmp" >"${SAMPLE_NAME}_Final_NoncoreTFs_file_In_Core_Enh.tmp"
wc -l "${SAMPLE_NAME}_Final_NoncoreTFs_file_In_Core_Enh.tmp"

########################
# Uncomment this chunk if you want to get interactions among non-core TFs as well (e.g. Neigbor1 -> Neighbor 3)
# Comment the chunk below in this case
# Chunk Start
awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' NC_Uniq_TF_Names_Bound_By_CoreTF.tmp "${SAMPLE_NAME}_All_TFs_In_NonCore_Enh.tmp" >"${SAMPLE_NAME}_Final_NoncoreTFs_file_In_NonCore_Enh.tmp"
wc -l "${SAMPLE_NAME}_Final_NoncoreTFs_file_In_NonCore_Enh.tmp"
# concatenating all into one
cat "${SAMPLE_NAME}_coreTFs_file_In_CoreEnh.tmp" "${SAMPLE_NAME}_Final_NoncoreTFs_file_In_Core_Enh.tmp" "${SAMPLE_NAME}_Final_coreTFs_file_In_NonCore_Enh.tmp" "${SAMPLE_NAME}_Final_NoncoreTFs_file_In_NonCore_Enh.tmp" >shell1_Enh.tmp
wc -l shell1_Enh.tmp
# Chunk End

################
# Uncomment this chunk if you DO NOT want to get interactions among non-core TFs (e.g. Neigbor1 -> Neighbor 3)
# comment the above chunk
# Chunk Start
# concatenating all into one
#cat "${SAMPLE_NAME}_coreTFs_file_In_CoreEnh.tmp" "${SAMPLE_NAME}_Final_NoncoreTFs_file_In_Core_Enh.tmp" "${SAMPLE_NAME}_Final_coreTFs_file_In_NonCore_Enh.tmp" > shell1_Enh.tmp
#wc -l shell1_Enh.tmp
# Chunk End

################

# removing TFs which are not expressed in the given cell type

printf "\n Removing TFs which are not expressed in the given cell type \n"

awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' $gene_ranking_file shell1_Enh.tmp >shell1_Enh_expTF.tmp
wc -l shell1_Enh_expTF.tmp
awk 'FNR==NR {a[$1]; next} ($(NF-7) in a)' $gene_ranking_file shell1_Enh_expTF.tmp >shell1_Enh_expTF_2.tmp
wc -l shell1_Enh_expTF_2.tmp

# processing the files to give output in the desired format (ChIP-peak, active enhancer region)

printf "\n Processing the files to give output in the desired format (ChIP-peak, active enhancer region) \n"

awk 'BEGIN {FS=OFS="\t"} {print $9,$10,$11,$12,$13,$1,$2,$3,$4,$5,$6,$7,$8,$14}' shell1_Enh_expTF_2.tmp >"${SAMPLE_NAME}_Shell1_Enh.bed"
wc -l "${SAMPLE_NAME}_Shell1_Enh.bed"

# getting interactions which are in active regulatory regions (accessibility not applied yet)
#cat "${SAMPLE_NAME}_Shell1_Enh.bed" | awk '{print $4"\t"$12}' | sort -n | uniq > "${SAMPLE_NAME}_Shell1_Enh_Int.txt"
#wc -l "${SAMPLE_NAME}_Shell1_Enh_Int.txt"

# Getting accessible bed files and interactions for shell1
# -f 0.5 means 50% of Enhancer Region (ChIP peak in Enh.) should be falling within the accessibility peak

printf "\n Checking if target_access_regions is provided \n"
# if $target_access_regions is provided
if [[ -n $target_access_regions ]]; then
	bedtools intersect -a "${SAMPLE_NAME}_Shell1_Enh.bed" -b $target_access_regions -f 0.50 -F 0.50 -e -wa >"${SAMPLE_NAME}_Shell1_Enh_Acc.bed"
else
	cat "${SAMPLE_NAME}_Shell1_Enh.bed" >"${SAMPLE_NAME}_Shell1_Enh_Acc.bed"
fi

# remove non-TF rows from the bed file

printf "\n Removing non-TF rows from the bed file \n"

awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' ${AnimalTFDB_TF} "${SAMPLE_NAME}_Shell1_Enh_Acc.bed" >"${SAMPLE_NAME}_Shell1_Enh_Acc_OnlyTF.tmp"
awk 'FNR==NR {a[$1]; next} ($(NF-10) in a)' ${AnimalTFDB_TF} "${SAMPLE_NAME}_Shell1_Enh_Acc_OnlyTF.tmp" >"${SAMPLE_NAME}_Shell1_Enh_Acc_OnlyTF.bed"
wc -l "${SAMPLE_NAME}_Shell1_Enh_Acc_OnlyTF.bed"
cat "${SAMPLE_NAME}_Shell1_Enh_Acc_OnlyTF.bed" | awk '{print $4"\t"$12}' | sort -n | uniq >"${SAMPLE_NAME}_Shell1_Enh_Int_Acc_OnlyTF.txt"
wc -l "${SAMPLE_NAME}_Shell1_Enh_Int_Acc_OnlyTF.txt"

# remove for only the TF genes (originally I also had CRM and co-factors in TF gene list)
#awk 'FNR==NR {a[$1]; next} ($(NF) in a)' ./AnimalTFDB_TF.txt "${SAMPLE_NAME}_Shell1_Enh_Int_Acc.tmp" > "${SAMPLE_NAME}_Shell1_Enh_Int_Acc_onlyTF.tmp"
#awk 'FNR==NR {a[$1]; next} ($(NF-1) in a)' ./AnimalTFDB_TF.txt "${SAMPLE_NAME}_Shell1_Enh_Int_Acc_onlyTF.tmp" > "${SAMPLE_NAME}_Shell1_Enh_Int_Acc.txt"
#wc -l "${SAMPLE_NAME}_Shell1_Enh_Int_Acc.txt"

# getting bed files and interactions for core TFs only

printf "\n Getting bed files and interactions for core TFs only \n"

: <<'end_long_comment'
awk 'BEGIN {FS=OFS="\t"} {print $9,$10,$11,$12,$13,$1,$2,$3,$4,$5,$6,$7,$8,$14}' "${SAMPLE_NAME}_coreTFs_file_In_CoreEnh.tmp" > "${SAMPLE_NAME}_core_Enh.bed"
wc -l "${SAMPLE_NAME}_core_Enh.bed"
cat "${SAMPLE_NAME}_core_Enh.bed" | awk '{print $4,$12}' | sort -n | uniq > "${SAMPLE_NAME}_core_Enh_Int.txt"
wc -l "${SAMPLE_NAME}_core_Enh_Int.txt"
end_long_comment

rm -f *.tmp

#######################
###### Promoters ######
#######################

#: <<'end_long_comment'
awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' $coreTFs_file ${GRCh38_promoter} >"${SAMPLE_NAME}_Promoters.tmp"
sort -k1,1 -k2,2n "${SAMPLE_NAME}_Promoters.tmp" >"sorted_${SAMPLE_NAME}_Promoters.tmp"
wc -l "sorted_${SAMPLE_NAME}_Promoters.tmp"

# consider the entire promoter for ChIP-scan if only a single K4me3 peak is falling within this promoter region (-u option)

printf "\n Checking if target_H3K4me3 is provided \n"
# if $target_H3K4me3 is provided
if [[ -n $target_H3K4me3 ]]; then
	bedtools intersect -a "sorted_${SAMPLE_NAME}_Promoters.tmp" -b $target_H3K4me3 -u >"${SAMPLE_NAME}_Act_Pro.tmp"
else
	cat "sorted_${SAMPLE_NAME}_Promoters.tmp" >"${SAMPLE_NAME}_Act_Pro.tmp"
fi

wc -l "${SAMPLE_NAME}_Act_Pro.tmp"

printf "\n Checking whether entire ChIP-peak should fall within Active promoter region \n"
# -F 1.0 means entire ChIP-peak should fall within Active promoter region
bedtools intersect -a "${SAMPLE_NAME}_Act_Pro.tmp" -b ${ChIPseq_for_networks} -F 1.0 -wo -sorted >"${SAMPLE_NAME}_All_TFs_In_Core_Pro.tmp"
wc -l "${SAMPLE_NAME}_All_TFs_In_Core_Pro.tmp"

awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' $coreTFs_file "${SAMPLE_NAME}_All_TFs_In_Core_Pro.tmp" >"${SAMPLE_NAME}_coreTFs_file_In_Core_Pro.tmp"
wc -l "${SAMPLE_NAME}_coreTFs_file_In_Core_Pro.tmp"

#### Getting first neighbours ####

# get unique non-core TFs bound in core promoters

printf "\n Get unique non-core TFs bound in core promoters \n"

awk 'FNR==NR {a[$1]; next} !($(NF-2) in a)' $coreTFs_file "${SAMPLE_NAME}_All_TFs_In_Core_Pro.tmp" >NC_TF_In_CorePro.tmp
wc -l NC_TF_In_CorePro.tmp
cat NC_TF_In_CorePro.tmp | awk '{print $10}' | sort -n | uniq >NC_Uniq_TF_Binding_CorePro.tmp
wc -l NC_Uniq_TF_Binding_CorePro.tmp

# get active promoter regions of non-core TFs

printf "\n Get active promoter regions of non-core TFs \n"

awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' NC_Uniq_TF_Binding_CorePro.tmp ${GRCh38_promoter} >NC_TFs_Promoters.tmp
wc -l NC_TFs_Promoters.tmp
sort -k1,1 -k2,2n NC_TFs_Promoters.tmp >sorted_NC_TFs_Promoters.tmp

printf "\n Checking if target_H3K4me3 is provided \n"
# if $target_H3K4me3 is provided
if [[ -n $target_H3K4me3 ]]; then
	bedtools intersect -a sorted_NC_TFs_Promoters.tmp -b $target_H3K4me3 -u >NC_TFs_Act_Pro.tmp
else
	cat sorted_NC_TFs_Promoters.tmp >NC_TFs_Act_Pro.tmp
fi

wc -l NC_TFs_Act_Pro.tmp

# scan non-core promoters for core-TF bindings and consider only those which are bound by core-TFs.

printf "\n Scan non-core promoters for core-TF bindings and consider only those which are bound by core-TFs \n"

bedtools intersect -a NC_TFs_Act_Pro.tmp -b ${ChIPseq_for_networks} -F 1.0 -wo -sorted >"${SAMPLE_NAME}_All_TFs_In_NonCore_Pro.tmp"
wc -l "${SAMPLE_NAME}_All_TFs_In_NonCore_Pro.tmp"
awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' $coreTFs_file "${SAMPLE_NAME}_All_TFs_In_NonCore_Pro.tmp" >"${SAMPLE_NAME}_Final_coreTFs_file_In_NonCore_Pro.tmp"
wc -l "${SAMPLE_NAME}_Final_coreTFs_file_In_NonCore_Pro.tmp"

# unique non-core TFs whose promoters are bound by core-TFs

printf "\n Get unique non-core TFs whose promoters are bound by core-TFs \n"

cat "${SAMPLE_NAME}_Final_coreTFs_file_In_NonCore_Pro.tmp" | awk '{print $4}' | sort -n | uniq >NC_Uniq_TF_Names_Bound_By_CoreTF.tmp
wc -l NC_Uniq_TF_Names_Bound_By_CoreTF.tmp

# get core-TFs promoters bound by non-core TFs (which in turn are regulated by the core TFs)

printf "\n Get core-TFs promoters bound by non-core TFs (which in turn are regulated by the core TFs) \n"

awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' NC_Uniq_TF_Names_Bound_By_CoreTF.tmp "${SAMPLE_NAME}_All_TFs_In_Core_Pro.tmp" >"${SAMPLE_NAME}_Final_NoncoreTFs_file_In_Core_Pro.tmp"
wc -l "${SAMPLE_NAME}_Final_NoncoreTFs_file_In_Core_Pro.tmp"

################

# Uncomment this chunk if you want to get interactions among non-core TFs as well (e.g. Neigbor1 -> Neighbor 3)
# Chunk Start
awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' NC_Uniq_TF_Names_Bound_By_CoreTF.tmp "${SAMPLE_NAME}_All_TFs_In_NonCore_Pro.tmp" >"${SAMPLE_NAME}_Final_NoncoreTFs_file_In_NonCore_Pro.tmp"
wc -l "${SAMPLE_NAME}_Final_NoncoreTFs_file_In_NonCore_Pro.tmp"
# concatenating all into one
cat "${SAMPLE_NAME}_coreTFs_file_In_Core_Pro.tmp" "${SAMPLE_NAME}_Final_NoncoreTFs_file_In_Core_Pro.tmp" "${SAMPLE_NAME}_Final_coreTFs_file_In_NonCore_Pro.tmp" "${SAMPLE_NAME}_Final_NoncoreTFs_file_In_NonCore_Pro.tmp" >shell1_Pro.tmp
wc -l shell1_Pro.tmp
# comment out the below chunk
# Chunk End

################

# Uncomment this chunk if you DO NOT want to get interactions among non-core TFs (e.g. Neigbor1 -> Neighbor 3)
# comment the above chunk
# Chunk Start
# concatenating all into one
#cat "${SAMPLE_NAME}_coreTFs_file_In_Core_Pro.tmp" "${SAMPLE_NAME}_Final_NoncoreTFs_file_In_Core_Pro.tmp" "${SAMPLE_NAME}_Final_coreTFs_file_In_NonCore_Pro.tmp" > shell1_Pro.tmp
#wc -l shell1_Pro.tmp
# Chunk End
################

# removing TFs which are not expressed in the given cell type

printf "\n Removing TFs which are not expressed in the given cell type \n"

awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' $gene_ranking_file shell1_Pro.tmp >shell1_pro_expTF.tmp
wc -l shell1_pro_expTF.tmp
awk 'FNR==NR {a[$1]; next} ($(NF-8) in a)' $gene_ranking_file shell1_pro_expTF.tmp >shell1_pro_expTF_2.tmp
wc -l shell1_pro_expTF_2.tmp

# processing the files to give output in the desired format (ChIP-peak, active enhancer region)

printf "\n Processing the files to give output in the desired format (ChIP-peak, active enhancer region) \n"

awk 'BEGIN {FS=OFS="\t"} {print $7,$8,$9,$10,$11,$1,$2,$3,$4,$5,$6,$12}' shell1_pro_expTF_2.tmp >"${SAMPLE_NAME}_Shell1_Pro.bed"
wc -l "${SAMPLE_NAME}_Shell1_Pro.bed"
# get interactions among active regulatory regions (accessibility not used yet)
#cat "${SAMPLE_NAME}_Shell1_Pro.bed" | awk '{print $4"\t"$9}' | sort -n | uniq > "${SAMPLE_NAME}_Shell1_Pro_Int.txt"
#wc -l "${SAMPLE_NAME}_Shell1_Pro_Int.txt"

# Remove inaccessible interactions

printf "\n Remove inaccessible interactions \n"

# if $target_access_regions is provided
if [[ -n $target_access_regions ]]; then
	bedtools intersect -a "${SAMPLE_NAME}_Shell1_Pro.bed" -b "$target_access_regions" -f 0.50 -F 0.50 -e -wa >"${SAMPLE_NAME}_Shell1_Pro_Acc.bed"
else
	cat "${SAMPLE_NAME}_Shell1_Pro.bed" >"${SAMPLE_NAME}_Shell1_Pro_Acc.bed"
fi

wc -l "${SAMPLE_NAME}_Shell1_Pro_Acc.bed"

# get output columns
#cat "${SAMPLE_NAME}_Shell1_Pro_Acc.bed" | awk '{print $4"\t"$9}' | sort -n | uniq > "${SAMPLE_NAME}_Shell1_Pro_Int_Acc.tmp"
#wc -l "${SAMPLE_NAME}_Shell1_Pro_Int_Acc.tmp"

# remove for only the TF genes (originally I also had CRM and co-factors in TF gene list)

printf "\n Remove for only the TF genes (originally I also had CRM and co-factors in TF gene list) \n"

awk 'FNR==NR {a[$1]; next} ($(NF-3) in a)' ${AnimalTFDB_TF} "${SAMPLE_NAME}_Shell1_Pro_Acc.bed" >"${SAMPLE_NAME}_Shell1_Pro_Acc_OnlyTF.tmp"
awk 'FNR==NR {a[$1]; next} ($(NF-8) in a)' ${AnimalTFDB_TF} "${SAMPLE_NAME}_Shell1_Pro_Acc_OnlyTF.tmp" >"${SAMPLE_NAME}_Shell1_Pro_Acc_OnlyTF.bed"
wc -l "${SAMPLE_NAME}_Shell1_Pro_Acc_OnlyTF.bed"
cat "${SAMPLE_NAME}_Shell1_Pro_Acc_OnlyTF.bed" | awk '{print $4"\t"$9}' | sort -n | uniq >"${SAMPLE_NAME}_Shell1_Pro_Int_Acc_OnlyTF.txt"
wc -l "${SAMPLE_NAME}_Shell1_Pro_Int_Acc_OnlyTF.txt"

# remove for only the TF genes (originally I also had CRM and co-factors in TF gene list)
#awk 'FNR==NR {a[$1]; next} ($(NF) in a)' ./AnimalTFDB_TF.txt "${SAMPLE_NAME}_Shell1_Pro_Int_Acc.tmp" > "${SAMPLE_NAME}_Shell1_Pro_Int_Acc_onlyTF.tmp"
#awk 'FNR==NR {a[$1]; next} ($(NF-1) in a)' ./AnimalTFDB_TF.txt "${SAMPLE_NAME}_Shell1_Pro_Int_Acc_onlyTF.tmp" > "${SAMPLE_NAME}_Shell1_Pro_Int_Acc.txt"
#wc -l "${SAMPLE_NAME}_Shell1_Pro_Int_Acc.txt"

# getting bed files and interactions for core TFs only

printf "\n Getting bed files and interactions for core TFs only \n"

: <<'end_long_comment'
awk 'BEGIN {FS=OFS="\t"} {print $7,$8,$9,$10,$11,$1,$2,$3,$4,$5,$6,$12}' "${SAMPLE_NAME}_coreTFs_file_In_Core_Pro.tmp" > "${SAMPLE_NAME}_core_Pro.bed"
wc -l "${SAMPLE_NAME}_core_Pro.bed"
# get output columns (Source Target)
cat "${SAMPLE_NAME}_core_Pro.bed" | awk '{print $4"\t"$9}' | sort -n | uniq > "${SAMPLE_NAME}_core_Pro_Int.txt"
wc -l "${SAMPLE_NAME}_core_Pro_Int.txt"
end_long_comment

rm -f *.tmp

echo "Interactions common between enhancers and promoters"
comm -12 "${SAMPLE_NAME}_Shell1_Enh_Int_Acc_OnlyTF.txt" "${SAMPLE_NAME}_Shell1_Pro_Int_Acc_OnlyTF.txt" | wc -l
echo "unique source in enh"
awk '{print $1}' "${SAMPLE_NAME}_Shell1_Enh_Int_Acc_OnlyTF.txt" | sort -n | uniq | wc -l
echo "unique source in pro"
awk '{print $1}' "${SAMPLE_NAME}_Shell1_Pro_Int_Acc_OnlyTF.txt" | sort -n | uniq | wc -l

rm -f *_Int_Acc_OnlyTF.txt
rm -f *_Acc.bed

# Check Output Files
check_file_content "${SAMPLE_NAME}_Shell1_Pro_Acc_OnlyTF.bed"
check_file_content "${SAMPLE_NAME}_Shell1_Enh_Acc_OnlyTF.bed"

cd $WORKDIR
