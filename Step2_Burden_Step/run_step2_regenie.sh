#!/bin/bash
pheno=${1}
catvars=${2}
output=${3}
chunk=${4}
ncon=${5}
masks=${6}
annotations=${7}
mafs=${8}
invnorm=${9}
conditional=${10}
binary=${11}
chunk_boundaries=${12}

# create output filename set
output_name="${output}/${pheno}_chr${chr}_${chunk_boundaries}_${masks}_burden"
if [ $conditional == "No" ];then
  output_name=${output_name}"_unconditional"
elif [ $conditional == "all" ];then
  output_name=${output_name}"_all_conditional"
fi

# build basic command
com="./regenie_v3.0.gz_x86_64_Linux  \
  --step 2 \
  --pgen chr${chr}_for_regenie \
  --phenoFile Phenotype_TSV \
  --covarFile Covariate_TSV \
  --bsize 400 \
  --minMAC=5 \
  --catCovarList ${catvars} \
  --maxCatLevels 30 \
  --par-region b38 \
  --vc-tests acato-full \
  --vc-maxAAF 0.001 \
  --joint acat \
  --extract variants_to_include \
  --anno-file ${annotations} \
  --set-list set_list_temp \
  --mask-def ${masks} \
  --threads 10 \
  --pred ${pheno}_Step1_pred.list \
  --aaf-bins ${mafs} \
  --check-burden-files \
  --out ${output_name}"

if [ $ncon -ne 0 ] && [ $conditional == "common" ]; then # if there are common variants to condition on
  com=${com}" \
    --condition-list conditional_list_local"
elif [ $ncon -ne 0 ] && [ $conditional == "all" ];then # if we're conditioning on all variants
  com=${com}" \
    --condition-list conditional_list_final"
fi

if [ "$invnorm" == "RINT" ];then # if inverse normalisation has been asked for
  com=${com}" \
    --apply-rint"
 elif [ "${invnorm}" == "MCC" ]; then
   com=${com}" \
    --mcc --mcc-skew"
fi

if [ $binary == "Binary" ];then
  com=${com}" \
  --bt \
  --firth --firth-se --approx \
  --pThresh 0.01"
fi

$com
