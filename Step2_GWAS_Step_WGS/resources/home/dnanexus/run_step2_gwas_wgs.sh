#!/bin/bash
pheno=${1}
catvars=${2}
output=${3}
chr=${4}
invnorm=${5}
ncon=${6}
output_prefix=${7}
mac=${8}
struct=${9}
binary=${10}


# build main command
command="
./regenie  \
--step 2 \
--bed chr${chr}_for_regenie \
--out ${output}/${output_prefix} \
--phenoFile Phenotype_TSV \
--covarFile Covariate_TSV \
--bsize 400 \
--catCovarList ${catvars} \
--minMAC 1 \
--par-region b38 \
--maxCatLevels 30 \
--pred ${pheno}_Step1_pred.list 
"
# and add in the optional parts if they're required
# inverse normalise
if [ "${invnorm}" == "RINT" ]; then
	command=${command}" \
	--apply-rint
	"

elif [ "${invnorm}" == "MCC" ]; then
	command=${command}" \
	--mcc --mcc-skew 
	"
fi
# conditional SNPs
if [ $ncon -gt 0 ]; then
command=${command}" \
--condition-list conditional_list_local 
"
fi


if [ "${binary}" == "Yes" ]; then
command=${command}" \
--bt \
--firth --approx \
--firth-se --pThresh 0.01 \
"
fi




# run command
${command}




