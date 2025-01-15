library(tidyverse)
library(magrittr)
library(vroom)
library(IRanges)
library(GenomicRanges)

args <- commandArgs(TRUE)
sim_start <- args[1]

nsim <- 100
freq_thresh <- 0.01
set.seed(20230912*as.numeric(sim_start))
# file locations
res_dir <- "/slade/projects/UKBB/DNA_Nexus_500k_WGS/Proteomics_50k_3k/"
ens_regions <- paste0("/slade/projects/UKBB/DNA_Nexus_WGS/Proteomics/ENSEMBL_REGULATORY/all_active_regions_by_cell")
noncoding_associations <- paste0(res_dir,"results_summaries/all_aggregates_including_cognate")
#protein_panels <- paste0("/slade/projects/UKBB/DNA_Nexus_500k_WGS/Proteomics_50k_3k/results_summaries/rob_working_area/protein_list_secreted")
protein_panels <- paste0("/slade/projects/UKBB/DNA_Nexus_500k_WGS/Proteomics_50k_3k/results_summaries/proteins_cov_995_max_diff_lt_90pctl")
burden_positions <- paste0("/slade/projects/UKBB/DNA_Nexus_500k_WGS/Proteomics_50k_3k/results_summaries/rob_working_area/mask_lengths_chr*")
pheno_file <- "/slade/projects/UKBB/DNA_Nexus_500k_WGS/Proteomics_50k_3k/phenotypes/proteins_50k_300k_cleaned.txt"

### read in data
# non-coding signals
signals_nc <- vroom(noncoding_associations,col_names=F) %>% dplyr::rename(
        protein_name=X20,
        CHR=X1,
        unit=X3,
    ) %>% 
    select(protein_name,CHR,unit) %>% distinct()
burden_positions_df <- vroom(Sys.glob(burden_positions)) %>% 
    dplyr::rename(CHROM=CHR) %>% 
    mutate(CHROM=ifelse(CHROM=="X",23,CHROM)) %>% 
    mutate(CHROM=as.numeric(CHROM))
burden_positions_df_tmp <- burden_positions_df %>% mutate(MASK=str_split_fixed(MASK,"[.]",2)[,1]) %>%
    group_by(MASK) %>%
    summarise(CHROM=max(CHROM),MIN_POS=min(MIN_POS),MAX_POS=max(MAX_POS)) %>%
    select(CHROM,MASK,MIN_POS,MAX_POS) %>%
    rbind(burden_positions_df)
signals_nc %<>% dplyr::rename(CHROM=CHR) %>%
    left_join(burden_positions_df_tmp,by=c("unit"="MASK","CHROM"="CHROM")) %>%
    filter(!is.na(MIN_POS))
# regions in tissues
active_regions <- vroom(ens_regions,col_names=F)
names(active_regions) <- c("CHROM","START","END","TYPE","CELL")
active_regions %<>% mutate(CHROM=as.numeric(ifelse(CHROM=="X",23,CHROM))) %>% 
    filter(!is.na(CHROM))

tissue_groups <- c(
    "A549"="lung",
    "A673"="muscle",
    "adrenal_gland"="adrenal_gland",
    "aorta"="vessel",
    "astrocyte"="neural",
    "B"="blood",
    "bipolar_neuron"="neural",
    "B_PB"="blood",
    "brain_1"="neural",
    "cardiac_muscle"="heart",
    "CD14_monocyte_1"="blood",
    "CD14_monocyte_PB"="blood",
    "CD4_ab_T"="blood",
    "CD4_ab_T_PB_1"="blood",
    "CD4_ab_T_PB_2"="blood",
    "CD4_ab_T_Th"="blood",
    "CD4_ab_T_VB"="blood",
    "CD4_CD25_ab_Treg_PB"="blood",
    "CD8_ab_T_CB"="blood",
    "CD8_ab_T_PB"="blood",
    "CM_CD4_ab_T_VB"="blood",
    "CMP_CD4_1"="blood",
    "CMP_CD4_2"="blood",
    "CMP_CD4_3"="blood",
    "dermal_fibroblast"="skin",
    "DND_41"="blood",
    "EB_CB"="embryo",
    "EM_CD4_ab_T_PB"="blood",
    "EM_CD8_ab_T_VB"="blood",
    "endocrine_pancreas"="pancreas",
    "endodermal"="endoderm",
    "eosinophil_VB"="blood",
    "EPC_VB"="blood",
    "esophagus"="esophagus",
    "foreskin_fibroblast_2"="skin",
    "foreskin_keratinocyte_1"="skin",
    "foreskin_keratinocyte_2"="skin",
    "foreskin_melanocyte_1"="skin",
    "foreskin_melanocyte_2"="skin",
    "germinal_matrix"="neural",
    "GM12878"="blood",
    "H1_hESC_2"="embryo",
    "H1_hESC_3"="embryo",
    "H9_1"="embryo",
    "HCT116"="intestine",
    "heart"="heart",
    "HeLa_S3"="cervix",
    "hepatocyte"="liver",
    "HepG2"="liver",
    "HSMM"="muscle",
    "HUES48"="embryo",
    "HUES6"="embryo",
    "HUES64"="embryo",
    "HUVEC"="vessel",
    "HUVEC_prol_CB"="vessel",
    "iPS_15b"="Stem_cells",
    "iPS_20b"="Stem_cells",
    "iPS_DF_19_11"="Stem_cells",
    "iPS_DF_6_9"="Stem_cells",
    "K562"="blood",
    "keratinocyte"="skin",
    "kidney"="kidney",
    "large_intestine"="intestine",
    "left_ventricle"="heart",
    "leg_muscle"="muscle",
    "lung_1"="lung",
    "lung_2"="lung",
    "M0_CB"="blood",
    "M0_VB"="blood",
    "M1_CB"="blood",
    "M1_VB"="blood",
    "M2_CB"="blood",
    "M2_VB"="blood",
    "mammary_epithelial_1"="breast",
    "mammary_epithelial_2"="breast",
    "mammary_myoepithelial"="breast",
    "MCF_7"="breast",
    "MM_1S"="blood",
    "monocyte_CB"="blood",
    "monocyte_VB"="blood",
    "mononuclear_PB"="blood",
    "MSC"="blood",
    "MSC_VB"="blood",
    "myotube"="muscle",
    "naive_B_VB"="blood",
    "neuron"="neural",
    "neurosphere_C"="neural",
    "neurosphere_GE"="neural",
    "neutro_myelocyte"="neural",
    "neutrophil_CB"="blood",
    "neutrophil_VB"="blood",
    "NHLF"="lung",
    "NK_PB"="blood",
    "NPC_1"="neural",
    "NPC_2"="neural",
    "NPC_3"="neural",
    "osteoblast"="bone",
    "ovary"="ovary",
    "pancreas"="pancreas",
    "PC_3"="prostate",
    "PC_9"="prostate",
    "placenta"="placenta",
    "psoas_muscle"="muscle",
    "right_atrium"="heart",
    "right_ventricle"="heart",
    "sigmoid_colon"="intestine",
    "SK_N_"="neural",
    "small_intestine_1"="intestine",
    "small_intestine_2"="intestine",
    "spleen"="spleen",
    "stomach_1"="stomach",
    "stomach_2"="stomach",
    "Th17"="blood",
    "thymus_1"="thymus",
    "thymus_2"="thymus",
    "T_PB"="blood",
    "trophoblast"="embryo",
    "trunk_muscle"="muscle",
    "UCSF_4"="embryo"
)
to_exclude <- data.frame(CELL=c("HSMM",
    "HUES48",
    "HUES6",
    "HUES64",
    "NHLF",
    "UCSF_4",
    "H1_hESC_2",
    "H1_hESC_3",
    "H9_1",
    "GM12878",
    "A549",
    "A673",
    "DND_41",
    "HCT116",
    "HeLa_S3",
    "HepG2",
    "K562",
    "MCF_7",
    "MM_1S",
    "PC_3",
    "PC_9",
    "SK_N_",
    "iPS_15b",
    "iPS_20b",
    "iPS_DF_19_11",
    "iPS_DF_6_9"))
active_regions %<>% filter(!(CELL %in% to_exclude$CELL)) %>% 
    mutate(CELL=tissue_groups[CELL])

# read in panel data
panel_names=c("Intracellular or membrane-bound"="membrane",                                                
    "Locally secreted"="secreted",
    "Locally secreted in brain"="secreted_brain",
    "Locally secreted in female tissues"="secreted",
    "Locally secreted in male tissues"="secreted",
    "Locally secreted to extracellular matrix"="secreted",                                   
    "Locally secreted to GI-tract"="secreted",
    "Secreted to blood"="secreted",
    "Secreted - unknown location"="secreted")
panel_data <- vroom(protein_panels,col_names=F,delim=" ") %>%
    rename("X1"="Assay") %>%
    mutate(Panel="covered") %>% 
    mutate(Assay=tolower(Assay)) #%>%
    #mutate(Panel=panel_names[Panel])

# add in panels with no secretion annotation
proteins_all <- vroom(pheno_file) %>%
    names() %>%
    tail(-1) %>%
    as_tibble() %>%
    rename(value="Assay")
panel_data %<>% full_join(proteins_all) %>%
    mutate(Panel=ifelse(is.na(Panel),"unannotated",Panel))


### read in all variant association stats
for(panel in c(levels(as.factor(panel_data$Panel)),"all")){
    # read in all single variant results
    single_variant_results <- data.frame()
    step <- 250
    if(panel=="all"){
        proteins <- panel_data %>% filter(Panel %in% c("unannotated","covered"))
    }else{
        proteins <- panel_data %>% filter(Panel==panel)
    }
    for(i in seq(1,nrow(proteins),step)){
        files <- c()
        lower <- i
        upper <- min(nrow(proteins),i+step-1)
        for(protein in proteins$Assay[lower:upper]){
            files <- c(files,Sys.glob(paste0(res_dir,"results/",protein,"_rint_burden_wgs/burden/",protein,"_*_masks_*_chr_*_burden_",protein,".regenie")))
        }
        single_variant_results <- vroom(files,altrep=F,skip=1) %>% 
            filter(!grepl("coding",ID) & !is.na(ALLELE0)) %>% 
            left_join(burden_positions_df,by=c("ID"="MASK","CHROM"="CHROM")) %>% 
            rbind(single_variant_results,.) %>% 
            filter(!is.na(MIN_POS)) %>% 
            as_tibble
        closeAllConnections()
    }
    # calculate nvar for the panel
    nvar <- signals_nc %>% 
        filter(protein_name %in% proteins$Assay) %>% 
        nrow()
    # subset all single variants to nvar random samples and run nsim times
    simulations_summary <- data.frame()
    for(sim in seq(sim_start,length.out=nsim)){
        # select nvar variants at random
        single_variant_results_simulated <- single_variant_results %>% 
            sample_n(nvar)
        single_variant_results_ranges <- single_variant_results_simulated %$% 
            GRanges(paste0("chr",CHROM),IRanges(start=MIN_POS,end=MAX_POS))
        # annotate them as which regions they're in for each tissue type and region type
        for(cell in unique(active_regions$CELL)){
            active_regions_ranges <- active_regions %>% 
                filter(CELL==cell) %$% 
                GRanges(paste0("chr",CHROM),IRanges(start=START,end=END))
            single_variant_results_simulated %<>% 
                mutate("{cell}":=ifelse(row_number() %in% queryHits(findOverlaps(single_variant_results_ranges,active_regions_ranges)),1,0))
            for(region in unique(active_regions$TYPE[which(active_regions$CELL==cell)])){
                active_regions_ranges <- active_regions %>% 
                    filter(CELL==cell & TYPE==region) %$% 
                    GRanges(paste0("chr",CHROM),IRanges(start=START,end=END))
                single_variant_results_simulated %<>% 
                    mutate("{cell}.{region}":=ifelse(row_number() %in% queryHits(findOverlaps(single_variant_results_ranges,active_regions_ranges)),1,0))
            }
        }
        # summarise
        simulations_summary <- single_variant_results_simulated %>% 
            summarise(sim_num=sim,across(names(single_variant_results_simulated[(ncol(single_variant_results)+1):ncol(single_variant_results_simulated)]),~sum(.x))) %>% 
            rbind(simulations_summary,.)
    }
    write.table(simulations_summary,paste0("burden_enrichment_simulations_list_",panel,"_",sim_start),quote=F,row.names=F)
}