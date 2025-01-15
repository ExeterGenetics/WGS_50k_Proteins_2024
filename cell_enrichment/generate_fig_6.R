library(tidyverse)
library(magrittr)
library(vroom)
library(IRanges)
library(GenomicRanges)
library(cowplot)
library(Cairo)
setEPS()

nsim <- 10000
freq_thresh <- 0.01
set.seed(20230912)
# file locations
res_dir <- "/slade/projects/UKBB/DNA_Nexus_500k_WGS/Proteomics_50k_3k/"

figure_panels <- list()
# read in enrichment stats
enrichment <- vroom("/slade/projects/UKBB/DNA_Nexus_500k_WGS/Proteomics_50k_3k/results_summaries/rob_working_area/enrichment/leigh_groups_covered/enrichment_p_values_all_joined.txt") %>% 
    mutate(analysis="single_vars") %>% select(-c(max_var_All,enrichment_p_All,n_variants_All)) 
enrichment <- vroom("/slade/projects/UKBB/DNA_Nexus_500k_WGS/Proteomics_50k_3k/results_summaries/rob_working_area/enrichment/leigh_groups_covered/enrichment_p_values_all_joined_burden.txt") %>% 
    mutate(analysis="burden") %>% select(-c(max_var_All,enrichment_p_All,n_variants_All)) %>% bind_rows(enrichment,.)
enrichment <- vroom("/slade/projects/UKBB/DNA_Nexus_500k_WGS/Proteomics_50k_3k/results_summaries/rob_working_area/enrichment/leigh_groups_covered/enrichment_p_values_all_joined_burden_sliding.txt") %>% 
    #mutate(analysis="sliding") %>% select(-c(max_var_All,enrichment_p_All,n_variants_All)) %>% rbind(enrichment,.)
    mutate(analysis="sliding") %>% bind_rows(enrichment,.)
enrichment <- vroom("/slade/projects/UKBB/DNA_Nexus_500k_WGS/Proteomics_50k_3k/results_summaries/rob_working_area/enrichment/leigh_groups_covered/enrichment_p_values_all_joined_burden_ensr.txt") %>% 
    mutate(analysis="ensr") %>% select(-c(max_var_All,enrichment_p_All,n_variants_All)) %>% bind_rows(enrichment,.)
# melt columns
enrichment %<>% select(-starts_with(c("n_variants","max_var"))) %>% 
    pivot_longer(cols=starts_with(c("enrichment_p")),values_to="enrichment_p",names_prefix="enrichment_p_",names_to="panel") %>% 
    filter(!(panel %in% c("All","secreted","signal","secreted_no_sig","signal_no_sec","membrane","unannotated"))) %>% 
    mutate(p_for_qq=ifelse(enrichment_p==0,1/20000,enrichment_p)) %>% 
    mutate(shape_col=as.factor(ifelse(cell_type %in% c("liver","vessel"),"liver/blood vessel","other"))) %>% 
    mutate(panel=ifelse(panel=="non_secreted","non secreted",ifelse(panel=="secreted_signal","secreted or signal",panel)))
# create qq
j<-1
#enrichment %<>% filter(panel !="secreted")
for(region_type in levels(as.factor(enrichment$region))){
    # for(analysis_type in levels(as.factor(enrichment$analysis))){
    for(analysis_type in c("single_vars","burden","sliding","ensr")){
        data_to_plot <- data.frame()
        index_all <- c()
        index_upper <- c()
        index_lower <- c()
        for(i in c(levels(as.factor(enrichment$panel)))){
            data_tmp <- enrichment %>% filter(panel==i & region==region_type & analysis==analysis_type)
            if(nrow(data_tmp)>0){
                data_tmp %<>% arrange(p_for_qq)
                index <- -log10(seq(1,nrow(data_tmp))/nrow(data_tmp))
                data_to_plot <- rbind(data_to_plot,cbind(data_tmp,index))
                if(length(index_all)<nrow(data_tmp)){
                    index_all <- -log10(seq(1,nrow(data_tmp))/nrow(data_tmp))
                    index_upper <- -log10(qbeta(0.975,seq(1,nrow(data_tmp)),(nrow(data_tmp)-seq(1,nrow(data_tmp))+1)))
                    index_lower <- -log10(qbeta(0.025,seq(1,nrow(data_tmp)),(nrow(data_tmp)-seq(1,nrow(data_tmp))+1)))
                }
            }
        }
        plot <- ggplot() + 
            geom_line(data=as_tibble(cbind(index_all,index_upper,index_lower)),aes(index_all,index_all),colour="black") +
            geom_line(data=as_tibble(cbind(index_all,index_upper,index_lower)),aes(index_all,index_upper),colour="black") +
            geom_line(data=as_tibble(cbind(index_all,index_upper,index_lower)),aes(index_all,index_lower),colour="black") +
            geom_point(data=data_to_plot,aes(x=index,y=-log10(p_for_qq),colour=panel,shape=shape_col)) + 
            xlab("Expected -log10p") + 
            ylab("Observed -log10p") + 
            theme_cowplot() + 
            theme(legend.position="top",legend.title=element_blank(),legend.justification="center")
        legend <- get_legend(plot)
        figure_panels[[j]]=plot + theme(legend.position="None")
        j<-j+1
    }
}

ncol <- 4
generate_command <- ""
for(i in seq(1,length(figure_panels))){
    generate_command %<>% paste0("figure_panels[[",i,"]],")
}
postscript("/slade/projects/UKBB/DNA_Nexus_500k_WGS/Proteomics_50k_3k/results_summaries/rob_working_area/enrichment/leigh_groups_covered/Supp_Fig6.eps",width=28,height=30,paper="special")
plot_main <- eval(parse(text=paste0("plot_grid(",generate_command,"ncol=",ncol,",labels=c(\"a)\",\"b)\",\"c)\",\"d)\",\"e)\",\"f)\",\"g)\",\"h)\",\"i)\",\"j)\",\"k)\",\"l)\",\"m)\",\"n)\",\"o)\",\"p)\",\"q)\",\"r)\",\"s)\",\"t)\",\"u)\",\"v)\",\"w)\",\"x)\"))")))
print(plot_grid(plot_main,legend,ncol=1,rel_heights=c(20,1)))
dev.off()
#CairoPNG("/slade/projects/UKBB/DNA_Nexus_500k_WGS/Proteomics_50k_3k/results_summaries/rob_working_area/enrichment/leigh_groups_covered/Supp_Fig6.png",width=1200,height=1400)
CairoTIFF("/slade/projects/UKBB/DNA_Nexus_500k_WGS/Proteomics_50k_3k/results_summaries/rob_working_area/enrichment/leigh_groups_covered/Supp_Fig6.tiff",width=1200,height=1400)
print(plot_grid(plot_main,legend,ncol=1,rel_heights=c(20,1)))
dev.off()
ncol <- 3
generate_command <- ""
#for(i in c(1,4,3,17,20,19,9,12,11)){
for(i in c(1,2,3,5,6,7,9,10,11)){
#for(i in c(1,4,3,5,8,7,9,12,11)){
    generate_command %<>% paste0("figure_panels[[",i,"]],")
}
postscript("/slade/projects/UKBB/DNA_Nexus_500k_WGS/Proteomics_50k_3k/results_summaries/rob_working_area/enrichment/leigh_groups_covered/Fig6.eps",width=28,height=30,paper="special")
plot_main <- eval(parse(text=paste0("plot_grid(",generate_command,"ncol=",ncol,",labels=c(\"a)\",\"b)\",\"c)\",\"d)\",\"e)\",\"f)\",\"g)\",\"h)\",\"i)\",\"j)\",\"k)\",\"l)\",\"m)\",\"n)\",\"o)\",\"p)\",\"q)\",\"r)\",\"s)\",\"t)\",\"u)\",\"v)\",\"w)\",\"x)\"))")))
print(plot_grid(plot_main,legend,ncol=1,rel_heights=c(20,1)))
dev.off()
#CairoPNG("/slade/projects/UKBB/DNA_Nexus_500k_WGS/Proteomics_50k_3k/results_summaries/rob_working_area/enrichment/leigh_groups_covered/Fig6.png",width=640,height=800)
CairoTIFF("/slade/projects/UKBB/DNA_Nexus_500k_WGS/Proteomics_50k_3k/results_summaries/rob_working_area/enrichment/leigh_groups_covered/Fig6.tiff",width=640,height=800)
print(plot_grid(plot_main,legend,ncol=1,rel_heights=c(20,1)))
dev.off()
