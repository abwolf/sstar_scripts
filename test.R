suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(tidyr))
suppressMessages(library(scales))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(stringr))
suppressMessages(library(optparse))

suppressMessages(require(bit64))
##########################

###########
###########
print('FUNCTION: Generate ECDFs')
generate.ecdf.region_ind.fn <- function(null.dt){
 print(' GENERATE ECDFS')
 for( i in sort(unique(as.numeric(null.dt$n_region_ind_snps)))){
   if(i>0){
   print(i)
   nam <<- paste0('null.f.region_ind.', i, '.ecdf')
   if(nrow(filter(null.dt, s_star>0, n_region_ind_snps==i))>0){
     print(nam)
     assign(nam, ecdf(filter(null.dt, s_star>0, n_region_ind_snps==i)$s_star), inherits = TRUE)
     max_snps_ecdf <<- i
   }}
 }
}
#############
#############
print('FUNCTION: Calculate S*-pvalue from ecdf')
estimate.pval.ecdf.region_ind.fn <- function(X, max_snps){
  s_star <- as.numeric(X[["s_star"]])
  n_snps <- as.numeric(X[["n_region_ind_snps"]])
  if (n_snps==0){
    X[["sstarpval_region_ind_snps"]] <- NA
    } else if (n_snps<=max_snps) {
    if( exists(paste0("null.f.region_ind.",n_snps,".ecdf")) ){
      ecdf.fn <- match.fun(paste0("null.f.region_ind.",n_snps,".ecdf"))
      s_star_pval <- 1-ecdf.fn(s_star)
      X[["sstarpval_region_ind_snps"]] <- round(x = s_star_pval, digits = 4)
      }
    } else if (n_snps>max_snps) {
    ecdf.fn <- match.fun(paste0("null.f.region_ind.",max_snps,".ecdf"))
    s_star_pval <- 1-ecdf.fn(s_star)
    X[["sstarpval_region_ind_snps"]] <- round(x = s_star_pval, digits = 4)
    }
  return(X[c("chrom","winstart","winend","ind_id","pop","s_star","n_region_ind_snps","sstarpval_region_ind_snps",
            "num_s_star_snps","hap_1_s_start","hap_1_s_end","hap_2_s_start","hap_2_s_end",
            "n_s_star_snps_hap1","n_s_star_snps_hap2","s_star_haps")])
}
#############
#############
print('FUNCTION: Write OUTPUT TABLES, FILTERED FOR S* AND MATCH PVALUES')
write.filtered.bed.fn <- function(dt, outputdir, mdl, admix, chrom, spval, matchpval){
 print(' Writing filtered .bed file')
 dat_1 <- dt %>%
     filter(s_star>0) %>%
     filter(sstarpval_region_ind_snps<=spval) %>%
     filter(match_pvalue<=matchpval) %>%
     filter(haplotype==0) %>%
     select(msp_ID, winstart, winend) %>%
     #select(msp_ID, hap_1_s_start, hap_1_s_end) %>%
     setnames(c('msp_ID','start','end')) %>%
     as.data.table()

 dat_2 <- dt %>%
     filter(s_star>0) %>%
     filter(sstarpval_region_ind_snps<=spval) %>%
     filter(match_pvalue<=matchpval) %>%
     filter(haplotype==1) %>%
     select(msp_ID, winstart, winend) %>%
     #select(msp_ID, hap_2_s_start, hap_2_s_end) %>%
     setnames(c('msp_ID','start','end')) %>%
     as.data.table()

 dat.bed <- rbind(dat_1,dat_2)

 options(scipen=10)
 dat.bed %<>% mutate(start=str_trim(as.character(start),side = "both")) %>% mutate(end=str_trim(as.character(end), side = "both"))
 write.table(x = dat.bed,
           file = paste0(outputdir,'/',mdl,'_',chrom,'_',admix,'.sstar_sig_',spval,'.match_sig_N_',matchpval,'.isc_0','.bed'),
           quote = FALSE,
           sep = '\t',
           row.names = FALSE,
           col.names = TRUE)
 options(scipen=0)
}
#############
#############
print('FUNCTION: Write OUTPUT TABLES, ALL S* AND MATCH PVALUES')
write.all.bed.fn <- function(dt, outputdir, mdl, admix, chrom){
 print(' Writing unfiltered .bed file')
 dat.bed <- rbind(dt %>%
     			filter(s_star>0) %>%
     			filter(haplotype==0) %>%
     			#select(msp_ID, winstart, winend, sstarpval_region_ind_snps) %>%
     			select(msp_ID, winstart, winend, sstarpval_region_ind_snps, match_pvalue) %>%
     			#select(msp_ID, hap_1_s_start, hap_1_s_end, sstarpval_region_ind_snps, match_pvalue) %>%
     			setnames(c('msp_ID','start','end', 'sstarpval_region_ind_snps', 'match_pvalue')) %>%
                #setnames(c('msp_ID','start','end', 'sstarpval_region_ind_snps')) %>%
     			as.data.table(),
		  dt %>%
     			filter(s_star>0) %>%
     			filter(haplotype==1) %>%
     			#select(msp_ID, winstart, winend, sstarpval_region_ind_snps) %>%
     			select(msp_ID, winstart, winend, sstarpval_region_ind_snps,match_pvalue) %>%
     			#select(msp_ID, hap_2_s_start, hap_2_s_end, sstarpval_region_ind_snps,match_pvalue) %>%
     			setnames(c('msp_ID','start','end', 'sstarpval_region_ind_snps', 'match_pvalue')) %>%
                #setnames(c('msp_ID','start','end', 'sstarpval_region_ind_snps')) %>%
     			as.data.table()
		)

 options(scipen=10)
 dat.bed %<>% mutate(start=str_trim(as.character(start),side = "both")) %>% mutate(end=str_trim(as.character(end), side = "both"))
 write.table(x = dat.bed,
           file = paste0(outputdir,'/',mdl,'_',chrom,'_',admix,'.sstar_sig_','ALL','.match_sig_N_MH_','ALL','.isc_0','.bed'),
           quote = FALSE,
           sep = '\t',
           row.names = FALSE,
           col.names = TRUE)
 options(scipen=0)
}


#####################################
#####################################

# Command line arguments
option_list = list(
    make_option(c("--inputdir"), action="store", default=NA, type='character', help="Directory containing input files"),
    make_option(c("--mdl"), action="store", default=NA, type='character', help="Model type, e.g. Tenn_nonAfr"),
    make_option(c("--null_dir"), action="store", default='null', type='character', help="Null directory name"),
    make_option(c("--null_tag"), action="store", default='n1_0.0_n2_0.0', type='character', help="Null model tag, e.g. n1_0.0_n2_0.0"),
    make_option(c("--ecdf"), action="store", default=NA, type='character', help='Specify stored ECDF RData set to use'),
    make_option(c("--admix_dir"), action="store", default=NA, type='character', help="Admix directory name"),
    make_option(c("--admix_tag"), action="store", default=NA, type='character', help="Admix model tag, e.g. n1_0.02_n2_0.0"),
    make_option(c("--max_chrm_admix"), action="store", default=NA, type='numeric', help="Number of chromosomes to test from admix data"),
    make_option(c("--max_chrm_null"), action="store", default=NA, type='numeric', help="Number of chromosomes to test from null data"),
    make_option(c("--sstarpval"), action="store", default=0.01, type='numeric', help="Sstar pvalue cutoff for significance"),
    make_option(c("--matchpval"), action="store", default=0.05, type='numeric', help="Match pvalue cutoff for significance"),
    make_option(c("--outputdir"), action="store", default=NA, type='character', help="Set output directory for bedfiles"),
    make_option(c("--long"), action="store_true", default=TRUE, help="Print complete output, w/o filtering [default]"),
    make_option(c("--short"), action="store_false", dest="long", help="Print the filtered output")
)

opt = parse_args(OptionParser(option_list=option_list))

opt

if(opt$long==FALSE){
    print('OKAY')
} else {
    print('NO')
}
