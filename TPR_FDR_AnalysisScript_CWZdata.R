###########################
# CWZ Data

outputdir <- '~/DATALab/SimulatedDemographic/Sstar/CZimmer/SstarSigFiles/'
mdl <- 'zimmer2'
admix <- 'CWZ'
dt <- data.table(NULL)

infile <- paste0(outputdir,'/',mdl,'_',admix,".sstar_sig_out")
print(infile)

dat <- read.table(infile, header=TRUE, na.strings=c("NA", "None",'.'), as.is = TRUE)
dat <- as.data.table(dat)
dt <- rbind(dt, dat)

CWZ.sstargsig.dt <- dt


dt <- CWZ.sstargsig.dt %>% filter(filter==FALSE) %>% filter(sig0.99==TRUE) %>% filter(as.numeric(hap_1_window_pval)<0.1 | as.numeric(hap_2_window_pval)<0.1) %>% select(chrom, winstart, winend) %>% as.data.table()

dat.bed <- dt

dat.bed <- out %>% filter(filter==FALSE) %>% filter(sstarpval_region_ind_snps<=0.01) %>% filter(0.8 <= as.numeric(callable_bases)/(winend-winstart)) %>% 
  filter(as.numeric(hap_1_window_pval)<0.1 | as.numeric(hap_2_window_pval)<0.1) %>% 
  select(msp_ID, winstart, winend, sstarpval_region_ind_snps, match_pvalue) %>% 
  setnames(c('msp_ID','start','end', 'sstarpval_region_ind_snps', 'match_pvalue')) %>% 
  as.data.table()

options(scipen=10)
write.table(x = dat.bed,
             file = paste0(outputdir,'zimmer2_CWZ.sstar_sig_','0.01','.matchpval_0.1.windows'),
             #file = paste0(outputdir,'test.rsf_',req.snp.frac,'.sstar_sig_','ALL','.windows'),
             #file = paste0(outputdir,mdl,'_','ALL','_',admix,'_msp_ALL','.sstar_sig_out_region_ind_snps.sstar_','ALL','.bed'),
             quote = FALSE,
             sep = '\t',
             row.names = FALSE,
             col.names = FALSE)
options(scipen=0)

########################################
########################################
# NA06984 Data

outputdir <- '~/DATALab/SimulatedDemographic/Sstar/CZimmer/SstarSigFiles/'
mdl <- 'zimmer2'
admix <- 'NA06984'
dt <- data.table(NULL)

#infile <- paste0(outputdir,'/',mdl,'_',admix,".windowcalc_out.gz")
#print(infile)

infile <- paste0(outputdir,'/',mdl,'_',admix,".sstar_sig_out.gz")
print(infile)

dat <- read.table(gzfile(infile), header=TRUE, na.strings=c("NA", "None",'.'), as.is = TRUE)
dat <- as.data.table(dat)
dt <- rbind(dt, dat)

NA06984.sstarsig.dt <- dt

dt <- NA06984.sstarsig.dt %>% filter(filter==FALSE) %>% filter(sig0.99==TRUE) %>% 
      filter(as.numeric(hap_1_window_pval)<0.1 | as.numeric(hap_2_window_pval)<0.1) %>% 
      filter(0.8 <= as.numeric(callable_bases)/(winend-winstart)) %>%
      mutate(chrom=paste0('NA06984_',chrom,'_NA')) %>%
      select(chrom, winstart, winend) %>% as.data.table()

dat.bed <- dt


## out file comes from ecdf estimation function...

dat.bed <- out %>% filter(filter==FALSE) %>% filter(sstarpval_region_ind_snps<=0.01) %>% filter(0.8 <= as.numeric(callable_bases)/(winend-winstart)) %>% 
  filter(as.numeric(hap_1_window_pval)<0.1 | as.numeric(hap_2_window_pval)<0.1) %>% 
  select(msp_ID, winstart, winend, sstarpval_region_ind_snps) %>% 
  setnames(c('msp_ID','start','end', 'sstarpval_region_ind_snps')) %>% 
  as.data.table()

options(scipen=10)
write.table(x = dat.bed,
            file = paste0(outputdir,'test.NA06984.mu_1.2.sstar_sig_','0.01','.matchpval_0.1.windows'),
            #file = paste0(outputdir,'test.rsf_',req.snp.frac,'.sstar_sig_','ALL','.windows'),
            #file = paste0(outputdir,mdl,'_','ALL','_',admix,'_msp_ALL','.sstar_sig_out_region_ind_snps.sstar_','ALL','.bed'),
            quote = FALSE,
            sep = '\t',
            row.names = FALSE,
            col.names = FALSE)
options(scipen=0)

########################################
########################################
# NA07000 Data

outputdir <- '~/DATALab/SimulatedDemographic/Sstar/CZimmer/SstarSigFiles//'
mdl <- 'zimmer2'
admix <- 'NA07000'
dt <- data.table(NULL)

#infile <- paste0(outputdir,'/',mdl,'_',admix,".windowcalc_out.gz")
#print(infile)

infile <- paste0(outputdir,'/',mdl,'_',admix,".sstar_sig_out.gz")
print(infile)

dat <- read.table(gzfile(infile), header=TRUE, na.strings=c("NA", "None",'.'), as.is = TRUE)
dat <- as.data.table(dat)
dt <- rbind(dt, dat)

#NA07000.dt <- dt

NA07000.sstarsig.dt <- dt

dt <- NA07000.sstarsig.dt %>% filter(filter==FALSE) %>% filter(sig0.99==TRUE) %>% 
  #filter(as.numeric(hap_1_window_pval)<0.1 | as.numeric(hap_2_window_pval)<0.1) %>% 
  filter(0.8 <= as.numeric(callable_bases)/(winend-winstart)) %>%
  mutate(chrom=paste0(admix,'_',chrom,'_NA')) %>%
  select(chrom, winstart, winend) %>% as.data.table()

dat.bed <- dt


## out comes from ecdf estimation function...

dat.bed <- out %>% filter(filter==FALSE) %>% filter(sstarpval_region_ind_snps<=0.01) %>% filter(0.8 <= as.numeric(callable_bases)/(winend-winstart)) %>% 
  #filter(as.numeric(hap_1_window_pval)<0.1 | as.numeric(hap_2_window_pval)<0.1) %>% 
  select(msp_ID, winstart, winend) %>% 
  setnames(c('msp_ID','start','end')) %>% 
  as.data.table()

options(scipen=10)
write.table(x = dat.bed,
            file = paste0(outputdir,'test.NA07000.ref_20.sstar_sig_','0.01','.windows'),
            #file = paste0(outputdir,'test.rsf_',req.snp.frac,'.sstar_sig_','ALL','.windows'),
            #file = paste0(outputdir,mdl,'_','ALL','_',admix,'_msp_ALL','.sstar_sig_out_region_ind_snps.sstar_','ALL','.bed'),
            quote = FALSE,
            sep = '\t',
            row.names = FALSE,
            col.names = FALSE)
options(scipen=0)


########################################
########################################
# NA06986 Data

outputdir <- '~/DATALab/SimulatedDemographic/Sstar/CZimmer/IndFiles/'
mdl <- 'zimmer2'
admix <- 'NA06986'
dt <- data.table(NULL)

infile <- paste0(outputdir,'/',mdl,'_',admix,".windowcalc_out.gz")
print(infile)

#infile <- paste0(outputdir,'/',mdl,'_',admix,".sstar_sig_out.gz")
#print(infile)

dat <- read.table(gzfile(infile), header=TRUE, na.strings=c("NA", "None",'.'), as.is = TRUE)
dat <- as.data.table(dat)
dt <- rbind(dt, dat)

NA06986.dt <- dt

#NA06986.sstarsig.dt <- dt

dt <- NA06986.sstarsig.dt %>% filter(filter==FALSE) %>% filter(sig0.99==TRUE) %>% 
  filter(as.numeric(hap_1_window_pval)<0.1 | as.numeric(hap_2_window_pval)<0.1) %>% 
  filter(0.8 <= as.numeric(callable_bases)/(winend-winstart)) %>%
  mutate(chrom=paste0(admix,'_',chrom,'_NA')) %>%
  select(chrom, winstart, winend) %>% as.data.table()

dat.bed <- dt
################
################
NA06986.sstarsig.dt <- left_join(NA06986.sstarsig.dt, select(NA06986.dt, chrom, winstart, winend, sstarpval_region_ind_snps)) %>% as.data.table()
NA06986.sstarsig.dt[,msp_ID:=paste0(ind_id,'_',chrom,'_NA')]

write.table(x = select(NA06986.sstarsig.dt, msp_ID, winstart,winend,filter,s_star,sstarpval_region_ind_snps,n_region_ind_snps,recomb, hap_1_window_pval, hap_2_window_pval), file = '~/DATALab/SimulatedDemographic/Sstar/CZimmer/SstarSigFiles/zimmer2_NA06986.sstar_sig_ALL.windows',quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

NA06986.sstarsig.LLcallset.dt <- fread('~/DATALab/SimulatedDemographic/Sstar/CZimmer/SstarSigFiles/zimmer2_NA06986.sstar_sig_ALL.windows.intersect_LL.callset', col.names = c('msp_ID', 'winstart','winend','filter','s_star','sstarpval_region_ind_snps','n_region_ind_snps','recomb', 'hap_1_window_pval', 'hap_2_window_pval', 'overlap_bases'))

NA06986.sstarsig.LLcallset.dt <- as.data.table(t(apply(X = NA06986.sstarsig.LLcallset.dt, MARGIN = 1, FUN = estimate.pval.ecdf.region_ind.fn, max_snps=235)))

NA06986.sstarsig.LLcallset.dt  = mutate_each(NA06986.sstarsig.LLcallset.dt, funs(as.numeric), winstart, winend, s_star,sstarpval_region_ind_snps,n_region_ind_snps,recomb,hap_1_window_pval,hap_2_window_pval, overlap_bases)

###############

NA06986.sstarsig.dt[,msp_ID:=paste0(ind_id,'_',chrom,'_NA')]
out <- NA06986.sstarsig.LLcallset.dt

out <- left_join(out, select(NA06986.sstarsig.dt, msp_ID, winstart, winend, num_s_star_snps, n_s_star_snps_hap1, n_s_star_snps_hap2, hap_1_s_start, hap_1_s_end,  hap_2_s_start, hap_2_s_end, s_start, s_end, min_callable), by=c('msp_ID', 'winstart', 'winend')) %>% as.data.table()



out.bed <- out %>% mutate(s_start=str_trim(as.character(s_start),side = "both")) %>% mutate(s_end=str_trim(as.character(s_end), side = "both")) %>%
  filter(filter==FALSE) %>%
  filter(recomb <= 1.1 & recomb >= 0.9) %>%
  filter(min_callable>=40000) %>% 
  filter(sstarpval_region_ind_snps <= 0.01) %>%
  filter(hap_1_window_pval <= 0.01 | hap_2_window_pval <= 0.01) %>%
  filter(!(s_star_hap_1 == TRUE & s_star_hap_2 == TRUE & haplotype==2))

out.bed[,msp_ID:=paste0(msp_ID,'_',haplotype)]

write.table(x = select(out.bed, msp_ID, s_start, s_end, haplotype), file = '~/DATALab/SimulatedDemographic/Sstar/CZimmer/test.NA06986.rsf_0.8.mincall_0.8.sstar_sig_0.01.matchpval_0.01.bed', quote = FALSE, sep = '\t',row.names =  FALSE, col.names = FALSE)

## out comes from ecdf estimation function...

dat.bed <- out %>% filter(filter==FALSE) %>% filter(sstarpval_region_ind_snps<=0.01) %>% filter(0.8 <= as.numeric(callable_bases)/(winend-winstart)) %>% 
  #filter(as.numeric(hap_1_window_pval)<0.1 | as.numeric(hap_2_window_pval)<0.1) %>% 
  select(msp_ID, winstart, winend) %>% 
  select(msp_ID, )
  setnames(c('msp_ID','start','end')) %>% 
  as.data.table()

options(scipen=10)
dat.bed %<>% mutate(start=str_trim(as.character(start),side = "both")) %>% mutate(end=str_trim(as.character(end), side = "both"))
write.table(x = dat.bed,
            file = paste0(outputdir,'test.NA07000.sstar_sig_','0.01','.bed'),
            #file = paste0(outputdir,'test.rsf_',req.snp.frac,'.sstar_sig_','ALL','.windows'),
            #file = paste0(outputdir,mdl,'_','ALL','_',admix,'_msp_ALL','.sstar_sig_out_region_ind_snps.sstar_','ALL','.bed'),
            quote = FALSE,
            sep = '\t',
            row.names = FALSE,
            col.names = FALSE)
options(scipen=0)

################

NA07000.sstarsig.dt <- fread('~/DATALab/SimulatedDemographic/Sstar/CZimmer/test.NA07000.sstar_sig_ALL.matchpval_ALL.windows.cM_filtered', header = TRUE)
NA07000.sstarsig.LLcallset.dt <- fread('~/DATALab/SimulatedDemographic/Sstar/CZimmer/SstarSigFiles/zimmer2_NA07000.sstar_sig_ALL.windows.intersect_LL.callset', 
                                       col.names = c('msp_ID','winstart','winend','filter','s_star','sstarpval_region_ind_snps','n_region_ind_snps','recomb','overlap'))

NA07000.sstarsig.LLcallset.dt <- NA07000.sstarsig.LLcallset.dt %>% left_join(select(NA07000.sstarsig.dt, msp_ID, winstart, winend, hap_1_window_pval, hap_2_window_pval), by = c('msp_ID', 'winstart', 'winend')) %>% as.data.table() %>% mutate_each(funs(as.numeric), winstart, winend, s_star,sstarpval_region_ind_snps,n_region_ind_snps,recomb, hap_1_window_pval, hap_2_window_pval, overlap)


ggplot() + theme_grey() +
  geom_histogram(data=NA07000.sstarsig.LLcallset.dt[filter==FALSE & s_star!=0 & recomb<=1.1 & recomb >=0.9], aes(x=s_star, fill='NA_07000.1cM_Mb'), alpha=0.5) + 
  geom_histogram(data=NA07000.sstarsig.LLcallset.dt[filter==FALSE & sstarpval_region_ind_snps<=0.01 & recomb<=1.1 & recomb >=0.9], aes(x=s_star, fill='NA_07000.1cM_Mb.S*pval'), alpha=0.5) +
  geom_histogram(data=NA07000.sstarsig.LLcallset.dt[filter==FALSE & sstarpval_region_ind_snps<=0.01 & (hap_1_window_pval<=0.05 | hap_2_window_pval<=0.05) & recomb<=1.1 & recomb >=0.9], aes(x=s_star, fill='NA_07000.1cM_Mb.S*pval.matchpval'), alpha=0.5) +
  #geom_histogram(data=NA07000.sstarsig.LLcallset.dt[filter==FALSE & s_star!=0 & recomb<=1.1 & recomb >=0.9 & overlap>0], aes(x=s_star, fill='NA_07000.1cM_Mb.LLcallset_intersect'), alpha=0.5)
  geom_histogram(data=NA07000.sstarsig.LLcallset.dt[filter==FALSE & recomb<=1.1 & recomb >=0.9 & overlap>0], aes(x=s_star, fill='NA_07000.1cM_Mb.LLcallset_intersect'), alpha=0.5)

##################
##################
dt <- as.data.table(NULL)
sstarsig.LLcallset.dt <- NA07000.sstarsig.LLcallset.dt
exists("sstarsig.LLcallset.dt")
for(sp in list(0.1, 0.05, 0.01)){
  #print(sp)
  for(mp in list(1, 0.1, 0.05, 0.01)){
    #print(mp)
    for(o in list(0.5, 0.2, 0.1, 0)){
      #print(o)
      TP <- sstarsig.LLcallset.dt %>% filter(filter==FALSE) %>% 
        filter(recomb<=1.1 & recomb>=0.9) %>% 
        filter(sstarpval_region_ind_snps<=sp) %>%
        filter(hap_1_window_pval<=mp | hap_2_window_pval<=mp) %>%
        filter(overlap_bases>o) %>% 
        nrow()
      FP <- sstarsig.LLcallset.dt %>% filter(filter==FALSE) %>% 
        filter(recomb<=1.1 & recomb>=0.9) %>% 
        filter(sstarpval_region_ind_snps<=sp) %>% 
        filter(hap_1_window_pval<=mp | hap_2_window_pval<=mp) %>% 
        filter(overlap_bases<=o) %>% 
        nrow()
      
      FDR = FP/(TP+FP)
      
      TN <- sstarsig.LLcallset.dt %>% filter(filter==FALSE) %>% 
        filter(recomb<=1.1 & recomb>=0.9) %>% 
        filter(sstarpval_region_ind_snps>sp) %>% 
        #filter(hap_1_window_pval>mp | hap_2_window_pval>mp) %>% 
        filter(overlap_bases<=o) %>% 
        nrow()
      FN <- sstarsig.LLcallset.dt %>% filter(filter==FALSE) %>% 
        filter(recomb<=1.1 & recomb>=0.9) %>% 
        filter(sstarpval_region_ind_snps>sp) %>% 
        #filter(hap_1_window_pval>mp | hap_2_window_pval>mp) %>% 
        filter(overlap_bases>o) %>% 
        nrow()
      
      TPR = TP/(TP+FN)
      
      d <- as.data.table(cbind(sp, mp, o, TP, FP, TN, FN, FDR, TPR))
      dt <- as.data.table(rbind(dt, d))
    }
  }
}

