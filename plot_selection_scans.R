require(ggplot2)
require(plyr)
require(reshape2)
require(VennDiagram)
require(tidyr)
require(data.table)
## functions

## Sven's scaler to scale values to a specified range
scaler <- function (x, r = c(0, 1), b = range(x, na.rm = TRUE)) {
  rl <- r[1]
  ru <- r[2]
  bl <- b[1]
  bu <- b[2]
  if (bl > min(x, na.rm = TRUE))
    stop("Lower border in b is greater than minimum of x.")
  if (bu < max(x, na.rm = TRUE))
    stop("Upper border in b is less than maximum of x.")
  (ru - rl) * (x - bl) / (bu - bl) + rl
}

## converting chromosome coordinates to circular for plotting with coord_polar
circular <- function(x){
  circ <- function(z){
    a <- z[1,which(colnames(z) == "chr")]
    b <- subset(circ_chrom, chr %in% a)[1,1]  
    z$pos_circ <- z$pos+b
    z
  }
  ddply(x, c("chr"), circ) 
}

## adjusting chr ends (+5 cM for 10 cM windows); !!! not always necessary, check min pos ####

adj.chrends <- function(y){
  a <- data.frame(
    chr=as.factor(rep(seq(chr.n))),
    pos=c(chr.max)
  )
  y$pos <- y$pos+windsize/2
  ddply(y,c("chr"),function(x){
    b=x[1,which(colnames(y) == "chr")]
    subset(x, pos < a[b,2])
  }
  )
}

## END functions ####

## INPUT parameters

windsize=10

nw_1=0
nw_2=.24
ng_1=.25
ng_2=.49
fw_1=.5
fw_2=0.74
#fst_1=0.61
#fst_2=0.8
sw_1=0.75
sw_2=1

# FayWu simulated confidence intervals p < 0.001

hdmincds <-  -4.23 # -2.99 # 
hdmaxcds <- 1.1
hwmincds <- -3.28
hwmaxcds <- 1.03

swmin <- quantile(dom1Kun$V2, .99)
swminw <- quantile(wild1K$V2, .995)


## END INPUT

## making circular chromsomes for plotting ####

chrom_lim <- read.delim(file="/biodata/dep_coupland/grp_korff/artem/scripts/vcf_to_pca/chr.ends.new.popseq", header = F, sep = "\t")
colnames(chrom_lim) <- c("Target", "Chr", "Pos")
chrom_lim$Chr <- as.factor(chrom_lim$Chr)
chr.n <- nlevels(chrom_lim$Chr)
# maximum chromosome positions
chr.max <- tapply(chrom_lim$Pos, chrom_lim$Chr, max, na.rm=T)

a <- as.vector(chr.max[1])
for (i in seq(2,12,by=2)){
  a[i] <- a[i-1] + 10
  a[i+1] <- a[i]+chr.max[i/2+1]
}
circ_chrom <- data.frame(pos=c(0,a),chr=sort(rep(seq(1:7),2)))

## calculated position of chr labels ####

chr_label <- ddply(circ_chrom, c("chr"), function(x){mean(x$pos)})

## calulating cM scale for drawing ####
ticks <- c()
for (i in 1:7){
  a <- subset(circ_chrom, chr %in% i)
  b <- floor(a[1,1])
  d <- ceiling(a[2,1])
  e <- seq(b,d,10)
  ticks <- c(ticks,e)
}
ticks <- as.data.frame(ticks)

temp_chr <- as.vector(ceiling(chr.max))
cM <- c()
for (i in 1:7){
  a <- temp_chr[i]
  b <- seq(0,a,10)
  cM <- c(cM, b)
}

ticks <- data.frame(label=cM, pos=ticks)
ticks_mod <- ticks[seq(1,nrow(ticks),by=2),]
#####

## END making circular chromosomes for plotting ####

## NeiPi MODULE - windows

# old popseq system("sed 's/wind//g' /netscratch/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/taj_fst/vcf_size10_step1_57samples_domest_homo_allconf_UG.miss0.5/*neipi.bywindow | sed 's/pop1/domest/g' > /biodata/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/taj_fst/neipi_window_domest_57samples_var_nonvar.txt")
# old popseq system("sed 's/wind//g' /netscratch/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/taj_fst/vcf_size10_step1_302samples_wild_homo_allconf_UG.miss0.5/*neipi.bywindow | sed 's/pop1/wild/g' > /biodata/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/taj_fst/neipi_window_wild_302samples_var_nonvar.txt")

system("bash /netscratch/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/taj_fst/vcf_size10_step1_302samples_wild_homo_allconf_UG.miss0.5.popseq17/mstats_to_neipi.sh > /biodata/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/taj_fst/neipi_window_wild_302samples_var_nonvar.popseq17.txt")
system("bash /netscratch/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/taj_fst/vcf_size10_step1_57samples_domest_homo_allconf_UG.miss0.5.popseq17/mstats_to_neipi.sh > /biodata/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/taj_fst/neipi_window_domest_57samples_var_nonvar.popseq17.txt")

## loading data, binding wild & domest

neipiwind101d <- read.delim(file="/biodata/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/taj_fst/neipi_window_domest_57samples_var_nonvar.popseq17.txt", header = F, sep = "\t")
neipiwind101w <- read.delim(file="/biodata/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/taj_fst/neipi_window_wild_302samples_var_nonvar.popseq17.txt", header = F, sep = "\t")
neipiwind101 <- rbind(neipiwind101w, neipiwind101d)
colnames(neipiwind101) <- c("test", "type", "val", "chr", "pos","count")

## filtering windows with > 10 SNPs

neipiwind101 <- subset(neipiwind101, count > 10)
neipiwind101 <- neipiwind101[,seq(1,5)]

## adjusting plotting positions to the middle of a window

neipiwind101  <- adj.chrends(neipiwind101)

## converting to wide format; calculating the ratio

neipiwind101_wide <- dcast(neipiwind101,test+chr+pos~type, value.var="val")
neipiwind101_wide$ratio <- neipiwind101_wide$wild/neipiwind101_wide$domest
neipiwind101_wide <- subset(neipiwind101_wide, ! ratio %in% NA)

## converting to circular coordinates

neipiwind101_wide <- circular(neipiwind101_wide)

## scaling the ratio for plotting

# sort(neipiwind101_wide$ratio, TRUE)[3] ## if only one value is very high, then check for the second highest

neipiwind101_wide$ratio[neipiwind101_wide$ratio > 6] <- 6

maxrat <- max(neipiwind101_wide$ratio)
neipiwind101_wide$ratio_sc <- scaler(neipiwind101_wide$ratio, c(nw_1,nw_2), b=c(0,maxrat))

## constructing scale/axis for NeiPi windows

neiscale <- seq(0,ceiling(maxrat),2)
plotneiscale <- data.frame(neival = neiscale, neiscale_sc = scaler(neiscale, c(nw_1,nw_2), b=c(min(neiscale),max(neiscale))))
plotneiscale$text <- plotneiscale[,1]
plotneiscale$text[length(plotneiscale$neival)] <- paste0(">",max(plotneiscale[,1]))


## nei's pi outliers, windows ## skip for the manuscript!!

#neipiwind101_wide$pval <- 2*pnorm(-abs(scale(neipiwind101_wide$ratio)))

#wind_outliers_nei <- subset(neipiwind101_wide, pval < 0.05)
#wind_outliers_nei <- wind_outliers_nei[order(wind_outliers_nei[,2],wind_outliers_nei[,3]),]
#nei_wind_thresh1 <- wind_outliers_nei[which.min(abs(wind_outliers_nei$pval-0.05)),which(colnames(wind_outliers_nei) %in% "ratio")]

## END nei's pi outliers, windows.

## END NeiPi MODULE - windows

## NeiPi MODULE - genes

# loading data

neipi_gene_d <- read.table("/biodata/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/taj_fst/NeiPi.bygene.57samples_domest_homo_allconf_UG.miss0.5.counts.mapping.newpopseq", sep="\t", header = F)
neipi_gene_w <- read.table("/biodata/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/taj_fst/NeiPi.bygene.302samples_wild_homo_allconf_UG.miss0.5.counts.mapping.newpopseq", sep="\t", header = F)
neipi_gene <- rbind(neipi_gene_d, neipi_gene_w)
colnames(neipi_gene) <- c("test","type","val","locus", "count", "chr_old","pos_old" , "chr","pos")
neipi_gene$chr <- as.integer(as.character(neipi_gene$chr))
neipi_gene$pos <- as.numeric(as.character(neipi_gene$pos))

# converting to wide format using data.table (multiple values)
neipi_gene_wide1 <- as.data.frame(dcast.data.table(setDT(neipi_gene), test+locus+chr+pos ~ type, value.var=c("val","count")))

neipi_gene_wide <- subset(neipi_gene_wide1, count_wild > 4)
neipi_gene_wide <- subset(neipi_gene_wide, count_domest > 4)

# calculating ratio wild / domest

neipi_gene_wide$ratio <- neipi_gene_wide$val_wild/neipi_gene_wide$val_domest

# filtering out NA and 0 diversity in domesticated

neipi_gene_widef_all <- subset(neipi_gene_wide, !is.na(ratio) & ! ratio %in% "Inf")

# extracting pvalues (mapped + unmapped) & outliers

neipi_gene_widef_all$pval <- 2*pnorm(-abs(scale(neipi_gene_widef_all$ratio)))

neipi_gene_widef <- subset(neipi_gene_widef_all, chr < 8 & pos < 777)

# scaling and circulating the genes

neipi_gene_widef$ratio[neipi_gene_widef$ratio > 100] <- 100
neipi_gene_widef$ratio_sc <- scaler(neipi_gene_widef$ratio, c(ng_1,ng_2), b = c(0,max(neipi_gene_widef$ratio)+1))
neipi_gene_widef <- circular(neipi_gene_widef)

# extracting outlier / nonoutliers

outliers_nei_gene <- subset(neipi_gene_widef, pval < 0.05)
nonoutliers_nei_gene <- subset(neipi_gene_widef, pval >= 0.05)

nei_gene_thresh <- outliers_nei_gene[which.max(abs(outliers_nei_gene$pval)),11]

## constructing scale for NeiPi genes ####

min_gen <- floor(min(neipi_gene_widef$ratio))
max_gen <- ceiling(max(neipi_gene_widef$ratio))

## testing even/odd ####
if(max_gen %% 2 == 0){
  step <- max_gen/4
} else {
  max_gen <- max_gen+1
  step <- (max_gen)/4
}
# ####
neiscale_genes <- c(min_gen,min_gen+step,min_gen+2*step,min_gen+3*step,max_gen)

plotneiscale_genes <- data.frame(neival = neiscale_genes, neiscale_sc = scaler(neiscale_genes, c(ng_1,ng_2), b=c(min(neiscale_genes),max(neiscale_genes))))
plotneiscale_genes$neiscale_sc <- round(plotneiscale_genes$neiscale_sc, digits = 2)

plotneiscale_genes$text <- plotneiscale_genes[,1]
plotneiscale_genes$text[length(plotneiscale_genes$neival)] <- paste0(">",max(plotneiscale_genes[,1]))

## END NeiPi MODULE - genes ####

## FayWU MODULE - windows ####

system("bash /netscratch/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/taj_fst/mstats_to_faywu_wild_domest.sh > /biodata/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/taj_fst/faywu_window_size10_step1_wild_dom_var_nonvar.popseq17.txt")

faywuwind101 <- read.delim(file="/biodata/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/taj_fst/faywu_window_size10_step1_wild_dom_var_nonvar.popseq17.txt", header = F, sep = "\t")
colnames(faywuwind101) <- c("test", "type", "val", "chr", "pos")
faywuwind101 <- adj.chrends(faywuwind101)
head(faywuwind101)
tmax <- max(faywuwind101$val)
tmin <- min(faywuwind101$val)

faywuwind101$val_sc <- scaler(faywuwind101$val, c(fw_1,fw_2), b=c(floor(tmin),ceiling(tmax)))
faywuwind101 <- circular(faywuwind101)

## END FayWU MODULE - windows ####

## FayWU MODULE - genes ####

#### loading by gene data, cds and all

faywu_genes_all <- read.delim(file="/netscratch/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/taj_fst/faywu.bygene.outgroup.nohead.478sampl.filtered.SNP.478samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo.with_sing.all.counts.cutoff5incl.newpopseq", header = F, sep = "\t")
colnames(faywu_genes_all) <- c("test","type","val","locus","count", "chr_old", "pos_old", "chr","pos")
faywu_genes_all <- subset(faywu_genes_all, count > 8 & type %in% "domest")
faywu_genes_all$chr <- as.integer(as.character(faywu_genes_all$chr))
faywu_genes_all$pos <- as.numeric(as.character(faywu_genes_all$pos))

faywu_genes_all$pval <- 2*pnorm(-abs(scale(faywu_genes_all$val)))
faywu_genes <- subset(faywu_genes_all, chr < 8 & pos < 777)

fwmin <- min(faywu_genes$val)
fwmax <- max(faywu_genes$val)

faywu_genes$val_sc <- scaler(faywu_genes$val, c(fw_1,fw_2),b=c(floor(fwmin),ceiling(fwmax)))

#faywu_gene_outliers <- subset(faywu_genes, pval <  0.05)
#faywu_gene_rest <- subset(faywu_genes, pval >=  0.05)

#faywu_genes_thresh <- faywu_gene_outliers[which.max(abs(faywu_gene_outliers$pval)),9]


faywu_gene_outliers <- subset(faywu_genes, val <=  hdmincds & type %in% "domest")
faywu_gene_rest <- subset(faywu_genes, (val > hdmincds & type %in% "domest") | (val > hwmincds & type %in% "wild"))

## converting to circular coordinates

faywu_gene_outliers <- circular(faywu_gene_outliers)
faywu_gene_rest <- circular(faywu_gene_rest)

## creating simulated borders

confidH <- data.frame(type = c(rep("domest",2), rep("wild",2)), group = c("hdmincds","hdmaxcds","hwmincds","hwmaxcds"), val = c(hdmincds,hdmaxcds,hwmincds,hwmaxcds))
confidH$val_sc <- scaler(confidH$val, c(fw_1,fw_2), b=c(floor(fwmin),ceiling(fwmax)))

## constructing scale/axis for faywu values

fayscale <- seq(floor(fwmin),ceiling(fwmax),2)
plotfayscale <- data.frame(fayjval = fayscale, fayscale_sc = scaler(fayscale, c(fw_1,fw_2), b=c(min(fayscale),max(fayscale))))

## END FayWU MODULE - genes ####

## Fst MODULE - genes ####

fst_genes <- read.delim(file="/netscratch/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/taj_fst/fst.bygene.478sampl.filtered.SNP.478samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo.with_sing.all.counts.position", header = F, sep = "\t")

fst_genes <- fst_genes[,c(2,5,7,8,9,10)]

colnames(fst_genes) <- c("val","p","locus", "count", "chr","pos")
fst_genes <- subset(fst_genes, p < 0.01)
fst_genes$val[fst_genes$val <0] <- 0
fst_genes <- subset(fst_genes, count >= 3)

fst_genes$val_sc <- scaler(fst_genes$val, c(fst_1,fst_2),b=c(0,1))

# extracting pvalues (mapped + unmapped) & outliers

fst_genes$pval <- 2*pnorm(-abs(scale(fst_genes$val_sc)))

# !!!! to extract outliers for the paper include also unmapped!!

fst_genesf <- subset(fst_genes, chr < 8 & pos < 777)

# circulating

fst_genesf <- circular(fst_genesf)

# extracting outlier / nonoutliers

outliers_fst_genes <- subset(fst_genesf, pval < 0.05)
nonoutliers_fst_genes <- subset(fst_genesf, pval >= 0.05)

fst_genes_thresh <- outliers_fst_genes[which.max(abs(outliers_fst_genes$pval)),6]


## constructing scale for Fst genes ####

fstscale_genes <- c(0,0.25,0.5,0.75,1)

plotfstscale_genes <- data.frame(fstval = fstscale_genes, fstscale_sc = scaler(fstscale_genes, c(fst_1,fst_2), b=c(0,1)))
plotfstscale_genes$text <- plotfstscale_genes[,1]

## END Fst MODULE - genes ####

## SweeD MODULE - genes ####

sweed_genes <- read.delim(file="/netscratch/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/sweed/57samples_domest.sweed.grid5.snp3.newpopseq", header = F, sep = "\t")

sweed_geneswild <- read.delim(file="/netscratch/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/sweed/302samples_wild.sweed.grid2.snp3.newpopseq", header = F, sep = "\t")

colnames(sweed_genes) <- c("locus","loc_pos","val","alpha", "chr_old", "pos_old", "chr","pos")
colnames(sweed_geneswild) <- c("locus","loc_pos","val","alpha", "chr_old", "pos_old", "chr","pos")

sweed_genes$chr <- as.integer(as.character(sweed_genes$chr))
sweed_geneswild$chr <- as.integer(as.character(sweed_geneswild$chr))
sweed_genes$pos <- as.numeric(as.character(sweed_genes$pos))
sweed_geneswild$pos <- as.numeric(as.character(sweed_geneswild$pos))

# extracting pvalues (mapped + unmapped) & outliers

sweed_genesmax <- ddply(sweed_genes, c("locus"),function(x){
  x[which.max(x$val),]  
})

sweed_genesmaxwild <- ddply(sweed_geneswild, c("locus"),function(x){
  x[which.max(x$val),]  
})

sweed_genesmax$pval <- 2*pnorm(-abs(scale(sweed_genesmax$val)))
sweed_genesmaxwild$pval <- 2*pnorm(-abs(scale(sweed_genesmaxwild$val)))

sweed_genesf <- subset(sweed_genesmax, chr < 8 & pos < 777)
sweed_geneswildf <- subset(sweed_genesmaxwild, chr < 8 & pos < 777)

# scaling and circulating the genes

sweed_genesf$val_sc <- scaler(sweed_genesf$val, c(sw_1,sw_2), b = c(0,max(sweed_genesf$val)))
sweed_genesf <- circular(sweed_genesf)

sweed_geneswildf$val_sc <- scaler(sweed_geneswildf$val, c(sw_1,sw_2), b = c(0,max(sweed_geneswildf$val)))
sweed_geneswildf <- circular(sweed_geneswildf)

# extracting outlier / nonoutliers

#outliers_sweed_genes <- subset(sweed_genesmax, val > 1.105207)
#nonoutliers_sweed_genes <- subset(sweed_genesmax, val <= 1.105207)

outliers_sweed_genes <- subset(sweed_genesf, val > swmin)
nonoutliers_sweed_genes <- subset(sweed_genesf, val <= swmin)

#outliers_sweed_genes <- subset(sweed_genesf, pval < 0.05)
#nonoutliers_sweed_genes <- subset(sweed_genesf, pval >= 0.05)

#outliers_sweed_genes_wild <- subset(sweed_geneswildf, pval < 0.05)
#nonoutliers_sweed_genes_wild <- subset(sweed_geneswildf, pval >= 0.05)
outliers_sweed_genes_wild <- subset(sweed_geneswildf, val > swminw)
#nonoutliers_sweed_genes_wild <- subset(sweed_geneswildf, val <= 2.505731)

outliers_sweed_genes <- outliers_sweed_genes[!outliers_sweed_genes$locus %in% outliers_sweed_genes_wild$locus,]

sweed_genes_thresh <- outliers_sweed_genes[which.max(abs(outliers_sweed_genes$pval)),10]

## constructing scale for SweeD genes ####

swmin_gen <- floor(min(sweed_genesf$val))
swmax_gen <- ceiling(max(sweed_genesf$val))

## testing even/odd ####
if(swmax_gen %% 2 == 0){
  step <- swmax_gen/4
} else {
  swmax_gen <- swmax_gen+1
  step <- (swmax_gen)/4
}
# ####

swscale_genes <- seq(swmin_gen,swmax_gen,length.out=5)

plotswscale_genes <- data.frame(swval = swscale_genes, swscale_sc = scaler(swscale_genes, c(sw_1,sw_2), b=c(min(swscale_genes),max(swscale_genes))))
plotswscale_genes$swscale_sc <- round(plotswscale_genes$swscale_sc, digits = 2)

plotswscale_genes$text <- plotswscale_genes[,1]
#plotneiscale_genes$text[length(plotneiscale_genes$neival)] <- paste0(">",max(plotneiscale_genes[,1]))


## circular plot, test values, wild and domest

plot_circ <- function(){ggplot() +
  geom_line(data = circ_chrom, aes(y = -.1, x = pos, group=chr), size = .5, colour = "black") +
  geom_line(data = circ_chrom, aes(x=pos+1, y=1.1, group=chr), alpha=0.5) +
  geom_text(data=chr_label,aes(x=V1,y=-.2,label=chr), size=3) + 
  #geom_segment(data = wind_outliers_nei, aes(x=pos_circ, xend=pos_circ, y=0, yend=1.2), colour=c("#d95f02"), alpha=0.2, size=2) + # window outliers
  geom_point(data = nonoutliers_nei_gene, aes(x=pos_circ,y=ratio_sc), colour=c("darkgrey"), size = .2, alpha = .6) + # gene nonoutlier dots
  geom_point(data = outliers_nei_gene, aes(x=pos_circ,y=ratio_sc), colour=c("#8254f2"),  size = .5, alpha = 1) + # gene outlier dots
  geom_line(data = neipiwind101_wide, aes(x=pos_circ, y=ratio_sc, group=chr), colour=c("#8254f2"), size=.5, alpha=.8) + #shape=19, 
  geom_segment(data=ticks,aes(x=ticks,xend=ticks,y=1.1,yend=1.14), alpha=.5) +
  geom_text(data = ticks_mod, aes(x=ticks, y=1.2, label=label), size = 2, alpha = .7) +
  #geom_segment(aes(x = 0,xend = max(circ_chrom$pos), y = nei_wind_thresh, yend = nei_wind_thresh), linetype="dashed", size=.1) +
  geom_segment(aes(x = 0,xend = max(circ_chrom$pos), y = nei_gene_thresh, yend = nei_gene_thresh), linetype="dashed", size=.1) +
  geom_point(data = circular(alloutplot), aes(x = pos_circ, y = 1.35), size = 2, colour=c("brown"), alpha=0.7) +
  #geom_text(data =  cand_genes, aes(x= pos_circ, y = 1.6, label = locus), size = 2) +
  geom_segment(aes(x=1,xend=1,y=0,yend=1.01), size = 0.2, alpha =.5, linetype="dashed") + # vertical line
  geom_segment(data = plotneiscale, aes(x=0,xend=1, y=neiscale_sc, yend = neiscale_sc), size= .4, alpha = .5) + # ticks
  geom_text(data = plotneiscale, aes(x=0,y=neiscale_sc,label=text), size=1, alpha=.7, hjust=1) + # text
  geom_segment(data = plotneiscale_genes, aes(x=0, xend=1, y=neiscale_sc, yend = neiscale_sc), size= .4, alpha = .5) + # ticks
  geom_text(data = plotneiscale_genes, aes(x=0,y=neiscale_sc,label=text), size=1, alpha=.7, hjust=1) + # text
  #geom_point(data = nonoutliers_fst_genes, aes(x=pos_circ,y=val_sc), colour=c("darkgrey"), size = .2, alpha = .6) + # gene nonoutlier dots
  #geom_point(data = outliers_fst_genes, aes(x=pos_circ,y=val_sc), colour=c("red"),  size = .5, alpha = 1) + # gene outlier dots
  #geom_segment(aes(x = 0,xend = max(circ_chrom$pos), y = fst_genes_thresh, yend = fst_genes_thresh), linetype="dashed", size=.1) +
  #geom_segment(data = plotfstscale_genes, aes(x=0,xend=1, y=fstscale_sc, yend = fstscale_sc), size= .4, alpha = .5) + # ticks
  #geom_text(data=plotfstscale_genes, aes(x=0,y=fstscale_sc,label=text), size=2, alpha=.5, hjust=1, angle = 0) + # text
  geom_point(data = nonoutliers_sweed_genes, aes(x=pos_circ,y=val_sc), colour=c("darkgrey"), size = .2, alpha = .6) + # gene nonoutlier dots
  geom_point(data = outliers_sweed_genes, aes(x=pos_circ,y=val_sc), colour=c("red"),  size = .5, alpha = 1) + # gene outlier dots
  geom_segment(aes(x = 0,xend = max(circ_chrom$pos), y = sweed_genes_thresh, yend = sweed_genes_thresh), linetype="dashed", size=.1) +
  geom_segment(data = plotswscale_genes, aes(x=0,xend=1, y=swscale_sc, yend = swscale_sc), size= .4, alpha = .5) + # ticks
  geom_text(data=plotswscale_genes, aes(x=0,y=swscale_sc,label=text), size=1, alpha=.7, hjust=1, angle = 0) + # text
  coord_polar(theta = "x")+
  scale_colour_manual(values=c("#1b9e77", "#d95f02")) +
  scale_x_continuous(limits=c(0,max(circ_chrom$pos)+25)) +
  scale_y_continuous(limits=c(-.6,1.6)) +
  geom_point(data = faywu_gene_outliers, aes(x=pos_circ,y=val_sc), colour = "#d95f02", size = .5, alpha = .8) + # gene outlier dots
  geom_point(data = faywu_gene_rest, aes(x=pos_circ,y=val_sc), colour="grey", size = .5, alpha = .5) + # gene outlier dots
  geom_line(data=subset(faywuwind101, type %in% "domest"), aes(x=pos_circ, y=val_sc, group=chr), colour=c("#d95f02"), size=.5, alpha=.7) +
  geom_line(data=subset(faywuwind101, type %in% "wild"), aes(x=pos_circ, y=val_sc, group=chr), colour=c("#1b9e77"), size=.5, alpha= .7) +
  geom_segment(aes(x = 0,xend = max(circ_chrom$pos), y = faywu_genes_thresh, yend = faywu_genes_thresh), linetype="dashed", size=.1) +
  #geom_segment(data = confidH, aes(x = 0,xend = max(circ_chrom$pos), y = val_sc, yend = val_sc, colour = type), linetype="dashed", size=.1) +
  geom_segment(data = plotfayscale, aes(x=0,xend=1, y=fayscale_sc, yend = fayscale_sc), size= .4, alpha = .5) + # ticks
  geom_text(data=plotfayscale, aes(x=0,y=fayscale_sc,label=fayscale), size=1, alpha=.7, hjust=1, angle = 0) + # text
  theme_bw() %+replace% theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks=element_blank(),axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x=element_blank(), axis.text.x = element_blank(), panel.border = element_blank(), axis.line = element_blank())
}

######## Extracting outliers (mapped + unmapped)

## mapped + unmapped
#outliers_sweed_genes_all <- subset(sweed_genesmax, pval < 0.05)
outliers_sweed_genes_all <- subset(sweed_genesmax, val > swmin)
outliers_sweed_genes_all_wild <- subset(sweed_genesmaxwild, val > swminw)
outliers_sweed_genes_all <- subset(outliers_sweed_genes_all, ! locus %in% outliers_sweed_genes_all_wild$locus)
#outliers_fst_genes_all <- subset(fst_genes, pval < 0.05)
outliers_nei_gene_all <- subset(neipi_gene_widef_all, pval < 0.05)
faywu_gene_outliers_all <- subset(faywu_genes_all, val <=  hdmincds & type %in% "domest")
#faywu_genes_all <- subset(faywu_genes_all, type %in% "domest")
#faywu_gene_outliers_all <- subset(faywu_genes_all, pval <=  0.05)


outliers_over <- list()
outliers_over$sweed <- outliers_sweed_genes_all$locus
#outliers_over$fst <- outliers_fst_genes_all$locus
outliers_over$pi_ratio <- outliers_nei_gene_all$loc
outliers_over$faywu <- faywu_gene_outliers_all$locus

threetests <- merge(merge(sweed_genesmax[,c(1,3)],
                  neipi_gene_wide[,c(2,5,6)], by=c("locus")),
                   subset(faywu_genes_all, type %in% "domest")[,c(4,3)], by=c("locus")) %>%
  `colnames<-`(c("locus","sweed","pi_domest","pi_wild","faywu"))

outliers_over1 <- lapply(outliers_over,function(x){x[x %in% threetests$locus]})

x <- outliers_over
A <- x[[1]]
B <- x[[2]]
C <- x[[3]]
dev.off()
draw.triple.venn(
  area1 = length(A),
  area2 = length(B),
  area3 = length(C),
  #area3 = length(D),
  n12 <- length(intersect(A, B)),
  n13 <- length(intersect(B, C)),
  #n14 <- length(intersect(A, D)),
  n23 <- length(intersect(A, C)),
  #n24 <- length(intersect(B, D)),
  #n34 <- length(intersect(C, D)),
  n123 <- length(intersect(intersect(A, B), C)),
  #n124 <- length(intersect(n12, D)),
  #n134 <- length(intersect(n13, D)),
  #n234 <- length(intersect(n23, D)),
  #n1234 <- length(intersect(n123, D)),
  category = c(names(x)[1], names(x)[2], names(x)[3]),
  fill = c("orange", "red", "green"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("orange", "red", "green")
)


# x <- outliers_over
# A <- x[[1]]
# B <- x[[2]]
# C <- x[[3]]
# D <- x[[4]]
# dev.off()
# venn.plot <- draw.quad.venn(
#   area1 = length(A),
#   area2 = length(B),
#   area3 = length(C),
#   area4 = length(D),
#   n12 <- length(intersect(A, B)),
#   n13 <- length(intersect(A, C)),
#   n14 <- length(intersect(A, D)),
#   n23 <- length(intersect(B, C)),
#   n24 <- length(intersect(B, D)),
#   n34 <- length(intersect(C, D)),
#   n123 <- length(intersect(intersect(A, B), C)),
#   n124 <- length(intersect(intersect(A, B), D)),
#   n134 <- length(intersect(intersect(A, C), D)),
#   n234 <- length(intersect(intersect(B, C), D)),
#   n1234 <- length(intersect(n123, D)),
#   category = c(names(x)[1], names(x)[2], names(x)[3], names(x)[4]),
#   fill = c("orange", "red", "green", "blue"),
#   lty = "dashed",
#   cex = 2,
#   cat.cex = 2,
#   cat.col = c("orange", "red", "green", "blue")
# )


## testing overlap 4 tests

hypotest <- data.frame(test=c("sweed_fst","sweed_pi","sweed_fay","fst_pi","fst_fay","pi_fay"),
                       area1=c(rep(length(A),3), rep(length(B),2),length(C)),
                       area2=c(length(B),length(C),length(D),length(C),length(D),length(D)),
                       overlap=c(length(intersect(A, B)),length(intersect(A, C)),length(intersect(A, D)),length(intersect(B, C)),length(intersect(B, D)),length(intersect(C, D))))

over_pval <- ddply(hypotest,c("test"),function(x){
  phyper(x[,4]-1,x[,2],nrow(fourtests)-x[,2],x[,3], lower.tail=F)
  #phyper(x[3,]-1,area1,nrow(fourtests)-area1,area2)
})

## END testing overlap 4 tests

## testing overlap pairwise

twotest <- merge(sweed_genesmax[,c(1,3)],
                          #neipi_gene_wide[,c(2,5,6)], #by=c("locus"))
                    subset(faywu_genes_all, type %in% "domest")[,c(4,3)], by=c("locus")) #%>%
#  `colnames<-`(c("locus","pi_domest","pi_wild","faywu"))

outliers_over2 <- lapply(outliers_over[-2],function(x){x[x %in% twotest$locus]})

A <- outliers_over2[[1]]
B <- outliers_over2[[2]]

length(intersect(A,B))
length(unique(A,B))


phyper(length(intersect(A,B))-1, length(B), nrow(twotest)-length(B), length(A), lower.tail=F)



## END testing overlap pairwise

## nei pi outliers that have 0 snps in domest but many in wild - not detected in any scans!!

subset(neipi_gene_wide1, count_domest == 0 & count_wild > 10 & val_wild > 0.01)

##

## extracting outliers table for the publication

alloutid <- unique(c(as.character(outliers_sweed_genes_all$locus), 
                      #as.character(outliers_fst_genes_all$locus), 
                        as.character(outliers_nei_gene_all$locus), 
                          as.character(faywu_gene_outliers_all$locus)))



swsub <- subset(sweed_genesmax, locus %in% alloutid)[,c(1,3,7,8)]
swsub$test <- c("sweed")
#fstsub <- subset(fst_genes, locus %in% alloutid)[,c(3,1,5,6)]
#fstsub$test <- c("fst")
fwsub <- subset(faywu_genes_all, locus %in% alloutid  & type %in% "domest")[,c(4,3,8,9,1)]
neisub <- subset(neipi_gene_widef_all, locus %in% alloutid)[,c(2,9,3,4,1)]
colnames(neisub)[2] <- c("val")

head(swsub)
head(fwsub)
head(neisub)
head(allout)
allout <- rbind(swsub,fwsub,neisub)
allout <- dcast(allout, locus+chr+pos~test, value.var=c("val"))

allout <- merge(allout, data.frame(locus=outliers_sweed_genes_all$locus, sweed_out5=c("sweed")), by.x=c("locus"), by.y=1, all.x=TRUE) %>%
  merge(., data.frame(locus=outliers_nei_gene_all$locus, nei_out5=c("pi")), by.x=c("locus"), by.y=1, all.x=TRUE) %>%
  merge(., data.frame(locus=faywu_gene_outliers_all$locus, fay_out5=c("faywu")), by.x=c("locus"), by.y=1, all.x=TRUE)
allout <-  merge(allout, data.frame(locus=subset(outliers_sweed_genes_all, pval < 0.01)$locus, sweed_out1=c("sweed")), by.x=c("locus"), by.y=1, all.x=TRUE) %>%
  merge(., data.frame(locus=subset(outliers_nei_gene_all, pval < 0.01)$locus, nei_out1=c("pi")), by.x=c("locus"), by.y=1, all.x=TRUE) %>%
  merge(., data.frame(locus=subset(faywu_gene_outliers_all, val < -4.23)$locus, fay_out1=c("faywu")), by.x=c("locus"), by.y=1, all.x=TRUE)

head(allout)

allout[,c(7:12)] <-apply(allout[,c(7:12)], 2, function(x) as.character(x))

allout[is.na(allout)]<-""
allout <- unite(allout,simul,c(sweed_out5,nei_out5,fay_out5), sep=",")
allout <- unite(allout,p0.01,c(sweed_out1,nei_out1,fay_out1), sep=",")

allout[,c(7:8)] <- apply(allout[,c(7:8)],2, function(x) gsub(",,",",",x, perl=T) %>%
  gsub("^,","",.) %>%
  gsub(",$","",.))

genes <- read.delim("/biodata/dep_coupland/grp_korff/artem/important_files_misc/selection_genes", header=F, sep="\t")

allout <- merge(allout, genes, by.x=c("locus"),by.y=c("V2"), all.x = T)
allout[,c(9:ncol(allout))] <- apply(allout[,c(9:ncol(allout))], 2, function(x) as.character(x))
allout[is.na(allout)]<-""
write.table(allout,file="/biodata/dep_coupland/grp_korff/artem/cgc_enrichment_98_wild/outliers.sweed_fst_fay_nei.out", quote=F, sep="\t", row.names = F)

## extracting outliers for the circular plot - colours - wild populations

alloutplot <- subset(allout, chr > 0)

alloutplot <- alloutplot[grep(",",alloutplot$simul),]
alloutplot$chr <- as.numeric(alloutplot$chr)
alloutplot$pos <- as.numeric(alloutplot$pos)

nrow(alloutplot)

outliers_pops_plot <- subset(mldist_pop_list[[8]], gene %in% alloutplot$locus)

ddply(unique(outliers_pops_plot), c("gene"), function(x){
  length(unique(x$pop))
})

plot_circ()

cbind(alloutplot[order(alloutplot$chr,alloutplot$pos),][,c(2,3)], seq(1,16))
