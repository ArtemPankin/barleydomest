## This script plots site frequency spectra (sfs) based on the allele frequencies ~
## ~ calculated using two different approaches - based on snp calling & on angsd analysis of raw mapping files. 
## It evaluates whether SNP calling & filtering biased the sfs estimates and thus any downstream PopGen analyses
## Part 1 : takes as an input individual allele freqs obtained from a vcf file (output of ...sh)
## Part 2 : takes .sfs output of angds software that estimates allele frequencies from raw .bam files - w/o subjectivity of SNP filters etc
## Part 3 : plots both sfs to compare the results

## required files: output of angsd_sfs.sh & vcf_to_maf.sh scripts

library(ggplot2)

## SET the number of MAF bins (e.g. 25, step 0.02)

nbr <- 25

## SET the missing data threshold for the MAF estimates

miss_data <- 0.3

## LOAD FUNCTIONS

# function to bin individual allele frequencies (Part 1)

hist_lines <- function(x){
  a <- hist(x$MAF, plot = F, breaks = nbr)
  a$density <- a$counts/sum(a$counts)
  data.frame(mids=a$mids, val=a$density)
}

# function to normalize allele counts (Part 2)
norm <- function(x) x/sum(x)

#### PART 1 ####

setwd("~/user/sfs_comparison")

## read in the output of vcf_to_maf.sh script

snp_freq <- read.delim(file="snp_freq.snp_stats.out", header=F, sep="\t")

colnames(snp_freq) <- c("Reference", "Position", "Missing", "Ref_allele", "Het_allele", "Alt_allele", "Het_Freq", "Missing_Freq", "HomHet_Freq", "Hom_Freq", "MAF")

## filter for polymorphic SNPs and missing data frequency (e.g.)

snp_freq_sub <- snp_freq[snp_freq$MAF>0 & snp_freq$Missing_Freq <= miss_data,]

snp_freq_sub_plot <- hist_lines(snp_freq_sub)

### END - PART 1 ####

### PART 2 ####

## read in output of angsd_sfs.sh scrips
sfs_d <- scan("~/user/output.sfs")

## remove non-variable sites
sfs_d <- sfs_d[-c(1,length(sfs_d))]

## fold the data, extract minor allele frequencies (MAF)

a <- length(sfs_d)+1

sfs_d <- data.frame(freq=seq(1:length(sfs_d))/a, count=sfs_d)

b <- c(1:(length(sfs_d[,1])/2))
c <- rev(c(length(sfs_d[,1])/2+1):length(sfs_d[,1]))
sfs_dfold <- data.frame()
for(i in b){
  x <- b[i]
  y <- c[i]
  f <- data.frame(freq=sfs_d[x,1],count=sfs_d[x,2]+sfs_d[y,2])
  sfs_dfold <- rbind(sfs_dfold,f)
}

## convert counts in every MAF bin to proportions  

sfs_dfold$countnorm <- norm(sfs_dfold$count)

## bin the data into the same number of bins as in PART 1

sfs_dfold1 <- transform(sfs_dfold, bin = cut(sfs_dfold$freq, breaks=(c(snp_freq_sub_plot$mids-(0.25/nbr), .5))))
sfs_dfold_bin <- data.frame(val= ddply(sfs_dfold1, "bin", summarize, freq=sum(countnorm))$freq, mids=snp_freq_sub_plot$mids)

### END - PART 2 ####

### PART 3 ####

ggplot() +
  geom_ribbon(data=snp_freq_sub_plot,aes(x=mids,ymin=0,ymax=val), fill="#b96a16", alpha= .7) +
  geom_line(data=sfs_dfold_bin, aes(x=mids,y=val), colour ="black",linetype = "dashed", alpha = 1) +
  geom_line(data=snp_freq_sub_plot, aes(x=mids,y=val), colour ="#ab671f",alpha =.7) +
  scale_y_continuous(label = percent, limits=c(0,1), breaks=c(0.00,.25,.50,.75,1)) +
  ylab("Percentage") +
  xlab("MAF") +
  theme_bw() %+replace% theme(panel.border = element_blank(), 
                              panel.grid.major = element_blank(), 
                              legend.position="", 
                              axis.line.x = element_line(colour = "black"), 
                              axis.line.y = element_line(colour = "black"))

### END - PART 3

