#### This script generates Fig. !!! 
####

## set number of populations K
k=9

## set threshold to define population
thresh=0.99

## Load and, if necessary, install the following packages
library(data.table)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)

## load ML distance data; output of !!!.sh script

mldist <- read.table("MLdist.data")

## add names to the data.frame columns

mldist <- mldist[,seq(3,6)]
colnames(mldist) <- c("idd","idw","dist","gene")

#### since the ML distance file is very large, 
#### the subsetting is done using memory-efficient data.table package

## MODULE select only the smallest ML distances

mldist$id <- paste0(mldist$gene,"_",mldist$idd,"_",mldist$dist)

mldist_min <- ddply(mldist, c("gene", "idd"), function(x){
  min(x$dist)
})

mldist_min$id <- paste0(mldist_min$gene,"_",mldist_min$idd,"_",mldist_min$V1)

mldist_dt <- data.table(mldist)
setkey(mldist_dt,id)
mldist_minid_dt <- data.table(mldist_min)

mldist_minall <- mldist_dt[mldist_minid_dt$id]
mldist_minall <- as.data.frame(mldist_minall)
mldist_minall <- mldist_minall[,seq(1,4)]

## END MODULE select only the smallest ML distances

# MODULE load fastSTRUCTURE data for all K

## fastSTRUCTURE files with *.out extension; output of !!!.sh script
## data format: 1st col - genotype name; 2nd - (K+1) columns are membership coefficients

dir=c("/path/to/dir/with/structure/files")

## read file names into the list; sort K as numeric, important if K >= 10
strfi <- list.files(path=dir, pattern="*.out", full.names=T)
z <- sub(".*(K[0-9]*).*", "\\1", strfi)
strfi <- strfi[order(as.numeric(sub('K','',z)))]

## load data into the list format, K number of elements

str_data <- lapply(strfi, function(x){read.delim(x, header=F)})

str_data <- lapply(seq_along(str_data), function(i){
  colnames(str_data[[i]]) <- c("a", paste0("V",seq(1:(i+1))))
  str_data[[i]]
})

### Rdistruct
### sort the different K dataframes similar to DISTRUCT software; 
### important to have the same colour for population comparing differen K
### only works if the K series start with K=2

##  initial step of Rdistruct: 
##  finds the highest correlation betweew two adjacent K dataframes and 
##  assigns it as "pop1" in the K=2 dataframe

correlations <- cor(str_data[[2]][,2:ncol(str_data[[2]])],str_data[[1]][,2:ncol(str_data[[1]])])
max1 <- as.vector(which.max(apply(correlations, MARGIN=2, max)))+1
colnames(str_data[[1]])[max1] <- c("pop1")
vcol <- grep("V",names(str_data[[1]]))
colnames(str_data[[1]])[vcol] <- c("pop2")
str_data[[1]] <- str_data[[1]][,order(names(str_data[[1]]))]

## continues with sorting the rest of the dataframes

for(x in 1:(length(str_data)-1)){
  y=x+1
  for(i in 2:ncol(str_data[[x]])){
    correlations <- cor(str_data[[y]][,2:ncol(str_data[[y]])],str_data[[x]][,i])
    name <- names(str_data[[x]])[i]
    max1 <- which.max(apply(correlations, MARGIN=1, max))+1
    names(str_data[[y]])[max1] <- name
  }
  vcol <- grep("V",names(str_data[[y]]))
  name_new <- paste0("pop",ncol(str_data[[x]]))
  names(str_data[[y]])[vcol] <- name_new
  str_data[[y]] <- str_data[[y]][,order(names(str_data[[y]]))]  
}

for (i in 1:(k-1)){
  str_data[[i]] <- cbind(colsplit(str_data[[i]][,1], "_", c("trait","id")), str_data[[i]][,seq(1,ncol(str_data[[i]]),1)])
}

## END Rdistruct

## END MODULE load fastSTRUCTURE data for all K

## MODULE assign population id for wild genotypes based on the threshold

str_datapop <- lapply(str_data, function(y){
  ddply(y, c("id"), function(x){
    if(max(x[,4:ncol(x)]) > thresh){
      cbind(x,pop=colnames(x)[which(x==max(x[,4:ncol(x)]))])
    } # closing if   
  }) # closing ddply
}) # closing lapply

## END MODULE assign population id for wild genotypes based on the threshold


## MODULE  select unique population hits; 
## domestication gene hits only 1 population

mldist_pop_list <- lapply(str_datapop, function(y){
  inner_join(mldist_minall, y[,c(2,ncol(y))], by=c("idw"="id"))
}) # closing lapply

mldist_pop_list_dt <- lapply(mldist_pop_list, function(y){
  data.table(y)
}) # closing lapply

mldist_pop_list_dt <- lapply(mldist_pop_list_dt, function(y){
  y[,id:=paste0(idd,"_",gene)]
  setkey(y, id)
}) # closing lapply

mldist_1pop_id_list <- lapply(mldist_pop_list, function(y){
  ddply(y, c("idd", "gene"), function(z){
    a <- length(unique(sort(z$pop)))
    if(a == 1){
      a
    } # closing if
  }) # closing ddply
}) # closing lapply

mldist_1pop_id_list_dt <- lapply(mldist_1pop_id_list, function(y){
  data.table(y)
})

mldist_1pop_id_list_dt <- lapply(mldist_1pop_id_list_dt, function(y){
  y[,id:=paste0(idd,"_",gene)]
})

mldist_pop_list_dt <- mapply(function(x,y){
  x[y$id]
}, x=mldist_pop_list_dt,y=mldist_1pop_id_list_dt, SIMPLIFY = FALSE)

mldist_pop_list_dt <- lapply(mldist_pop_list_dt, function(y){
  setkey(y, id)
})

mldist_pop_list_dt <- lapply(mldist_pop_list_dt, function(y){
  y[, .SD[1,], by=id]
})

mldist_pop_list <- lapply(mldist_pop_list_dt, function(y){
  as.data.frame(y)
}) # closing lapply

## END MODULE  select unique population hits

## MODULE: plots

## set factor levels so that the populations remain in the same order

mldist_pop_list <- lapply(seq_along(mldist_pop_list), function(i){
mldist_pop_list[[i]][,6] <- factor(mldist_pop_list[[i]][,6], levels = sort(unique(as.character(mldist_pop_list[[i]][,6]))))
mldist_pop_list[[i]]
}) # closing lapply

##  plots for each value of K

plots <- lapply(seq_along(mldist_pop_list), function(i){
  ggplot(mldist_pop_list[[i]][,c(1,6)],aes(x=1)) +
    geom_bar(aes(y=(..count..)/sum(..count..),  fill=pop)) +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() %+replace% theme(panel.background=element_rect(fill="transparent", colour=NA),legend.position="right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x=element_blank(),axis.title.y = element_blank(), axis.text.y = element_text(), axis.title.x=element_blank(), axis.text.x = element_blank(), panel.border = element_blank(), axis.line.x = element_blank()) 
})

lapply(seq_along(plots), function(x){
  ggsave(filename=paste0("pop_K",(x+1),".pdf"),
         plot=plots[[x]],
         width = 5, height = 15, unit="cm")
})

## END plot barplots; unique population hits


## Fig_S!!! histogram - number of hits wild populations / domesticated gene

gene_pop_hits <- ddply(mldist_pop_list[[8]], c("gene"), function(x){
  length(unique(x$pop))
})

ggplot(gene_pop_hits) +
  geom_bar(aes(x=V1, y=(..count..)/sum(..count..)), colour= "black", fill = NA) +
  scale_y_continuous(labels= percent) +
  theme_bw()

### END MODULE: histogram - number of hits wild populations / domesticated gene

## Fig_S!!!! plot for each genotype separately, only K9

ggplot(mldist_pop_list[[8]][,c(1,2,6)],aes(x=1)) +
  geom_bar(aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..], fill=pop)) +
  scale_fill_brewer(palette = "Set1") +
  facet_grid(~idd) +
  scale_y_continuous(labels= percent) +
  theme_bw() %+replace% theme(strip.text.x = element_text(size=10,angle=90),strip.background = element_blank(), panel.background=element_rect(fill="transparent", colour=NA),legend.position="right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x=element_blank(),axis.title.y = element_blank(), axis.text.y = element_text(), axis.title.x=element_blank(), axis.text.x = element_blank(), panel.border = element_blank(), axis.line.x = element_blank())
