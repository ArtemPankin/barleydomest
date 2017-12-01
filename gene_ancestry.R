#### This script generates Fig. 2de2de 
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
library(maptools)
library(foreign)
library(grid)
library(broom)
library(RColorBrewer)
library(gplots)
library(stringr)

## blank theme

blk_theme <- function(){theme(
  line = element_blank(),
  rect = element_blank(),
  text = element_blank(),
  axis.ticks.length = unit(0, "cm"),
  #axis.ticks.margin = unit(0.01, "cm"),
  legend.position = "none",
  panel.spacing = unit(0, "lines"),
  plot.margin = unit(c(0, 0, -.5, 0), "lines"),
  complete = TRUE)
}

## load geo coordinates of wild accessions

#wild_coord <- read.table("/netscratch/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/wild_no_adm90_coord_adj")
wild_coord <- read.table("~/wild_coord") ## download file from GitHub repo
colnames(wild_coord) <- c("id","lat","long")

## loading population colours
popcol <- brewer.pal(9,name="Set1")
attributes(popcol) <- list(names=paste0("pop",c(1:9)))

## loading ML distance data; output of !!!.sh script

mldist <- read.table("MLdist.data")

#mldist <- read.table("/netscratch/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/phylogeny_all/gene_by_gene/distances/all.distances.coord")

## adding names to the data.frame columns

mldist <- mldist[,seq(3,6)]
colnames(mldist) <- c("idd","idw","dist","gene")

#### since the ML distance file is very large, 
#### subsetting is done using memory-efficient data.table package

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
# dir=c("/netscratch/dep_coupland/grp_korff/artem/EIG5.0.2/final_struct_pca_ld/CLUMPAK_fig_wild_no_adm90_478sampl.filtered.SNP.478samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo.wo_sing_0.05_0.5_0.99")

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
## selecting genes from domesticated barley that hit population-specific wild allele

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
  y[,id1:=paste0(gene,"_",pop,"_",idd)]
})

mldist_pop_list_dt <- lapply(mldist_pop_list_dt, function(y){
  setkey(y, id1)
})

mldist_pop_list_dt <- lapply(mldist_pop_list_dt, function(y){
  y[, .SD[1,], by=id1]
})

mldist_pop_list <- lapply(mldist_pop_list_dt, function(y){
  as.data.frame(y)
}) # closing lapply

## END MODULE  select unique population hits

## set factor levels so that the populations remain in the same order

mldist_pop_list <- lapply(seq_along(mldist_pop_list), function(i){
  mldist_pop_list[[i]][,6] <- factor(mldist_pop_list[[i]][,6], levels = sort(unique(as.character(mldist_pop_list[[i]][,6]))))
  mldist_pop_list[[i]]
}) 

### MODULE: Figure 4a; plot distribution of ancestral wild haplotypes on a map of the Fertile Crescent

# loading map data from NaturalEarth

pol <- tidy(readShapePoly("~/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp", IDvar="ADMIN"))
ids <- read.dbf("~/ne_10m_admin_0_countries/ne_10m_admin_0_countries.dbf")
pol$id <- as.factor(pol$id)
dat <- merge(pol, ids, all.x=T, by.x="id", by.y="ADMIN")

dat.sub_rev <- subset(dat, id=="Turkmenistan" | id=="Uzbekistan" | id=="Tajikistan" |SUBREGION=="Western Asia" | id=="Egypt") # | id=="China" | id=="Sri Lanka" | SUBREGION=="Southern Asia")

rivers <- tidy(readShapeSpatial("~/natural_earth/ne_10m_rivers_lake_centerlines_scale_rank.shp"))
lakes <- tidy(readShapePoly("~/natural_earth/ne_10m_lakes.shp"))
bas3000 <- tidy(readShapePoly("~/natural_earth/ne_10m_bathymetry_H_3000.shp"))
bas2000 <- tidy(readShapePoly("~/natural_earth/ne_10m_bathymetry_I_2000.shp"))
bas1000 <- tidy(readShapePoly("~/natural_earth/ne_10m_bathymetry_J_1000.shp"))
bas200 <- tidy(readShapePoly("~/natural_earth/ne_10m_bathymetry_K_200.shp"))
bas0 <- tidy(readShapePoly("~/natural_earth/ne_10m_bathymetry_L_0.shp"))
ids_pol <- read.dbf("~/natural_earth/ne_10m_geography_regions_polys.dbf")

rivers_sub <- subset(rivers, lat < max(dat.sub_rev$lat) & lat > min(dat.sub_rev$lat) & long > min(dat.sub_rev$long) & long < max(dat.sub_rev$long))
lakes_sub <- subset(lakes, lat < max(dat.sub_rev$lat) & lat > min(dat.sub_rev$lat) & long > min(dat.sub_rev$long) & long < max(dat.sub_rev$long))
bas3000_sub <- subset(bas3000, lat < max(dat.sub_rev$lat) & lat > min(dat.sub_rev$lat) & long > min(dat.sub_rev$long) & long < max(dat.sub_rev$long))
bas2000_sub <- subset(bas2000, lat < max(dat.sub_rev$lat) & lat > min(dat.sub_rev$lat) & long > min(dat.sub_rev$long) & long < max(dat.sub_rev$long))
bas200_sub <- subset(bas2000, lat < max(dat.sub_rev$lat) & lat > min(dat.sub_rev$lat) & long > min(dat.sub_rev$long) & long < max(dat.sub_rev$long))
bas1000_sub <- subset(bas1000, lat < max(dat.sub_rev$lat) & lat > min(dat.sub_rev$lat) & long > min(dat.sub_rev$long) & long < max(dat.sub_rev$long))
bas0_sub <- subset(bas0, lat < max(dat.sub_rev$lat) & lat > min(dat.sub_rev$lat) & long > min(dat.sub_rev$long) & long < max(dat.sub_rev$long))

## extract dataframe for K=9

mldist_pop_8_all <- mldist_pop_list[[8]]

## add geo coordinates

mldist_pop_8_all <- inner_join(mldist_pop_8_all, wild_coord, by=c("idw"="id"))

## add location id, latitude_longitude

mldist_pop_8_all$lat_long <- paste0(mldist_pop_8_all$lat,"_", mldist_pop_8_all$long)

## extract ancestral haplotypes specific to a single location

mldist_pop_8_all <- ddply(mldist_pop_8_all, c("gene"), function(x){
  a <- length(unique(x$lat_long))
  if(a==1){
    x
  }
})

## extract number of wild haplotypes / location - size of a circle

mldist_toplot <- ddply(mldist_pop_8_all, c("lat_long"), nrow)

## link locations with population colours  - colour of a circle

mldist_toplot <- join(mldist_toplot,mldist_pop_8_all, match="first")[,c(1,2,8,10,11)]

## sort factors to match population colours

mldist_toplot$pop <- factor(mldist_toplot$po, levels = sort(unique(as.character(mldist_toplot$pop))))

## plot

ggplot(data= dat.sub_rev, aes(x=long,y=lat,group=group)) +
  geom_polygon(data=lakes_sub, size=.05, colour= "#A3C7CF",fill ="#ADD8E6") +
  geom_polygon(data=bas0_sub, size=.05, colour= "#acd6f3",fill ="#acd6f3") +
  geom_polygon(data=bas200_sub, size=.05, colour= "#8ec3e6",fill ="#8ec3e6", alpha=0.8) +
  geom_polygon(data=bas1000_sub, size=.05, colour= "#7cbae0",fill ="#7cbae0", alpha=0.8) +
  geom_polygon(data=bas2000_sub, size=.05, colour= "#63afd9",fill ="#63afd9", alpha=0.8) +
  geom_polygon(data=bas3000_sub, size=.05, colour= "#4aa3d2",fill ="#4aa3d2", alpha=0.8) +
  geom_polygon(size=.1, colour= "darkgrey",fill ="#eee79f") +
  geom_point(data=mldist_toplot, aes(y=lat, x=long, group="", size=V1, colour=pop), stroke = 0, alpha=.5) +
  scale_color_brewer(palette = "Set1") +
  geom_path(data=rivers_sub, aes(x=long, y=lat, group=id), colour="#A3C7CF", size=0.07) +
  theme_bw() %+replace% theme(panel.background=element_rect(fill="#eee79f", colour=NA),legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks=element_blank(),axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x=element_blank(), 
                              axis.text.x = element_blank(), panel.border = element_blank(), axis.line = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))

### END Fig. 4a; plot distribution of ancestral wild haplotypes on a map of the Fertile Crescent

## MODULE: Fig. 4b - longitudinal distribution of wild ancestral haplotypes

## load selection scan outliers; see plot_selection_scans.R script

alloutid <- read.table(file = "/biodata/dep_coupland/grp_korff/artem/scripts_git/korffgroup/alloutid")$x

## extract ancestral haplotypes for domestication loci

mldist_domest_loci <- subset(mldist_pop_8_all, gene %in% alloutid)

mldist_domest_loci <- ddply(mldist_domest_loci, c("gene"), function(x){
  a <- length(unique(x$lat_long))
  if(a==1){
    x
  }
})

## extract ancestral haplotypes for neutral loci

mldist_pop_8_all_neut <- subset(mldist_pop_8_all, ! gene %in% alloutid)

mldist_pop_8_all_neut <- ddply(mldist_pop_8_all_neut, c("gene"), function(x){
  a <- length(unique(x$lat_long))
  if(a==1){
    x
  }
})

## plot

ggplot() + 
  geom_density(data= mldist_pop_8_all_neut, aes(x=long, y=..scaled..),colour="#1B9E77", fill = "#1B9E77", alpha=.1, size=.5, adjust=.1) + 
  geom_density(data= mldist_domest_loci, aes(x=long, y=..scaled..), colour="#D95F02", fill = "#D95F02",  alpha=.1, size=.5, adjust=.1) + 
  theme_bw()

### END Fig. 4b - longitudinal distribution of wild ancestral haplotypes

## MODULE: Figs. S13, S14 - plot unsorted and sorted ancestry palletes for all domesticated genotypes

## plot unsorted palletes Fig. S13

ggplot(subset(mldist_pop_list[[8]],gene %in% alloutid)[,c(2,6)],aes(x=1)) +
  geom_bar(aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..], fill=pop)) +
  facet_grid(~idd) +
  scale_fill_manual(values=popcol) +
  scale_y_continuous(labels= percent) +
  theme_bw() %+replace% theme(strip.text.x = element_text(size=9,angle=90),strip.background = element_blank(), panel.background=element_rect(fill="transparent", colour=NA),legend.position="right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x=element_blank(),axis.title.y = element_blank(), axis.text.y = element_text(), axis.title.x=element_blank(), axis.text.x = element_blank(), panel.border = element_blank(), axis.line.x = element_blank())

## plot sorted palettes Fig. S14 (only the loci with with < 20% of unknown ancestry across the accessions)

c <- dcast(mldist_pop_list[[8]][,c(2,5,6)], gene ~ idd)

## sort palettes based on 1st genotype

q <- c[order(c[,2]),]

## calculate % of uknown ancestry across the accessions, select with < 20%

d <- as.data.frame(apply(q,2, function(x){
  x[is.na(x)] <- "unknown"
  x
}))

d$unknown <- apply(d, 1, function(x){
  sum(str_count(unlist(x),"unknown")) / length(x)
})

d <- subset(d, unknown < 0.2)

## plot (takes long time!)

dev.off()  
pushViewport(viewport(layout = grid.layout(nrow(d)+2, ncol(d)-1)))

invisible(lapply(seq(2,ncol(d)), function(y){
  invisible(lapply(seq(1,nrow(d)), function(x){
    f <-  d[x,c(1,y)]
    colnames(f)[2] <- "pop"
    e <- ggplot(data=f) + geom_rect(aes(xmin=1,xmax=2,ymin=1,ymax=2, fill=factor(pop))) +
      scale_fill_manual(values=popcol) + blk_theme()
    print(e, vp = viewport(layout.pos.row = c(x), layout.pos.col = c(y-1)))
  }))
}))

## END Figs. S13, S14 - plot unsorted and sorted ancestry palletes for all domesticated genotypes



## caclulating jaccard distances between domesticated genotypes & clustering

excl <- c("FT-477","FT-494","FT-492","FT-489","FT-500","FT-478","FT-498","FT-502","FT-483","FT-493","FT-486","FT-488")

## neutral loci

jacc_matrn <- subset(mldist_pop_list[[8]], ! gene %in% alloutid)
jacc_matrn$id2 <- paste0(jacc_matrn$gene,"_",jacc_matrn$pop)
jacc_matrn <- dcast(jacc_matrn[,c(2,5,8)], gene ~ idd) 
jacc_matrn <- jacc_matrn[, ! colnames(jacc_matrn) %in% excl]

corr_matrn <- sapply(colnames(jacc_matrn[,c(2:length(jacc_matrn))]), function(y){
  Y <- apply(jacc_matrn[,c(2:length(jacc_matrn))], 2, function(x){
    a <- gsub("_pop[0-9]", "",na.omit(jacc_matrn[,y]))
    b <- gsub("_pop[0-9]", "",na.omit(x))
    c <- intersect(a,b)
    I <- length(intersect(x,jacc_matrn[,y]))-1  
    J <- I/(length(pmatch(c,na.omit(x)))+length(pmatch(c,na.omit(jacc_matrn[,y])))-I)
    J
  })
})


# domestication loci

jacc_matrd <- subset(mldist_pop_list[[8]], gene %in% alloutid)
jacc_matrd$id2 <- paste0(jacc_matrd$gene,"_",jacc_matrd$pop)
jacc_matrd <- dcast(jacc_matrd[,c(2,5,8)], gene ~ idd) 
jacc_matrd <- jacc_matrd[, ! colnames(jacc_matrd) %in% excl]

corr_matrd <- sapply(colnames(jacc_matrd[,c(2:length(jacc_matrd))]), function(y){
  Y <- apply(jacc_matrd[,c(2:length(jacc_matrd))], 2, function(x){
    a <- gsub("_pop[0-9]", "",na.omit(jacc_matrd[,y]))
    b <- gsub("_pop[0-9]", "",na.omit(x))
    c <- intersect(a,b)
    I <- length(intersect(x,jacc_matrd[,y]))-1  
    J <- I/(length(pmatch(c,na.omit(x)))+length(pmatch(c,na.omit(jacc_matrd[,y])))-I)
    J
  })
})


a <- as.vector(corr_matrd)
b <- as.vector(corr_matrn)

as.numeric(corr_matr)

c <- data.frame(val=as.numeric(as.vector(corr_matrd)), name=rep(c("dom"),length(as.vector(corr_matrd))))
d <- data.frame(val=as.numeric(as.vector(corr_matrn)), name=rep(c("neutral"),length(as.vector(corr_matrn))))


e <- rbind(c,d)     

ggplot(e, aes(x=val)) + geom_density(aes(colour=name, fill=name, y=..scaled..), alpha=.4, size=1, adjust=2) + 
  geom_vline(xintercept=med_d, colour=c("#D95F02"), linetype="dashed", size=1) +
  geom_vline(xintercept=med_wd, colour=c("#1B9E77"), linetype="dashed", size=1) +
  scale_fill_manual(values=c("#D95F02","#1B9E77")) +
  scale_colour_manual(values=c("#D95F02","#1B9E77")) +
  theme_bw()




### from hereafter - draft



corr_matr_dist <- as.dist(corr_matr)

a <- hclust(corr_matr_dist)
b <- cmdscale(corr_matr_dist)

plot(a)
plot(b)


corr_matr[c("M-45"),c("FT-495")]



test <- heatmap.2(corr_matr, trace="none")
corr_matr[rev(test$rowInd),test$colInd]

corr_df <- as.data.frame(corr_matr)
corr_df$id <- row.names(corr_matr)

# END all



y <- colnames(jacc_matrd[,2:3])[1]
x <- jacc_matrd[,3]







ggplot(data = 
         as.data.frame(prcomp(corr_matrd, scale. =T)$x), 
       aes(x=PC1, y= PC2)) + geom_point()


med_d <- median(as.vector(as.dist(corr_matrd)))

heatmap.2(corr_matrd, trace="none")

boxplot(as.vector(as.dist(corr_matrd)))

corr_matr_distd <- as.dist(corr_matrd)
a <- hclust(corr_matr_distd)
b <- cmdscale(corr_matr_distd)

plot(a)
plot(b)


heatmap(corr_matrd)

corr_df <- as.data.frame(corr_matr)
corr_df$id <- row.names(corr_matr)

ggplot(data= melt(corr_df), aes(x=id, y = variable, fill=value)) + geom_tile()

# END domestication loci

# without domestication loci
jacc_matrwd <- subset(mldist_pop_list[[8]], ! gene %in% alloutid)
jacc_matrwd$id2 <- paste0(jacc_matrwd$gene,"_",jacc_matrwd$pop)
#jacc_matrwd <- jacc_matrwd[sample(1:20224,2389),]
jacc_matrwd <- dcast(jacc_matrwd[,c(2,5,8)], gene ~ idd) 
#jacc_matrwd <- jacc_matrwd[sample(1:1141,91),]
jacc_matrwd <- jacc_matrwd[, ! colnames(jacc_matrwd) %in% excl]

head(jacc_matrwd)

nrow(jacc_matrwd)




heatmap.2(corr_matrwd, trace="none",Rowv=T)

med_wd <- median(as.vector(as.dist(corr_matrwd)))

corr_medians <- sapply(seq(1:100),function(z){
  jacc_matrwd <- subset(mldist_pop_list[[8]], ! gene %in% alloutid)
  jacc_matrwd$id2 <- paste0(jacc_matrwd$gene,"_",jacc_matrwd$pop)
  #jacc_matrwd <- jacc_matrwd[sample(1:20224,2389),]
  jacc_matrwd <- dcast(jacc_matrwd[,c(2,5,8)], gene ~ idd) 
  jacc_matrwd <- jacc_matrwd[sample(1:1141,91),]
  jacc_matrwd <- jacc_matrwd[, ! colnames(jacc_matrwd) %in% excl]
  corr_matrwd <- sapply(colnames(jacc_matrwd[,c(2:length(jacc_matrwd))]), function(y){
    Y <- apply(jacc_matrwd[,c(2:length(jacc_matrwd))], 2, function(x){
      a <- gsub("_pop[0-9]", "",na.omit(jacc_matrwd[,y]))
      b <- gsub("_pop[0-9]", "",na.omit(x))
      c <- intersect(a,b)
      I <- length(intersect(x,jacc_matrwd[,y]))-1  
      J <- I/(length(pmatch(c,na.omit(x)))+length(pmatch(c,na.omit(jacc_matrwd[,y])))-I)
      J
    })
  })
  med <- median(as.vector(as.dist(corr_matrwd)))
  med
}, simplify = T)


ggplot(data=as.data.frame(corr_medians), aes(x=1.75,y=corr_medians)) + geom_boxplot(width=.5) + geom_point(position=position_jitter(width=.1)) +
  geom_segment(aes(x=1, xend=1.5,y=med_wd,yend=med_wd)) + theme_bw() +
  geom_segment(aes(x=0.5, xend=1, y=med_d, yend=med_d)) +
  xlab(label="Domesticated loci (91) / W/o domesticated loci (1141) / 100 x 91 W/o domesticated") +
  ylab(label="Ancestry similarity") +
  ylim(0.5,1)

corr_matr_distwd <- as.dist(corr_matrwd)
a <- hclust(corr_matr_distwd)
b <- cmdscale(corr_matr_distwd, k=2)

plot(a)
plot(b)

corr_dfwd <- as.data.frame(corr_matrwd)
corr_dfwd$id <- row.names(corr_matrwd)

ggplot(data= melt(corr_dfd), aes(x=id, y = variable, fill=value)) + geom_tile()

# END without domestication loci

a <- heatmap.2(corr_matr, trace="none")
b <- heatmap.2(corr_matrd, trace="none")
c <- heatmap.2(corr_matrwd, trace="none")

densityplot(subset (mldist_pop_list[[8]],  ! gene %in% alloutid)$dist)




## set factor levels so that the populations remain in the same order


mldist9pop <- mldist_pop_list[[8]]

new_popseq <- read.table("/biodata/dep_coupland/grp_korff/artem/scripts/ref_to_chrom_new_popseq")
mldist9pop_chr <- merge(mldist9pop, new_popseq[,c(1,4)], by.x = "gene", by.y = "V1")

## to plot sorted palettes 



## END
