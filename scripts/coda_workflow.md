Before we start, we need a working directory.

Make a directory (folder) in a place where you can output files. For example, in your Documents folder


library(compositions)
library(zCompositions)


#d<-read.table("test_data/otu_table.txt", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, row.names=1, skip=1, comment.char="")
d<-read.table("OTU_table_brazil.txt", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, row.names=1, skip=1, comment.char="")
d<-read.table("OTU_td_bm.txt", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, row.names=1)


head(d)
dim(d)
nrow(d)
colnames(d)


#remove refseqs with mean read count <=1
count <- 10
d.1 <- data.frame(d[which(apply(d[,1:ncol(d)-1], 1, function(x){mean(x)}) > count),], check.names=F)

#OR remove taxonomy column for now
rownames(d)<-paste(rownames(d), d$taxonomy, sep="_")
d$taxonomy <- NULL

d.1 <- data.frame(d[which(apply(d, 1, function(x){mean(x)}) > count),], check.names=F)


#replace zeros with estimate
#samples must be ROWS (so transpose)
d.czm <- cmultRepl(t(d.1),  label=0, method="CZM")

#calculate the CLR
#need to transpose because of cmultRepl
d.clr <- t(apply(d.czm, 1, function(x){log(x) - mean(log(x))}))


#????
#Need this to calculate %variance explained
d.mvar <- mvar(d.clr)


#calculate principal components
#features are COLUMNS
d.pcx <- prcomp(d.clr)

# Make the number of points equal to the number of features (for labels)
#use: "o" or "."
points <- c(rep(".", length(dimnames(d.pcx$rotation)[[1]])))

#color and text size for labels and points (vector of 2)
col=c("black",rgb(1,0,0,0.2))
c=c(0.6, 2)

pdf("test.pdf")
biplot(d.pcx, cex=c, col=col, var.axes=F,
    xlab=paste("PC1: ", round(sum(d.pcx$sdev[1]^2)/d.mvar, 3)),
    ylab=paste("PC2: ", round(sum(d.pcx$sdev[2]^2)/d.mvar, 3)),
    scale=0, ylabs=points
)
dev.off()


library(ALDEx2)

conds<-c(rep("td", length(grep("td", colnames(d.1)))), rep("bm", length(grep("bm", colnames(d.1)))))

#or
td<-rep("td", length(grep("td", colnames(d.1))))
bm<-rep("bm", length(grep("bm", colnames(d.1))))
conds<-c(td,bm)


x <- aldex.clr(d.1, mc.samples=128, verbose=TRUE)


x.tt <- aldex.ttest(x, conds, paired.test=FALSE)
x.effect <- aldex.effect(x, conds, include.sample.summary=FALSE, verbose=TRUE)

x.all <- data.frame(x.tt, x.effect)


psig <- which(x.all$we.eBH < 0.05 & abs(x.all$effect) > 1)

subset(x.all, x.all$we.eBH < 0.05 & abs(x.all$effect) > 1)

# A LOT OF DATA! What about control vs BV

#d.f<-as.matrix(d[,grep("cont", colnames(d))])

con<-colnames(d)[grep("cont", colnames(d))]
bv<-colnames(d)[grep("0_bv", colnames(d))]

d.f<-d[,c(con, bv)]

count <- 10
d.1 <- data.frame(d.f[which(apply(d, 1, function(x){mean(x)}) > count),], check.names=F)

d.czm <- cmultRepl(t(d.1),  label=0, method="CZM")

d.clr <- t(apply(d.czm, 1, function(x){log(x) - mean(log(x))}))
d.mvar <- mvar(d.clr)
d.pcx <- prcomp(d.clr)

l<-c("2_0_cont_NA", "8_0_cont_NA", "20_0_cont_NA", "34_0_cont_NA", "42_0_cont_NA", "48_0_cont_NA", "50_0_cont_NA", "51_0_cont_NA", "56_0_cont_NA", "64_0_cont_NA", "153_0_bv_y", "153_1_bv_y", "158_0_bv_n", "158_1_bv_n", "159_0_bv_n", "159_1_bv_n", "160_0_bv_y", "160_1_bv_y", "171_0_bv_n", "171_1_bv_n", "174_0_bv_y", "174_1_bv_y", "177_0_bv_y", "177_1_bv_y")

biplot(d.pcx, cex=c, col=col, var.axes=F,
    xlab=paste("PC1: ", round(sum(d.pcx$sdev[1]^2)/d.mvar, 3)),
    ylab=paste("PC2: ", round(sum(d.pcx$sdev[2]^2)/d.mvar, 3)),
    scale=0, ylabs=points
)
