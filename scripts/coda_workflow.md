#### Setup

Before we start, we need a working directory.

Make a directory (folder) in a place where you can output files. For example, in your Documents folder.

Next, download the otu_table from data to your working directory

Now open R (or RStudio) and set your working directory by typing the commands in your Console:
````
setwd("...")
````

in which, the "..." is the specific pathway, e.g.,

in Windows: 
````
setwd("C:/Users/User Name/Documents/FOLDER")
````
in Macs:
````
setwd("/Users/User Name/Documents/FOLDER")
````

#### Load packages
We need to load the packages we'll be using. You will have had to previously installed the packages.
````
library(compositions)
````
And then
````
library(zCompositions)
````
### Part 1: Exploratory Compositional PCA biplot

#### Read in the data table
In R, you load your OTU table into a `data frame`. You need yo have some idea of what your table looks like before attempting to import. You may want to open the table in Excel to see what format it has and what headers.

````
d<-read.table("otu_table.txt", skip=1, header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, row.names=1, comment.char="")
````
- `d` is the name we are using for our dataframe
- `skip=1` Becasue we are using QIIME-style format, the table has an extra header line that must be skipped 
- `header=T` our table has a header row (with sample names)
- `sep="\t"` The columns are tab-separated
- `stringsAsFactors=F, quote = "", check.names=F` these are here for parsing
- `row.names=1` we are using the OTU IDs as rownames
- `comment.char=""` this ignores the comment # in the headers

Whenever you import a data table, you'll want to inspect it to ensure the data looks correct. Try some of the following commands:

````
head(d)
dim(d)
nrow(d)
colnames(d)
````

We are going to filter our OTU table to remove very low count reads and zeros.

Count the number of zeros in the table: 
`sum(d == 0)`

Since there are lots of zeros in this table, we should filter out low abundace OTUs. There are a number of ways you could filter [some examples here](https://github.com/mmacklai/16S/blob/master/manipulating_counts_table.md)

````
# Remove OTUs <= mean read count
# But first, we need to remove the "taxonomy" column

rownames(d)<-paste(rownames(d), d$taxonomy, sep="_")
d$taxonomy <- NULL

count <- 10
d.1 <- data.frame(d[which(apply(d, 1, function(x){mean(x)}) > count),], check.names=F)
````
Now we will transform the data. We are using the the centred log-ratio, or CLR (Aitchison). This will turn our read counts into a relative abundace (abundance relative to the geometric mean of all OTUs per sample). This will allow us to retain the relationships between the individual components, but also puts the data in a geometric space where we can perform familiar statistics (and not violate CoDa).

We are taking a logm so the zeros must be replaced with an estimate value
````
#samples must be ROWS (so transpose)
d.czm <- cmultRepl(t(d.1),  label=0, method="CZM")
```
Calculate the CLR
```
#need to transpose because of cmultRepl
d.clr <- t(apply(d.czm, 1, function(x){log(x) - mean(log(x))}))
```
Principal component analysis (PCA)
````
#features are COLUMNS
d.pcx <- prcomp(d.clr)
````

Calculate the total variance and percent variance explained for our PCA axes
```r
d.mvar <- mvar(d.clr)
```

Let's plot the CoDa PCA biplot

````
# These are some hacks to get the plot to look OK
# Make the number of points equal to the number of features (for labels)
#   Use: "o" or "."

points <- c(rep(".", length(dimnames(d.pcx$rotation)[[1]])))

# We can also make the samples as points instead of labels

samples <- c(rep("o", length(dimnames(d.pcx$x)[[1]])))

# Color and text size for labels and points (vector of 2)
#   The first is the sample lables, the second is the points (OTUs). 
col=c("black",rgb(1,0,0,0.2))
c=c(0.5, 2) #Relative scale, 1 is 100%
````
````
biplot(d.pcx, cex=c, col=col, var.axes=F,
    xlab=paste("PC1: ", round(sum(d.pcx$sdev[1]^2)/d.mvar, 3)),
    ylab=paste("PC2: ", round(sum(d.pcx$sdev[2]^2)/d.mvar, 3)),
    scale=0, ylabs=points, xlabs=samples
)
````
If you would like a PDF of the plot, do the following:
````
pdf("PCA_plot.pdf")
# The plotting code goes here
dev.off()
````

### Part 2: ALDEx2 for differential expression analysis

Load the ALDEx library (needs to have been previously installed)
````
library(ALDEx2)
````
We will use the same table from before (d) and make a copy for ALDEx
````
aldex.in<-d
````
Make a vector of conditions. This must be in the same order and the same number as the columns (samples) of the input table (aldex.in)
````
conds<-c(rep("td", length(grep("td", colnames(d.1)))), rep("bm", length(grep("bm", colnames(d.1)))))

td<-rep("td", length(grep("td", colnames(d.1))))
bm<-rep("bm", length(grep("bm", colnames(d.1))))

conds<-c(td,bm)

#shortform
#conds<-c(rep("td", length(grep("td", colnames(d.1)))), rep("bm", length(grep("bm", colnames(d.1)))))
````

Get the clr values. This is using a Dirichlet estimate for zeros, rather than a constant point estimate we used in the biplots (czm function)
````
# this is the main ALDEx function for all downstream analyses
# mc.samples=128 is often sufficient. This is the number of Monte-Carlo samples from the Dirichlet distribution
x <- aldex.clr(d.1, mc.samples=128, verbose=TRUE)
````
Perform t-test between conditions
````
# both Welches and Wilcoxon t-test, plus a Benjamini-Hochberg multiple test correction
x.tt <- aldex.ttest(x, conds, paired.test=FALSE)x.effect <- aldex.effect(x, conds, include.sample.summary=FALSE, verbose=TRUE)
````
Estimate effect size and the within and between condition values
````
#include indiv. samples or not with include.sample.summary
x.effect <- aldex.effect(x, conds, include.sample.summary=TRUE, verbose=TRUE)

# Combine the results in one table
x.all <- data.frame(x.tt, x.effect)
````
Write a .txt table with your results
````
write.table(x.all, file="aldex_output.txt", sep="\t", quote=F, col.names=NA)
````

MA and MV (Effect) Plots
````
pdf("MA.pdf")
aldex.plot(x.all, type="MA", test="welch")
dev.off()

pdf("MW.pdf")
aldex.plot(x.all, type="MW", test="welch")
dev.off()
````
Sidebar: Find significant OTUs by p-value, and effect size cutoff
````
psig <- which(x.all$we.eBH < 0.05 & abs(x.all$effect) > 1)
# or
subset(x.all, x.all$we.eBH < 0.05 & abs(x.all$effect) > 1)
````


