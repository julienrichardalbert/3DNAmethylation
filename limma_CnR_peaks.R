
library(edgeR)
counts <- read.delim("/Users/jra/CTCF/CTCF_peaks_indexed_colnames.txt", row.names = 1)
head(counts)


d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
d0
cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(d0) # number of peaks we began with
dim(d) # number of peaks left
colnames(counts)
dnmt <- c('TKO', 'TKO', 'WT', 'WT', 'TKO', 'TKO', 'WT', 'WT' )
days  <- c('d0', 'd0', 'd0', 'd0', 'd4', 'd4', 'd4', 'd4')
group <- interaction(days, dnmt)
group
head(counts, 0)

#snames <- colnames(counts) # Sample names
#cultivar <- substr(snames, 1, nchar(snames) - 2)
#time <- substr(snames, nchar(snames) - 1, nchar(snames) - 1)
#group <- interaction(cultivar, time)
plotMDS(d, col = as.numeric(group), cex=0.5)

mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
tmp <- voom(d0, mm, plot = T)

fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupd4.TKO - groupd4.WT, levels = colnames(coef(fit)))
contr

tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.table(top.table, file = "/Users/jra/Desktop/limma_results.txt", row.names = T, sep = "\t", quote = F)


# awk '{OFS=FS="\t"} { if ($2>=1 && $6<0.05) print $0,"UP_TKO_LIMMA"; else print $0,"."}' limma_results.txt > limma_results_thresholded.txt

