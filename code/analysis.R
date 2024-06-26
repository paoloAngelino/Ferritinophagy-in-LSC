## Find differences from cycling to quiescent, make prototypes
## Derive classifier of cycling / quiescent
## Apply to public data
library(DESeq2)
library(edgeR)
library(tidyverse)
library(openxlsx)
library(nclust)
library(broom)
library(glmnet)
library(biomaRt) ## get packages below
library(EDASeq)
library(parallel)
library(oligo)
library(affycoretools)
library(arrayQualityMetrics)
## library(clariomshumantranscriptcluster.db) ## standard
## Custom from MBNI
## library(clariomshumanhsrefseq.db)
## library(clariomshumanhsrefseqprobe)
## library(pd.clariomshuman.hs.refseq)
## library(clariomshumanhsrefseqcdf)
library(future.apply)
library(survival)
library(survminer)
library(umap)

## Write outputs (possibly overwriting?)
output <- FALSE

plan(multicore, workers=4)


## Convenience function, esp. for os x (cairo)
tiffopen <- function(fname, w=20, h=12, tres=150) {
    tiff(fname,width=w, height=h, units="cm", res=tres, compression="lzw", type="cairo")
}

## If we want to make several tiff with same sizes/res, we can simplify by
## generating a function that remembers these params
tiffparam <- function(w, h, tres)
{
    function(fname) { tiffopen(fname, w, h, tres) }
}

## small function to return a continuous vector in quantiles with labels
## essentially a wrapper for cut() and quantiles()

quantizevector <- function(v, p=seq(0,1,by=0.25))
{
    vq <- quantile(v, probs=p)
    cut(v, breaks=vq, labels=paste0("Q",p[-1]*100), include.lowest=TRUE)
}

## load preprocessed data: LSC17 signature, TCGA and OHSU
lsc17 <- readRDS("~/Documents/People/JTamb/StemCell/rds/lsc17.rds")
tcga <- readRDS("~/Documents/R.lib/datasets/TCGA/AML/tcga.aml.rds")
ohsu <- readRDS("~/Documents/R.lib/datasets/aml_ohsu_2018/ohsu.aml.rds")

######## Read Jérôme data
qcells <- read.xlsx("Quiescent vs Cycling 220204.xlsx")

jmat <- as.matrix(qcells[,2:11])[!duplicated(qcells$Gene.Symbol),]
rownames(jmat) <- qcells$Gene.Symbol[!duplicated(qcells$Gene.Symbol)]

qmat <- jmat
snames <- colnames(qmat)
snames <- gsub("quiescent", "G0", snames)
snames <- gsub("cycling", "G1S", snames)
sampledf <- data.frame(sample=snames,
                       patient=gsub("[.]G0|[.]G1S","", snames),
                       state=gsub("IM[0-9]+[.]","",snames))
colnames(qmat) <- sampledf$sample
rownames(sampledf) <- sampledf$sample
sdf <- sampledf


## Do PCA
tmp <- prcomp(t(qmat))
plotdf <- data.frame(sample=rownames(tmp$x), PC1=tmp$x[,"PC1"], PC2=tmp$x[,"PC2"])
plotdf$state <- sdf$state ; plotdf$patient <- sdf$patient

if (output) {
    tiffopen("scatter-pca-patient.tif", w=20, h=18, tres=300)
    ggplot(plotdf, mapping=aes(x=PC1, y=PC2, color=patient)) + geom_point()
    dev.off()
}

## Just remove 83, clear outlier
qmat <- qmat[,-grep("83",sdf$sample)]
sdf <- sdf[colnames(qmat),]

## Rescale per patient, to remove patient specific effects
## adjmat <- t(apply(qmat, 1, function(x) rep(tapply(x, gsub("[.]quiescent|[.]cycling","", colnames(qmat)), mean), each=2)))
adjmat <- t(apply(qmat, 1, function(x) rep(tapply(x, sdf$patient, mean), each=2)))
rmat <- qmat - adjmat


ccl <- nclust::coldmap(t(scale(t(qmat))), clab=colnames(qmat))
if (output) {
    tiffopen("coldmap-qmat.tif", w=20, h=20, tres=300)
    coldmap(t(scale(t(qmat))), clust=ccl, clab=colnames(qmat), rlab=c("MKI67", "TOP2A", "BUB1", "CENPA", "CD3", "CD33", "CD5"))
    dev.off()
}

mvdf <- data.frame(m=rowMeans(qmat), v=rowVars(qmat))
keep.genes <- rep(TRUE, nrow(qmat))

tmpl <- future_apply(qmat[keep.genes,], 1, function(x) tidy(lm(x ~ patient + state, sdf)))
## tmpl <- apply(qmat, 1, function(x) tidy(lm(x ~ state, sdf)))
genedf <- bind_rows(Map(function(x) x[grep("state",x$term),], tmpl))
genedf$term <- gsub("state", "", genedf$term)
genedf$gene <- rownames(qmat)[keep.genes]
genedf$p.adj <- p.adjust(genedf$p.value, method="fdr")
genedf$patient.stable <- unlist(Map(function(x) sum(x$p.value[grep("patient",x$term)]>0.1), tmpl))

## Also try limma
library(limma)
design <- cbind(rep(1, ncol(qmat)), rep(1:(ncol(qmat)/2), each=2), as.numeric(grepl("G1S", colnames(qmat))))
colnames(design) <- c("intercept" ,"patient", "state")
rownames(design) <- colnames(qmat)

model <- lmFit(qmat, design=design)
plotSA(model)
## Filter at about 5 expression? This will definitely improve the power!
## model <- eBayes(model[model$Amean>5,], trend=TRUE, robust=TRUE)
model <- eBayes(model, trend=TRUE, robust=TRUE)
## model <- eBayes(model, trend=FALSE, robust=TRUE)

if (output) {
    tiffopen("limma-sa.tif",w=16, h=12, tres=300)
    plotSA(model)
    dev.off()
}

limmadf <- topTable(model, coef="state", n=1000, sort.by="P", lfc=1, p.value=0.1)
sgenes <- rownames(limmadf)
sgenes <- sgenes[!grepl("HIST", sgenes)]
hgenes <- sgenes[grepl("HIST", sgenes)]
limmadf$gene <- rownames(limmadf)

if (output) {
    saveRDS(limmadf, "rds/limmadf.rds")
    saveRDS(sgenes, "rds/sgenes.rds")
    write.xlsx(limmadf, "limma-output.xlsx", overwrite=TRUE)
    tiffopen("volcanoplot.tif", w=20, h=16, tres=300)
    volcanoplot(model, coef="state", highlight=50, names=rownames(qmat), pch=16)
    dev.off()
}

## define our signature from limma output
cgenes <- limmadf %>% filter(logFC>0, gene %in% sgenes)
qgenes <- limmadf %>% filter(logFC<0, gene %in% sgenes)

if (output) {
    sink("cycling-sig.txt")
    dput(cgenes$gene)
    sink()
    saveRDS(cgenes,"rds/cycling-genes.rds")

    sink("quiescent-sig.txt")
    dput(qgenes$gene)
    sink()
    saveRDS(qgenes,"rds/quiescent-genes.rds")

    tiffopen("coldmap-diffex.tif", w=20, h=18, tres=300)
    nclust::coldmap(t(scale(t(qmat)))[sgenes,], clab=colnames(qmat))
    dev.off()
}


cnames <- colnames(qmat)
cnames <- gsub("IM06","PDX1",cnames) ; cnames <- gsub("IM10","PDX2",cnames)
cnames <- gsub("IM84","PDX4",cnames) ; cnames <- gsub("IM93","PDX5",cnames)
cnames <- gsub("[.]","/", cnames)
colnames(qmat) <- cnames

## Output heatmap in PDF/SVG
if (output) {
    pdf("plots/coldmap-diffex.pdf", width=8, height=7, pointsize=12)
    nclust::coldmap(t(scale(t(qmat)))[sgenes,], clab=cnames)
    dev.off()
}

## With gene names
gnames.all <- list(High=list(limmadf %>% filter(logFC>0 & !grepl("HIST",gene)) %>% pull(gene)),
               Low=list(limmadf %>% filter(logFC<0 & !grepl("HIST",gene)) %>% pull(gene)))
gnames <- list(High=list("CDK.*","CENP.*","CCN.*", "BUB1*", "MKI67", "EZH2", "CHEK2",
                      "BRCA2", "AURK.*","APOBEC.*", "H2AFΧ", "MAD2L1", "BIRC5"),
               Low=gnames.all$Low)

if (output) {
    pdf("plots/coldmap-diffex-genes.pdf", width=8, height=9, pointsize=10)
    nclust::coldmap(t(scale(t(qmat)))[sgenes,], clab=cnames, rlab=gnames)
    dev.off()

    tiffopen("coldmap-rmat.tif", w=24, h=24, tres=300)
    set.seed(1234)
    nclust::coldmap(rmat, clab=colnames(rmat), rlab=gnames)
    dev.off()
}

## Now translate to "bulk"
set.seed(1234)
isg <- intersect(sgenes, rownames(tcga$expmat.hgnc.tpm))
ccl <- coldmap(cor(t(tcga$expmat.hgnc.tpm[isg,])), rlab=list(list(qgenes$gene),list(c("MKI67","BUB1","TOP2A","CENPA")))) ## , rdend.col=list(list(c(20:105), col="red"), list(c(118:137), col="blue")))

if (output) {
    tiffopen("coldmap-cor-tcga-sgenes.tif", w=28, h=28, tres=300)
    coldmap(cor(t(tcga$expmat.hgnc.tpm[isg,])), clust=ccl,
            rlab=list(list(qgenes$gene),list(c("MKI67","BUB1","TOP2A","CENPA"))))
    ## , rdend.col=list(list(c(20:105), col="red"), list(c(118:137), col="blue")))
    dev.off()
}
curated.qsig <- isg[ccl$rclust$order][106:137]


## Make 2 templates (cycling/quiescent)
cprot <- rowMeans(qmat[sgenes, grep("G1",colnames(qmat))], na.rm=TRUE)
qprot <- rowMeans(qmat[sgenes, grep("G0", colnames(qmat))], na.rm=TRUE)
## cprot <- rowMeans(t(scale(t(qmat)))[sgenes, grep("cycling",colnames(qmat))], na.rm=TRUE)
## qprot <- rowMeans(t(scale(t(qmat)))[sgenes, grep("quiescent", colnames(qmat))], na.rm=TRUE)
## cprot <- rowMeans(scale(qmat)[sgenes, grep("cycling",colnames(qmat))], na.rm=TRUE)
## qprot <- rowMeans(scale(qmat)[sgenes, grep("quiescent", colnames(qmat))], na.rm=TRUE)

tmat <- as.matrix(cbind(cprot, qprot))
colnames(tmat) <- c("G1","G0")
tmat <- tmat[!apply(tmat,1,function(x) any(is.na(x))),]
tmat <- tmat[!apply(tmat,1,function(x) any(is.nan(x))),]

## Verify that histone genes are not very well expressed in RNAseq
tcga.hgenes <- intersect(rownames(tcga$expmat.hgnc.tpm), hgenes)
tcga.mvdf <- do.call(rbind, apply(tcga$expmat.hgnc.tpm,1,function(x) return(data.frame(m=mean(x), v=var(x)))))

######## Now apply to tcga and ohsu
## TEMPLATE-based is slow, also not significant
## tmp <- calcsig(tcga$expmat.hgnc.tpm, bs.iter=20, metric="pearson")
tclindf <- tcga$clindf
tclindf$qg <- apply(tcga$expmat.hgnc.tpm[intersect(rownames(tcga$expmat.hgnc.tpm), qgenes$gene),],2,mean)
tclindf$cg <- apply(tcga$expmat.hgnc.tpm[intersect(rownames(tcga$expmat.hgnc.tpm), cgenes$gene),],2,mean)
## tclindf$qsig <- tmp$G0
## tclindf$csig <- tmp$G1
tclindf$cqcall <- tmp$confident.icms

## Also add the LSC17 signature
tclindf$lsc17 <- as.vector(lsc17$coef %*% tcga$expmat.hgnc.tpm[lsc17$genes,])

## new signature from scRNAseq
qs <- readRDS("~/Documents/People/JTamb/StemCell/rds/qcluster.sig.rds")
tqs <- intersect(qs, rownames(tcga$expmat.hgnc.tpm))
tclindf$qs <- colMeans(tcga$expmat.hgnc.tpm[tqs,])
tclindf$qs.call <- ifelse(tclindf$qs>median(tclindf$qs), "high", "low")

## quantized
tclindf$qg.quant <- quantizevector(tclindf$qg)
tclindf$qs.quant <- quantizevector(tclindf$qs)


## curated
qs2 <- c("ADGRE2", "BIN2", "CD37", "FRMD4A", "HOXA7", "MAPRE2", "RAB13",
         "SAT2", "TESPA1", "TFPI")
qs2 <- intersect(qs2, rownames(tcga$expmat.hgnc.tpm))
tclindf$qs2 <- colMeans(tcga$expmat.hgnc.tpm[qs2,])

isg <- intersect(sgenes, rownames(tcga$expmat.hgnc.tpm))
tmpmat <- cor(tcga$expmat.hgnc.tpm[isg,], qmat[isg,])
tclindf$nearest <- apply(tmpmat,1,function(x) sdf$state[which.max(x)])
tclindf$maxg0 <- apply(tmpmat,1,function(x) max(x[sdf$state=="G0"]))


model <- coxph(Surv(os.time, os.event) ~ qg, tclindf) ; summary(model)
model <- coxph(Surv(os.time, os.event) ~ qs, tclindf) ; summary(model) ## scRNA sig
model <- coxph(Surv(os.time, os.event) ~ qs2, tclindf) ; summary(model) ## scRNA sig
model <- coxph(Surv(os.time, os.event) ~ lsc17, tclindf) ; summary(model)
model <- coxph(Surv(os.time, os.event) ~ lsc17 + qg, tclindf) ; summary(model)
model <- coxph(Surv(os.time, os.event) ~ maxg0, tclindf) ; summary(model)
model <- coxph(Surv(os.time, os.event) ~ CALGB.risk + maxg0, tclindf) ; summary(model)

model <- coxph(Surv(os.time, os.event) ~ CALGB.risk + maxg0, tclindf) ; summary(model)
model <- coxph(Surv(os.time, os.event) ~ CALGB.risk + qg.call, tclindf) ; summary(model)
model <- coxph(Surv(os.time, os.event) ~ CALGB.risk + qg, tclindf) ; summary(model)
model <- coxph(Surv(os.time, os.event) ~ CALGB.risk + lsc17, tclindf) ; summary(model)

## tiffopen("surv-tcga-qg.tif", w=20, h=18, tres=300)
pdf("plots/surv-tcga-qg.pdf", w=8, h=7, onefile=FALSE)
tclindf$qg.call <- ifelse(tclindf$qg<5.7, "low", "high")
model <- coxph(Surv(os.time, os.event) ~ qg.call, tclindf) ; summary(model)
model <- survfit(Surv(os.time/30, os.event) ~ qg.call, tclindf)
ggsurvplot(model, tclindf, conf.int=FALSE, risk.table=TRUE, xlim=c(0,60), legend.labs=c("high", "low"), xlab="Months", break.time.by=6, ylab="Overall survival", legend.title="Quiescence\nsignature", palette=c("red","blue"), pval=TRUE)
dev.off()

pdf("plots/surv-tcga-qs.pdf", w=8, h=7, onefile=FALSE)
model <- coxph(Surv(os.time, os.event) ~ qs.call, tclindf) ; summary(model)
pv <- paste0("P=",formatC(tidy(model)$p.value, digits=3))
model <- survfit(Surv(os.time/30, os.event) ~ qs.call, tclindf)
ggsurvplot(model, tclindf, conf.int=FALSE, risk.table=TRUE, xlim=c(0,60), legend.labs=c("high", "low"), xlab="Months", break.time.by=6, ylab="Overall survival", legend.title="Filtered\nsignature", palette=c("red","blue"), pval=TRUE)
dev.off()


## tiffopen("surv-tcga-lsc17.tif", w=20, h=18, tres=300)
pdf("plots/surv-tcga-lsc17.pdf", w=8, h=7, onefile=FALSE)
tclindf$lsc17.call <- ifelse(tclindf$lsc17<0.63, "low", "high")
model <- coxph(Surv(os.time, os.event) ~ lsc17.call, tclindf) ; summary(model)
pv <- paste0("P=",formatC(tidy(model)$p.value, digits=3))
model <- survfit(Surv(os.time/30, os.event) ~ lsc17.call, tclindf)
ggsurvplot(model, tclindf, conf.int=FALSE, risk.table=TRUE, xlim=c(0,60), legend.labs=c("high", "low"), xlab="Months", break.time.by=6, ylab="Overall survival", legend.title="LSC17\nsignature", palette=c("red","blue"), pval=TRUE)
dev.off()


## tiffopen("surv-tcga-lsc17-qg.tif", w=20, h=18, tres=300)
pdf("plots/surv-tcga-lsc17-qg.pdf", w=10, h=8)
model <- coxph(Surv(os.time, os.event) ~ lsc17.call + qg.call, tclindf) ; summary(model)
model <- survfit(Surv(os.time/30, os.event) ~ lsc17.call + qg.call, tclindf)
ggsurvplot(model, legend.labs=c("H/H", "H/L", "L/H", "L/L"), tclindf, conf.int=FALSE, risk.table=TRUE, xlim=c(0,60), xlab="Months", break.time.by=6, ylab="Overall survival", legend.title="LSC17 / Quiescence", palette=c("red4","red2", "orange", "blue2"))
dev.off()

######## Automatic p-value output

sink("tcga-survival.txt")
model <- coxph(Surv(os.time, os.event) ~ qg, tclindf) ; print("QG --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ qg.call, tclindf) ; print("QG call") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ qg.quant, tclindf) ; print("QS --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ qs, tclindf) ; print("QS --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ qs.call, tclindf) ; print("QS call --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ qs.quant, tclindf) ; print("QS --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ lsc17, tclindf) ; print("LSC17 --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ lsc17 + qg, tclindf) ; print("LSC17 + QG --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ lsc17.call + qg.call, tclindf) ; print("LSC17 call + QG call --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ lsc17 + qs, tclindf) ; print("LSC17 + QS --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ lsc17.call + qs.call, tclindf) ; print("LSC17 call + QS call --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
sink()

## Make some exploratory plots of QG and QS with more quantiles
## for QG the top quantile seems most interesting
pdf("plots/surv-tcga-qg-quant.pdf", w=8, h=7, onefile=FALSE)

model <- coxph(Surv(os.time, os.event) ~ qg.quant, tclindf) ; tidy(model) ; anova(model)
model <- survfit(Surv(os.time/30, os.event) ~ qg.quant, tclindf)
ggsurvplot(model, tclindf, conf.int=FALSE, risk.table=TRUE, xlim=c(0,60), xlab="Months", break.time.by=6, ylab="Overall survival", legend.title="Quiescence\nsignature", palette=c("red1", "red4", "blue1", "blue4"), pval=TRUE)
## dev.off()

pdf("plots/surv-tcga-qs-quant.pdf", w=8, h=7, onefile=FALSE)
model <- coxph(Surv(os.time, os.event) ~ qs.quant, tclindf) ; summary(model)
pv <- paste0("P=",formatC(tidy(model)$p.value, digits=3))
model <- survfit(Surv(os.time/30, os.event) ~ qs.quant, tclindf)
ggsurvplot(model, tclindf, conf.int=FALSE, risk.table=TRUE, xlim=c(0,60), xlab="Months", break.time.by=6, ylab="Overall survival", legend.title="Filtered\nsignature", palette=c("red1", "red4", "blue1", "blue4"), pval=TRUE)
dev.off()


## Must add Verhaak et Metzeler, 2 signatures (qg et LSC17)
vr <- readRDS("~/Documents/People/JTamb/StemCell/rds/verhaak.rds")
vclindf <- vr$clindf

## export for Jerome to add eln?
## write.xlsx(vclindf, "verhaak-petros.xlsx")

vclindf$qg <- apply(vr$expmat[intersect(qgenes$gene, rownames(vr$expmat)),],2, mean)
vclindf$qg.call <- ifelse(vclindf$qg>median(vclindf$qg), "high", "low")
vclindf$qs <- apply(vr$expmat[intersect(qs, rownames(vr$expmat)),],2, mean)
vclindf$qs.call <- ifelse(vclindf$qg>median(vclindf$qs), "high", "low")
vlsc17 <- lsc17[lsc17$genes %in% rownames(vr$expmat),]
vclindf$lsc17 <- as.vector(vlsc17$coef %*% vr$expmat[vlsc17$genes,])
vclindf$lsc17.call <- ifelse(vclindf$lsc17>median(vclindf$lsc17),"high", "low")
vclindf$qg.quant <- quantizevector(vclindf$qg)
vclindf$qs.quant <- quantizevector(vclindf$qs)


sink("verhaak-survival.txt")
model <- coxph(Surv(os.time, os.event) ~ qg, vclindf) ; print("QG --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ qg.call, vclindf) ; print("QG call") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ qg.quant, vclindf) ; print("QS --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ qs, vclindf) ; print("QS --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ qs.call, vclindf) ; print("QS call --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ qs.quant, vclindf) ; print("QS --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ lsc17, vclindf) ; print("LSC17 --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ lsc17 + qg, vclindf) ; print("LSC17 + QG --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ lsc17.call + qg.call, vclindf) ; print("LSC17 call + QG call --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ lsc17 + qs, vclindf) ; print("LSC17 + QS --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
model <- coxph(Surv(os.time, os.event) ~ lsc17.call + qs.call, vclindf) ; print("LSC17 call + QS call --- ") ; tidy(model, exp=TRUE, conf.int=TRUE)
sink()

######## Verhaak plots

pdf("plots/surv-verhaak-qg.pdf", w=10, h=8)
model <- coxph(Surv(os.time, os.event) ~ qg.call, vclindf) ; summary(model)
pv <- paste0("P=",formatC(tidy(model)$p.value, digits=3))
model <- survfit(Surv(os.time/30, os.event) ~ qg.call, vclindf)
ggsurvplot(model, vclindf, conf.int=FALSE, risk.table=TRUE, xlim=c(0,60), legend.labs=c("high", "low"), xlab="Months", break.time.by=6, ylab="Overall survival", legend.title="Quiescence\nsignature", palette=c("red","blue"), pval=TRUE)
dev.off()

pdf("plots/surv-verhaak-qs.pdf", w=10, h=8)
model <- coxph(Surv(os.time, os.event) ~ qs.call, vclindf) ; summary(model)
model <- survfit(Surv(os.time/30, os.event) ~ qs.call, vclindf)
ggsurvplot(model, vclindf, conf.int=FALSE, risk.table=TRUE, xlim=c(0,60), legend.labs=c("high", "low"), xlab="Months", break.time.by=6, ylab="Overall survival", legend.title="Filtered\nsignature", palette=c("red","blue"), pval=TRUE)
dev.off()


pdf("plots/surv-verhaak-lsc17.pdf", w=8, h=7, onefile=FALSE)
model <- coxph(Surv(os.time, os.event) ~ lsc17.call, vclindf) ; summary(model)
pv <- paste0("P=",formatC(tidy(model)$p.value, digits=3))
model <- survfit(Surv(os.time/30, os.event) ~ lsc17.call, vclindf)
ggsurvplot(model, vclindf, conf.int=FALSE, risk.table=TRUE, xlim=c(0,60), legend.labs=c("high", "low"), xlab="Months", break.time.by=6, ylab="Overall survival", legend.title="LSC17\nsignature", palette=c("red","blue"), pval=TRUE)
dev.off()


######## Generate signature with glmnet
library(glmnet)
scgenes <- readRDS("~/Documents/People/JTamb/StemCell/rds/qcluster.genes.rds")
cgenes <- Reduce(union, scgenes)
cgenes <- intersect(cgenes, rownames(vr$expmat))
cgenes <- intersect(cgenes, rownames(tcga$expmat.hgnc.tpm))

keep.pt <- !is.na(vclindf$os.event)
keep.tc <- !is.na(tclindf$os.time)

## propagate ELN
##### I AM ASSUMING HERE THAT RISK.CYTOGEN IS ELN!!
#### NOT TRUE --> is different thing
vclindf$ELN <- NA
vclindf$age <- vclindf$Age
vclindf$WBC <- NA
vclindf$cyto.risk <- factor(c("1"="fav","2"="int","3"="adv")[vclindf$Risk.Cytogen],levels=c("fav","int","adv"))

tclindf$ELN <- tclindf$ELN2017
tclindf$ELN <- c("fav"="fav","int"="int","adverse"="adv")[tclindf$ELN]
tclindf$ELN <- factor(tclindf$ELN, levels=c("fav","int","adv"))
tclindf$age <- tclindf$age.diag
tclindf$cyto.risk <- NA ## should calculate from cytogenetic, but this is for Jerome

traindf <- rbind(vclindf[keep.pt,c("ID", "os.time","os.event", "age", "ELN", "WBC", "cyto.risk")], tclindf[keep.tc,c("ID", "os.time","os.event", "age", "ELN", "WBC", "cyto.risk")])
traindf$origin <- c(rep("v",sum(keep.pt)), rep("t", sum(keep.tc)))
traindf$age <- as.numeric(traindf$age)

traindf$set <- "test"

## original
## trainmat <- cbind(t(scale(t(vr$expmat[cgenes,keep.pt]),scale=FALSE)), tcga$expmat.hgnc.tpm[cgenes,keep.tc])
## mean center both
trainmat <- cbind(t(scale(t(vr$expmat[cgenes,keep.pt]),scale=FALSE)),
                  t(scale(t(tcga$expmat.hgnc.tpm[cgenes,keep.tc]), scale=FALSE)))
## full scale both
## trainmat <- cbind(t(scale(t(vr$expmat[cgenes,keep.pt]))),
##                   t(scale(t(tcga$expmat.hgnc.tpm[cgenes,keep.tc]))))
stopifnot(all(colnames(trainmat)==traindf$ID))

set.seed(1337) ; trainpt <- sample(1:nrow(traindf), size=450)
traindf[trainpt,"set"] <- "train"
traindf$os.time[traindf$os.time==0] <- 1 ### glmnet will not tolerate
testpt <- setdiff(1:nrow(traindf), trainpt)
stopifnot(all(1:nrow(traindf) %in% union(trainpt , testpt)))

tmpl <- apply(trainmat[cgenes,trainpt], 1, function(x) tidy(coxph(Surv(os.time, os.event) ~ x, traindf[trainpt,])))
tmpdf <- bind_rows(tmpl) ; tmpdf$gene <- cgenes ; tmpdf$p.adj <- p.adjust(tmpdf$p.value, method="fdr")
trainmat.new <- trainmat[tmpdf %>% filter(p.adj<0.1, estimate>0) %>% pull(gene),]

set.seed(1337)
model.glm <- cv.glmnet(t(trainmat[,trainpt]), with(traindf[trainpt,], Surv(os.time, os.event)), family="cox")

## new sig
tmp <- coef(model.glm, s="lambda.min")
rownames(tmp)[tmp[,1]!=0]
glmsig <- tmp[tmp[,1]!=0,]
## also translate lsc17
vlsc17 <- lsc17$coef
names(vlsc17) <- lsc17$genes

## export to text for Jerome
sink("glmsig.txt")
dput(glmsig)
sink()

saveRDS(glmsig, "rds/qs-glm.rds")

calcsig.glm <- function(x, sig)
{
    cg <- intersect(names(sig), rownames(x))
    csig <- sig[cg]

    as.vector(csig %*% x[names(csig),])
}

traindf$lsc17 <- calcsig.glm(trainmat, vlsc17)
traindf$qs35 <- calcsig.glm(trainmat, glmsig)
traindf$qs35.call <- ifelse(traindf$qs35>median(traindf$qs35), "high", "low")

## output for Clement
write.xlsx(traindf, "tcga-ver-qs35.xlsx")


#### Must add multivariable with forest plot
model <- coxph(Surv(os.time, os.event) ~ age + qs35, traindf[traindf$set=="train",]) ; summary(model) ; anova(model)
model <- coxph(Surv(os.time, os.event) ~ I(age>60) + qs35, traindf[traindf$set=="train",]) ; summary(model) ; anova(model)
model <- coxph(Surv(os.time, os.event) ~ ELN + qs35, traindf[traindf$set=="train",]) ; summary(model) ; anova(model)

library(forestplot)

## Multivariable with ELN in TCGA data (test+train)
model <- coxph(Surv(os.time, os.event) ~ (ELN) + qs35, traindf[traindf$origin=="t",]) ; summary(model) ; anova(model)
fordf <- tidy(model, conf.int=TRUE, exp=TRUE)
fordf$pred <- gsub("ELN","ELN ",fordf$term)
fordf$pred <- gsub("qs35","QS35",fordf$pred)
fordf$HR <- sprintf("%3.2f (%1.2f–%1.2f)",fordf$estimate, fordf$conf.low, fordf$conf.high)
fordf$PV <- sprintf("%3.4f", fordf$p.value)

pdf("plots/forest-tcga-ELN-QS35.pdf", width=5, height=3, onefile=FALSE)
fordf |> forestplot(labeltext=c("pred", "HR", "PV"), mean="estimate", lower="conf.low", upper="conf.high", zero=1) |> fp_add_header(pred="Variable", HR="HR", PV="P-value") |>
    fp_set_style(txt_gp=fpTxtGp(ticks=gpar(fontfamily="", cex=0.8)))
dev.off()

## Multivariable with age in test data
model <- coxph(Surv(os.time, os.event) ~ I(age>60) + qs35, traindf[traindf$set=="test",]) ; summary(model) ; anova(model)
fordf <- tidy(model, conf.int=TRUE, exp=TRUE)
fordf$pred <- fordf$term
## fordf$pred <- gsub("ELN","ELN ",fordf$term)
fordf$pred[1] <- "Age > 60"
fordf$pred <- gsub("qs35","QS35",fordf$pred)
fordf$HR <- sprintf("%3.2f (%1.2f–%1.2f)",fordf$estimate, fordf$conf.low, fordf$conf.high)
fordf$PV <- sprintf("%3.4f", fordf$p.value)

pdf("plots/forest-test-age-QS35.pdf", width=5, height=3, onefile=FALSE)
fordf |> forestplot(labeltext=c("pred", "HR", "PV"), mean="estimate", lower="conf.low", upper="conf.high", zero=1) |> fp_add_header(pred="Variable", HR="HR", PV="P-value") |>
    fp_set_style(txt_gp=fpTxtGp(ticks=gpar(fontfamily="", cex=0.8)))
dev.off()


## Multivariable with age + WBV + ELN in TCGA
model <- coxph(Surv(os.time, os.event) ~ I(age>60) + I(WBC>10) + ELN + qs35, traindf[traindf$origin=="t",]) ; summary(model) ; anova(model)
fordf <- tidy(model, conf.int=TRUE, exp=TRUE)
fordf$pred <- gsub("ELN","ELN ",fordf$term)
fordf$pred <- gsub("qs35","QS35",fordf$pred)
fordf$pred[1] <- "Age > 60"
fordf$pred[2] <- "WBC > 10"
fordf$HR <- sprintf("%3.2f (%1.2f–%1.2f)",fordf$estimate, fordf$conf.low, fordf$conf.high)
fordf$PV <- sprintf("%3.4f", fordf$p.value)

pdf("plots/forest-tcga-age-WBC-ELN-QS35.pdf", width=5, height=4, onefile=FALSE)
fordf |> forestplot(labeltext=c("pred", "HR", "PV"), mean="estimate", lower="conf.low", upper="conf.high", zero=1) |> fp_add_header(pred="Variable", HR="HR", PV="P-value") |>
    fp_set_style(txt_gp=fpTxtGp(ticks=gpar(fontfamily="", cex=0.8)))
dev.off()


## Multivariable with age + WBV + ELN==adverse in TCGA
model <- coxph(Surv(os.time, os.event) ~ I(age>60) + I(WBC>10) + I(ELN=="adv") + qs35, traindf[traindf$origin=="t",]) ; summary(model) ; anova(model)
fordf <- tidy(model, conf.int=TRUE, exp=TRUE)
fordf$pred <- gsub("ELN","ELN ",fordf$term)
fordf$pred <- gsub("qs35","QS35",fordf$pred)
fordf$pred[1] <- "Age > 60"
fordf$pred[2] <- "WBC > 10"
fordf$pred[3] <- "ELN adverse"
## fordf$pred[4] <- "ELN intermediate"
fordf$HR <- sprintf("%3.2f (%1.2f–%1.2f)",fordf$estimate, fordf$conf.low, fordf$conf.high)
fordf$PV <- sprintf("%3.4f", fordf$p.value)

pdf("plots/forest-tcga-age-WBC-ELNadv-QS35.pdf", width=5.5, height=4.5, onefile=FALSE)
fordf |> forestplot(labeltext=c("pred", "HR", "PV"), mean="estimate", lower="conf.low", upper="conf.high", zero=1) |> fp_add_header(pred="Variable", HR="HR", PV="P-value") |>
    fp_set_style(txt_gp=fpTxtGp(ticks=gpar(fontfamily="", cex=0.8)))
dev.off()


## Multivariable with cytogenetics + QS35 in Verhaak
model <- coxph(Surv(os.time, os.event) ~ cyto.risk + qs35, traindf[traindf$origin=="v",]) ; summary(model) ; anova(model)
fordf <- tidy(model, conf.int=TRUE, exp=TRUE)
fordf$pred <- fordf$term
fordf$pred <- gsub("cyto.risk","Cyto ",fordf$term)
fordf$pred <- gsub("qs35","QS35",fordf$pred)
fordf$HR <- sprintf("%3.2f (%1.2f–%1.2f)",fordf$estimate, fordf$conf.low, fordf$conf.high)
fordf$PV <- sprintf("%3.4f", fordf$p.value)

pdf("plots/forest-verh-cyto-QS35.pdf", width=5, height=3, onefile=FALSE)
fordf |> forestplot(labeltext=c("pred", "HR", "PV"), mean="estimate", lower="conf.low", upper="conf.high", zero=1) |> fp_add_header(pred="Variable", HR="HR", PV="P-value") |>
    fp_set_style(txt_gp=fpTxtGp(ticks=gpar(fontfamily="", cex=0.8)))
dev.off()


## Multivariable with cytogenetics + QS35 in Verhaak
model <- coxph(Surv(os.time, os.event) ~ I(cyto.risk=="adv") + qs35, traindf[traindf$origin=="v",]) ; summary(model) ; anova(model)
fordf <- tidy(model, conf.int=TRUE, exp=TRUE)
fordf$pred <- fordf$term
fordf$pred[1] <- "Cyto adverse"
fordf$pred <- gsub("qs35","QS35",fordf$pred)
fordf$HR <- sprintf("%3.2f (%1.2f–%1.2f)",fordf$estimate, fordf$conf.low, fordf$conf.high)
fordf$PV <- sprintf("%3.4f", fordf$p.value)

pdf("plots/forest-verh-cytoadv-QS35.pdf", width=5, height=3, onefile=FALSE)
fordf |> forestplot(labeltext=c("pred", "HR", "PV"), mean="estimate", lower="conf.low", upper="conf.high", zero=1) |> fp_add_header(pred="Variable", HR="HR", PV="P-value") |>
    fp_set_style(txt_gp=fpTxtGp(ticks=gpar(fontfamily="", cex=0.8)))
dev.off()

## Multivariable with cytogenetics + QS35 in Verhaak
model <- coxph(Surv(os.time, os.event) ~ I(age>50) + I(cyto.risk=="adv") + qs35, traindf[traindf$origin=="v",]) ; summary(model) ; anova(model)
fordf <- tidy(model, conf.int=TRUE, exp=TRUE)
fordf$pred <- fordf$term
fordf$pred[1] <- "Age > 50"
fordf$pred[2] <- "Cyto adverse"
fordf$pred <- gsub("qs35","QS35",fordf$pred)
fordf$HR <- sprintf("%3.2f (%1.2f–%1.2f)",fordf$estimate, fordf$conf.low, fordf$conf.high)
fordf$PV <- sprintf("%3.4f", fordf$p.value)

pdf("plots/forest-verh-age50-cytoadv-QS35.pdf", width=5, height=3, onefile=FALSE)
fordf |> forestplot(labeltext=c("pred", "HR", "PV"), mean="estimate", lower="conf.low", upper="conf.high", zero=1) |> fp_add_header(pred="Variable", HR="HR", PV="P-value") |>
    fp_set_style(txt_gp=fpTxtGp(ticks=gpar(fontfamily="", cex=0.8)))
dev.off()

######## BEAT
beat <- readRDS("~/Documents/R.lib/datasets/AML/BEAT/beat.rds")
beat$clin$ELN <- beat$clin$ELN2017clean
beat$clin$age <- as.numeric(beat$clin$ageAtDiagnosis)
levels(beat$clin$ELN) <- c("fav","int","adv")
beat$clin$qs35 <- calcsig.glm(beat$expmat,  glmsig)

sink("beat-survival.txt")
model <- coxph(Surv(os.time, as.numeric(os.event)) ~ qs35, beat$clin) ; summary(model)
model <- coxph(Surv(os.time, os.event) ~ ELN + qs35, beat$clin) ; summary(model)
model <- coxph(Surv(os.time, os.event) ~ I(WBC>10) + qs35, beat$clin) ; summary(model) ; anova(model)
model <- coxph(Surv(os.time, os.event) ~ I(age>60) + ELN + qs35, beat$clin) ; summary(model) ; anova(model)
model <- coxph(Surv(os.time, os.event) ~ I(age>60) + I(WBC>10) + I(ELN=="adv") +  qs35, beat$clin) ; summary(model) ; anova(model)
sink()

## Multivariable with ELN and QS35 in BEAT
model <- coxph(Surv(os.time, os.event) ~ ELN + qs35, beat$clin) ; summary(model) ; anova(model)
fordf <- tidy(model, conf.int=TRUE, exp=TRUE)
fordf$pred <- gsub("ELN2017clean","ELN ",fordf$term)
fordf$pred <- gsub("qs35","QS35",fordf$pred)
fordf$HR <- sprintf("%3.2f (%1.2f–%1.2f)",fordf$estimate, fordf$conf.low, fordf$conf.high)
fordf$PV <- sprintf("%3.4f", fordf$p.value)

pdf("plots/forest-beat-ELN-QS35.pdf", width=5, height=3, onefile=FALSE)
fordf |> forestplot(labeltext=c("pred", "HR", "PV"), mean="estimate", lower="conf.low", upper="conf.high", zero=1) |> fp_add_header(pred="Variable", HR="HR", PV="P-value") |>
    fp_set_style(txt_gp=fpTxtGp(ticks=gpar(fontfamily="", cex=0.8)))
dev.off()


## Multivariable with age + ELN + QS35 in BEAT
model <- coxph(Surv(os.time, os.event) ~ I(age>60) + ELN + qs35, beat$clin) ; summary(model) ; anova(model)
fordf <- tidy(model, conf.int=TRUE, exp=TRUE)
fordf$pred <- gsub("ELN","ELN ",fordf$term)
fordf$pred <- gsub("qs35","QS35",fordf$pred)
fordf$pred[1] <- "Age > 60"
fordf$HR <- sprintf("%3.2f (%1.2f–%1.2f)",fordf$estimate, fordf$conf.low, fordf$conf.high)
fordf$PV <- sprintf("%3.4f", fordf$p.value)

pdf("plots/forest-beat-age-ELN-QS35.pdf", width=5, height=4, onefile=FALSE)
fordf |> forestplot(labeltext=c("pred", "HR", "PV"), mean="estimate", lower="conf.low", upper="conf.high", zero=1) |> fp_add_header(pred="Variable", HR="HR", PV="P-value") |>
    fp_set_style(txt_gp=fpTxtGp(ticks=gpar(fontfamily="", cex=0.8)))
dev.off()

## Multivariable with ELN=="adv"
model <- coxph(Surv(os.time, os.event) ~ I(age>60) + I(ELN=="adv") + qs35, beat$clin) ; summary(model) ; anova(model)
fordf <- tidy(model, conf.int=TRUE, exp=TRUE)
fordf$pred <- fordf$term
fordf$pred[1] <- "Age > 60"
fordf$pred[2] <- "ELN adv"
fordf$pred <- gsub("qs35","QS35",fordf$pred)
fordf$HR <- sprintf("%3.2f (%1.2f–%1.2f)",fordf$estimate, fordf$conf.low, fordf$conf.high)
fordf$PV <- sprintf("%3.4f", fordf$p.value)

pdf("plots/forest-beat-age60-ELNadv-QS35.pdf", width=5, height=3, onefile=FALSE)
fordf |> forestplot(labeltext=c("pred", "HR", "PV"), mean="estimate", lower="conf.low", upper="conf.high", zero=1) |> fp_add_header(pred="Variable", HR="HR", PV="P-value") |>
    fp_set_style(txt_gp=fpTxtGp(ticks=gpar(fontfamily="", cex=0.8)))
dev.off()

## Also add WBC
beat$clin$WBC <- beat$clin$wbcCount
model <- coxph(Surv(os.time, os.event) ~ I(age>60) + I(WBC>10) + I(ELN=="adv") +  qs35, beat$clin) ; summary(model) ; anova(model)
fordf <- tidy(model, conf.int=TRUE, exp=TRUE)
fordf$pred <- fordf$term
fordf$pred[1] <- "Age > 60"
fordf$pred[2] <- "WBC > 10"
fordf$pred[3] <- "ELN adv"
fordf$pred <- gsub("qs35","QS35",fordf$pred)
fordf$HR <- sprintf("%3.2f (%1.2f–%1.2f)",fordf$estimate, fordf$conf.low, fordf$conf.high)
fordf$PV <- sprintf("%3.4f", fordf$p.value)

pdf("plots/forest-beat-age60-WBC-ELNadv-QS35.pdf", width=5, height=3, onefile=FALSE)
fordf |> forestplot(labeltext=c("pred", "HR", "PV"), mean="estimate", lower="conf.low", upper="conf.high", zero=1) |> fp_add_header(pred="Variable", HR="HR", PV="P-value") |>
    fp_set_style(txt_gp=fpTxtGp(ticks=gpar(fontfamily="", cex=0.8)))
dev.off()




## Plot seperate training and validation plots
pdf("plots/surv-train-qs35.pdf", w=8, h=7, onefile=FALSE)
sfit <- survfit(Surv(os.time/30, as.numeric(os.event)) ~ qs35.call, traindf[trainpt,])
ggsurvplot(sfit, risk.table=TRUE,
           legend.labs=c("high", "low"), xlab="Months", break.time.by=6, ylab="Overall survival",
           legend.title="QS35 signature\n Training", palette=c("red","blue"), pval=TRUE, xlim=c(0,72))
dev.off()

pdf("plots/surv-test-qs35.pdf", w=8, h=7, onefile=FALSE)
sfit <- survfit(Surv(os.time/30, as.numeric(os.event)) ~ qs35.call, traindf[testpt,])
ggsurvplot(sfit, risk.table=TRUE,
           legend.labs=c("high", "low"), xlab="Months", break.time.by=6, ylab="Overall survival",
           legend.title="QS35 signature\n Test", palette=c("red","blue"), pval=TRUE, xlim=c(0,72))
dev.off()

## Also ohsu
ohsu$clindf$lsc17 <- calcsig.glm(ohsu$expmat, vlsc17)
ohsu$clindf$lsc17.call <- ifelse(ohsu$clindf$lsc17>median(ohsu$clindf$lsc17), "high", "low")
ohsu$clindf$qs35 <- calcsig.glm(ohsu$expmat,  glmsig)
ohsu$clindf$os.event <- as.numeric(ohsu$clindf$os.event)
model <- coxph(Surv(os.time, as.numeric(os.event)) ~ qs35, ohsu$clindf) ; summary(model)
pv <- paste0("P=",formatC(tidy(model)$p.value, digits=3))

ohsu$clindf$qs35.call <- ifelse(ohsu$clindf$qs35>median(ohsu$clindf$qs35), "high", "low")

## automatic cox output
sink("signature-survival.txt")
print("Training set (used for the model)")
model <- coxph(Surv(os.time, os.event) ~ qs35, traindf[trainpt,]) ; tidy(model, exp=TRUE, conf.int=TRUE)
print("Test set (not used for the model but from same datasets)")
model <- coxph(Surv(os.time, os.event) ~ qs35, traindf[testpt,]) ; tidy(model, exp=TRUE, conf.int=TRUE)
print("Ohsu data (external validation)")
model <- coxph(Surv(os.time, os.event) ~ qs35, ohsu$clindf) ; tidy(model, exp=TRUE, conf.int=TRUE)
## Correlation?
cor.test(traindf$lsc17, traindf$qs35)
cor.test(ohsu$clindf$lsc17, ohsu$clindf$qs35)
model <- coxph(Surv(os.time, os.event) ~ lsc17 + qs35, ohsu$clindf) ; tidy(model, exp=TRUE, conf.int=TRUE)
sink()


pdf("plots/cor-ohsu-qs35-lsc17.pdf", w=8, h=8, onefile=FALSE)
ggplot(traindf, mapping=aes(x=lsc17, y=qs35)) + geom_point() + theme_bw()
dev.off()


pdf("plots/surv-ohsu-qs35.pdf", w=8, h=7, onefile=FALSE)
sfit <- survfit(Surv(os.time, as.numeric(os.event)) ~ qs35.call, ohsu$clindf)
ggsurvplot(sfit, risk.table=TRUE,
           legend.labs=c("high", "low"), xlab="Months", break.time.by=6, ylab="Overall survival",
           legend.title="QS35 signature", palette=c("red","blue"), pval=TRUE)
dev.off()


pdf("plots/surv-ohsu-lsc17-qs35.pdf", w=8, h=6, onefile=FALSE)
sfit <- survfit(Surv(os.time, as.numeric(os.event)) ~ lsc17.call + qs35.call, ohsu$clindf)
model <- coxph(Surv(os.time, as.numeric(os.event)) ~ lsc17.call + qs35.call, ohsu$clindf)
pv <- paste0("P(LSC17)=",formatC(tidy(model)$p.value[1], digits=2),"\n",
             "P(QS35)=",formatC(tidy(model)$p.value[2], digits=2))
## model <- coxph(Surv(os.time, as.numeric(os.event)) ~ lsc17 + qs35, ohsu$clindf)
## pv <- paste0("P(LSC17)=",formatC(tidy(model)$p.value[1], digits=2),"\n",
##             "P(QS35)=",formatC(tidy(model)$p.value[2], digits=2))

ggsurvplot(sfit, risk.table=FALSE,
           legend.labs=c("h/h", "h/l", "l/h", "l/l"), xlab="Months",
           break.time.by=6, ylab="Overall survival",
           legend.title="LSC17/QS35 signature", palette=c("red1","red4", "blue1", "blue4"),
           pval=pv)
dev.off()

## Also ROC curves
library(pROC)
roc.qs35 <- with(ohsu$clindf, roc(os.event, qs35))
roc.lsc17 <- with(ohsu$clindf, roc(os.event, lsc17))

pdf("plots/roc-ohsu-lsc17-qs35.pdf", w=6, h=6, onefile=FALSE)
plot(roc.lsc17)
plot(roc.qs35, add=TRUE, col="red")
dev.off()
roc.test(roc1=roc.lsc17, roc2=roc.yqs35)

## show genes in verhaak
coldmap(t(scale(t(vr$expmat))), rlab=intersect(names(glmsig), rownames(vr$expmat)))


stop();
