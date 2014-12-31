### R code from vignette source 'DiffBind.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()
savewd = getwd()


###################################################
### code chunk number 2: DiffBind.Rnw:95-99
###################################################
tmp =  tempfile(as.character(Sys.getpid()))
pdf(tmp)
savewarn = options("warn")
options(warn=-1)


###################################################
### code chunk number 3: DiffBind.Rnw:103-105
###################################################
library(DiffBind)
setwd(system.file("extra", package="DiffBind"))


###################################################
### code chunk number 4: DiffBind.Rnw:110-116 (eval = FALSE)
###################################################
## print(savewd)
## tamoxifen = dba(sampleSheet="tamoxifen.csv")
## tamoxifen = dba.count(tamoxifen)
## tamoxifen = dba.contrast(tamoxifen, categories=DBA_CONDITION)
## tamoxifen = dba.analyze(tamoxifen)
## tamoxifen.DB = dba.report(tamoxifen)


###################################################
### code chunk number 5: sampSheet
###################################################
samples = read.csv(file.path(system.file("extra", package="DiffBind"),
                             "tamoxifen.csv"))
samples


###################################################
### code chunk number 6: dbaContruct
###################################################
tamoxifen = dba(sampleSheet="tamoxifen.csv")


###################################################
### code chunk number 7: DiffBind.Rnw:140-141
###################################################
tamoxifen


###################################################
### code chunk number 8: tamox_occ_corhm
###################################################
plot(tamoxifen)


###################################################
### code chunk number 9: DiffBind.Rnw:164-165 (eval = FALSE)
###################################################
## tamoxifen = dba.count(tamoxifen, minOverlap=3) 


###################################################
### code chunk number 10: DiffBind.Rnw:168-169
###################################################
data(tamoxifen_counts)


###################################################
### code chunk number 11: tamox_aff_corhm
###################################################
plot(tamoxifen)


###################################################
### code chunk number 12: DiffBind.Rnw:183-184
###################################################
tamoxifen = dba.contrast(tamoxifen, categories=DBA_CONDITION)


###################################################
### code chunk number 13: tamox_sdb_corhm
###################################################
tamoxifen = dba.analyze(tamoxifen) 


###################################################
### code chunk number 14: DiffBind.Rnw:205-206
###################################################
tamoxifen    


###################################################
### code chunk number 15: DiffBind.Rnw:219-220
###################################################
tamoxifen.DB = dba.report(tamoxifen)


###################################################
### code chunk number 16: DiffBind.Rnw:225-226
###################################################
tamoxifen.DB


###################################################
### code chunk number 17: tamox_sdb_ma
###################################################
data(tamoxifen_analysis)
dba.plotMA(tamoxifen)


###################################################
### code chunk number 18: tamox_aff_pca
###################################################
dba.plotPCA(tamoxifen,DBA_TISSUE,label=DBA_CONDITION)


###################################################
### code chunk number 19: tamox_sdb_pca
###################################################
dba.plotPCA(tamoxifen, contrast=1,th=.05,label=DBA_TISSUE)


###################################################
### code chunk number 20: DiffBind.Rnw:293-295
###################################################
sum(tamoxifen.DB$Fold<0)
sum(tamoxifen.DB$Fold>0)


###################################################
### code chunk number 21: tamox_sdb_box
###################################################
pvals = dba.plotBox(tamoxifen)


###################################################
### code chunk number 22: DiffBind.Rnw:311-312
###################################################
pvals


###################################################
### code chunk number 23: DiffBind.Rnw:322-323
###################################################
corvals = dba.plotHeatmap(tamoxifen)


###################################################
### code chunk number 24: tamox_sdb_hm
###################################################
corvals = dba.plotHeatmap(tamoxifen, contrast=1, correlations=FALSE)


###################################################
### code chunk number 25: DiffBind.Rnw:346-349
###################################################
data(tamoxifen_counts)
tamoxifen = dba.contrast(tamoxifen,categories=DBA_CONDITION,
                         block=tamoxifen$masks$MCF7)


###################################################
### code chunk number 26: DiffBind.Rnw:354-356
###################################################
tamoxifen = dba.analyze(tamoxifen)
tamoxifen


###################################################
### code chunk number 27: tamox_block_ma
###################################################
dba.plotMA(tamoxifen,method=DBA_EDGER_BLOCK)


###################################################
### code chunk number 28: tamox_block_corhm
###################################################
dba.plotHeatmap(tamoxifen,contrast=1,method=DBA_EDGER_BLOCK,
                attributes=c(DBA_TISSUE,DBA_CONDITION,DBA_REPLICATE))


###################################################
### code chunk number 29: tamox_block_pca
###################################################
dba.plotPCA(tamoxifen,contrast=1,method=DBA_EDGER_BLOCK,
            attributes=DBA_CONDITION,label=DBA_TISSUE)


###################################################
### code chunk number 30: DiffBind.Rnw:397-399
###################################################
tamoxifen = dba.analyze(tamoxifen,method=DBA_ALL_METHODS)
tamoxifen


###################################################
### code chunk number 31: tamox_block_venn
###################################################
tam.block = dba.report(tamoxifen,method=DBA_ALL_BLOCK,bDB=TRUE,bAll=TRUE)
tam.block
dba.plotVenn(tam.block,1:3,label1="edgeR",label2="DESeq",label3="DESeq2")


###################################################
### code chunk number 32: DiffBind.Rnw:425-426
###################################################
data(tamoxifen_peaks)


###################################################
### code chunk number 33: DiffBind.Rnw:442-444
###################################################
olap.rate = dba.overlap(tamoxifen,mode=DBA_OLAP_RATE)
olap.rate


###################################################
### code chunk number 34: tamox_rate
###################################################
plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets') 


###################################################
### code chunk number 35: DiffBind.Rnw:465-466
###################################################
names(tamoxifen$masks)


###################################################
### code chunk number 36: DiffBind.Rnw:471-473
###################################################
dba.overlap(tamoxifen,tamoxifen$masks$MCF7 & tamoxifen$masks$Responsive,
            mode=DBA_OLAP_RATE)


###################################################
### code chunk number 37: tamox_mcf7_venn
###################################################
dba.plotVenn(tamoxifen, tamoxifen$masks$MCF7 & tamoxifen$masks$Responsive)


###################################################
### code chunk number 38: DiffBind.Rnw:495-498
###################################################
tamoxifen = dba.peakset(tamoxifen, consensus = c(DBA_TISSUE,DBA_CONDITION), 
                        minOverlap=0.66)
tamoxifen


###################################################
### code chunk number 39: DiffBind.Rnw:506-508
###################################################
tamoxifen_consensus = dba(tamoxifen, mask = tamoxifen$masks$Consensus)
tamoxifen_consensus


###################################################
### code chunk number 40: tamox_lines_venn
###################################################
data(tamoxifen_peaks)
tamoxifen = dba.peakset(tamoxifen, consensus = DBA_TISSUE, minOverlap=0.66)
dba.plotVenn(tamoxifen, tamoxifen$masks$Consensus)


###################################################
### code chunk number 41: DiffBind.Rnw:530-531
###################################################
data(tamoxifen_peaks)


###################################################
### code chunk number 42: DiffBind.Rnw:536-538
###################################################
dba.overlap(tamoxifen,tamoxifen$masks$Resistant,mode=DBA_OLAP_RATE)
dba.overlap(tamoxifen,tamoxifen$masks$Responsive,mode=DBA_OLAP_RATE)


###################################################
### code chunk number 43: tamox_cons_venn
###################################################
tamoxifen = dba.peakset(tamoxifen, consensus =  DBA_CONDITION, minOverlap = 0.33)
dba.plotVenn(tamoxifen,tamoxifen$masks$Consensus)


###################################################
### code chunk number 44: DiffBind.Rnw:567-568
###################################################
tamoxifen.OL = dba.overlap(tamoxifen, tamoxifen$masks$Consensus)


###################################################
### code chunk number 45: DiffBind.Rnw:573-575
###################################################
tamoxifen.OL$onlyA
tamoxifen.OL$onlyB


###################################################
### code chunk number 46: tamox_compare_venn
###################################################
tamoxifen = dba.peakset(tamoxifen,tamoxifen$masks$Consensus, 
                        minOverlap=1,sampID="OL Consensus")
tamoxifen = dba.peakset(tamoxifen,!tamoxifen$masks$Consensus,
                        minOverlap=3,sampID="Consensus_3") 
dba.plotVenn(tamoxifen,14:15)


###################################################
### code chunk number 47: DiffBind.Rnw:602-603
###################################################
data(tamoxifen_analysis)


###################################################
### code chunk number 48: DiffBind.Rnw:608-609
###################################################
tamoxifen.rep = dba.report(tamoxifen,bCalled=T,th=1)


###################################################
### code chunk number 49: DiffBind.Rnw:614-620
###################################################
onlyResistant = tamoxifen.rep$Called1>=2 & tamoxifen.rep$Called2<3
sum(onlyResistant )
onlyResponsive = tamoxifen.rep$Called2>=3 &  tamoxifen.rep$Called1<2
sum(onlyResponsive)
bothGroups = tamoxifen.rep$Called1>= 2 & tamoxifen.rep$Called2>=3
sum(bothGroups)


###################################################
### code chunk number 50: DiffBind.Rnw:632-639
###################################################
tamoxifen.DB = dba.report(tamoxifen,bCalled=T,th=.1)
onlyResistant.DB = tamoxifen.DB$Called1>=2 & tamoxifen.DB$Called2<3
sum(onlyResistant.DB)
onlyResponsive.DB = tamoxifen.DB$Called2>=3 & tamoxifen.DB$Called1<2
sum(onlyResponsive.DB)
bothGroups.DB = tamoxifen.DB$Called1>=2 & tamoxifen.DB$Called2>=3
sum(bothGroups.DB)


###################################################
### code chunk number 51: DiffBind.Rnw:780-781 (eval = FALSE)
###################################################
## file.path(system.file("extra", package="DiffBind"),"tamoxifen_GEO.csv")


###################################################
### code chunk number 52: sessionInfo
###################################################
toLatex(sessionInfo())


###################################################
### code chunk number 53: DiffBind.Rnw:811-812
###################################################
setwd(savewd)
