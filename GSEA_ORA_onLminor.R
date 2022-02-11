#######################################################
###    GSEA & ORA for Lemna minor DESeq2 results    ###
#######################################################

HOME = getwd()
setwd(HOME)

# go2gene database for GSEA / ORA via clusterProfiler
load(url("https://zenodo.org/record/6045874/files/Lminor.GO2gen4clusterProfiler.RData?download=1"))


### PACKAGES ### 
#require(stats)
require(clusterProfiler)
require(ggplot2)
require(AnnotationDbi)
require(org.Lminor.eg.db) #this needs to manually installed first!!!
# Follow instructions from here to install the org.Lminor.eg.db package:
# https://zenodo.org/record/6045874 


## Input for GSEA & ORA  ---------------------------------------------------
INPUT <- list()
for(i in list.files(pattern = "_reslfs_.+[.]csv", recursive = T)[1:4]){
  message("Reading in:\t",i)
  s = sub(".csv","",sub("_reslfs_",".",gsub("^.+/","",i)))
  tmp = read.csv2(i, header = T)
  colnames(tmp)[1] <- "geneID"
  INPUT[[s]] = tmp[order(-tmp$log2FoldChange),] # Sort by decreasing LFC
}

# Reduce all dfs to the common set
all = unique(Reduce(union, lapply(INPUT, function(x) x$geneID))) # 15435
univ = unique(Reduce(intersect, lapply(INPUT, function(x) x$geneID))) # 15278


## get DEGs in separate list 
DEG = list()
for(i in list.files(pattern = "_lfs_pcut_LFcut_.+[.]csv", recursive = T)[1:4]){
  message("Reading in:\t",i)
  s = sub(".csv","",sub("_lfs_pcut_LFcut_",".",gsub("^.+/","",i)))
  tmp = read.csv2(i, header = T)
  colnames(tmp)[1] <- "geneID"
  DEG[[s]] = tmp[order(-tmp$log2FoldChange),] # Sort by decreasing LFC
}
rm(i,s,tmp)

# subset DEG to the common univ gene set
sampleOrd = c(2,1,4,3) # set custom sample order here!
DEG = lapply(DEG, function(x) x[x$geneID %in% univ,] )[sampleOrd]
COMMON = lapply(INPUT, function(x) x[x$geneID %in% univ,] )[sampleOrd]

# Get the right format for GSEA - and ORA
GSE = lapply(COMMON, function(x) na.omit(setNames(x$log2FoldChange, as.character(x$geneID))) )
ORA = lapply(DEG, function(x) na.omit(unique(x$geneID)) )


dir.create("dataQC") # output for dataQC plots
## minor data check ## -----------------------------
df = as.data.frame(unlist(lapply(INPUT, function(x) length(x$geneID))) )
x = setNames(round((length(univ) / df[,1]) * 100,2), row.names(df))
barplot(x[c(1,3)], las = 1, ylab = "% remaining ProtIDs for common set", main = "ENSEMBL IDs", ylim = c(0,100))
# Range: 99.56 99.42 %
seed = 42

pdf("dataQC/geneIDcoverage_venn.pdf", width = 10, height = 7)
plot(eulerr::euler(list(Atorvastatin = INPUT$Atorvastatin.High$geneID,
                        Bentazon = INPUT$Bentazon.High$geneID)),
     fills = list(fill = c("lightblue","green4"), alpha = 0.7),
     labels = list(col = "black", font = 4),
     legend = list(col = "black", font = 4),
     main = "Shared Lemna queryIDs",
     quantities = TRUE, shape = "ellipse", lty = 0)
dev.off()


## Check how many GO annotated Lemna IDs we can cover with this 
all = unique(Reduce(union, lapply(TERM2GENE, function(x) x$gene)))
# Nice! In total we have a set of 18094 Lemna genes assigned to at least one GO term
# Let's see how this looks like for each ontology term
ls = lapply(TERM2GENE, function(x) unique(x$gene))
ls[["univ"]]  <- univ
ls[["all.ont"]] <- all

# venn plot 
require(RColorBrewer)
#display.brewer.all(n = NULL, type = "all", select = NULL, colorblindFriendly = T)

# venn type 1 ------------------------
pdf("dataQC/GOcoverage_venn.pdf", width = 10, height = 7)
venn::venn(ls[-5], ilab=TRUE, zcolor =  brewer.pal(n = length(ls), name = "Dark2"), lty = 3, 
           ilcs = 1, sncs = 1, box = F)
venn::venn(ls, ilab=TRUE, zcolor =  brewer.pal(n = length(ls), name = "Dark2"), lty = 3, 
           ilcs = 1, sncs = 1, box = F)
dev.off()

# venn type 2 -----------------------
gg = list()
gg[["v"]] = plot(eulerr::euler(ls[-5]),
     fills = list(fill = brewer.pal(n = length(ls), name = "Dark2"), alpha = 0.4),
     labels = list(col = "black", font = 4),
     legend = list(col = "black", font = 4),
     main = "Shared Lemna queryIDs",
     quantities = TRUE, shape = "ellipse", lty = 3)
# Bar plot with % coverage 
df = as.data.frame(unlist(lapply(ls, length))[c(1,3,2,4,5)])
x = as.data.frame(lapply(ls, function(x) (length(intersect(univ,x))/length(univ))*100))
df = merge(df, t(x), by=0)
colnames(df) = c("type","count","perc")
df$type <- as.factor(df$type)
gg[["bar"]] = ggplot(df, aes(x=type, y=count)) + geom_bar(stat="identity", alpha=.9, width = .8) +
  labs(subtitle = "Number of Lemna gene IDs for each group", x = "Group", y = 'Counts') +
  geom_text(aes(x=type, y=max(count)*1.05),
            #label = round(df$perc,1),
            label= paste0(df$count,"\n(",round(df$perc,1)," % shared)"),
            angle = "0") + theme_light()

pdf("dataQC/GOcoverage_bar.venn.pdf", width = 10, height = 10)
ggpubr::ggarrange(plotlist = gg, ncol = 1, nrow = 2)
dev.off()
# Looks good :)
######################


### prepare GO DATA for measuring semantic similarity ### ---------
GOdis = list()
for(ont in c("BP","CC","MF")){
  GOdis[[ont]] = GOSemSim::godata("org.Lminor.eg.db", ont=ont,
                                  computeIC=F, keytype = "GID")
}
#########################################################


### GSEA ### -------------------------------------------------------------------
gseGO.res = list() # Output
for(ont in c("BP","CC","MF")){
  message("\nRunning GSEA for GO terms associated with:\t", ont)
  for(i in names(GSE)){
    message("\nAnalyzing ",i)
    res = clusterProfiler::GSEA(GSE[[i]], TERM2GENE=TERM2GENE[[ont]],
                                TERM2NAME=TERM2NAME[[ont]],
                                pvalueCutoff = 1, pAdjustMethod = "BH", seed = T)#@result
    res@setType  <- ont
    res@organism <- "Lemna minor"
    res@keytype  <- "GID"
    
    # Better to use the org.Lminor.db ... somehow this is not working ... 
    #res = gseGO(GSE[[i]], ont = ont, OrgDb = "org.Lminor.eg.db", keyType = "GID", 
    #            pvalueCutoff = 1, pAdjustMethod = "BH", seed = T)
    # Error in (function (classes, fdef, mtable)  : 
    # unable to find an inherited method for function 'species' for signature '"character"'
    
    gseGO.res[[paste0(ont,"_",i)]] = res
    message("Done!\n")
  }
}
rm(ont,i,res)
gc(verbose = F)


## Simplify with custom distance matrix object ## ---------------
# See: https://github.com/YuLab-SMU/enrichplot/issues/79 

# skip this ----------------
x = gseGO.res$BP_Bentazon.High
# rmv everything which is "response to" GO term
df= x@result
df= df[grep("response to",df$Description, invert = T),]
df= df[grep("light",df$Description, invert = F),]
x@result = df

x@result = subset(x@result, p.adjust <= .1)
#d = godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)
#x = enrichplot::pairwise_termsim(x, method = "Wang", semData = d)
x = enrichplot::pairwise_termsim(x, method = "JC")
emapplot(x, showCategory = 10)
cnetplot(x, showCategory = 30, node_label = "category")

x2 = simplify(x, cutoff = 0.7, measure = "Wang")
x2 = enrichplot::pairwise_termsim(x2, method = "JC")
emapplot(x2, showCategory = 10)
cnetplot(x2, showCategory = 30, node_label = "category")

x = gseGO.res$BP_Atorvastatin.Low
x@result = subset(x@result, p.adjust <= 0.05)
emapplot(gseGO.res$BP_Atorvastatin.High, showCategory = 10)
###########
## Filter out non-sign results (pvalue <= 0.05) --------------------
pcut = .1
gseGO.cl = lapply(gseGO.res, function(x){
  res = x
  res@result <- subset(x@result, p.adjust <= pcut)
  return(res)
})

# Simplify each of the filtered result tables
gseGO = list()
for(ont in c("BP","CC","MF")){
  message("Simplifying and merging GSEA results for:\t",ont)
  tmp = gseGO.cl[grep(pattern = paste0("^",ont,"_"), names(gseGO.res))]
  names(tmp) = sub("^.._","",names(tmp))
  # simplify results in tmp list object
  n = .7 #simCut
  tmp = lapply(tmp, function(x){
    if(nrow(x) > 0){
      message("Computing pairwise_termsim ...")
      X = enrichplot::pairwise_termsim(x, method = "Wang", semData = GOdis[[ont]])
      message("Filtering redudant ",ont,"-GO terms (sim.cutoff: ",n,") ...")
      res = simplify(X, cutoff = n, measure = "Wang")
      message("Done!\n")
      return(res)
    } else {
      return(NA)
      message("No enriched terms provided!")
    }
  })
  gseGO[[ont]] = tmp
}
gseGO = lapply(gseGO, function(x){ x[!is.na(x)] }) # rmv empty df!!!
rm(tmp,ont,gseGO.cl)

## skip this ---------------------------------------------
# merge together the simplified GO results 
gseGO.cl = list() #output for df
ggGO = list()     #output for plots
pcut2 = .05
top = 15
filt = "p.adjust"

for(ont in c("BP","CC","MF")){
  tmp = gseGO.res[grep(pattern = paste0("^",ont,"_"), names(gseGO.res))]
  names(tmp) = sub("^.._","",names(tmp))
  gseGO.cl[[ont]] = mergeGSEAres(tmp, filter = filt, pcut = pcut2)
  ggGO[[ont]] = multiGSEAplot(tmp, top = top, title=paste("GO tems:",ont), filter = filt, pcut = pcut2)
}

dir.create("GSEA", showWarnings = F) #store output there
for(i in names(gseGO.cl)){
  write.csv2(gseGO.cl[[i]],file = paste0("GSEA/GSEAres_",i,"_pcut",pcut2,".csv"), row.names = F)
}
pdf(paste0("GSEA/GSEAres_pcut",pcut2,"top",top,".pdf"), width = 10, height = 10)
ggGO$BP
ggGO$MF
ggGO$CC
dev.off()
rm(ont, i, tmp, s, x, ls, df, gg)
############
## Merge individual GSEA results to long df 2.0 ## ----------------------------------
# I am not too hapy with the above function. Let's try again simpler ... 
mergeGOgseaRes = function(ls){
  df = rlist::list.rbind(ls)
  df$Count = apply(df,1,function(x){ length(stringr::str_split(x['core_enrichment'], pattern = "/")[[1]]) })
  df$perc.enriched = round(( df$Count / df$setSize )*100, 4)
  df$Cluster = sub("[.]GO:[0-9].+$","",row.names(df))
  row.names(df) = NULL
  return(df)
}

# Simplified
gseResS = lapply(gseGO, function(x){
  ls = lapply(x, function(x2){ x2@result })
  mergeGOgseaRes(ls)
})
# All - non simplified!
gseRes = list()
for(ont in c("BP","CC","MF")){
  tmp = gseGO.res[grep(pattern = paste0("^",ont,"_"), names(gseGO.res))]
  names(tmp) = sub("^.._","",names(tmp))
  tmp2 = lapply(tmp, function(x){ x@result })
  gseRes[[ont]] = mergeGOgseaRes(tmp2)
}
rm(tmp, tmp2, ont)

## Export results ## -----------------------------------
dir.create("GSEA", showWarnings = F) #store output there
# ALL
for(i in names(gseRes)){
  message("Exporting GSEA results for:\t",i)
  write.csv2(gseRes[[i]],file = paste0("GSEA/GSEAres_",i,"_all.csv"), row.names = F)
}
# Simplified and filtered
for(i in names(gseResS)){
  message("Exporting GSEA results for:\t",i)
  write.csv2(gseResS[[i]],file = paste0("GSEA/GSEAres.Simpl_",i,"_pcut",pcut,".csv"), row.names = F)
}


## Plot results ## ------------------------------------
dotPlot.GSE = function(res, top = 10, filter = "p.adjust", topBy="pvalue", pcut = .05, title = ""){
  if(filter != "p.adjust" & filter != "pvalue" & topBy != "qvalues"){
    stop("'filter' must be specified as 'p.adjust', 'qvalues' or 'pvalue'!") }
  if(topBy != "p.adjust" & topBy != "pvalue"  & topBy != "qvalues"){
    stop("'filter' must be specified as 'p.adjust', 'qvalues' or 'pvalue'!") }
  
  # subset df based on filter features 
  df = subset(res, res[,filter] <= pcut)
  
  # subset df back into cluster groups and pick the top N GO terms from each Cluster
  ls = list()
  for(i in unique(df$Cluster)){
    x = df[df$Cluster %in% i,]
    # select topN from x
    if(nrow(x) > 0){
      x = x[order(x[,topBy]),]
      if(top > nrow(x)) {
        x = x[1:nrow(x),]
      } else {
        x = x[1:top,]
      }
    }
    ls[[i]] = x
  }
  df1 = rlist::list.rbind(ls)
  df1$Cluster = factor(df1$Cluster, levels = unique(df1$Cluster))
  
  ## PLOT ##
  #g = ggplot(df1, aes(x=perc.enriched, y=Description, color=NES, size= -log10(df1[,filter]) )) + geom_point() +
  #  scale_colour_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0) +
  #  labs(size= paste0("-log10(",filter,")")) +
  #  facet_grid(~Cluster) + theme_light() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  #  aes(x=perc.enriched, reorder(stringr::str_wrap(Description, 80), -p.adjust)) + ylab(NULL) +
  #  labs(title = paste("GSEA -",title,":",length(unique(df1$ID)),"terms"),
  #       subtitle = paste(filter,"<",pcut,"/ Top:",top,"/ Top select by:",topBy))
  
  g = ggplot(df1, aes(x=-log10(df1[,filter]), y=Description, color=NES, size=perc.enriched )) + geom_point() +
    scale_colour_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0) +
    facet_grid(~Cluster) + theme_light() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    aes(x=-log10(df1[,filter]), reorder(stringr::str_wrap(Description, 80), -df1[,filter])) +
    labs(size= "% enriched") + ylab(NULL) + xlab(paste0("-log10(",filter,")")) +
    labs(title = paste("GSEA -",title,":",length(unique(df1$ID)),"terms"),
         subtitle = paste(filter,"<",pcut,"/ Top:",top,"/ Top select by:",topBy))
  
  return(g)
}

n = 25
pdf(paste0("GSEA/GSEAresSimpl.LFC_top",n,".pdf"), width = 10, height = 10)
dotPlot.GSE(gseResS$BP, top = n, title = "BP")
dotPlot.GSE(gseResS$CC, top = n, title = "CC")
dotPlot.GSE(gseResS$MF, top = n, title = "MF")
dev.off()

pdf(paste0("GSEA/GSEAres.LFC_top",n,".pdf"), width = 10, height = 10)
dotPlot.GSE(gseRes$BP, top = n, title = "BP")
dotPlot.GSE(gseRes$CC, top = n, title = "CC")
dotPlot.GSE(gseRes$MF, top = n, title = "MF")
dev.off()

p.gse = .1
pdf(paste0("GSEA/GSEAresSimpl.LFC_TermSearch.pval",p.gse,".pdf"), width = 12, height = 10)
toi = c("steroid","pathway","light","lipid","metabolic process","catabolic process") #terms of interest 
for(n in toi){
  for(ont in c("BP","CC","MF")){
    message(n," - ",ont)
    x = gseResS[[ont]]
    x = x[grep(n, x$Description),]
    if(nrow(x)>0 & min(x$pvalue) <= p.gse){ 
      if(any(grep("bolic process",n))){
        print(dotPlot.GSE(x, filter = "p.adjust", title = paste(n,ont), 
                          pcut = 0.05, top = nrow(x)))
      } else {
        print(dotPlot.GSE(x, filter = "pvalue", title = paste(n,ont), 
                          pcut = p.gse, top = nrow(x)))
      }
    }
  }
}
dev.off()


### RANK METRIC based GSEA ### -------------------------------------------------
# Ok this is a little bit messy ... let's try to play around with integrating pval & lfc
# Computing a ranked score metric based on lfc & p.adj using the following equation:
#n = log10((abs(lfc)+1) / padj)
#if(lfc < 0){n*-1}else{n}
# The idea is to keep the initial up or down regulation info by considering the input lfc.
# if lfc < 0 multiply n by -1 to return a negativ score.
# We are adding +1 to lfc to avoid values < 1 to be log transformed. This would result in
# negative n which compromises the initial idea to keep up/down regulation info in the score.
# It might be possible that p values need to be -log transformed as well if numbers tend to
# spread to far appart. But let's see for now.

df = COMMON$Atorvastatin.Low
my_rankMetric = function(df, lfc = "log2FoldChange", pval = "padj",
                         penalize = F, penaltyLfc = .1, penaltyFactor = 10){
  stopifnot(length(df[,lfc])==length(df[,pval]))
  stopifnot( any(colnames(df) == lfc) & any(colnames(df) == pval) )
  # Compute the rank score
  apply(df, 1, function(x){
    L = as.numeric(x[grep(lfc, colnames(df))])
    P = as.numeric(x[grep(pval, colnames(df))])
    n = log10( (abs(L)+1)/P )
    if(L < 0){S = n*(-1)}else{S = n}
    # penalize low lfc values! 
    if(penalize == T){
      if( abs(L) <= penaltyLfc ){
        return( S/penaltyFactor )
      }else{return( S )}
    }else{
      return( S )
    }
  })
}

df$rankScore = my_rankMetric(df)
df$rankScoreP = my_rankMetric(df, penalize = T)
h = 0.1 # => It would be a good idea to penalize very low lfc values as some of them tend to have an inflated rankScore

my_rankScorePlots = function(df, title = "", penaltyLfc = .1) {
  h = penaltyLfc
  par(mfrow = c(3,3))
  
  hist(df$log2FoldChange, breaks =100, ylim= c(0,75), xlab = "log2FC", main = title)
  abline(v = c(-h,h), lty = 3, col = "red2")
  text(h, 60, paste("+/-",h), pos = 4, col = "red2")
  
  hist(df$rankScore, breaks = 100, ylim = c(0,75), xlab = "Rank score", main = title)
  hist(df$rankScoreP, breaks = 100, ylim = c(0,75), xlab = "Rank score - penalized", main = title)
  
  # Check for correlation
  plot(abs(df$log2FoldChange), -log10(df$padj))
  abline(v = h, lty = 3, col = "red2")
  text(h, max(-log10(df$padj))*.9, h, pos = 4, col = "red2")
  
  plot(abs(df$rankScore ), -log10(df$padj))
  plot(abs(df$rankScoreP ), -log10(df$padj), main = "Penalized rank scores")
  
  plot(abs(df$log2FoldChange), -log10(df$padj))
  pcc = round(cor(abs(df$log2FoldChange), -log10(df$padj), method = "pearson"),4)
  text(0, max(-log10(df$padj))*.9, pos=4, labels = paste("PCC:",pcc))
  
  plot(df$rankScore, df$log2FoldChange)
  abline(h = c(-h,h), lty = 3, col = "red2")
  pcc = round(cor(df$log2FoldChange, df$rankScore, method = "pearson"),4)
  text(min(df$rankScore), max(df$log2FoldChange)*.9, pos=4, labels = paste("PCC:",pcc))
  text(min(df$rankScore)*.9, h, paste("+/-",h), pos = 3, col = "red2")
  
  plot(df$rankScoreP, df$log2FoldChange, main = "Penalized rank scores")
  abline(h = c(-h,h), lty = 3, col = "red2")
  pcc = round(cor(df$log2FoldChange, df$rankScoreP, method = "pearson"),4)
  text(min(df$rankScoreP), max(df$log2FoldChange)*.9, pos=4, labels = paste("PCC:",pcc))
  text(min(df$rankScoreP)*.9, h, paste("+/-",h), pos = 3, col = "red2")
}

# Compute the rank score for all datasets
COMMON = lapply(COMMON, function(x){
  df = x
  df$rankScore = my_rankMetric(df)
  df$rankScoreP = my_rankMetric(df, penalize = T)
  return(df)
})

for(i in names(COMMON)){ 
  png(paste0("GSEA/RankScores_",i,".png"), units = "mm", width = 225, height = 250, res = 350)
  print( my_rankScorePlots(COMMON[[i]], title = i))
  dev.off()
  }

# Now build input list for GSEA based on rank score
GSE.rankS = lapply(COMMON, function(x) {
  X = na.omit(setNames(x$rankScoreP, as.character(x$geneID)))
  X = sort(X, decreasing = T)
  })

# Run GSEA with rank scores
gse = list() # Output
for(ont in c("BP","CC","MF")){
  message("\nRunning GSEA for GO terms associated with:\t", ont)
  for(i in names(GSE.rankS)){
    message("\nAnalyzing ",i)
    res = clusterProfiler::GSEA(GSE[[i]], TERM2GENE=TERM2GENE[[ont]],
                                TERM2NAME=TERM2NAME[[ont]],
                                pvalueCutoff = .1, pAdjustMethod = "BH", seed = T)#@result
    res@setType  <- ont
    res@organism <- "Lemna minor"
    res@keytype  <- "GID"
    gse[[paste0(ont,"_",i)]] = res
    message("Done!\n")
  }
}

# Simplify each of the filtered result tables
gseS = list()
for(ont in c("BP","CC","MF")){
  message("Simplifying and merging GSEA results for:\t",ont)
  tmp = gse[grep(pattern = paste0("^",ont,"_"), names(gse))]
  names(tmp) = sub("^.._","",names(tmp))
  # simplify results in tmp list object
  n = .7 #simCut
  tmp = lapply(tmp, function(x){
    if(nrow(x) > 0){
      message("Computing pairwise_termsim ...")
      X = enrichplot::pairwise_termsim(x, method = "Wang", semData = GOdis[[ont]])
      message("Filtering redudant ",ont,"-GO terms (sim.cutoff: ",n,") ...")
      res = simplify(X, cutoff = n, measure = "Wang")
      message("Done!\n")
      return(res)
    } else {
      return(NA)
      message("No enriched terms provided!")
    }
  })
  gseS[[ont]] = tmp
}
gseS = lapply(gseS, function(x){ x[!is.na(x)] }) # rmv empty df!!!
rm(tmp,ont)

# merge gse results by ONT
gseRes2 = lapply(gseS, function(x){
  ls = lapply(x, function(x2){ x2@result })
  mergeGOgseaRes(ls)
})

# Export Simplified and filtered on rankS
for(i in names(gseRes2)){
  write.csv2(gseRes2[[i]],file = paste0("GSEA/GSEAresRankS.Simpl_",i,"_pcut",pcut,".csv"), row.names = F)
}

n = 25
pdf(paste0("GSEA/GSEAresSimpl.rankS_top",n,".pdf"), width = 10, height = 10)
dotPlot.GSE(gseRes2$BP, top = n, title = "BP")
dotPlot.GSE(gseRes2$CC, top = n, title = "CC")
dotPlot.GSE(gseRes2$MF, top = n, title = "MF")
dev.off()

p.gse = .1
pdf(paste0("GSEA/GSEAresSimpl.rankS_TermSearch.pval",p.gse,".pdf"), width = 12, height = 10)
toi = c("steroid","pathway","light","lipid","metabolic process","catabolic process") #terms of interest 
for(n in toi){
  for(ont in c("BP","CC","MF")){
    message(n," - ",ont)
    x = gseRes2[[ont]]
    x = x[grep(n, x$Description),]
    if(nrow(x)>0 & min(x$pvalue) <= p.gse){ 
      if(any(grep("bolic process",n))){
        print(dotPlot.GSE(x, filter = "p.adjust", title = paste(n,ont), 
                          pcut = 0.05, top = nrow(x)))
      } else {
        print(dotPlot.GSE(x, filter = "pvalue", title = paste(n,ont), 
                          pcut = p.gse, top = nrow(x)))
      }
    }
  }
}
dev.off()

##############################


### ORA ### --------------------------------------------------------------------
#?enricher()
### ORA plot function ###
dotPlot.ORA = function(res, top = 10, filter = "p.adjust", topBy="pvalue", pcut = .05, title = ""){
  if(filter != "p.adjust" & filter != "pvalue"){ stop("'filter' must be specified as 'p.adjust' or 'pvalue'!") }
  if(topBy != "p.adjust" & topBy != "pvalue"){ stop("'filter' must be specified as 'p.adjust' or 'pvalue'!") }
  
  # subset df based on filter features 
  df = subset(res, res[,filter] <= pcut)
  
  # subset df back into cluster groups and pick the top N GO terms from each Cluster
  ls = list()
  for(i in unique(df$Cluster)){
    x = df[df$Cluster %in% i,]
    # select topN from x
    if(nrow(x) > 0){
      x = x[order(x[,topBy]),]
      if(top > nrow(x)) {
        x = x[1:nrow(x),]
      } else {
        x = x[1:top,]
      }
    }
    ls[[i]] = x
  }
  df1 = rlist::list.rbind(ls)
  
  # compute "Bg ratio"
  z = as.integer(sub("/.+$","",df1$BgRatio))
  n = as.integer(sub("^.+/","",df1$BgRatio))
  stopifnot(length(z) == length(n))
  df1$bgR = z/n
  # compute "gene ratio"
  z = df1$Count
  N = as.integer(sub("^.+/","",df1$GeneRatio))
  stopifnot(length(z) == length(N))
  df1$geneR = z/N
  # get Nbr of all DEGs identified for a GO term set in a cluster
  df1$Cluster = paste0(df1$Cluster," (",N,")")
  df1$Cluster = factor(df1$Cluster, levels = unique(df1$Cluster)) #for correct input order; can be adjusted here
  row.names(df1) = NULL
  # compute geneR over bg-ratio
  df1$genR.bgR <- df1$geneR / df1$bgR
  
  #col = viridis::viridis(13)[c(1,7,13)]
  #col = viridis::plasma(13)[c(1,7,13)]
  #g = ggplot(df1, aes(x= geneR, y= Description, color=p.adjust, size= log2(Count) )) + geom_point() +
    #scale_colour_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0) +
    #scale_colour_gradient2(low = "firebrick2", mid = "yellow", high = "#0D0887FF", midpoint = .05) +
    ###scale_colour_gradient2(low = col[1], mid = col[2], high = col[3]) +
    ###scale_colour_gradient(low = "firebrick2", high = "yellow") +
    #facet_grid(~Cluster) + theme_light() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    ###aes(geneR, reorder(stringr::str_wrap(Description, 80), geneR)) +
    ###aes(Count, reorder(stringr::str_wrap(Description, 80), geneR)) +
    #aes(geneR, reorder(stringr::str_wrap(Description, 80), Count)) + ylab(NULL) + xlab("Gene ratio")# +
    #labs(title = paste("ORA -",title,":",length(unique(df1$ID)),"terms"), subtitle = paste(filter,"<",pcut,"/ Top:",top,"/ Top select by:",topBy))
  # aes(genR.bgR, ...)
  g = ggplot(df1, aes(x= geneR, y= Description, color=df1[,filter], size= log2(Count) )) + geom_point() +
    scale_colour_gradient2(low = "firebrick2", mid = "yellow", high = "#0D0887FF", midpoint = .05) +
    #labs(col=if(filter=="p.adjust"){"p.adjust"}else{"pvalue"}) +
    labs(col= paste0(filter)) +
    facet_grid(~Cluster) + theme_light() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    aes(geneR, reorder(stringr::str_wrap(Description, 80), -p.adjust)) + ylab(NULL) +
    labs(title = paste("ORA -",title,":",length(unique(df1$ID)),"terms"),
         subtitle = paste(filter,"<",pcut,"/ Top:",top,"/ Top select by:",topBy))
  return(g)
}

## ORA for each treatment (ALL DEGs) ------------------------------------------------------
oraGO.res = list() # output
for(ont in c("BP","CC","MF")){
  message("\nRunning ORA for GO terms associated with:\t", ont)
  for(i in names(ORA)){
    message("\nAnalyzing ",i)
    oraGO.res[[paste0(ont,"_",i)]] = clusterProfiler::enricher(
      ORA[[i]], pvalueCutoff = 1, pAdjustMethod = "BH", universe = univ,
      TERM2GENE= TERM2GENE[[ont]], TERM2NAME = TERM2NAME[[ont]])@result
    message("Done!\n")
  }
}
# Let's filter out non.sign
#p.ora = 1 #remove everything with a p.adjust <= p.ora
#tmpORA = lapply(oraGO.res, function(x) subset(x, p.adjust <= p.ora))

# merge ORA results into single df
oraRes = list()
for(ont in c("BP","CC","MF")){
  tmp = oraGO.res[grep(pattern = paste0("^",ont,"_"), names(oraGO.res))]
  names(tmp) = sub("^.._","",names(tmp))
  # merge
  df = rlist::list.rbind(tmp)
  df$Cluster = sub("[.]GO:[0-9].+$","",row.names(df))
  row.names(df) = NULL
  oraRes[[ont]] <- df
}

# Export ORA#1 results
dir.create("ORA", showWarnings = F)
for(i in names(oraRes)){
  write.csv2(oraRes[[i]],file = paste0("ORA/ORAres_allDEGs_",i,".csv"), row.names = F)
}
# Plot results 
p.ora = 0.1
n = 20
pdf(paste0("ORA/ORAres_allDEGs_top",n,"pcut",p.ora,".pdf"), width = 12, height = n/2.5)
dotPlot.ORA(oraRes$BP, title = "BP", top=n, pcut = p.ora)
dotPlot.ORA(oraRes$CC, title = "CC", top=n, pcut = p.ora)
dotPlot.ORA(oraRes$MF, title = "MF", top=n, pcut = p.ora)
dev.off()


# Subset results for Lipid an light related GO terms and plot that
pdf(paste0("ORA/ORAres_allDEGs_TermSearch.pval",p.ora,".pdf"), width = 12, height = 10)
toi = c("steroid","pathway","light","lipid","metabolic process","catabolic process") #terms of interest 
for(n in toi){
  for(ont in c("BP","CC","MF")){
    message(n," - ",ont)
    x = oraRes[[ont]]
    x = x[grep(n, x$Description),]
    if(nrow(x)>0 & min(x$pvalue) <= p.ora){ 
      if(any(grep("bolic process",n))){
        print(dotPlot.ORA(x, filter = "p.adjust", title = paste(n,ont), 
                          pcut = .1, top = nrow(x)))
      } else {
        print(dotPlot.ORA(x, filter = "pvalue", title = paste(n,ont), 
                          pcut = p.ora, top = nrow(x)))
      }
    }
  }
}
dev.off()



## ORA on common set of DEGs for each substance (core DEGs) --------------------------------
ORA2 = list()
for(s in unique(sub("[.].+$","",names(ORA))) ){
  ORA2[[s]] = Reduce(intersect, lapply(ORA[grep(s, names(ORA))], function(x){x}) )
}

oraGO.res2 = list() # output
p.ora = 0.1
for(ont in c("BP","CC","MF")){
  message("\nRunning ORA for GO terms associated with:\t", ont)
  for(i in names(ORA2)){
    message("\nAnalyzing ",i)
    
    #oraGO.res2[[paste0(ont,"_",i)]] = clusterProfiler::enricher(
    #  ORA2[[i]], pvalueCutoff = 1, pAdjustMethod = "BH", universe = univ,
    #  TERM2GENE= TERM2GENE[[ont]], TERM2NAME = TERM2NAME[[ont]])@result
    res = enrichGO(ORA2[[i]], OrgDb = "org.Lminor.eg.db", keyType = "GID",
                   pAdjustMethod = "none", universe = univ, ont = ont, pvalueCutoff = 1)
    oraGO.res2[[paste0(ont,"_",i)]] <- res@result
    
    #if(nrow(subset(res@result, p.adjust <= p.ora)) > 0){
    #  res = simplify(res, .6, measure = "Wang", by = "pvalue")
    #  oraGO.res2[[paste0(ont,"_",i)]] <- res@result
    #} else { message("No sign. enriched terms found in ",i," for ",ont) }
    message("Done!\n")
    rm(res)
  }
}

# merge ORA results into single df
oraRes2 = list()
for(ont in c("BP","CC","MF")){
  tmp = oraGO.res2[grep(pattern = paste0("^",ont,"_"), names(oraGO.res2))]
  names(tmp) = sub("^.._","",names(tmp))
  # merge
  df = rlist::list.rbind(tmp)
  df$Cluster = sub("[.]GO:[0-9].+$","",row.names(df))
  row.names(df) = NULL
  oraRes2[[ont]] <- df
}

# Export ORA#2 results
dir.create("ORA", showWarnings = F)
for(i in names(oraRes2)){
  write.csv2(oraRes2[[i]],file = paste0("ORA/ORAres_coreDEGs_",i,".csv"), row.names = F)
}

# Plot results 
n = 20
pdf(paste0("ORA/ORAres_coreDEGs_top",n,"pcut",p.ora,".pdf"), width = 10, height =  n/2.5)
dotPlot.ORA(oraRes2$BP, title = "BP", top=n, pcut = p.ora) #+aes(log2(Count), reorder(stringr::str_wrap(Description, 80), genR.bgR))
dotPlot.ORA(oraRes2$CC, title = "CC", top=n, pcut = p.ora)
dotPlot.ORA(oraRes2$MF, title = "MF", top=n, pcut = p.ora)
dev.off()


# Subset results for Lipid an light related GO terms and plot that
pdf(paste0("ORA/ORAres_coreDEGs_TermSearch.pval",p.ora,".pdf"), width = 9, height = 6)
toi = c("steroid","pathway","light","lipid",#"metabolic process","
        "catabolic process") #terms of interest 
for(n in toi){
  for(ont in c("BP","CC","MF")){
    message(n," - ",ont)
    x = oraRes2[[ont]]
    x = x[grep(n, x$Description),]
    if(nrow(x)>0 & min(x$pvalue) <= .1){ 
      if(any(grep("bolic process",n))){
        print(dotPlot.ORA(x, filter = "p.adjust", title = paste(n,ont), 
                          pcut = p.ora, top = nrow(x)))
      } else {
        print(dotPlot.ORA(x, filter = "pvalue", title = paste(n,ont), 
                          pcut = p.ora, top = nrow(x)))
      }
    }
  }
}
dev.off()



pdf(paste0("ORA/ORAres_coreDEGs_TermSearch.pval",p.ora,".metabolic.pdf"), width = 9, height = 6)
tmp = lapply(oraRes2, function(x){ subset(x, p.adjust < .02)}) # subset data first, otherwise too large
for(n in "metabolic process"){
  for(ont in c("BP")){
    message(n," - ",ont)
    x = tmp[[ont]]
    x = x[grep(n, x$Description),]
    if(nrow(x)>0 & min(x$pvalue) <= .1){ 
      if(any(grep("bolic process",n))){
        print(dotPlot.ORA(x, filter = "p.adjust", title = paste(n,ont), 
                          pcut = p.ora, top = nrow(x)))
      } else {
        print(dotPlot.ORA(x, filter = "pvalue", title = paste(n,ont), 
                          pcut = p.ora, top = nrow(x)))
      }
    }
  }
}
dev.off()


### custom term selection ###
# all from lipid and light terms (BP) with pvalue <= .1
df = oraRes2[["BP"]]

pdf(paste0("ORA/ORAres_coreDEGs_customTerms.pdf"), width = 8, height = 6)
dotPlot.ORA(df[grep("light|lipid|pathway", df$Description),], title = "lipid|pathway|light, p sorted", 
            filter="p.adjust", pcut = .05, top=99) +
  aes(geneR, reorder(stringr::str_wrap(Description, 80), -pvalue ))+ylab(NULL)

#dotPlot.ORA(df[grep("bolic", df$Description),], title = "cata- & metabolic terms, padj.sorted", 
 #           filter="p.adjust", pcut = .05, top=99)
dev.off()

customTerms = df$Description[grep("light|lipid", df$Description)]
############



### compareCluster - ORA ### ---------------------------------------------------
# Input: ORA & ORA2
oraClust = list()
for(ont in c("BP","CC","MF")){
  res = compareCluster(ORA, fun = "enrichGO", OrgDb = "org.Lminor.eg.db", 
                       keyType = "GID", ont = ont, universe = univ,
                       pvalueCutoff = .1, pAdjustMethod = "BH")
  res = enrichplot::pairwise_termsim(res, method = "Wang", semData = GOdis[[ont]])
  #res = simplify(res, .7, measure = "Wang")
  res@compareClusterResult <- res@compareClusterResult[order(res@compareClusterResult$pvalue),]
  oraClust[[ont]] = res
}

oraClust2 = list()
for(ont in c("BP","CC","MF")){
  res = compareCluster(ORA2, fun = "enrichGO", OrgDb = "org.Lminor.eg.db", 
                       keyType = "GID", ont = ont, universe = univ,
                       pvalueCutoff = .1, pAdjustMethod = "BH")
  res = enrichplot::pairwise_termsim(res, method = "Wang", semData = GOdis[[ont]])
  #res = simplify(res, .7, measure = "Wang")
  res@compareClusterResult <- res@compareClusterResult[order(res@compareClusterResult$pvalue),]
  oraClust2[[ont]] = res
}
rm(res,ont)

x = oraClust2
for(k in names(x)){
  pdf(paste0("ORA/networkPlot2.coreDEG_",k,".pdf"), width = 9, height = 7)
  print(emapplot(x[[k]], showCategory = 20))
  print(emapplot(x[[k]], showCategory = 30))
  #print(emapplot(x[[k]], showCategory = 40))
  print(cnetplot(x[[k]], showCategory = 20, node_label = "category"))
  print(cnetplot(x[[k]], showCategory = 10, layout = "kk"))
  dev.off()
}

## subset BP to the processes shwon in dotplot with custom terms (see line 787)
tmp = oraClust2$BP
df = oraClust2$BP@compareClusterResult

pdf(paste0("ORA/networkPlot2.coreDEG_BP_CustomSet.pdf"), width = 9, height = 9)
# custom set with "light|lipid|pathway"
tmp@compareClusterResult <- df[grep("light|lipid|pathway", df$Description),]
emapplot(tmp, showCategory = 99)
cnetplot(tmp, showCategory = 10, node_label = "all")
df.light_lip_pw = df[grep("light|lipid|pathway", df$Description),]


# everything but "response to"
tmp@compareClusterResult <- df[grep("response to", df$Description, invert = T),]
emapplot(tmp, showCategory = 10)
cnetplot(tmp, showCategory = 10, node_label = "all")
emapplot(tmp, showCategory = 20)
emapplot(tmp, showCategory = 30)
emapplot(tmp, showCategory = nrow(tmp@compareClusterResult))
dev.off()


## rm "response to" then add everything with "light|lipid|pathway" again.
df1 = df[grep("light|lipid|pathway", df$Description),]
df2 = df[grep("response to", df$Description, invert = T),]
df2 = rbind(df1,df2)
df2 = df2[!duplicated(df2$ID),]
df2 = df2[order(df2$pvalue),]
tmp@compareClusterResult <- rbind(df1,df2)

pdf(paste0("ORA/networkPlot2.coreDEG_BP_CustomSet2.pdf"), width = 10, height = 8)
emapplot(tmp, showCategory = 10)
cnetplot(tmp, showCategory = 10, node_label = "all")
emapplot(tmp, showCategory = 15)
cnetplot(tmp, showCategory = 15, node_label = "all")
emapplot(tmp, showCategory = 20)
emapplot(tmp, showCategory = 25)
cnetplot(tmp, showCategory = 25, node_label = "all")
emapplot(tmp, showCategory = nrow(df2))
dev.off()
dev.off()


# get all gene IDs from this list:
lemnaAnno = read.csv2("~/BLASTdb/Lemna.minor_ref.org.blast.db/uniprotFASTA/Lminor.blastp.eNOG.GOanno.csv")

getIDs = function(df){
  s = unlist(strsplit(paste(df$geneID, sep = "/"), split = "/"))
  return(unique(s))
}
lemna.idOI = getIDs(df2) #163

for(i in c("light","lipid","pathway")){
  x = lemnaAnno[lemnaAnno$query %in% getIDs(df2[grep(i, df2$Description),]), ]
  write.csv2(x, paste0("Lminor_GOI.based.on.clustProf_",i,".csv"), row.names = F)
}

getIDs(df2[grep("light", df2$Description),])
getIDs(df2[grep("lipid", df2$Description),])
getIDs(df2[grep("pathway", df2$Description),])

df = lemnaAnno[lemnaAnno$query %in% lemna.idOI,]
write.csv2(df, "Lminor_GOI.based.on.clusterProfiler.enriched.sets.csv", row.names = F)





x = oraClust
for(k in names(x)){
  pdf(paste0("ORA/networkPlot2.allDEG_",k,".pdf"), width = 9, height = 7)
  print(emapplot(x[[k]], showCategory = 18))
  print(cnetplot(x[[k]], showCategory = 7, node_label = "category"))
  print(cnetplot(x[[k]], showCategory = 7, layout = "kk"))
  dev.off()
}

### Export compare cluster Res
for(k in names(oraClust)){
  write.csv2(oraClust[[k]]@compareClusterResult, row.names = F,
             paste0("ORA/compClust.allDEG_",k,".csv"))
  write.csv2(oraClust2[[k]]@compareClusterResult, row.names = F,
             paste0("ORA/compClust.coreDEG_",k,".csv"))
}
############################