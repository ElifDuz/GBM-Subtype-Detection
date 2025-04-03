#                      
#             ****Subtype detection analysis for GBM datasets****
#                       ****from normalized data****

# Elif Duz
# 03.04.25

# libraries
library(dplyr)
library(limma)
library(readxl)
library(ggplot2)
library(GSVA)
library(pheatmap)
library(RColorBrewer)

gseid="GSE4271" #------> GSEID

# load the normalized unique data
load("GSE4271-normalizeduniquedata.RData")#-------> upload the normalized expression values

#you can use the metadata or arranged metadata
remove(metadata)
#-----------> upload the metadata 
metadata <- read_excel("C:/Users/elif.duz/Desktop/GBM data analysis/toplantisonrasimakale/toplantisonrasianaliz/celllines.xlsx", sheet = "GSE62153paired")

# always order the data based on metadata (Disease vs Control)
# in this analysis we used only paired samples.
# check the order of metadata (primary recurrent or vice versa)

paireddata= df_unique[, c(metadata[metadata$type2=="recurrent",]$Accession,metadata[metadata$type2=="primary",]$Accession )]
rownames(paireddata)= rownames(df_unique)

# Gene lists for subtypes
# in this analysis we have 3 gene lists for 3 GBM related subtypes

genlists <- read_excel("C:/Users/elif.duz/Desktop/GBM data analysis/toplantisonrasimakale/toplantisonrasianaliz/paireddataanalysis/githubGSE222515.xlsx")
proneuralgenes= na.omit(genlists$preneural)
mesenchymalgenes=na.omit(genlists$mesenchymal)
classicalgenes=na.omit(genlists$classical)
genelists=list(classicalgenes,mesenchymalgenes,proneuralgenes )
names(genelists)=c("Classical","Mesenchymal", "Proneural")

# Subtype detection with GSVA package

gsva_ssgsva <- ssgseaParam(as.matrix(paireddata),genelists) # perform the analysis for each samples seperately
output_ssgsea <- gsva(gsva_ssgsva) 
heatmap(output_ssgsea)


custom_colors <- colorRampPalette(c("white", "pink", "red"))(100)


png(paste0(gseid, "_subtypes.png"), width = 3000, height = 1500, res = 300)
pheatmap(output_ssgsea, 
         color = custom_colors, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = round(output_ssgsea, 2),  
         main = "ssGSEA Heatmap")
dev.off()

molec_subtype <- names(genelists)[apply(output_ssgsea,2,function(z){which.max(z)})]
names(molec_subtype) <- colnames(output_ssgsea)
molec_subtype=data.frame(sampID=names(molec_subtype),grup=molec_subtype )

#merge the metadata with subtype information
metadat= merge(metadata, molec_subtype, by.x="Accession", by.y= "sampID", all.x=TRUE)
write.table(metadat, file= paste0(gseid, "_subtypes.xls"))

# metadata have order (paired smples) and subtype information
# detect the samples that have same order and same subtypes

matched_samples <- metadat %>% 
  group_by(order, grup) %>%  
  filter(n() == 2)   

groupes= c("Proneural","Mesenchymal", "Classical")
groupes= intersect(matched_samples$grup, groupes) # to check if all the groupes have paired samples
for (i in 1: length(groupes)){
  filtered_res <- matched_samples %>% filter(type2 == "recurrent" & grup == groupes[i])
  filtered_sens <- matched_samples %>% filter(type2 == "primary" & grup == groupes[i])
  sdata= data.frame(df_unique[, filtered_res$Accession], df_unique[, filtered_sens$Accession]) # disease vs control
  rownames(sdata)= rownames(df_unique)
  meta= rbind(filtered_res, filtered_sens)
  meta$index <- ave(meta$order, meta$order, FUN = seq_along) 
  meta <- meta[order(meta$index, meta$order), ]  
  meta$index <- NULL 
  #---------------> you can check the specific genes 
  oxgenes <-read_excel("C:/Users/elif.duz/Desktop/GBM data analysis/toplantisonrasimakale/toplantisonrasianaliz/celllines.xlsx", sheet = "oxidativestressgenes")
  
  if (length(filtered_res$Accession)>2){ #-----------> if we have more than 3 samples we can perform limma paired DEG analysis
  samples <- data.frame(
    patient <- factor(c(1:length(filtered_res$Accession) ,1:length(filtered_res$Accession))),  
    group <- factor(meta$type2))
  design= model.matrix(~ patient + group, data=samples)
  # Lineer model fit
  fit <- lmFit(sdata, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef="grouprecurrent", adjust="BH", number=Inf) #-----> recurrent vs primary
  
  DEGSpaired= results[results$P.Value<=0.05,] #------------> set the p-value threshold
  results$FC= 2^(results$logFC)


  results$genes <- rownames(results) 
  df_merged <- merge(results, oxgenes, by = "genes", all.x = TRUE)  
  
  upgenes <- df_merged %>% filter(P.Value < 0.05 & FC>1.5)
  upgenes$dataset= paste0(gseid,"_",groupes[i],"_upregulated")
  downgenes <- df_merged %>% filter(P.Value < 0.05 & FC<0.666)
  downgenes$dataset=paste0(gseid,"_",groupes[i],"_downregulated")
  
  res=rbind(upgenes, downgenes)
  
  save(res, file= paste0(gseid,"_",groupes[i],"paired.RData"))
  write.table(res, file= paste0(gseid,"_",groupes[i],"paired.xls"))
  assign(paste0(gseid,"_",groupes[i]),res)
  }
  
  if (length(filtered_res$Accession)==2){ #-------> if we have only two samples for subtype
    
    res=data.frame(genes= rownames(sdata), FC=rowMeans(sdata[,1:2])/rowMeans(sdata[,3:4])) 
    df_merged <- merge(res, oxgenes, by = "genes", all.x = TRUE) 
    upgenes= df_merged %>% filter(FC>1.5)
    upgenes$dataset= paste0(gseid,"_",groupes[i],"_upregulated")
    downgenes= df_merged %>% filter(FC<0.6666)
    downgenes$dataset=paste0(gseid,"_",groupes[i],"_downregulated")
    res=rbind(upgenes, downgenes)
    save(res, file= paste0(gseid,"_",groupes[i],"paired.RData"))
    write.table(res, file= paste0(gseid,"_",groupes[i],"paired.RData"))
    assign(paste0(gseid,"_",groupes[i]),res)}
  
  if (length(filtered_res$Accession)==1){ #-------> if we have only one samples for subtype
    
    res=data.frame(genes= rownames(sdata), FC=sdata[,1]/sdata[,2]) 
    df_merged <- merge(res, oxgenes, by = "genes", all.x = TRUE) 
    upgenes= df_merged %>% filter(FC>1.5)
    upgenes$dataset= paste0(gseid,"_",groupes[i],"_upregulated")
    downgenes= df_merged %>% filter(FC<0.6666)
    downgenes$dataset=paste0(gseid,"_",groupes[i],"_downregulated")
    res=rbind(upgenes, downgenes)
    save(res, file= paste0(gseid,"_",groupes[i],"paired.RData"))
    write.table(res, file= paste0(gseid,"_",groupes[i],"paired.xls"))
    assign(paste0(gseid,"_",groupes[i]),res)}
    
}

gseid
#---------> save the results in one file
save(GSE4271_Classical,GSE4271_Proneural,GSE4271_Mesenchymal, file="GSE4271pairedsubtypesDEG.RData")
