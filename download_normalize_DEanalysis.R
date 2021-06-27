set.seed(1991)

dir_data <- "/home/rapiduser/mount/General/ASR_Gay"
site <- c("BRCA")

library(dplyr)
library(tibble)
library(tidyverse)
library(limma) #version 3.40.2
library(TCGAbiolinks) # version 2.13.3
library(SummarizedExperiment) 
library(edgeR)
library(biomaRt)

# query <- GDCquery(project = paste0("TCGA-",site),
#                   data.category = "Transcriptome Profiling",
#                   data.type = "Gene Expression Quantification", 
#                   workflow.type = "HTSeq - Counts")
# GDCdownload(query, 
#             directory = dir_data)
# data <- GDCprepare(query, 
#                    directory = dir_data,
#                    summarizedExperiment = T)
# 
# assign(paste0("rse_", site), data) 
# dim(rse_BRCA) 
# # 56512  1221
#   save(list = paste0("rse_", site),
#        file = file.path(dir_data, paste0("TCGA-", site, "-queriedv2.RData")))

####################

#normalization
# mart <- useMart("ensembl", 
#                 dataset = "hsapiens_gene_ensembl")
# load(file.path(dir_data, paste0("TCGA-", site, "-queriedv2.RData")))
# rse <- eval(parse(text = paste0("rse_", site)))
# rse <- rse[,-which((is.na(rse$paper_BRCA_Subtype_PAM50))&(rse$definition != "Solid Tissue Normal"))]
# rse <- rse[, !is.na(rse$race)]
# rse <- rse[, !is.na(rse$definition)]
# dim(rse)
# # 56512  1210
# pdat_tmp <- as.data.frame(colData(rse)) %>%
#   rownames_to_column(var = "rname") %>%
#   transmute(rname,
#             patient, 
#             barcode, 
#             definition = definition,
#             race = race, 
#             subtype_BRCA_Subtype_PAM50= paper_BRCA_Subtype_PAM50,
#             status = factor(vital_status,
#                             levels = c("Alive", "Dead"),
#                             labels = c("0", "1")),
#             month = ifelse(status == "0", 
#                            days_to_last_follow_up/365.25*12, 
#                            days_to_death/365.25*12)) 
# 
# # if there are duplicates, keep the samples with higher read depth
# if(length(pdat_tmp$patient) != length(unique(pdat_tmp$patient))){
#   pdat <- pdat_tmp %>%
#     mutate(depth = colSums(assay(rse)[, rname])) %>%
#     group_by(patient, definition) %>%
#     filter(depth == max(depth)) %>%
#     column_to_rownames(var = "rname")
# } else {
#   pdat <- pdat_tmp %>%
#     column_to_rownames(var = "rname")
# }
# subls <- pdat$subtype_BRCA_Subtype_PAM50
# subls[which(is.na(subls))] <- "NormAdj"
# pdat$subtype_BRCA_Subtype_PAM50 <- subls
# 
# # Assay (filter low-count genes, normalize and transform counts using TMM+Voom method)
# rse_sub <- rse[,colData(rse)$barcode %in% pdat$barcode]
# adat_count <- assay(rse_sub)
# designS <- model.matrix(~ subtype_BRCA_Subtype_PAM50 -1 , data = pdat)
# designR <- model.matrix(~ race-1 , data = pdat)
# designD <- model.matrix(~ definition-1 , data = pdat)
# design <- cbind(designS, designR, designD)
# rids_kept <- filterByExpr(adat_count, design =design)
# adat_count_sub <- adat_count[rids_kept, ]
# dim(adat_count_sub)
# # 46077  1210
# norm_factor <- calcNormFactors(adat_count)
# lib_size <- colSums(adat_count) * norm_factor
# adat_expr <- voom(adat_count, 
#                   lib.size = lib_size, 
#                   design = design)$E
# # Row data of the assay (gene annotation)
# rdat_expr <- as.data.frame(rowData(rse)) %>%
#   left_join(getBM(values = .$ensembl_gene_id,
#                   filters = "ensembl_gene_id",
#                   attributes = c("ensembl_gene_id",
#                                  "gene_biotype"),
#                   mart = mart)) %>%
#   dplyr::rename(Ens = ensembl_gene_id,
#                 Symbol = external_gene_name,
#                 Ens_ver = original_ensembl_gene_id,
#                 Biotype = gene_biotype) %>%
#   distinct() %>%
#   column_to_rownames(var = "Ens")
# rdat_expr$Enz <- row.names(rdat_expr)
# 
# # Column data of the assay (survival, race)
# cdat_expr <- pdat
# 
# # RSE object
# cids <- rownames(cdat_expr)
# rids <- rownames(adat_expr)
# rse_expr <- SummarizedExperiment(assays = adat_expr[, cids],
#                                  colData = cdat_expr[cids, ],
#                                  rowData = rdat_expr[rids,])
# dim(rse_expr)
# # 46077  1210
# 
# assign(paste0("rse_expr_", site), rse_expr)
# save(list = paste0("rse_expr_", site), 
#      file = file.path(dir_data, paste0("TCGA-", site, "-expressionDifRaceSubType.RData")), row)

######################
#differential analysis Subtype vs Normal (DE_SN)

load(file.path(dir_data, paste0("TCGA-", site, "-expressionDifRaceSubType.RData")))

rse_expr <- eval(parse(text = paste0("rse_expr_", site)))
adat_expr <- assay(rse_expr)
cdat_expr <- as.data.frame(colData(rse_expr))
rdat_expr <- as.data.frame(rowData(rse_expr))
rm(list = paste0("rse_expr_",site))
ASRgenes <- read.csv(file.path(dir_data,"ASRgenes.csv"), header = TRUE, stringsAsFactors = FALSE)

subnamelst <- c("Basal", "Her2", "LumB", "LumA", "Normal")

subtypNall <- list()
subtypN <- list()
enzN <- list()
for (typ in subnamelst){
  adat_exprt <- adat_expr[,(cdat_expr$subtype_BRCA_Subtype_PAM50 %in% c(typ,"NormAdj"))]
  cdat_exprt <- cdat_expr[(cdat_expr$subtype_BRCA_Subtype_PAM50 %in% c(typ,"NormAdj")),]
  cdat_exprt$subtype_BRCA_Subtype_PAM50 <- factor(cdat_exprt$subtype_BRCA_Subtype_PAM50, levels = c(typ,"NormAdj"))
  design <- model.matrix(~ subtype_BRCA_Subtype_PAM50 , data = cdat_exprt)
  fit <- lmFit(adat_exprt, design)
  fit <- eBayes(fit)
  res_de <- topTable(fit, n=Inf) %>%
    merge(rdat_expr, ., by = 0) %>%
    filter(Biotype == "protein_coding") %>%
    transmute(Symbol, logFC = -logFC, AveExpr, 
              p_value = P.Value,
              adj_p_value = p.adjust(p_value, metho = "BH")) %>%
    arrange(p_value) 
  
  ASRdiff <- res_de[res_de$Symbol %in% ASRgenes$Gene.Symbol,]
  ASRdifftyp <- ASRdiff[which((abs(ASRdiff$logFC)>1) &(ASRdiff$adj_p_value<0.05)),]
  print(paste0(typ, " vs Normal_adj: ", nrow(ASRdifftyp)))
  if(typ == "Normal"){
    typ <- "Normal-like"
  }
  write.csv(ASRdiff, file.path(dir_data, paste0(typ, "vsNormalall.csv")), row.names = FALSE)
  
  subtypNall[paste0(typ, "vsNormal")] <- list(ASRdiff) 
  subtypN[paste0(typ, "vsNormal")] <- list(ASRdifftyp) 
  enzN[paste0(typ, "vsNormal")] <- list(row.names(rse_expr)[which(rowData(rse_expr)$Symbol %in% ASRdifftyp$Symbol)])
}
save(subtypNall, subtypN,enzN, file = "/home/rapiduser/mount/General/ASR_Gay/typevsNormallist.RData")

#########################
#differential analysis Subtype vs Subtype (DE_SS)

subnamelst <- c("Basal", "Her2", "LumB", "LumA", "Normal")
subtypall <- list()
subtypS <- list()
EnzS <- list()
alc <- c()
for (i in 1:5){
  typ1 <- subnamelst[i]
  for(ii in 1:5){
    if (i == ii){next}
    typ2 <- subnamelst[ii]
    if(paste0(typ1, "vs", typ2) %in% alc){next}
    adat_exprt <- adat_expr[,(cdat_expr$subtype_BRCA_Subtype_PAM50 %in% c(typ1, typ2)) ]
    cdat_exprt <- cdat_expr[(cdat_expr$subtype_BRCA_Subtype_PAM50 %in% c(typ1, typ2)),]
    cdat_exprt$subtype_BRCA_Subtype_PAM50 <- factor(cdat_exprt$subtype_BRCA_Subtype_PAM50, levels = c(typ1, typ2))
    design <- model.matrix(~ subtype_BRCA_Subtype_PAM50 , data = cdat_exprt)
    fit <- lmFit(adat_exprt, design)
    fit <- eBayes(fit)
    res_de <- topTable(fit, n=Inf) %>%
      merge(rdat_expr, ., by = 0) %>%
      filter(Biotype == "protein_coding") %>%
      transmute(Symbol, logFC = -logFC, AveExpr, 
                p_value = P.Value,
                adj_p_value = p.adjust(p_value, metho = "BH")) %>%
      arrange(p_value) 
    
    ASRdiff <- res_de[res_de$Symbol %in% ASRgenes$Gene.Symbol,]
    ASRdifftyp <- ASRdiff[which((abs(ASRdiff$logFC)>1) &(ASRdiff$adj_p_value<0.05)),]
    print(paste0(typ1, "vs", typ2, " : ", nrow(ASRdifftyp) ))
    
    alc <- append(alc, paste0(typ2, "vs", typ1))
    
    if(typ1 == "Normal"){
      typ3 <- "Normal-like"
    } else {typ3 <- typ1}    
    if(typ2 == "Normal"){
      typ2 <- "Normal-like"
    }
    
    write.csv(ASRdiff, file.path(dir_data, paste0(typ3, "vs", typ2, "all.csv")), row.names = FALSE)
    
    subtypall[paste0(typ3, "vs", typ2)] <- list(ASRdiff)
    subtypS[paste0(typ3, "vs", typ2)] <- list(ASRdifftyp)
    EnzS[paste0(typ3, "vs", typ2)] <- list(row.names(rse_expr)[which(rowData(rse_expr)$Symbol %in% ASRdifftyp$Symbol)])
  }
}
save(subtypall, subtypS,EnzS, file = "/home/rapiduser/mount/General/ASR_Gay/typevstypelist.RData")

##### 
##### 
#differential analysis within W (W-DE_SN)

rse_exprW <- rse_expr
rse_exprW <- rse_exprW[,(colData(rse_exprW)$race == "white")]

adat_expr <- assay(rse_exprW)
cdat_expr <- as.data.frame(colData(rse_exprW))
rdat_expr <- as.data.frame(rowData(rse_exprW))
row.names(rdat_expr) <- row.names(adat_expr)

unique(cdat_expr$subtype_BRCA_Subtype_PAM50)[!is.na(unique(cdat_expr$subtype_BRCA_Subtype_PAM50))]
#[1] "LumB"     "Normal"   "NormAdj" "LumA"     "Basal"    "Her2" 

subnamelst <- c("Basal", "Her2", "LumB", "LumA", "Normal")

subtypNallRaceW <- list()
subtypNRaceW <- list()
enzNRaceW <- list()
for (typ in subnamelst){
  adat_exprt <- adat_expr[,(cdat_expr$subtype_BRCA_Subtype_PAM50 %in% c(typ,"NormAdj"))]
  cdat_exprt <- cdat_expr[(cdat_expr$subtype_BRCA_Subtype_PAM50 %in% c(typ,"NormAdj")),]
  cdat_exprt$subtype_BRCA_Subtype_PAM50 <- factor(cdat_exprt$subtype_BRCA_Subtype_PAM50, levels = c(typ,"NormAdj"))
  design <- model.matrix(~ subtype_BRCA_Subtype_PAM50 , data = cdat_exprt)
  fit <- lmFit(adat_exprt, design)
  fit <- eBayes(fit)
  res_de <- topTable(fit, n=Inf) %>%
    merge(rdat_expr, ., by = 0) %>%
    filter(Biotype == "protein_coding") %>%
    transmute(Symbol, logFC = -logFC, AveExpr, 
              p_value = P.Value,
              adj_p_value = p.adjust(p_value, metho = "BH")) %>%
    arrange(p_value) 
  
  ASRdiffRace <- res_de[res_de$Symbol %in% ASRgenes$Gene.Symbol,]
  ASRdifftypRace <- ASRdiffRace[which((abs(ASRdiffRace$logFC)>1) &(ASRdiffRace$adj_p_value<0.05)),]
  print(nrow(ASRdifftypRace))
  if(typ == "Normal"){
    typ <- "Normal-like"
  }
  write.csv(ASRdiffRace, file.path(dir_data, paste0(typ, "vsNormalallRaceW.csv")), row.names = FALSE) #ASRdifftyp
  
  subtypNallRaceW[paste0(typ, "vsNormal")] <- list(ASRdiffRace)
  subtypNRaceW[paste0(typ, "vsNormal")] <- list(ASRdifftypRace)
  enzNRaceW[paste0(typ, "vsNormal")] <- list(row.names(rse_exprW)[which(rowData(rse_exprW)$Symbol %in% ASRdifftypRace$Symbol)])
  
}

#differential analysis within W (W-DE_SS)

unique(cdat_expr$subtype_BRCA_Subtype_PAM50)[!is.na(unique(cdat_expr$subtype_BRCA_Subtype_PAM50))]
#"LumB"     "Normal"   "NormAdj" "LumA"     "Basal"    "Her2" 
subnamelst <- c("Basal", "Her2", "LumB", "LumA", "Normal")
subtypallRaceW <- list()
subtypSRaceW <- list()
EnzSRaceW <- list()
for (i in 1:5){
  typ1 <- subnamelst[i]
  for(ii in 1:5){
    if (i == ii){next}
    typ2 <- subnamelst[ii]
    if(paste0(typ1, "vs", typ2) %in% alc){next}
    adat_exprt <- adat_expr[,(cdat_expr$subtype_BRCA_Subtype_PAM50 %in% c(typ1, typ2)) ]
    cdat_exprt <- cdat_expr[(cdat_expr$subtype_BRCA_Subtype_PAM50 %in% c(typ1, typ2)),]
    cdat_exprt$subtype_BRCA_Subtype_PAM50 <- factor(cdat_exprt$subtype_BRCA_Subtype_PAM50, levels = c(typ1, typ2))
    design <- model.matrix(~ subtype_BRCA_Subtype_PAM50 , data = cdat_exprt)
    fit <- lmFit(adat_exprt, design)
    fit <- eBayes(fit)
    res_de <- topTable(fit, n=Inf) %>%
      merge(rdat_expr, ., by = 0) %>%
      filter(Biotype == "protein_coding") %>%
      transmute(Symbol, logFC = -logFC, AveExpr, 
                p_value = P.Value,
                adj_p_value = p.adjust(p_value, metho = "BH")) %>%
      arrange(p_value) 
    
    ASRdiffRace <- res_de[res_de$Symbol %in% ASRgenes$Gene.Symbol,]
    ASRdifftypRace <- ASRdiffRace[which((abs(ASRdiffRace$logFC)>1) &(ASRdiffRace$adj_p_value<0.05)),]
    print(paste0(typ1, "vs", typ2, " : ", nrow(ASRdifftypRace) ))
    
    
    if(typ1 == "Normal"){
      typ3 <- "Normal-like"
    } else {typ3 <- typ1}    
    if(typ2 == "Normal"){
      typ2 <- "Normal-like"
    }
    
    write.csv(ASRdiffRace, file.path(dir_data, paste0(typ3, "vs", typ2, "RaceWall.csv")), row.names = FALSE) #ASRdifftyp
    
    subtypallRaceW[paste0(typ3, "vs", typ2)] <- list(ASRdiffRace)
    subtypSRaceW[paste0(typ3, "vs", typ2)] <- list(ASRdifftypRace)
    EnzSRaceW[paste0(typ3, "vs", typ2)] <- list(row.names(rse_exprW)[which(rowData(rse_exprW)$Symbol %in% ASRdifftypRace$Symbol)])
  }
}

##### 
#differential analysis within AA (AA-DE_SN)

rse_exprAA <- rse_expr
rse_exprAA <- rse_exprAA[,(colData(rse_exprAA)$race == "black or african american")]

adat_expr <- assay(rse_exprAA)
cdat_expr <- as.data.frame(colData(rse_exprAA))
rdat_expr <- as.data.frame(rowData(rse_exprAA))
row.names(rdat_expr) <- row.names(adat_expr)

unique(cdat_expr$subtype_BRCA_Subtype_PAM50)[!is.na(unique(cdat_expr$subtype_BRCA_Subtype_PAM50))]
#[1] "LumB"     "Normal"   "NormAdj" "LumA"     "Basal"    "Her2" 

subnamelst <- c("Basal", "Her2", "LumB", "LumA", "Normal")

subtypNallRaceAA <- list()
subtypNRaceAA <- list()
enzNRaceAA <- list()
for (typ in subnamelst){
  adat_exprt <- adat_expr[,(cdat_expr$subtype_BRCA_Subtype_PAM50 %in% c(typ,"NormAdj"))]
  cdat_exprt <- cdat_expr[(cdat_expr$subtype_BRCA_Subtype_PAM50 %in% c(typ,"NormAdj")),]
  cdat_exprt$subtype_BRCA_Subtype_PAM50 <- factor(cdat_exprt$subtype_BRCA_Subtype_PAM50, levels = c(typ,"NormAdj"))
  design <- model.matrix(~ subtype_BRCA_Subtype_PAM50 , data = cdat_exprt)
  fit <- lmFit(adat_exprt, design)
  fit <- eBayes(fit)
  res_de <- topTable(fit, n=Inf) %>%
    merge(rdat_expr, ., by = 0) %>%
    filter(Biotype == "protein_coding") %>%
    transmute(Symbol, logFC = -logFC, AveExpr, 
              p_value = P.Value,
              adj_p_value = p.adjust(p_value, metho = "BH")) %>%
    arrange(p_value) 
  
  ASRdiffRace <- res_de[res_de$Symbol %in% ASRgenes$Gene.Symbol,]
  ASRdifftypRace <- ASRdiffRace[which((abs(ASRdiffRace$logFC)>1) &(ASRdiffRace$adj_p_value<0.05)),]
  print(nrow(ASRdifftypRace))
  if(typ == "Normal"){
    typ <- "Normal-like"
  }
  write.csv(ASRdiffRace, file.path(dir_data, paste0(typ, "vsNormalallRaceAA.csv")), row.names = FALSE)
  
  subtypNallRaceAA[paste0(typ, "vsNormal")] <- list(ASRdiffRace)
  subtypNRaceAA[paste0(typ, "vsNormal")] <- list(ASRdifftypRace)
  enzNRaceAA[paste0(typ, "vsNormal")] <- list(row.names(rse_exprAA)[which(rowData(rse_exprAA)$Symbol %in% ASRdifftypRace$Symbol)])
}
#differential analysis within AA (AA-DE_SS)

unique(cdat_expr$subtype_BRCA_Subtype_PAM50)[!is.na(unique(cdat_expr$subtype_BRCA_Subtype_PAM50))]
#"LumB"     "Normal"   "NormAdj" "LumA"     "Basal"    "Her2" 
subnamelst <- c("Basal", "Her2", "LumB", "LumA", "Normal")
subtypallRaceAA <- list()
subtypSRaceAA <- list()
EnzSRaceAA <- list()
alc <- c()
for (i in 1:5){
  typ1 <- subnamelst[i]
  for(ii in 1:5){
    if (i == ii){next}
    typ2 <- subnamelst[ii]
    if(paste0(typ1, "vs", typ2) %in% alc){next}
    adat_exprt <- adat_expr[,(cdat_expr$subtype_BRCA_Subtype_PAM50 %in% c(typ1, typ2)) ]
    cdat_exprt <- cdat_expr[(cdat_expr$subtype_BRCA_Subtype_PAM50 %in% c(typ1, typ2)),]
    cdat_exprt$subtype_BRCA_Subtype_PAM50 <- factor(cdat_exprt$subtype_BRCA_Subtype_PAM50, levels = c(typ1, typ2))
    design <- model.matrix(~ subtype_BRCA_Subtype_PAM50 , data = cdat_exprt)
    fit <- lmFit(adat_exprt, design)
    fit <- eBayes(fit)
    res_de <- topTable(fit, n=Inf) %>%
      merge(rdat_expr, ., by = 0) %>%
      filter(Biotype == "protein_coding") %>%
      transmute(Symbol, logFC = -logFC, AveExpr, 
                p_value = P.Value,
                adj_p_value = p.adjust(p_value, metho = "BH")) %>%
      arrange(p_value) 
    
    ASRdiffRace <- res_de[res_de$Symbol %in% ASRgenes$Gene.Symbol,]
    ASRdifftypRace <- ASRdiffRace[which((abs(ASRdiffRace$logFC)>1) &(ASRdiffRace$adj_p_value<0.05)),]
    print(paste0(typ1, "vs", typ2, " : ", nrow(ASRdifftypRace) ))
    
    alc <- append(alc, paste0(typ2, "vs", typ1))
    
    if(typ1 == "Normal"){
      typ3 <- "Normal-like"
    } else {typ3 <- typ1}    
    if(typ2 == "Normal"){
      typ2 <- "Normal-like"
    }
    
    write.csv(ASRdiffRace, file.path(dir_data, paste0(typ3, "vs", typ2, "RaceAAall.csv")), row.names = FALSE)
    
    subtypallRaceAA[paste0(typ3, "vs", typ2)] <- list(ASRdiffRace)
    subtypSRaceAA[paste0(typ3, "vs", typ2)] <- list(ASRdifftypRace)
    EnzSRaceAA[paste0(typ3, "vs", typ2)] <- list(row.names(rse_exprAA)[which(rowData(rse_exprAA)$Symbol %in% ASRdifftypRace$Symbol)])
  }
} 

save(subtypNallRaceW, subtypNRaceW, enzNRaceW, subtypallRaceW, subtypSRaceW,EnzSRaceW, 
     subtypNallRaceAA, subtypNRaceAA, enzNRaceAA, subtypallRaceAA, subtypSRaceAA,EnzSRaceAA,
     file = "/home/rapiduser/mount/General/ASR_Gay/typevsNormalandtypevstypeRacelist.RData")

###############
# differential race-related ASR
print(1)
adat_expr <- assay(rse_expr)
cdat_expr <- as.data.frame(colData(rse_expr))
rdat_expr <- as.data.frame(rowData(rse_expr))
row.names(rdat_expr) <- row.names(adat_expr)
print(2) 
adat_expr <- adat_expr[,(cdat_expr$race %in% c("white", "black or african american")) & (cdat_expr$definition != "Solid Tissue Normal")]
cdat_expr <- cdat_expr[(cdat_expr$race %in% c("white", "black or african american")) & (cdat_expr$definition != "Solid Tissue Normal"),]
cdat_expr$race <- factor(cdat_expr$race, levels = c("white", "black or african american"))
levels(cdat_expr$race) <- c("W", "AA")
print(3)
cdat_expr$racesubtype <- factor(paste0(cdat_expr$race, "_", cdat_expr$subtype_BRCA_Subtype_PAM50))
design <- model.matrix(~ 0+racesubtype, data = cdat_expr)
print("design")
fit <- lmFit(adat_expr, design)
head(design,1)
cont.matrix <- makeContrasts(Basal="racesubtypeAA_Basal-racesubtypeW_Basal", 
                             Her2="racesubtypeAA_Her2-racesubtypeW_Her2", 
                             LumB="racesubtypeAA_LumB-racesubtypeW_LumB", 
                             LumA="racesubtypeAA_LumA-racesubtypeW_LumA", 
                             Normal="racesubtypeAA_Normal-racesubtypeW_Normal",
                             levels=design)
fit2  <- contrasts.fit(fit, cont.matrix)
fit3  <- eBayes(fit2, trend = TRUE)
res_de <- topTable(fit3, n=Inf) %>%
  merge(rdat_expr, ., by = 0) %>%
  filter(Biotype == "protein_coding")
print("res_de")
ASRdiffRaceDE <- res_de[res_de$Symbol %in% ASRgenes$Gene.Symbol,]
print(nrow(ASRdiffRaceDE))
ASRdifftypRaceDEbasal <- ASRdiffRaceDE[which((abs(ASRdiffRaceDE$Basal)>1) &(ASRdiffRaceDE$P.Value<0.05)),]
ASRdifftypRaceDEHer2 <- ASRdiffRaceDE[which((abs(ASRdiffRaceDE$Her2)>1) &(ASRdiffRaceDE$P.Value<0.05)),]
ASRdifftypRaceDELumB <- ASRdiffRaceDE[which((abs(ASRdiffRaceDE$LumB)>1) &(ASRdiffRaceDE$P.Value<0.05)),]
ASRdifftypRaceDELumA <- ASRdiffRaceDE[which((abs(ASRdiffRaceDE$LumA)>1) &(ASRdiffRaceDE$P.Value<0.05)),]
ASRdifftypRaceDENormallike <- ASRdiffRaceDE[which((abs(ASRdiffRaceDE$Normal)>1) &(ASRdiffRaceDE$P.Value<0.05)),]

print(nrow(ASRdifftypRaceDEbasal))
print(nrow(ASRdifftypRaceDEHer2))
print(nrow(ASRdifftypRaceDELumB))
print(nrow(ASRdifftypRaceDELumA))
print(nrow(ASRdifftypRaceDENormallike))

ASRdifftypRaceDE <- rbind(ASRdifftypRaceDEbasal,ASRdifftypRaceDEHer2, ASRdifftypRaceDELumB, ASRdifftypRaceDELumA, ASRdifftypRaceDENormallike)

write.csv(ASRdifftypRaceDE, file = file.path(dir_data, "AAvsW.csv"), row.names = FALSE)


sessionInfo()
