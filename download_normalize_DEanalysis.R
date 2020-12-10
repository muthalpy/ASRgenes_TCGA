root <- "/home/guest/Muthal"
dir_data <- file.path(root, "Data")

site <- c("BRCA")

library(dplyr)
library(tibble)
library(tidyverse)
library(SummarizedExperiment)
library(limma)
library(edgeR)
library(biomaRt)

query <- GDCquery(project = paste0("TCGA-",site),
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, 
            directory = dir_gdc)
data <- GDCprepare(query, 
                   directory = dir_gdc,
                   summarizedExperiment = T)

assign(paste0("rse_", site), data)
save(list = paste0("rse_", site),
     file = paste0("TCGA-", site, "-queriedv2.RData"))
####################
#normalization
mart <- useMart("ensembl", 
                dataset = "hsapiens_gene_ensembl")#,
site <- c("BRCA")
load(file.path(dir_data, paste0("TCGA-", site, "-queriedv2.RData")))
rse <- eval(parse(text = paste0("rse_")))
rm(rse_)
rse <- rse[, !is.na(rse$race)]
rse <- rse[, !is.na(rse$definition)]

pdat_tmp <- as.data.frame(colData(rse)) %>%
  rownames_to_column(var = "rname") %>%
  transmute(rname,
            patient, 
            barcode, 
            definition = definition,
            race = race, 
            subtype_BRCA_Subtype_PAM50= subtype_BRCA_Subtype_PAM50,
            status = factor(vital_status,
                            levels = c("Alive", "Dead"),
                            labels = c("0", "1")),
            month = ifelse(status == "0", 
                           days_to_last_follow_up/365.25*12, 
                           days_to_death/365.25*12)) 

# if there are duplicates, keep the samples with higher read depth
if(length(pdat_tmp$patient) != length(unique(pdat_tmp$patient))){
  pdat <- pdat_tmp %>%
    mutate(depth = colSums(assay(rse)[, rname])) %>%
    group_by(patient, definition) %>%
    filter(depth == max(depth)) %>%
    column_to_rownames(var = "rname")
} else {
  pdat <- pdat_tmp %>%
    column_to_rownames(var = "rname")
}
subls <- pdat$subtype_BRCA_Subtype_PAM50
subls[which(is.na(subls))] <- "atissues"
pdat$subtype_BRCA_Subtype_PAM50 <- subls
# Assay (filter low-count genes, normalize and transform counts using TMM+Voom method)
rse_sub <- rse[,colData(rse)$barcode %in% pdat$barcode]
adat_count <- assay(rse_sub)
designS <- model.matrix(~ subtype_BRCA_Subtype_PAM50 -1 , data = pdat)
designR <- model.matrix(~ race-1 , data = pdat)
designD <- model.matrix(~ definition-1 , data = pdat)
design <- cbind(designS, designR, designD)
rids_kept <- filterByExpr(adat_count, group =design) # use `group` instead of `design` due to a bug of filterByExpr
adat_count_sub <- adat_count[rids_kept, ]
norm_factor <- calcNormFactors(adat_count_sub)
lib_size <- colSums(adat_count_sub) * norm_factor
adat_expr <- voom(adat_count_sub, 
                  lib.size = lib_size, 
                  design = design)$E
# Row data of the assay (gene annotation)
rdat_expr <- as.data.frame(rowData(rse)) %>%
  left_join(getBM(values = .$ensembl_gene_id,
                  filters = "ensembl_gene_id",
                  attributes = c("ensembl_gene_id",
                                 "gene_biotype"),
                  mart = mart)) %>%
  dplyr::rename(Ens = ensembl_gene_id,
                Symbol = external_gene_name,
                Ens_ver = original_ensembl_gene_id,
                Biotype = gene_biotype) %>%
  distinct() %>%
  column_to_rownames(var = "Ens")
rdat_expr$Enz <- row.names(rdat_expr)
# Column data of the assay (survival, race)
cdat_expr <- pdat
# RSE object
cids <- rownames(cdat_expr)
rids <- rownames(adat_expr)
rse_expr <- SummarizedExperiment(assays = adat_expr[, cids],
                                 colData = cdat_expr[cids, ],
                                 rowData = rdat_expr[rids,]) #which(row.names(adat_expr) %in% rids)

assign(paste0("rse_expr_", site), rse_expr)
save(list = paste0("rse_expr_", site), 
     file = file.path(dir_data, paste0("TCGA-", site, "-expressionDifRaceSubType.RData")), row)

######################
#differential analysis Subtype vs Normal (DE_SN)

IBCgenes <- read.csv("./Muthal/Adaptive Stress Response GenesetMA.csv", header = TRUE, stringsAsFactors = FALSE)

unique(cdat_expr$subtype_BRCA_Subtype_PAM50)[!is.na(unique(cdat_expr$subtype_BRCA_Subtype_PAM50))]
#[1] "LumB"     "Normal"   "atissues" "LumA"     "Basal"    "Her2" 

subnamelst <- c("Basal", "Her2", "LumB", "LumA", "Normal")

subtypNall <- list()
subtypN <- list()
enzN <- list()
for (typ in subnamelst){
  adat_exprt <- adat_expr[,(cdat_expr$subtype_BRCA_Subtype_PAM50 %in% c(typ,"atissues"))]
  cdat_exprt <- cdat_expr[(cdat_expr$subtype_BRCA_Subtype_PAM50 %in% c(typ,"atissues")),]
  cdat_exprt$subtype_BRCA_Subtype_PAM50 <- factor(cdat_exprt$subtype_BRCA_Subtype_PAM50, levels = c(typ,"atissues"))
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

  save(res_de, 
       file = file.path(dir_data, paste0("TCGA-", site, "-expressionDifRaceDiff.RData")), row)
  
  IBCdiff <- res_de[res_de$Symbol %in% IBCgenes$Gene.Symbol,]
  IBCdifftyp <- IBCdiff[which((abs(IBCdiff$logFC)>1) &(IBCdiff$adj_p_value<0.05)),]
  print(nrow(IBCdifftyp))
  if(typ == "Normal"){
    typ <- "Normal-like"
  }
  write.csv(IBCdiff, paste0(typ, "vsNormalall.csv"), row.names = FALSE) #IBCdifftyp
  write.csv(IBCdifftyp, paste0(typ, "vsNormal.csv"), row.names = FALSE) #IBCdifftyp
  
  subtypNall[paste0(typ, "vsNormal")] <- list(IBCdiff) 
  subtypN[paste0(typ, "vsNormal")] <- list(IBCdifftyp) 
  enzN[paste0(typ, "vsNormal")] <- list(row.names(rse_expr)[which(rowData(rse_expr)$Symbol %in% IBCdifftyp$Symbol)])
  
}
save(subtypNall, subtypN,enzN, file = "./Muthal/Data/typevsNormallist.RData")
#########################
#differential analysis Subtype vs Subtype (DE_SS)

unique(cdat_expr$subtype_BRCA_Subtype_PAM50)[!is.na(unique(cdat_expr$subtype_BRCA_Subtype_PAM50))]
#"LumB"     "Normal"   "atissues" "LumA"     "Basal"    "Her2" 
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
    
    IBCdiff <- res_de[res_de$Symbol %in% IBCgenes$Gene.Symbol,]
    IBCdifftyp <- IBCdiff[which((abs(IBCdiff$logFC)>1) &(IBCdiff$adj_p_value<0.05)),]
    print(paste0(typ1, "vs", typ2, " : ", nrow(IBCdifftyp) ))
    
    alc <- append(alc, paste0(typ2, "vs", typ1))
    
    if(typ1 == "Normal"){
      typ3 <- "Normal-like"
    } else {typ3 <- typ1}    
    if(typ2 == "Normal"){
      typ2 <- "Normal-like"
    }
    
    write.csv(IBCdiff, paste0(typ3, "vs", typ2, "all.csv"), row.names = FALSE) #IBCdifftyp
    write.csv(IBCdifftyp, paste0(typ3, "vs", typ2, ".csv"), row.names = FALSE) 
    
    subtypall[paste0(typ3, "vs", typ2)] <- list(IBCdiff)
    subtypS[paste0(typ3, "vs", typ2)] <- list(IBCdifftyp)
    EnzS[paste0(typ3, "vs", typ2)] <- list(row.names(rse_expr)[which(rowData(rse_expr)$Symbol %in% IBCdifftyp$Symbol)])
  }
}
save(subtypall, subtypS,EnzS, file = "./Muthal/Data/typevstypelist.RData")

##### 
#differential analysis within AA or white (AA- or W-DE_SN)
#the code below will run twice, once for W-DE_SN and once for AA-DE_SN

load(file.path(dir_data, paste0("TCGA-", site, "-expressionDifRaceSubType.RData")))
rse_expr <- eval(parse(text = paste0("rse_expr_", site)))
rse_expr <- rse_expr[, -which((rse_expr$subtype_BRCA_Subtype_PAM50 == "atissues")&(rse_expr$definition != "Solid Tissue Normal"))]

rse_expr <- rse_expr[,(colData(rse_expr)$race == "white")]
#rse_expr <- rse_expr[,(colData(rse_expr)$race == "black or african american")]

adat_expr <- assay(rse_expr)
cdat_expr <- as.data.frame(colData(rse_expr))
rdat_expr <- as.data.frame(rowData(rse_expr))
row.names(rdat_expr) <- row.names(adat_expr)

unique(cdat_expr$subtype_BRCA_Subtype_PAM50)[!is.na(unique(cdat_expr$subtype_BRCA_Subtype_PAM50))]
#[1] "LumB"     "Normal"   "atissues" "LumA"     "Basal"    "Her2" 

subnamelst <- c("Basal", "Her2", "LumB", "LumA", "Normal")

subtypNallRaceW <- list()
subtypNRaceW <- list()
enzNRaceW <- list()
for (typ in subnamelst){
  adat_exprt <- adat_expr[,(cdat_expr$subtype_BRCA_Subtype_PAM50 %in% c(typ,"atissues"))]
  cdat_exprt <- cdat_expr[(cdat_expr$subtype_BRCA_Subtype_PAM50 %in% c(typ,"atissues")),]
  cdat_exprt$subtype_BRCA_Subtype_PAM50 <- factor(cdat_exprt$subtype_BRCA_Subtype_PAM50, levels = c(typ,"atissues"))
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

  IBCdiffRace <- res_de[res_de$Symbol %in% IBCgenes$Gene.Symbol,]
  IBCdifftypRace <- IBCdiffRace[which((abs(IBCdiffRace$logFC)>1) &(IBCdiffRace$adj_p_value<0.05)),]
  print(nrow(IBCdifftypRace))
  if(typ == "Normal"){
    typ <- "Normal-like"
  }
  write.csv(IBCdiffRace, paste0(typ, "vsNormalallRaceW.csv"), row.names = FALSE) #IBCdifftyp
  write.csv(IBCdifftypRace, paste0(typ, "vsNormalRaceW.csv"), row.names = FALSE) #IBCdifftyp
  
  subtypNallRaceW[paste0(typ, "vsNormal")] <- list(IBCdiffRace)
  subtypNRaceW[paste0(typ, "vsNormal")] <- list(IBCdifftypRace)
  enzNRaceW[paste0(typ, "vsNormal")] <- list(row.names(rse_expr)[which(rowData(rse_expr)$Symbol %in% IBCdifftypRace$Symbol)])
  
}
#differential analysis within AA or white (AA- or W-DE_SS)
#the code below will run twice, once for W-DE_SS and once for AA-DE_SS

unique(cdat_expr$subtype_BRCA_Subtype_PAM50)[!is.na(unique(cdat_expr$subtype_BRCA_Subtype_PAM50))]
#"LumB"     "Normal"   "atissues" "LumA"     "Basal"    "Her2" 
subnamelst <- c("Basal", "Her2", "LumB", "LumA", "Normal")
subtypallRaceW <- list()
subtypSRaceW <- list()
EnzSRaceW <- list()
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
    
    IBCdiffRace <- res_de[res_de$Symbol %in% IBCgenes$Gene.Symbol,]
    IBCdifftypRace <- IBCdiffRace[which((abs(IBCdiffRace$logFC)>1) &(IBCdiffRace$adj_p_value<0.05)),]
    print(paste0(typ1, "vs", typ2, " : ", nrow(IBCdifftypRace) ))
    
    alc <- append(alc, paste0(typ2, "vs", typ1))
    
    if(typ1 == "Normal"){
      typ3 <- "Normal-like"
    } else {typ3 <- typ1}    
    if(typ2 == "Normal"){
      typ2 <- "Normal-like"
    }
    
    write.csv(IBCdiffRace, paste0(typ3, "vs", typ2, "RaceWall.csv"), row.names = FALSE) #IBCdifftyp
    write.csv(IBCdifftypRace, paste0(typ3, "vs", typ2, "RaceW.csv"), row.names = FALSE) 
    
    subtypallRaceW[paste0(typ3, "vs", typ2)] <- list(IBCdiffRace)
    subtypSRaceW[paste0(typ3, "vs", typ2)] <- list(IBCdifftypRace)
    EnzSRaceW[paste0(typ3, "vs", typ2)] <- list(row.names(rse_expr)[which(rowData(rse_expr)$Symbol %in% IBCdifftypRace$Symbol)])
  }
}

save(subtypNallRaceW, subtypNRaceW, enzNRaceW, subtypallRaceW, subtypSRaceW,EnzSRaceW, 
     subtypNallRaceW, subtypNRaceW, enzNRaceW, subtypallRaceW, subtypSRaceW,EnzSRaceW,
     file = "./Muthal/Data/typevsNormalandtypevstypeRacelist.RData")


##############################################################################
##############################################################################
#differntial analysis for race-related DEs
#this was modified on Oct 22nd as Steven suggested

root <- "/home/guest/Muthal"
dir_data <- file.path(root, "Data")

site <- "BRCA"
load(file.path(dir_data, paste0("TCGA-", site, "-expressionDifRaceSubType.RData")))
rse_expr <- eval(parse(text = paste0("rse_expr_", site)))
rse_expr <- rse_expr[, -which((rse_expr$subtype_BRCA_Subtype_PAM50 == "atissues")&(rse_expr$definition != "Solid Tissue Normal"))]
adat_expr <- assay(rse_expr)
cdat_expr <- as.data.frame(colData(rse_expr))
rdat_expr <- as.data.frame(rowData(rse_expr))
row.names(rdat_expr) <- row.names(adat_expr)

adat_expr <- adat_expr[,(cdat_expr$race %in% c("white", "black or african american")) & (cdat_expr$definition != "Solid Tissue Normal")]
cdat_expr <- cdat_expr[(cdat_expr$race %in% c("white", "black or african american")) & (cdat_expr$definition != "Solid Tissue Normal"),]
cdat_expr$race <- factor(cdat_expr$race, levels = c("white", "black or african american"))
levels(cdat_expr$race) <- c("W", "AA")

cdat_expr$racesubtype <- factor(paste0(cdat_expr$race, "_", cdat_expr$subtype_BRCA_Subtype_PAM50))
design <- model.matrix(~ 0+racesubtype, data = cdat_expr)
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

write.csv(res_de, "res_de.csv", row.names = FALSE)

IBCdiffRaceDE <- res_de[res_de$Symbol %in% IBCgenes$Gene.Symbol,]
IBCdifftypRaceDE <- IBCdiffRaceDE[which((abs(IBCdiffRaceDE$Basal)>1) &(IBCdiffRace$adj_p_value<0.05)),6]
