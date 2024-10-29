---
title: "FMD data alysis"
author: "Phong BH Nguyen"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    code_folding: hide
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, fig.path='figures/fig-', warning = F, message = F)
```

```{r, warning=FALSE}
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(tidyr))
suppressPackageStartupMessages(require(tibble))
suppressPackageStartupMessages(require(ComplexHeatmap))
suppressPackageStartupMessages(require(ConsensusClusterPlus))
suppressPackageStartupMessages(require(ggrepel))
suppressPackageStartupMessages(require(limma))
suppressPackageStartupMessages(require(affy)) # ExpressionSet
suppressPackageStartupMessages(require(fgsea))
suppressPackageStartupMessages(require(circlize))
suppressPackageStartupMessages(require(RColorBrewer))
suppressPackageStartupMessages(require(ggpubr))
suppressPackageStartupMessages(require(caret))
suppressPackageStartupMessages(require(pROC))
suppressPackageStartupMessages(require(impute))
suppressPackageStartupMessages(require(ggConvexHull))
suppressPackageStartupMessages(require(fmsb))
suppressPackageStartupMessages(require(ggalluvial))
suppressPackageStartupMessages(require(limma))
suppressPackageStartupMessages(require(readxl))
suppressPackageStartupMessages(require(openxlsx))
suppressPackageStartupMessages(require(XLConnect))
suppressPackageStartupMessages(require(xlsx))
suppressPackageStartupMessages(require(globaltest))
suppressPackageStartupMessages(require(igraph))
suppressPackageStartupMessages(require(ggraph))
suppressPackageStartupMessages(require(scatterpie))
suppressPackageStartupMessages(require(ggstats))
suppressPackageStartupMessages(require(ggsignif))
```

# Processing clinical variables

```{r}
clin <- read_excel("FMD_Master_V0.xlsx", sheet = 1)
#phys <- read_excel("FMD_Master_V0.xlsx", sheet = 2, range = "A1:CT41")
lab <- read_excel("FMD_Master_V0.xlsx", sheet = 3)

clin_proc <- clin %>% as.data.frame()
colnames(clin_proc) <- gsub(" ","_", colnames(clin_proc)) %>%
  gsub("-","_",.)
clin_proc$Gender <- ifelse(clin_proc$Gender == "M",1,0)

lab_proc <- lab %>%
  dplyr::select(-c(Chlorid...5, Magnesium...6, `B-Ketonkörper`, `U-Ketonkörper`))
colnames(lab_proc) <- gsub(" ","_", colnames(lab_proc)) %>%
  gsub("-","_",.)

comb_proc <- left_join(clin_proc, lab_proc, by = "Patient_ID")
keep_col <- apply(comb_proc, 2, function(x) sum(is.na(x)) <= length(x)*30/100)
keep_row <- apply(comb_proc, 1, function(x) sum(is.na(x)) <= length(x)*30/100)
comb_proc <- comb_proc[keep_row, keep_col]
pat_ID <- comb_proc$Patient_ID
comb_proc <- apply(dplyr::select(comb_proc, -Patient_ID), 2, function(x) as.numeric(x)) %>% as.data.frame()
comb_proc <- impute.knn(t(comb_proc))$data %>% t() %>% as.data.frame()
comb_proc$Patient_ID <- pat_ID
```


```{r}
lab_f0 <- read_excel("FMD_Master_V0.xlsx", sheet = 3)
lab_f0_proc <- lab_f0 %>%
  dplyr::select(-c(Chlorid...5, Magnesium...6, `U-Ketonkörper`,`IL-6`, `Lp(a)`,`HOMA-IR`,Insulin))
lab_f0_proc$`B-Ketonkörper` <- ifelse(lab_f0_proc$`B-Ketonkörper` == "<0,1", '0.09', lab_f0_proc$`B-Ketonkörper`) %>%
  as.numeric()
phys_f0 <- read_excel("FMD_Master_V0.xlsx", sheet = 2, range = "A1:CT41") %>%
  dplyr::select(`Patient-ID`, BMI, `Systole OA re`)
lab_f0_proc <- inner_join(phys_f0, lab_f0_proc, by = "Patient-ID")
colnames(lab_f0_proc) <- gsub(" ","_", colnames(lab_f0_proc)) %>%
  gsub("-","_",.) %>%
  paste0(.,"f0")
colnames(lab_f0_proc)[1] <- "Patient_ID"
keep_col <- apply(lab_f0_proc, 2, function(x) sum(is.na(x)) <= length(x)*30/100)
keep_row <- apply(lab_f0_proc, 1, function(x) sum(is.na(x)) <= length(x)*30/100)
lab_f0_proc <- lab_f0_proc[keep_row, keep_col]

lab_f3 <- read_excel("FMD_Master_V3.xlsx", sheet = 2)
lab_f3_proc <- lab_f3 %>%
  dplyr::select(-c(Chlorid...5, Magnesium...6, `U-Ketonkörper`,`IL-6`, `Lp(a)`,`HOMA-IR`,Insulin))
lab_f3_proc$`B-Ketonkörper` <- ifelse(lab_f3_proc$`B-Ketonkörper` == "<0,1", '0.09', lab_f3_proc$`B-Ketonkörper`) %>%
  as.numeric()
phys_f3 <- read_excel("FMD_Master_V3.xlsx", sheet = 1, range = "A1:CT41") %>%
  dplyr::select(`Patient-ID`, BMI, `Systole OA re`)
lab_f3_proc <- inner_join(phys_f3, lab_f3_proc, by = "Patient-ID")
colnames(lab_f3_proc) <- gsub(" ","_", colnames(lab_f3_proc)) %>%
  gsub("-","_",.) %>%
  paste0(.,"f3")
colnames(lab_f3_proc)[1] <- "Patient_ID"
keep_col <- apply(lab_f3_proc, 2, function(x) sum(is.na(x)) <= length(x)*30/100)
keep_row <- apply(lab_f3_proc, 1, function(x) sum(is.na(x)) <= length(x)*30/100)
lab_f3_proc <- lab_f3_proc[keep_row, keep_col]

lab_f6 <- read_excel("FMD_Master_V6.xlsx", sheet = 2)
lab_f6_proc <- lab_f6 %>%
  dplyr::select(-c(Chlorid...5, Magnesium...6, `U-Ketonkörper`, `IL-6`, `Lp(a)`,`HOMA-IR`,Insulin))
lab_f6_proc$`B-Ketonkörper` <- ifelse(lab_f6_proc$`B-Ketonkörper` == "<0,1", '0.09', lab_f6_proc$`B-Ketonkörper`) %>%
  as.numeric()
phys_f6 <- read_excel("FMD_Master_V6.xlsx", sheet = 1, range = "A1:CT41") %>%
  dplyr::select(`Patient-ID`, BMI, `Systole OA re`)
lab_f6_proc <- inner_join(phys_f6, lab_f6_proc, by = c("Patient-ID" = "FMD Nr."))
colnames(lab_f6_proc) <- gsub(" ","_", colnames(lab_f6_proc)) %>%
  gsub("-","_",.) %>%
  paste0(.,"f6")
colnames(lab_f6_proc)[1] <- "Patient_ID"
keep_col <- apply(lab_f6_proc, 2, function(x) sum(is.na(x)) <= length(x)*30/100)
keep_row <- apply(lab_f6_proc, 1, function(x) sum(is.na(x)) <= length(x)*30/100)
lab_f6_proc <- lab_f6_proc[keep_row, keep_col]

lab_fu <- read_excel("FMD_Master_FU.xlsx", sheet = 2)
lab_fu_proc <- lab_fu %>%
  dplyr::select(-c(Chlorid...5, Magnesium...6, `U-Ketonkörper`,`IL-6`, `Lp(a)`,`HOMA-IR`,Insulin))
lab_fu_proc$`B-Ketonkörper` <- ifelse(lab_fu_proc$`B-Ketonkörper` == "<0,1", '0.09', lab_fu_proc$`B-Ketonkörper`) %>%
  as.numeric()
phys_fu <- read_excel("FMD_Master_FU.xlsx", sheet = 1, range = "A1:CT41") %>%
  dplyr::select(`Patient-ID`, BMI, `Systole OA re`)
lab_fu_proc <- inner_join(phys_fu, lab_fu_proc, by = "Patient-ID")
colnames(lab_fu_proc) <- gsub(" ","_", colnames(lab_fu_proc)) %>%
  gsub("-","_",.) %>%
  paste0(.,"fu")
colnames(lab_fu_proc)[1] <- "Patient_ID"
keep_col <- apply(lab_fu_proc, 2, function(x) sum(is.na(x)) <= length(x)*30/100)
keep_row <- apply(lab_fu_proc, 1, function(x) sum(is.na(x)) <= length(x)*30/100)
lab_fu_proc <- lab_fu_proc[keep_row, keep_col]

#Filter for variables that exist in all timepoints
n_f0 <- gsub("f0","",colnames(lab_f0_proc)[-1])
n_f3 <- gsub("f3","",colnames(lab_f3_proc)[-1])
n_f6 <- gsub("f6","",colnames(lab_f6_proc)[-1])
n_fu <- gsub("fu","",colnames(lab_fu_proc)[-1])
n_int <- intersect(n_f0,intersect(n_f3, intersect(n_f6,n_fu)))
lab_f0_proc <- lab_f0_proc[,c("Patient_ID",paste0(n_int,"f0"))]
lab_f3_proc <- lab_f3_proc[,c("Patient_ID",paste0(n_int,"f3"))]
lab_f6_proc <- lab_f6_proc[,c("Patient_ID",paste0(n_int,"f6"))]
lab_fu_proc <- lab_fu_proc[,c("Patient_ID",paste0(n_int,"fu"))]

#Filter for patients that exist in all time points
lab_comb <- inner_join(lab_f0_proc, 
                       inner_join(lab_f3_proc, 
                                  inner_join(lab_f6_proc, 
                                             lab_fu_proc, by = "Patient_ID"), by = "Patient_ID"), by = "Patient_ID")
pat_ID <- lab_comb$Patient_ID
lab_comb <- apply(dplyr::select(lab_comb, -Patient_ID), 2, function(x) as.numeric(x)) %>% as.data.frame()

#Impute missing data
lab_comb <- impute.knn(t(lab_comb))$data %>% t() %>% as.data.frame()
lab_comb$Patient_ID <- pat_ID

#Make final variable list
lab_f0_proc <- lab_comb[,intersect(colnames(lab_f0_proc), colnames(lab_comb))] %>%
  column_to_rownames("Patient_ID")
colnames(lab_f0_proc) <- gsub("f0","",colnames(lab_f0_proc))
lab_f3_proc <- lab_comb[,intersect(colnames(lab_f3_proc), colnames(lab_comb))] %>%
  column_to_rownames("Patient_ID")
colnames(lab_f3_proc) <- gsub("f3","",colnames(lab_f3_proc))
lab_f6_proc <- lab_comb[,intersect(colnames(lab_f6_proc), colnames(lab_comb))] %>%
  column_to_rownames("Patient_ID")
colnames(lab_f6_proc) <- gsub("f6","",colnames(lab_f6_proc))
lab_fu_proc <- lab_comb[,intersect(colnames(lab_fu_proc), colnames(lab_comb))] %>%
  column_to_rownames("Patient_ID")
colnames(lab_fu_proc) <- gsub("fu","",colnames(lab_fu_proc))
lab_list <- list(lab_f0_proc, lab_f3_proc, lab_f6_proc, lab_fu_proc)

#Make confounder dfs
comb_proc <- comb_proc %>% column_to_rownames("Patient_ID")
cat_conf <- comb_proc[rownames(lab_f0_proc), c("Gender"
                                               #,"shortacting_Insulin","long_acting_Insulin"
                                               ), drop = F]
con_conf <- comb_proc[rownames(lab_f0_proc), c("Age","Diabetes_duration"), drop = F]
```

# Process metabolomics data

```{r, eval=FALSE}
met <- read_excel("Metabolites_ASulaj_Plasma_All_values_µM.xlsx", sheet = 1, na = c("", "NA")) 

wb <- xlsx::loadWorkbook("Metabolites_ASulaj_Plasma_All_values_µM.xlsx")
sheet1 <- getSheets(wb)[[1]]

rows  <- getRows(sheet1)
cells <- getCells(rows)

styles <- sapply(cells, getCellStyle)

cellColor <- function(style) 
   {
    fg  <- style$getFillForegroundXSSFColor()
    rgb <- tryCatch(fg$getRgb(), error = function(e) NULL)
    rgb <- paste(rgb, collapse = "")
    return(rgb)
   }

c <- sapply(styles, cellColor)

na_colors <- c("6a5acd", "87ceeb")
c_na <- c[c %in% na_colors]

for (i in 1:length(c_na)){
  loc <- names(c_na)[i] %>% strsplit(.,split = "\\.") %>% unlist() %>% as.numeric()
  met[loc[1]-1,loc[2]] <- NA
}

met <- met %>% dplyr::select(-c(1,4,5))
colnames(met) <- gsub(" ","_", colnames(met)) %>%
  gsub("-","_",.)
met$Patient_ID <- ifelse(nchar(met$Patient_ID) == 1, paste0("00",met$Patient_ID), paste0("0", met$Patient_ID))

met_annotation <- read.csv("metabolite_annotation.csv", sep = ";", nrows = 630)
met_annotation$Metabolites <- gsub(" ","_", met_annotation$Metabolites) %>%
  gsub("-","_",.)

to_rm <- setdiff(setdiff(colnames(met), met_annotation$Metabolites),c("Patient_ID","Timepoint"))
met <- dplyr::select(met, -to_rm)

met_v0 <- dplyr::filter(met, Timepoint == "V0") %>%
  dplyr::select(-Timepoint) 
colnames(met_v0) <- paste0(colnames(met_v0),"v0")
colnames(met_v0)[1] <- "Patient_ID"
keep_col <- apply(met_v0, 2, function(x) sum(is.na(x)) <= length(x)*30/100)
keep_row <- apply(met_v0, 1, function(x) sum(is.na(x)) <= length(x)*30/100)
met_v0 <- met_v0[keep_row, keep_col]

met_v3 <- dplyr::filter(met, Timepoint == "V3") %>%
  dplyr::select(-Timepoint) 
colnames(met_v3) <- paste0(colnames(met_v3),"v3")
colnames(met_v3)[1] <- "Patient_ID"
keep_col <- apply(met_v3, 2, function(x) sum(is.na(x)) <= length(x)*30/100)
keep_row <- apply(met_v3, 1, function(x) sum(is.na(x)) <= length(x)*30/100)
met_v3 <- met_v3[keep_row, keep_col]

met_v6 <- dplyr::filter(met, Timepoint == "V6") %>%
  dplyr::select(-Timepoint) 
colnames(met_v6) <- paste0(colnames(met_v6),"v6")
colnames(met_v6)[1] <- "Patient_ID"
keep_col <- apply(met_v6, 2, function(x) sum(is.na(x)) <= length(x)*30/100)
keep_row <- apply(met_v6, 1, function(x) sum(is.na(x)) <= length(x)*30/100)
met_v6 <- met_v6[keep_row, keep_col]

met_fu <- dplyr::filter(met, Timepoint == "FU") %>%
  dplyr::select(-Timepoint) 
colnames(met_fu) <- paste0(colnames(met_fu),"fu")
colnames(met_fu)[1] <- "Patient_ID"
keep_col <- apply(met_fu, 2, function(x) sum(is.na(x)) <= length(x)*30/100)
keep_row <- apply(met_fu, 1, function(x) sum(is.na(x)) <= length(x)*30/100)
met_fu <- met_fu[keep_row, keep_col]

#Filter for variables that exist in all timepoints
n_v0 <- gsub("v0","",colnames(met_v0))[-1]
n_v3 <- gsub("v3","",colnames(met_v3))[-1]
n_v6 <- gsub("v6","",colnames(met_v6))[-1]
n_fu <- gsub("fu","",colnames(met_fu))[-1]
n_int <- intersect(n_v0,intersect(n_v3, intersect(n_v6,n_fu)))
met_v0 <- met_v0[c("Patient_ID", paste0(n_int,"v0"))]
met_v3 <- met_v3[c("Patient_ID", paste0(n_int,"v3"))]
met_v6 <- met_v6[c("Patient_ID", paste0(n_int,"v6"))]
met_fu <- met_fu[c("Patient_ID", paste0(n_int,"fu"))]

#Filter for patients that exist in all time points
met_comb <- inner_join(met_v0, 
                       inner_join(met_v3, 
                                  inner_join(met_v6, 
                                             met_fu, by = "Patient_ID"), by = "Patient_ID"), by = "Patient_ID")
pat_ID <- met_comb$Patient_ID
met_comb <- apply(dplyr::select(met_comb, -Patient_ID), 2, function(x) as.numeric(x)) %>% as.data.frame()

#Impute missing data
met_comb <- impute.knn(t(met_comb))$data %>% t() %>% as.data.frame()
met_comb$Patient_ID <- pat_ID

#Make final variable list
met_v0 <- met_comb[,intersect(colnames(met_v0), colnames(met_comb))] %>%
  column_to_rownames("Patient_ID")
colnames(met_v0) <- gsub("v0","",colnames(met_v0))
met_v3 <- met_comb[,intersect(colnames(met_v3), colnames(met_comb))] %>%
  column_to_rownames("Patient_ID")
colnames(met_v3) <- gsub("v3","",colnames(met_v3))
met_v6 <- met_comb[,intersect(colnames(met_v6), colnames(met_comb))] %>%
  column_to_rownames("Patient_ID")
colnames(met_v6) <- gsub("v6","",colnames(met_v6))
met_fu <- met_comb[,intersect(colnames(met_fu), colnames(met_comb))] %>%
  column_to_rownames("Patient_ID")
colnames(met_fu) <- gsub("fu","",colnames(met_fu))
met_list <- list(met_v0, met_v3, met_v6, met_fu)
saveRDS(met_list, "met_list.rds")
saveRDS(met_v0, "met_v0.rds")
saveRDS(met_v3, "met_v3.rds")
saveRDS(met_v6, "met_v6.rds")
saveRDS(met_fu, "met_fu.rds")
saveRDS(met_annotation, "met_annotation.rds")
```

```{r}
met_annotation <- readRDS("met_annotation.rds")
met_list <- readRDS("met_list.rds")
met_v0 <- readRDS("met_v0.rds")
met_v3 <- readRDS("met_v3.rds")
met_v6 <- readRDS("met_v6.rds")
met_fu <- readRDS("met_fu.rds")
#Make confounder dfs
cat_conf <- comb_proc[rownames(met_v0), c("Gender"
                                          #,"shortacting_Insulin","long_acting_Insulin"
                                          ), drop = F]
con_conf <- comb_proc[rownames(met_v0), c("Age","Diabetes_duration"), drop = F]
```

# Define functions

```{r}
DMA <- function(
  #Funtion to do differential metabolomic analysis
  met_list = NULL, #list of metabolomic profiles in chronological order
  groups = NULL, #named vector of groupping info
  cat_confounders = NULL, #Df of categorical confounders, max. 2 confounders recommended 
  con_confounders = NULL, #Df of continous confounders
  comparisons = c("baseline","longitudinal"), #type of comparisons
  baseline_modes = c("one_vs_rest","pairwise") #type of baseline comparisons
){
  DMARes <- list()
  
  if ("baseline" %in% comparisons){
    met_bl <- met_list[[1]]
    ind <- intersect(rownames(con_confounders), intersect(rownames(cat_confounders), 
                                                          intersect(rownames(met_bl), names(groups))))
    met_bl <- met_bl[ind,]
    groups <- groups[ind]
    cat_confounders_b <- cat_confounders[ind,, drop = F]
    con_confounders_b <- con_confounders[ind,, drop = F]
    #met_bl <- apply(met_bl,2,function(x) (x - mean(x))/sd(x)) %>% t()
    inf <- data.frame(cat_confounders_b, con_confounders_b, cluster = as.character(groups))
    
    if ("pairwise" %in% baseline_modes){
      mod <- model.matrix(~ 0 + ., data = inf)
      
      fit = lmFit(log(t(met_bl)), mod)
      contrast = makeContrasts(cluster1 - cluster2,
                               levels = colnames(coef(fit)))
      tmp <- contrasts.fit(fit, contrast)
      tmp <- eBayes(tmp)
      DMARes[["baseline"]][["pairwise"]] <- topTable(tmp, coef = colnames(contrast), sort.by = "P", n = Inf) 
    }
    
    if ("one_vs_rest" %in% baseline_modes){
      for (i in c("1","2","3")){
        inf_tmp <- inf
        inf_tmp$cluster <- ifelse(inf_tmp$cluster == i, "1", "0")
        mod <- model.matrix(~ 0 + ., data = inf_tmp)
  
        fit = lmFit(log2(t(met_bl)), mod)
        contrast = makeContrasts(cluster1 - cluster0, levels = colnames(coef(fit)))
        tmp <- contrasts.fit(fit, contrast)
        tmp <- eBayes(tmp)
        DMARes[["baseline"]][["one_vs_rest"]][[as.character(i)]] <- topTable(tmp, sort.by = "P", n = Inf) 
      }
    }
  }
  
  if ("longitudinal" %in% comparisons){
    for (i in unique(groups)){
      sp <- names(groups[groups == i])
      ind <- intersect(rownames(con_confounders), intersect(rownames(cat_confounders), 
                                                          intersect(rownames(met_list[[1]]), sp)))
      met_list_i <- lapply(met_list, function(x) x[ind, ]  %>% t())
      cat_confounders_l <- cat_confounders[ind,, drop = F]
      con_confounders_l <- con_confounders[ind,, drop = F]
      inf <- data.frame(cat_confounders_l, con_confounders_l, patients = rownames(cat_confounders_l))
      
      for (j in 2:length(met_list)){
        inf_j <- rbind(inf,inf)
        inf_j$time <- rep(c(0,j), each = nrow(inf)) 
        inf_j$time <- ifelse(inf_j$time == j, "1", "0")
        inf_j <- dplyr::select(inf_j, time, everything())
        mod <- model.matrix(~ 0 + ., data = inf_j)
        met_j <- cbind(met_list_i[[1]], met_list_i[[j]])
        
        fit = lmFit(log2(met_j), mod)
        contrast = makeContrasts(time1 - time0, levels = colnames(coef(fit)))
        tmp <- contrasts.fit(fit, contrast)
        tmp <- eBayes(tmp)
        tt <- topTable(tmp, sort.by = "P", n = Inf)
        DMARes[["longitudinal"]][[as.character(i)]][[as.character(j)]] <- tt
      }
    }
  }
  
  return(DMARes)
}

integrated_clustering <- function(
    variable_list = NULL, #List of df in chronological order, with matched rows & columns  
    cat_confounders = NULL, #Df of categorical confounders, max. 2 confounders recommended 
    con_confounders = NULL, #Df of continous confounders
    distance = "euclidean" #Distance metric
    ##OUTPUT: A pairwise distance matrix showing relationship between samples
){
  
  #Normalize measurements by baseline measurements
  for (i in 1:length(variable_list)){
    if (i == 1){
      next
    }
    variable_list[[i]] <- variable_list[[i]]/variable_list[[1]]
  }
 
  #Remove baseline measurements
  variable_list[[1]] <- NULL
  
  #Clustering
  dist_list <- list()
  for (v in colnames(variable_list[[1]])){
    
    #Extract specific variable and make a df of samples as columns and time points as rows
    v_list <- lapply(variable_list, function(x) return(x[,v]))
    v_df <- do.call(rbind, v_list)
    colnames(v_df) <- rownames(variable_list[[1]])
    
    #Remove counfounding effects
    v_df <- removeBatchEffect(log(v_df), 
                              batch = cat_confounders[,1],
                              #batch2 = cat_confounders[,2],
                              covariates = con_confounders)
    # v_df <- removeBatchEffect(v_df, 
    #                           batch = cat_confounders[,3])
    
    #Compute distance matrix
    dist_list[[v]] <- dist(t(v_df), method = distance, diag=TRUE, upper=TRUE) %>% as.matrix()
    #dimnames(dist_list[[v]]) <- list(rownames(variable_list[[1]]), rownames(variable_list[[1]]))
  }
  
  #Compute an aggregated distance matrix 
  distance_matrix <- Reduce(`+`, dist_list)/length(dist_list)
  
  return(distance_matrix)
}

significant_vars <- function(
    variable_list = NULL, #List of df in chronological order, with matched rows & columns
    cluster = NULL #Samples' cluster information
    ###OUTPUT: new variable list of significant variables
){
  
  #Normalize measurements by baseline measurements
  for (i in 1:length(variable_list)){
    if (i == 1){
      next
    }
    variable_list[[i]] <- variable_list[[i]]/variable_list[[1]]
  }
 
  #Remove baseline measurements
  variable_list[[1]] <- NULL
  
  #Do ANOVA
  for (i in 1:length(variable_list)){
    keep <- c()
    for(v in colnames(variable_list[[i]])){
      df <- data.frame(var = variable_list[[i]][,v], cluster = cluster[rownames(variable_list[[i]])])
      aov <- aov(var ~ cluster, data = df)
      p_val <- summary(aov)[[1]]["cluster","Pr(>F)"]
      
      if (p_val < 0.05){
        keep <- c(keep, v)
      }
    }
    
    variable_list[[i]] <- variable_list[[i]][,keep, drop = F]
  }
  return(variable_list)
}

make_boxplots <- function(
    variable_list = NULL, #List of df in chronological order, with matched rows
    cluster = NULL, #Samples' cluster information
    color_code = NULL 
    ###OUTPUT: print boxplots to output
){
  
  my_comparisons <- list( c("1", "2"), c("1", "3"), c("2", "3") )
  for(i in 1:length(variable_list)){
    for (v in colnames(variable_list[[i]])){
      df <- data.frame(var = variable_list[[i]][,v],cluster = 
                         as.character(cluster[rownames(variable_list[[i]])]))
      p <- ggplot(df, aes(x = cluster, y = var, fill = cluster)) +
        geom_boxplot(width = 0.5) +
        geom_jitter(position = position_jitterdodge(jitter.width = 0.1)) +
        scale_fill_manual(values = color_code)+
        geom_hline(yintercept = 1) +
        ylab(v) +
        ggtitle(paste(v)) +
        theme_bw() +
        stat_compare_means(method = "wilcox.test", comparisons = my_comparisons) +
        theme(legend.position = "none",
              axis.title = element_blank(),
              axis.text.x = element_blank(),
              plot.title = element_text(size = 8, face = "bold"))
      print(p)
    }
  }
  
}
```

# Metabolic profiles of diet groups

```{r}
info <- read_excel("Metabolites_ASulaj_Plasma_All_values_µM.xlsx", sheet = 2)
colnames(info) <- gsub(" ","_", colnames(info)) %>%
  gsub("-","_",.)
info$Patient_ID <- ifelse(nchar(info$Patient_ID) == 1, paste0("00",info$Patient_ID), paste0("0", info$Patient_ID))
participants <- unique(info$Patient_ID)
ind <- match(participants, info$Patient_ID)
diet_groups <- info$Diet_group[ind]
names(diet_groups) <- participants
```

```{r}
DMRes_dietGr <- DMA(met_list = met_list, groups = diet_groups,
             cat_confounders = cat_conf, con_confounders = con_conf,
             comparisons = "longitudinal")
```

## M-Diet

### 3 vs 0

```{r}
dm_m_w3w0 <- DMRes_dietGr$longitudinal$`M-Diet`$`2` %>% rownames_to_column("Name")
ggplot(dm_m_w3w0, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.5, col = "black") +
  geom_point(data = dplyr::filter(dm_m_w3w0, adj.P.Val < 0.05), col = "red") +
  geom_label_repel(data = dplyr::filter(dm_m_w3w0, adj.P.Val < 0.05 & logFC < 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  geom_label_repel(data = dplyr::filter(dm_m_w3w0, adj.P.Val < 0.05 & logFC > 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  theme_classic()
```

### 6 vs 0

```{r}
dm_m_w6w0 <- DMRes_dietGr$longitudinal$`M-Diet`$`3` %>% rownames_to_column("Name")
ggplot(dm_m_w6w0, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.5, col = "black") +
  geom_point(data = dplyr::filter(dm_m_w6w0, adj.P.Val < 0.05), col = "red") +
  geom_label_repel(data = dplyr::filter(dm_m_w6w0, adj.P.Val < 0.05 & logFC < 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  geom_label_repel(data = dplyr::filter(dm_m_w6w0, adj.P.Val < 0.05 & logFC > 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  theme_classic()
```

## FMD

### 3 vs 0

```{r, fig.width=9, fig.height=5}
dm_fmd_w3w0 <- DMRes_dietGr$longitudinal$FMD$`2` %>% rownames_to_column("Name")
dm_fmd_w3w0$Type <- met_annotation$Categories[match(dm_fmd_w3w0$Name,met_annotation$Metabolites)]
met_col <- c(brewer.pal(9,"Set1"),brewer.pal(8, "Set2"),brewer.pal(6,"Set3"))
names(met_col) <- levels(as.factor(dm_fmd_w3w0$Type))
ggplot(dm_fmd_w3w0, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = dplyr::filter(dm_fmd_w3w0, adj.P.Val < 0.05), aes(col = Type), alpha = 0.7) +
  scale_color_manual(values = met_col) +
  geom_point(data = dplyr::filter(dm_fmd_w3w0, adj.P.Val >= 0.05), col = "grey60", alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(dm_fmd_w3w0, adj.P.Val < 0.05 & logFC < 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  geom_label_repel(data = dplyr::filter(dm_fmd_w3w0, adj.P.Val < 0.05 & logFC > 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  theme_classic()
```

### 6 vs 0

```{r, fig.width=9, fig.height=5}
dm_fmd_w6w0 <- DMRes_dietGr$longitudinal$FMD$`3` %>% rownames_to_column("Name")
dm_fmd_w6w0$Type <- met_annotation$Categories[match(dm_fmd_w6w0$Name,met_annotation$Metabolites)]
ggplot(dm_fmd_w6w0, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = dplyr::filter(dm_fmd_w6w0, adj.P.Val < 0.05), aes(col = Type), alpha = 0.7) +
  scale_color_manual(values = met_col) +
  geom_point(data = dplyr::filter(dm_fmd_w6w0, adj.P.Val >= 0.05), col = "grey60", alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(dm_fmd_w6w0, adj.P.Val < 0.05 & logFC < 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  geom_label_repel(data = dplyr::filter(dm_fmd_w6w0, adj.P.Val < 0.05 & logFC > 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  theme_classic()
```

### FU vs 0

```{r, fig.width=9, fig.height=5}
dm_fmd_fuw0 <- DMRes_dietGr$longitudinal$FMD$`4` %>% rownames_to_column("Name")
dm_fmd_fuw0$Type <- met_annotation$Categories[match(dm_fmd_fuw0$Name,met_annotation$Metabolites)]
ggplot(dm_fmd_fuw0, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = dplyr::filter(dm_fmd_fuw0, adj.P.Val < 0.05), aes(col = Type), alpha = 0.7) +
  scale_color_manual(values = met_col) +
  geom_point(data = dplyr::filter(dm_fmd_fuw0, adj.P.Val >= 0.05), col = "grey60", alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(dm_fmd_fuw0, adj.P.Val < 0.05 & logFC < 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  geom_label_repel(data = dplyr::filter(dm_fmd_fuw0, adj.P.Val < 0.05 & logFC > 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  theme_classic()
```

### Summary

```{r, fig.width=4, fig.height=6}
sig_type_fmd <- rbind(data.frame(Type = dm_fmd_w3w0$Type[dm_fmd_w3w0$adj.P.Val<0.05],
                                 t = dm_fmd_w3w0$t[dm_fmd_w3w0$adj.P.Val<0.05],
                                 Month = "Three"),
                      data.frame(Type = dm_fmd_w6w0$Type[dm_fmd_w6w0$adj.P.Val<0.05],
                                 t = dm_fmd_w6w0$t[dm_fmd_w6w0$adj.P.Val<0.05],
                                 Month = "Six"))
sig_type_fmd <- mutate(sig_type_fmd, Direction = ifelse(sig_type_fmd$t < 0, "down","up")) %>%
  dplyr::select(-t)%>%
  group_by(Type, Month, Direction) %>%
  summarise(Count = n())
sig_type_fmd$Count <- ifelse(sig_type_fmd$Direction == "down", sig_type_fmd$Count*-1, sig_type_fmd$Count)
sig_type_fmd$Month <- factor(sig_type_fmd$Month, levels = c("Three","Six"))
# sig_type_fmd$Type <- factor(sig_type_fmd$Type, 
#                             levels = table(sig_type_fmd$Type) %>% sort(decreasing = T) %>% names())
ggplot() +
  geom_bar(data = filter(sig_type_fmd, Direction == "up"), aes(x = Month, y = Count, fill = Type),
           stat = "identity", position = "stack") +
  geom_bar(data = filter(sig_type_fmd, Direction == "down"), aes(x = Month, y = Count, fill = Type),
           stat = "identity", position = "stack") +
  scale_fill_manual(values = met_col) +
  theme_bw() +
  coord_cartesian(ylim = c(-300,50))
```

## FMD - Responders

```{r}
ind <- match(participants, info$Patient_ID)
diet_groups_strat <- info$Subgroup[ind]
names(diet_groups_strat) <- participants
DMRes_dietGr_strat <- DMA(met_list = met_list, groups = diet_groups_strat,
             cat_confounders = cat_conf, con_confounders = con_conf,
             comparisons = "longitudinal")
```

### 3 vs 0

```{r, fig.width=9, fig.height=5}
dm_fmdr_w3w0 <- DMRes_dietGr_strat$longitudinal$FMD_Resp$`2` %>% rownames_to_column("Name")
dm_fmdr_w3w0$Type <- met_annotation$Categories[match(dm_fmdr_w3w0$Name,met_annotation$Metabolites)]
ggplot(dm_fmdr_w3w0, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = dplyr::filter(dm_fmdr_w3w0, adj.P.Val < 0.05), aes(col = Type), alpha = 0.7) +
  scale_color_manual(values = met_col) +
  geom_point(data = dplyr::filter(dm_fmdr_w3w0, adj.P.Val >= 0.05), col = "grey60", alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(dm_fmdr_w3w0, adj.P.Val < 0.05 & logFC < 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  geom_label_repel(data = dplyr::filter(dm_fmdr_w3w0, adj.P.Val < 0.05 & logFC > 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  theme_classic()
```

### 6 vs 0

```{r, fig.width=9, fig.height=5}
dm_fmdr_w6w0 <- DMRes_dietGr_strat$longitudinal$FMD_Resp$`3` %>% rownames_to_column("Name")
dm_fmdr_w6w0$Type <- met_annotation$Categories[match(dm_fmdr_w6w0$Name,met_annotation$Metabolites)]
ggplot(dm_fmdr_w6w0, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = dplyr::filter(dm_fmdr_w6w0, adj.P.Val < 0.05), aes(col = Type), alpha = 0.7) +
  scale_color_manual(values = met_col) +
  geom_point(data = dplyr::filter(dm_fmdr_w6w0, adj.P.Val >= 0.05), col = "grey60", alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(dm_fmdr_w6w0, adj.P.Val < 0.05 & logFC < 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  geom_label_repel(data = dplyr::filter(dm_fmdr_w6w0, adj.P.Val < 0.05 & logFC > 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  theme_classic()
```

### FU vs 0

```{r, fig.width=9, fig.height=5}
dm_fmdr_fuw0 <- DMRes_dietGr_strat$longitudinal$FMD_Resp$`4` %>% rownames_to_column("Name")
dm_fmdr_fuw0$Type <- met_annotation$Categories[match(dm_fmdr_fuw0$Name,met_annotation$Metabolites)]
ggplot(dm_fmdr_fuw0, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = dplyr::filter(dm_fmdr_fuw0, adj.P.Val < 0.05), aes(col = Type), alpha = 0.7) +
  scale_color_manual(values = met_col) +
  geom_point(data = dplyr::filter(dm_fmdr_fuw0, adj.P.Val >= 0.05), col = "grey60", alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(dm_fmdr_fuw0, adj.P.Val < 0.05 & logFC < 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  geom_label_repel(data = dplyr::filter(dm_fmdr_fuw0, adj.P.Val < 0.05 & logFC > 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  theme_classic()
```

### Summary

```{r, fig.width=4, fig.height=6}
sig_type_fmdr <- rbind(data.frame(Type = dm_fmdr_w3w0$Type[dm_fmdr_w3w0$adj.P.Val<0.05],
                                 t = dm_fmdr_w3w0$t[dm_fmdr_w3w0$adj.P.Val<0.05],
                                 Month = "Three"),
                      data.frame(Type = dm_fmdr_w6w0$Type[dm_fmdr_w6w0$adj.P.Val<0.05],
                                 t = dm_fmdr_w6w0$t[dm_fmdr_w6w0$adj.P.Val<0.05],
                                 Month = "Six"))
sig_type_fmdr <- mutate(sig_type_fmdr, Direction = ifelse(sig_type_fmdr$t < 0, "down","up")) %>%
  dplyr::select(-t)%>%
  group_by(Type, Month, Direction) %>%
  summarise(Count = n())
sig_type_fmdr$Count <- ifelse(sig_type_fmdr$Direction == "down", sig_type_fmdr$Count*-1, sig_type_fmdr$Count)
sig_type_fmdr$Month <- factor(sig_type_fmdr$Month, levels = c("Three","Six"))
# sig_type_fmdr$Type <- factor(sig_type_fmdr$Type, 
#                             levels = table(sig_type_fmdr$Type) %>% sort(decreasing = T) %>% names())
ggplot() +
  geom_bar(data = filter(sig_type_fmdr, Direction == "up"), aes(x = Month, y = Count, fill = Type),
           stat = "identity", position = "stack") +
  geom_bar(data = filter(sig_type_fmdr, Direction == "down"), aes(x = Month, y = Count, fill = Type),
           stat = "identity", position = "stack") +
  scale_fill_manual(values = met_col) +
  theme_bw() +
  coord_cartesian(ylim = c(-300,50))
```

## FMD - Non responders

### 3 vs 0

```{r, fig.width=9, fig.height=5}
dm_fmdn_w3w0 <- DMRes_dietGr_strat$longitudinal$FMD_Non_Resp$`2` %>% rownames_to_column("Name")
dm_fmdn_w3w0$Type <- met_annotation$Categories[match(dm_fmdn_w3w0$Name,met_annotation$Metabolites)]
ggplot(dm_fmdn_w3w0, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = dplyr::filter(dm_fmdn_w3w0, adj.P.Val < 0.05), aes(col = Type), alpha = 0.7) +
  scale_color_manual(values = met_col) +
  geom_point(data = dplyr::filter(dm_fmdn_w3w0, adj.P.Val >= 0.05), col = "grey60", alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(dm_fmdn_w3w0, adj.P.Val < 0.05 & logFC < 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  geom_label_repel(data = dplyr::filter(dm_fmdn_w3w0, adj.P.Val < 0.05 & logFC > 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  theme_classic()
```

### 6 vs 0

```{r, fig.width=9, fig.height=5}
dm_fmdn_w6w0 <- DMRes_dietGr_strat$longitudinal$FMD_Non_Resp$`3` %>% rownames_to_column("Name")
dm_fmdn_w6w0$Type <- met_annotation$Categories[match(dm_fmdn_w6w0$Name,met_annotation$Metabolites)]
ggplot(dm_fmdn_w6w0, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = dplyr::filter(dm_fmdn_w6w0, adj.P.Val < 0.05), aes(col = Type), alpha = 0.7) +
  scale_color_manual(values = met_col) +
  geom_point(data = dplyr::filter(dm_fmdn_w6w0, adj.P.Val >= 0.05), col = "grey60", alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(dm_fmdn_w6w0, adj.P.Val < 0.05 & logFC < 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  geom_label_repel(data = dplyr::filter(dm_fmdn_w6w0, adj.P.Val < 0.05 & logFC > 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  theme_classic()
```

### FU vs 0

```{r, fig.width=9, fig.height=5}
dm_fmdn_fuw0 <- DMRes_dietGr_strat$longitudinal$FMD_Non_Resp$`4` %>% rownames_to_column("Name")
dm_fmdn_fuw0$Type <- met_annotation$Categories[match(dm_fmdn_fuw0$Name,met_annotation$Metabolites)]
ggplot(dm_fmdn_fuw0, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = dplyr::filter(dm_fmdn_fuw0, adj.P.Val < 0.05), aes(col = Type), alpha = 0.7) +
  scale_color_manual(values = met_col) +
  geom_point(data = dplyr::filter(dm_fmdn_fuw0, adj.P.Val >= 0.05), col = "grey60", alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(dm_fmdn_fuw0, adj.P.Val < 0.05 & logFC < 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  geom_label_repel(data = dplyr::filter(dm_fmdn_fuw0, adj.P.Val < 0.05 & logFC > 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  theme_classic()
```

### Summary

```{r, fig.width=4, fig.height=6}
sig_type_fmdn <- rbind(data.frame(Type = dm_fmdn_w3w0$Type[dm_fmdn_w3w0$adj.P.Val<0.05],
                                 t = dm_fmdn_w3w0$t[dm_fmdn_w3w0$adj.P.Val<0.05],
                                 Month = "Three"),
                      data.frame(Type = dm_fmdn_w6w0$Type[dm_fmdn_w6w0$adj.P.Val<0.05],
                                 t = dm_fmdn_w6w0$t[dm_fmdn_w6w0$adj.P.Val<0.05],
                                 Month = "Six"))
sig_type_fmdn <- mutate(sig_type_fmdn, Direction = ifelse(sig_type_fmdn$t < 0, "down","up")) %>%
  dplyr::select(-t)%>%
  group_by(Type, Month, Direction) %>%
  summarise(Count = n())
sig_type_fmdn$Count <- ifelse(sig_type_fmdn$Direction == "down", sig_type_fmdn$Count*-1, sig_type_fmdn$Count)
sig_type_fmdn$Month <- factor(sig_type_fmdn$Month, levels = c("Three","Six"))
# sig_type_fmdn$Type <- factor(sig_type_fmdn$Type, 
#                             levels = table(sig_type_fmdn$Type) %>% sort(decreasing = T) %>% names())
ggplot() +
  geom_bar(data = filter(sig_type_fmdn, Direction == "up"), aes(x = Month, y = Count, fill = Type),
           stat = "identity", position = "stack") +
  geom_bar(data = filter(sig_type_fmdn, Direction == "down"), aes(x = Month, y = Count, fill = Type),
           stat = "identity", position = "stack") +
  scale_fill_manual(values = met_col) +
  theme_bw() +
  coord_cartesian(ylim = c(-300,50))
```

# Unsupervised clustering using clinical variables

```{r}
dist_clin <- integrated_clustering(variable_list = lab_list, cat_confounders = cat_conf, 
                                   con_confounders = con_conf, distance = "manhattan")
clRes_c <- ConsensusClusterPlus(d = as.dist(dist_clin), maxK = 6, reps = 500, 
                              pItem = 0.8, seed = 993, clusterAlg = "km")
```

## Boxplots

```{r, fig.width=1.5, fig.height=2.3}
lab_list_sig <- significant_vars(variable_list = lab_list, cluster = clRes_c[[2]]$consensusClass)
colnames(lab_list_sig[[1]]) <- c("Chloride month 3", "CHE month 3", "Tot bili month 3", "Iron month 3",
                      "HbA1c(%) month 3", "HbA1c(mmol) month 3", "TFS month 3", "IGF1 month 3")
colnames(lab_list_sig[[2]]) <- c("Sodium month 6", "Chloride month 6", "Glucose month 6", "Tot bili month 6",
                      "Iron month 6", "TG month 6", "HbA1c(%) month 6", "HbA1c(mmol) month 6",
                      "TFS month 6", "C-peptide month 6")
colnames(lab_list_sig[[3]]) <- c("CHE follow-up", "Tot bili follow-up", 
                      "Iron follow-up", "TFS follow-up", "hsCRP follow-up")
make_boxplots(lab_list_sig, cluster = clRes_c[[2]]$consensusClass, color_code = c("1" = "aquamarine4", 
                                                                                "2" = "cornflowerblue", 
                                                                                "3" = "brown"))
```

## Heatmaps

```{r, fig.width=10, fig.height=8}
for (i in 1:length(lab_list_sig)){
  colnames(lab_list_sig[[i]]) <- paste0(colnames(lab_list_sig[[i]]), "_",i)
}
sig_df <- do.call(cbind, lab_list_sig)
colnames(sig_df) <- c("Chloride month 3", "CHE month 3", "Tot bili month 3", "Iron month 3",
                      "HbA1c(%) month 3", "HbA1c(mmol) month 3", "TFS month 3", "IGF1 month 3", 
                      "Sodium month 6", "Chloride month 6", "Glucose month 6", "Tot bili month 6",
                      "Iron month 6", "TG month 6", "HbA1c(%) month 6", "HbA1c(mmol) month 6",
                      "TFS month 6", "C-peptide month 6", "CHE follow-up", "Tot bili follow-up", 
                      "Iron follow-up", "TFS follow-up", "hsCRP follow-up")
sig_df <- sig_df[,-c(6,16)]

info_fil <- dplyr::filter(info, Patient_ID %in% rownames(sig_df))
info_fil$Cluster <- clRes_c[[2]]$consensusClass[match(info_fil$Patient_ID, names(clRes_c[[2]]$consensusClass))]

ind <- match(rownames(sig_df), info_fil$Patient_ID)
column_an = HeatmapAnnotation(Cluster = info_fil$Cluster[ind], 
                              `Diet groups` = info_fil$Diet_group[ind],
                              `Response groups` = info_fil$Subgroup[ind],
                              col = list(Cluster = c("1" = "aquamarine4", "2" = "cornflowerblue", 
                                                     "3" = "brown"),
                                         `Diet groups` = c("M-Diet" = "burlywood4", "FMD" = "brown4"),
                                         `Response groups` = c("M-Diet" = "burlywood4", "FMD_Resp" = "brown1",
                                                      "FMD_Non_Resp" = "darksalmon")))
cn <- colnames(sig_df)
rsplt <- ifelse(grepl("*month 3",cn), "Month 3", ifelse(grepl("*month 6", cn), "Month 6", "Follow-up"))
col_fun = colorRamp2(c(0, 1, 2.5), c("darkviolet", "white", "orange"))
Heatmap(t(sig_df), top_annotation = column_an, col = col_fun,
            column_split = factor(info_fil$Cluster[ind], levels = c("1","2")), 
            row_split = factor(rsplt, levels = c("Month 3", "Month 6", "Follow-up")),
        row_gap = unit(5, "mm"),
            clustering_distance_columns = "manhattan",
            clustering_method_columns = "ward.D2",
            cluster_rows = F, row_title_rot = 0,
            cluster_columns = F, name = "Baseline-normalized value")
```

# Metabolic profile of the clusters

```{r}
DMRes_cluster <- DMA(met_list = met_list, groups = clRes_c[[2]]$consensusClass,
             cat_confounders = cat_conf, con_confounders = con_conf,
             comparisons = "longitudinal")
```

## Cluster 1

### 3 vs 0

```{r, fig.width=9, fig.height=5}
dm_1_w3w0 <- DMRes_cluster$longitudinal$`1`$`2` %>% rownames_to_column("Name")
dm_1_w3w0$Type <- met_annotation$Categories[match(dm_1_w3w0$Name,met_annotation$Metabolites)]
ggplot(dm_1_w3w0, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = dplyr::filter(dm_1_w3w0, adj.P.Val < 0.05), aes(col = Type), alpha = 0.7) +
  scale_color_manual(values = met_col) +
  geom_point(data = dplyr::filter(dm_1_w3w0, adj.P.Val >= 0.05), col = "grey60", alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(dm_1_w3w0, adj.P.Val < 0.05 & logFC < 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  geom_label_repel(data = dplyr::filter(dm_1_w3w0, adj.P.Val < 0.05 & logFC > 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  theme_classic()
```

### 6 vs 0

```{r, fig.width=9, fig.height=5}
dm_1_w6w0 <- DMRes_cluster$longitudinal$`1`$`3` %>% rownames_to_column("Name")
dm_1_w6w0$Type <- met_annotation$Categories[match(dm_1_w6w0$Name,met_annotation$Metabolites)]
ggplot(dm_1_w6w0, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = dplyr::filter(dm_1_w6w0, adj.P.Val < 0.05), aes(col = Type), alpha = 0.7) +
  scale_color_manual(values = met_col) +
  geom_point(data = dplyr::filter(dm_1_w6w0, adj.P.Val >= 0.05), col = "grey60", alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(dm_1_w6w0, adj.P.Val < 0.05 & logFC < 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  geom_label_repel(data = dplyr::filter(dm_1_w6w0, adj.P.Val < 0.05 & logFC > 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  theme_classic()
```

### FU vs 0

```{r, fig.width=9, fig.height=5}
dm_1_fuw0 <- DMRes_cluster$longitudinal$`1`$`4` %>% rownames_to_column("Name")
dm_1_fuw0$Type <- met_annotation$Categories[match(dm_1_fuw0$Name,met_annotation$Metabolites)]
ggplot(dm_1_fuw0, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = dplyr::filter(dm_1_fuw0, adj.P.Val < 0.05), aes(col = Type), alpha = 0.7) +
  scale_color_manual(values = met_col) +
  geom_point(data = dplyr::filter(dm_1_fuw0, adj.P.Val >= 0.05), col = "grey60", alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(dm_1_fuw0, adj.P.Val < 0.05 & logFC < 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  geom_label_repel(data = dplyr::filter(dm_1_fuw0, adj.P.Val < 0.05 & logFC > 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  theme_classic()
```

### Summary

```{r, fig.width=4, fig.height=6}
sig_type_1 <- rbind(data.frame(Type = dm_1_w3w0$Type[dm_1_w3w0$adj.P.Val<0.05],
                                 t = dm_1_w3w0$t[dm_1_w3w0$adj.P.Val<0.05],
                                 Month = "Three"),
                      data.frame(Type = dm_1_w6w0$Type[dm_1_w6w0$adj.P.Val<0.05],
                                 t = dm_1_w6w0$t[dm_1_w6w0$adj.P.Val<0.05],
                                 Month = "Six"))
sig_type_1 <- mutate(sig_type_1, Direction = ifelse(sig_type_1$t < 0, "down","up")) %>%
  dplyr::select(-t)%>%
  group_by(Type, Month, Direction) %>%
  summarise(Count = n())
sig_type_1$Count <- ifelse(sig_type_1$Direction == "down", sig_type_1$Count*-1, sig_type_1$Count)
sig_type_1$Month <- factor(sig_type_1$Month, levels = c("Three","Six"))
# sig_type_1$Type <- factor(sig_type_1$Type, 
#                             levels = table(sig_type_1$Type) %>% sort(decreasing = T) %>% names())
ggplot() +
  geom_bar(data = filter(sig_type_1, Direction == "up"), aes(x = Month, y = Count, fill = Type),
           stat = "identity", position = "stack") +
  geom_bar(data = filter(sig_type_1, Direction == "down"), aes(x = Month, y = Count, fill = Type),
           stat = "identity", position = "stack") +
  scale_fill_manual(values = met_col) +
  theme_bw() +
  coord_cartesian(ylim = c(-300,50))
```

## Cluster 2

### 3 vs 0

```{r, fig.width=9, fig.height=5}
dm_2_w3w0 <- DMRes_cluster$longitudinal$`2`$`2` %>% rownames_to_column("Name")
dm_2_w3w0$Type <- met_annotation$Categories[match(dm_2_w3w0$Name,met_annotation$Metabolites)]
ggplot(dm_2_w3w0, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = dplyr::filter(dm_2_w3w0, adj.P.Val < 0.05), aes(col = Type), alpha = 0.7) +
  scale_color_manual(values = met_col) +
  geom_point(data = dplyr::filter(dm_2_w3w0, adj.P.Val >= 0.05), col = "grey60", alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(dm_2_w3w0, adj.P.Val < 0.05 & logFC < 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  geom_label_repel(data = dplyr::filter(dm_2_w3w0, adj.P.Val < 0.05 & logFC > 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  theme_classic()
```

### 6 vs 0

```{r, fig.width=9, fig.height=5}
dm_2_w6w0 <- DMRes_cluster$longitudinal$`2`$`3` %>% rownames_to_column("Name")
dm_2_w6w0$Type <- met_annotation$Categories[match(dm_2_w6w0$Name,met_annotation$Metabolites)]
ggplot(dm_2_w6w0, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = dplyr::filter(dm_2_w6w0, adj.P.Val < 0.05), aes(col = Type), alpha = 0.7) +
  scale_color_manual(values = met_col) +
  geom_point(data = dplyr::filter(dm_2_w6w0, adj.P.Val >= 0.05), col = "grey60", alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(dm_2_w6w0, adj.P.Val < 0.05 & logFC < 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  geom_label_repel(data = dplyr::filter(dm_2_w6w0, adj.P.Val < 0.05 & logFC > 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  theme_classic()
```

### FU vs 0

```{r, fig.width=9, fig.height=5}
dm_2_fuw0 <- DMRes_cluster$longitudinal$`2`$`4` %>% rownames_to_column("Name")
dm_2_fuw0$Type <- met_annotation$Categories[match(dm_2_fuw0$Name,met_annotation$Metabolites)]
ggplot(dm_2_fuw0, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = dplyr::filter(dm_2_fuw0, adj.P.Val < 0.05), aes(col = Type), alpha = 0.7) +
  scale_color_manual(values = met_col) +
  geom_point(data = dplyr::filter(dm_2_fuw0, adj.P.Val >= 0.05), col = "grey60", alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(dm_2_fuw0, adj.P.Val < 0.05 & logFC < 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  geom_label_repel(data = dplyr::filter(dm_2_fuw0, adj.P.Val < 0.05 & logFC > 0) %>% head(n=10),
                  aes(label = Name), max.overlaps = 30, alpha = 0.6) +
  theme_classic()
```

### Summary

```{r, fig.width=3.5, fig.height=6}
sig_type_2 <- data.frame(Type = dm_2_w6w0$Type[dm_2_w6w0$adj.P.Val<0.05],
                                 t = dm_2_w6w0$t[dm_2_w6w0$adj.P.Val<0.05],
                                 Month = "Six")
sig_type_2 <- mutate(sig_type_2, Direction = ifelse(sig_type_2$t < 0, "down","up")) %>%
  dplyr::select(-t)%>%
  group_by(Type, Month, Direction) %>%
  summarise(Count = n())
sig_type_2$Count <- ifelse(sig_type_2$Direction == "down", sig_type_2$Count*-1, sig_type_2$Count)
sig_type_2$Month <- factor(sig_type_2$Month, levels = c("Three","Six"))
# sig_type_2$Type <- factor(sig_type_2$Type, 
#                             levels = table(sig_type_2$Type) %>% sort(decreasing = T) %>% names())
ggplot() +
  geom_bar(data = filter(sig_type_2, Direction == "up"), aes(x = Month, y = Count, fill = Type),
           stat = "identity", position = "stack") +
  geom_bar(data = filter(sig_type_2, Direction == "down"), aes(x = Month, y = Count, fill = Type),
           stat = "identity", position = "stack") +
  scale_fill_manual(values = met_col) +
  theme_bw() +
  coord_cartesian(ylim = c(-300,50))
```

# Metabolite set enrichment analysis

## FMD - Responders

### 3 vs 0

```{r}
pathway_list <- readRDS("pathway_list.rds")
pathway_list <- pathway_list[grep("reactome", names(pathway_list))]
ranklist <- dm_fmdr_w3w0$t
names(ranklist) <- dm_fmdr_w3w0$Name
ranklist <- sort(ranklist)

set.seed(993)
fgseaRes <- fgsea(pathways = pathway_list,
                stats    = ranklist,
                minSize  = 5,
                maxSize  = 200) %>% arrange(padj)
fgseaRes$leadingEdge <- lapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ",")) %>% unlist()
fgseaRes$pathway <- gsub("reactome:","",fgseaRes$pathway)
saveRDS(fgseaRes, "fgsea_fmdr_m3_v_m0.rds")
write.csv(fgseaRes,"fgsea_fmdr_m3_v_m0.csv", quote = F, row.names = F)
```

```{r, fig.width=7.5, fig.height=2.5}
sigpw <- fgseaRes$pathway[fgseaRes$padj < 0.05]
fgseaRes <- fgseaRes[-c(4,5,11,24),]
ggplot(head(fgseaRes[fgseaRes$pathway %in% sigpw,],10), aes(x = NES, y = reorder(pathway, NES), col = padj))+
  geom_point(aes(size = size)) +
  scale_size_continuous(breaks = c(25,50,75,100), limits = c(0,100)) +
  scale_color_continuous(breaks = c(0.01,0.02,0.03,0.04,0.05), limits = c(0,0.05)) +
  xlab("Enrichment score")+
  ylab("Pathways") +
  theme_bw()
```

### 6 vs 0

```{r}
ranklist <- dm_fmdr_w6w0$t
names(ranklist) <- dm_fmdr_w6w0$Name
ranklist <- sort(ranklist)

set.seed(993)
fgseaRes <- fgsea(pathways = pathway_list,
                stats    = ranklist,
                minSize  = 5,
                maxSize  = 200) %>% arrange(padj)
fgseaRes$leadingEdge <- lapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ",")) %>% unlist()
fgseaRes$pathway <- gsub("reactome:","",fgseaRes$pathway)
saveRDS(fgseaRes, "fgsea_fmdr_m6_v_m0.rds")
write.csv(fgseaRes,"fgsea_fmdr_m6_v_m0.csv", quote = F, row.names = F)
```

```{r, fig.width=7.5, fig.height=2.5}
sigpw <- fgseaRes$pathway[fgseaRes$padj < 0.05]
fgseaRes <- fgseaRes[-c(10,43,57,63,72,77,79),]
ggplot(head(fgseaRes[fgseaRes$pathway %in% sigpw,],10), aes(x = NES, y = reorder(pathway, NES), col = padj))+
  geom_point(aes(size = size)) +
  scale_size_continuous(breaks = c(25,50,75,100), limits = c(0,100)) +
  scale_color_continuous(breaks = c(0.01,0.02,0.03,0.04,0.05), limits = c(0,0.05)) +
  xlab("Enrichment score")+
  ylab("Pathways") +
  theme_bw()
```

## FMD - Non Responders

### 3 vs 0

```{r}
ranklist <- dm_fmdn_w3w0$t
names(ranklist) <- dm_fmdn_w3w0$Name
ranklist <- sort(ranklist)

set.seed(993)
fgseaRes <- fgsea(pathways = pathway_list,
                stats    = ranklist,
                minSize  = 5,
                maxSize  = 200) %>% arrange(padj)
fgseaRes$leadingEdge <- lapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ",")) %>% unlist()
fgseaRes$pathway <- gsub("reactome:","",fgseaRes$pathway)
saveRDS(fgseaRes, "fgsea_fmdn_m3_v_m0.rds")
write.csv(fgseaRes,"fgsea_fmdn_m3_v_m0.csv", quote = F, row.names = F)
```

```{r, fig.width=6, fig.height=2.5}
sigpw <- fgseaRes$pathway[fgseaRes$padj < 0.05]
fgseaRes <- fgseaRes[-c(6,11,26),]
ggplot(head(fgseaRes[fgseaRes$pathway %in% sigpw,],10), aes(x = NES, y = reorder(pathway, NES), col = padj))+
  geom_point(aes(size = size)) +
  scale_size_continuous(breaks = c(25,50,75,100), limits = c(0,100)) +
  scale_color_continuous(breaks = c(0.01,0.02,0.03,0.04,0.05), limits = c(0,0.05)) +
  xlab("Enrichment score")+
  ylab("Pathways") +
  theme_bw()
```

### 6 vs 0

```{r}
ranklist <- dm_fmdn_w6w0$t
names(ranklist) <- dm_fmdn_w6w0$Name
ranklist <- sort(ranklist)

set.seed(993)
fgseaRes <- fgsea(pathways = pathway_list,
                stats    = ranklist,
                minSize  = 5,
                maxSize  = 200) %>% arrange(padj)
fgseaRes$leadingEdge <- lapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ",")) %>% unlist()
fgseaRes$pathway <- gsub("reactome:","",fgseaRes$pathway)
saveRDS(fgseaRes, "fgsea_fmdn_m6_v_m0.rds")
write.csv(fgseaRes,"fgsea_fmdn_m6_v_m0.csv", quote = F, row.names = F)
```

```{r, fig.width=7.5, fig.height=2.5}
sigpw <- fgseaRes$pathway[fgseaRes$padj < 0.05]
fgseaRes <- fgseaRes[-c(1,2,26,30,3165,79,86,97,98),]
ggplot(head(fgseaRes[fgseaRes$pathway %in% sigpw,],10), aes(x = NES, y = reorder(pathway, NES), col = padj))+
  geom_point(aes(size = size)) +
  scale_size_continuous(breaks = c(25,50,75,100), limits = c(0,100)) +
  scale_color_continuous(breaks = c(0.01,0.02,0.03,0.04,0.05), limits = c(0,0.05)) +
  xlab("Enrichment score")+
  ylab("Pathways") +
  theme_bw()
```

## Cluster 1

### 3 vs 0

```{r}
ranklist <- dm_1_w3w0$t
names(ranklist) <- dm_1_w3w0$Name
ranklist <- sort(ranklist)

set.seed(993)
fgseaRes <- fgsea(pathways = pathway_list,
                stats    = ranklist,
                minSize  = 5,
                maxSize  = 200) %>% arrange(padj)
fgseaRes$leadingEdge <- lapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ",")) %>% unlist()
fgseaRes$pathway <- gsub("reactome:","",fgseaRes$pathway)
saveRDS(fgseaRes, "fgsea_cl1_m3_v_m0.rds")
write.csv(fgseaRes,"fgsea_cl1_m3_v_m0.csv", quote = F, row.names = F)
```

```{r, fig.width=7.5, fig.height=15}
sigpw <- fgseaRes$pathway[fgseaRes$padj < 0.05]
fgseaRes <- fgseaRes[-c(2,13,24,37,51,65,72),]
ggplot(fgseaRes[fgseaRes$pathway %in% sigpw,], aes(x = NES, y = reorder(pathway, NES), col = padj))+
  geom_point(aes(size = size)) +
  scale_size_continuous(breaks = c(25,50,75,100), limits = c(0,100)) +
  scale_color_continuous(breaks = c(0.01,0.02,0.03,0.04,0.05), limits = c(0,0.05)) +
  xlab("Enrichment score")+
  ylab("Pathways") +
  theme_bw()
```

### 6 vs 0

```{r}
ranklist <- dm_1_w6w0$t
names(ranklist) <- dm_1_w6w0$Name
ranklist <- sort(ranklist)

set.seed(993)
fgseaRes <- fgsea(pathways = pathway_list,
                stats    = ranklist,
                minSize  = 5,
                maxSize  = 200) %>% arrange(padj)
fgseaRes$leadingEdge <- lapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ",")) %>% unlist()
fgseaRes$pathway <- gsub("reactome:","",fgseaRes$pathway)
saveRDS(fgseaRes, "fgsea_cl1_m6_v_m0.rds")
write.csv(fgseaRes,"fgsea_cl1_m6_v_m0.csv", quote = F, row.names = F)
```

```{r, fig.width=7.5, fig.height=18}
sigpw <- fgseaRes$pathway[fgseaRes$padj < 0.05]
fgseaRes <- fgseaRes[-c(1,2,26,30,3165,79,86,97,98),]
ggplot(fgseaRes[fgseaRes$pathway %in% sigpw,], aes(x = NES, y = reorder(pathway, NES), col = padj))+
  geom_point(aes(size = size)) +
  scale_size_continuous(breaks = c(25,50,75,100), limits = c(0,100)) +
  scale_color_continuous(breaks = c(0.01,0.02,0.03,0.04,0.05), limits = c(0,0.05)) +
  xlab("Enrichment score")+
  ylab("Pathways") +
  theme_bw()
```

## Cluster 2

### 3 vs 0

```{r}
ranklist <- dm_2_w3w0$t
names(ranklist) <- dm_2_w3w0$Name
ranklist <- sort(ranklist)

set.seed(993)
fgseaRes <- fgsea(pathways = pathway_list,
                stats    = ranklist,
                minSize  = 5,
                maxSize  = 200) %>% arrange(padj)
fgseaRes$leadingEdge <- lapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ",")) %>% unlist()
fgseaRes$pathway <- gsub("reactome:","",fgseaRes$pathway)
saveRDS(fgseaRes, "fgsea_cl2_m3_v_m0.rds")
write.csv(fgseaRes,"fgsea_cl2_m3_v_m0.csv", quote = F, row.names = F)
```

```{r, fig.width=7.5, fig.height=5}
sigpw <- fgseaRes$pathway[fgseaRes$padj < 0.05]
fgseaRes <- fgseaRes[-c(6,11,26),]
ggplot(fgseaRes[fgseaRes$pathway %in% sigpw,], aes(x = NES, y = reorder(pathway, NES), col = padj))+
  geom_point(aes(size = size)) +
  scale_size_continuous(breaks = c(25,50,75,100), limits = c(0,100)) +
  scale_color_continuous(breaks = c(0.01,0.02,0.03,0.04,0.05), limits = c(0,0.05)) +
  xlab("Enrichment score")+
  ylab("Pathways") +
  theme_bw()
```

### 6 vs 0

```{r}
ranklist <- dm_2_w6w0$t
names(ranklist) <- dm_2_w6w0$Name
ranklist <- sort(ranklist)

set.seed(993)
fgseaRes <- fgsea(pathways = pathway_list,
                stats    = ranklist,
                minSize  = 5,
                maxSize  = 200) %>% arrange(padj)
fgseaRes$leadingEdge <- lapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ",")) %>% unlist()
fgseaRes$pathway <- gsub("reactome:","",fgseaRes$pathway)
saveRDS(fgseaRes, "fgsea_cl2_m6_v_m0.rds")
write.csv(fgseaRes,"fgsea_cl2_m6_v_m0.csv", quote = F, row.names = F)
```

```{r, fig.width=7.5, fig.height=18}
sigpw <- fgseaRes$pathway[fgseaRes$padj < 0.05]
fgseaRes <- fgseaRes[-c(1,2,26,30,3165,79,86,97,98),]
ggplot(fgseaRes[fgseaRes$pathway %in% sigpw,], aes(x = NES, y = reorder(pathway, NES), col = padj))+
  geom_point(aes(size = size)) +
  scale_size_continuous(breaks = c(25,50,75,100), limits = c(0,100)) +
  scale_color_continuous(breaks = c(0.01,0.02,0.03,0.04,0.05), limits = c(0,0.05)) +
  xlab("Enrichment score")+
  ylab("Pathways") +
  theme_bw()
```

# Trajectory

## MD

```{r, fig.width=3, fig.height=2}
dm_m_fuw0 <- DMRes_dietGr$longitudinal$`M-Diet`$`4` %>% rownames_to_column("Name")
sig_met_m <- "TrpBetaine"
ind1 <- match(sig_met_m, dm_m_w3w0$Name)
ind2 <- match(sig_met_m, dm_m_w6w0$Name)
ind3 <- match(sig_met_m, dm_m_fuw0$Name)

sig_met_m_df <- data.frame(Name = rep(sig_met_m, times = 4), Time = rep(c("0","3","6","9"), 
                                                                        each = length(sig_met_m)),
                           
                           logFC = c(rep(0,length(sig_met_m)),
                                     dm_m_w3w0$logFC[ind1],
                                     dm_m_w6w0$logFC[ind2],
                                     dm_m_fuw0$logFC[ind3]),
                           adjP = c(rep(1,length(sig_met_m)),
                                  dm_m_w3w0$adj.P.Val[ind1],
                                  dm_m_w6w0$adj.P.Val[ind2],
                                  dm_m_fuw0$adj.P.Val[ind3])) %>%
  mutate(Significance = ifelse(.$adjP < 0.05, "Yes", "No"))

ggplot(sig_met_m_df, aes(x = Time, y = logFC)) + 
  geom_point(alpha = 0.3) + 
  geom_line(aes(group = Name), alpha = 0.3) +
  geom_point(data = dplyr::filter(sig_met_m_df, Significance == "Yes"), alpha = 0.3, col = "red") +
  theme_classic()
```

## FMD

```{r, fig.width=3, fig.height=2}
dm_fmd_fuw0 <- DMRes_dietGr$longitudinal$FMD$`4` %>% rownames_to_column("Name")
sig_met_fmd <- unique(c(dm_fmd_w3w0$Name[dm_fmd_w3w0$adj.P.Val < 0.05],
                      dm_fmd_w6w0$Name[dm_fmd_w6w0$adj.P.Val < 0.05],
                      dm_fmd_fuw0$Name[dm_fmd_fuw0$adj.P.Val < 0.05]))
ind1 <- match(sig_met_fmd, dm_fmd_w3w0$Name)
ind2 <- match(sig_met_fmd, dm_fmd_w6w0$Name)
ind3 <- match(sig_met_fmd, dm_fmd_fuw0$Name)

sig_met_fmd_df <- data.frame(Name = rep(sig_met_fmd, times = 4), Time = rep(c("0","3","6","9"), 
                                                                        each = length(sig_met_fmd)),
                           
                           logFC = c(rep(0,length(sig_met_fmd)),
                                     dm_fmd_w3w0$logFC[ind1],
                                     dm_fmd_w6w0$logFC[ind2],
                                     dm_fmd_fuw0$logFC[ind3]),
                           adjP = c(rep(1,length(sig_met_fmd)),
                                  dm_fmd_w3w0$adj.P.Val[ind1],
                                  dm_fmd_w6w0$adj.P.Val[ind2],
                                  dm_fmd_fuw0$adj.P.Val[ind3])) %>%
  mutate(Significance = ifelse(.$adjP < 0.05, "Yes", "No"))

ggplot(sig_met_fmd_df, aes(x = Time, y = logFC)) + 
  geom_point(alpha = 0.3) + 
  geom_line(aes(group = Name), alpha = 0.3) +
  geom_point(data = dplyr::filter(sig_met_fmd_df, Significance == "Yes"), alpha = 0.3, col = "red") +
  theme_classic()
```

## FMD - Responders

```{r, fig.width=3, fig.height=2}
dm_fmdr_fuw0 <- DMRes_dietGr_strat$longitudinal$FMD_Resp$`4` %>% rownames_to_column("Name")
sig_met_fmdr <- unique(c(dm_fmdr_w3w0$Name[dm_fmdr_w3w0$adj.P.Val < 0.05],
                      dm_fmdr_w6w0$Name[dm_fmdr_w6w0$adj.P.Val < 0.05],
                      dm_fmdr_fuw0$Name[dm_fmdr_fuw0$adj.P.Val < 0.05]))
ind1 <- match(sig_met_fmdr, dm_fmdr_w3w0$Name)
ind2 <- match(sig_met_fmdr, dm_fmdr_w6w0$Name)
ind3 <- match(sig_met_fmdr, dm_fmdr_fuw0$Name)

sig_met_fmdr_df <- data.frame(Name = rep(sig_met_fmdr, times = 4), Time = rep(c("0","3","6","9"), 
                                                                        each = length(sig_met_fmdr)),
                           
                           logFC = c(rep(0,length(sig_met_fmdr)),
                                     dm_fmdr_w3w0$logFC[ind1],
                                     dm_fmdr_w6w0$logFC[ind2],
                                     dm_fmdr_fuw0$logFC[ind3]),
                           adjP = c(rep(1,length(sig_met_fmdr)),
                                  dm_fmdr_w3w0$adj.P.Val[ind1],
                                  dm_fmdr_w6w0$adj.P.Val[ind2],
                                  dm_fmdr_fuw0$adj.P.Val[ind3])) %>%
  mutate(Significance = ifelse(.$adjP < 0.05, "Yes", "No"))

ggplot(sig_met_fmdr_df, aes(x = Time, y = logFC)) + 
  geom_point(alpha = 0.3) + 
  geom_line(aes(group = Name), alpha = 0.3) +
  geom_point(data = dplyr::filter(sig_met_fmdr_df, Significance == "Yes"), alpha = 0.3, col = "red") +
  theme_classic()
```

## FMD - Non Responders

```{r, fig.width=3, fig.height=2}
dm_fmdn_fuw0 <- DMRes_dietGr_strat$longitudinal$FMD_Resp$`4` %>% rownames_to_column("Name")
sig_met_fmdn <- unique(c(dm_fmdn_w3w0$Name[dm_fmdn_w3w0$adj.P.Val < 0.05],
                      dm_fmdn_w6w0$Name[dm_fmdn_w6w0$adj.P.Val < 0.05],
                      dm_fmdn_fuw0$Name[dm_fmdn_fuw0$adj.P.Val < 0.05]))
ind1 <- match(sig_met_fmdn, dm_fmdn_w3w0$Name)
ind2 <- match(sig_met_fmdn, dm_fmdn_w6w0$Name)
ind3 <- match(sig_met_fmdn, dm_fmdn_fuw0$Name)

sig_met_fmdn_df <- data.frame(Name = rep(sig_met_fmdn, times = 4), Time = rep(c("0","3","6","9"), 
                                                                        each = length(sig_met_fmdn)),
                           
                           logFC = c(rep(0,length(sig_met_fmdn)),
                                     dm_fmdn_w3w0$logFC[ind1],
                                     dm_fmdn_w6w0$logFC[ind2],
                                     dm_fmdn_fuw0$logFC[ind3]),
                           adjP = c(rep(1,length(sig_met_fmdn)),
                                  dm_fmdn_w3w0$adj.P.Val[ind1],
                                  dm_fmdn_w6w0$adj.P.Val[ind2],
                                  dm_fmdn_fuw0$adj.P.Val[ind3])) %>%
  mutate(Significance = ifelse(.$adjP < 0.05, "Yes", "No"))

ggplot(sig_met_fmdn_df, aes(x = Time, y = logFC)) + 
  geom_point(alpha = 0.3) + 
  geom_line(aes(group = Name), alpha = 0.3) +
  geom_point(data = dplyr::filter(sig_met_fmdn_df, Significance == "Yes"), alpha = 0.3, col = "red") +
  theme_classic()
```

# Extra plots

```{r, fig.width=8, fig.height=3}
#met_v0_norm <- apply(met_v0, 2, function(x) (x - mean(x))/sd(x))
met_v0_md <- met_v0[intersect(names(diet_groups[diet_groups == "M-Diet"]), rownames(met_v0)),
                    intersect(met_annotation$Metabolites[met_annotation$Categories == "Amino acids"],
                              colnames(met_v0))] %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Metabolites", values_to = "Normalized abundance") %>%
  mutate(Time = "Baseline", Diet = "MD", 
         Category = met_annotation$Subcategories[match(.$Metabolites, met_annotation$Metabolites)])

met_v3_md <- met_v3[intersect(names(diet_groups[diet_groups == "M-Diet"]), rownames(met_v3)),
                    intersect(met_annotation$Metabolites[met_annotation$Categories == "Amino acids"],
                              colnames(met_v3))] %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Metabolites", values_to = "Normalized abundance") %>%
  mutate(Time = "Month 3", Diet = "MD", 
         Category = met_annotation$Subcategories[match(.$Metabolites, met_annotation$Metabolites)])
met_v3_md$`Normalized abundance` <- met_v3_md$`Normalized abundance`/met_v0_md$`Normalized abundance`

met_v6_md <- met_v6[intersect(names(diet_groups[diet_groups == "M-Diet"]), rownames(met_v6)),
                    intersect(met_annotation$Metabolites[met_annotation$Categories == "Amino acids"],
                              colnames(met_v6))] %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Metabolites", values_to = "Normalized abundance") %>%
  mutate(Time = "Month 6", Diet = "MD", 
         Category = met_annotation$Subcategories[match(.$Metabolites, met_annotation$Metabolites)])
met_v6_md$`Normalized abundance` <- met_v6_md$`Normalized abundance`/met_v0_md$`Normalized abundance`

met_v0_fmd <- met_v0[intersect(names(diet_groups[diet_groups == "FMD"]), rownames(met_v0)),
                    intersect(met_annotation$Metabolites[met_annotation$Categories == "Amino acids"],
                              colnames(met_v0))] %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Metabolites", values_to = "Normalized abundance") %>%
  mutate(Time = "Baseline", Diet = "FMD", 
         Category = met_annotation$Subcategories[match(.$Metabolites, met_annotation$Metabolites)])

met_v3_fmd <- met_v3[intersect(names(diet_groups[diet_groups == "FMD"]), rownames(met_v3)),
                    intersect(met_annotation$Metabolites[met_annotation$Categories == "Amino acids"],
                              colnames(met_v3))] %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Metabolites", values_to = "Normalized abundance") %>%
  mutate(Time = "Month 3", Diet = "FMD", 
         Category = met_annotation$Subcategories[match(.$Metabolites, met_annotation$Metabolites)])
met_v3_fmd$`Normalized abundance` <- met_v3_fmd$`Normalized abundance`/met_v0_fmd$`Normalized abundance`

met_v6_fmd <- met_v6[intersect(names(diet_groups[diet_groups == "FMD"]), rownames(met_v6)),
                    intersect(met_annotation$Metabolites[met_annotation$Categories == "Amino acids"],
                              colnames(met_v6))] %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Metabolites", values_to = "Normalized abundance") %>%
  mutate(Time = "Month 6", Diet = "FMD", 
         Category = met_annotation$Subcategories[match(.$Metabolites, met_annotation$Metabolites)])
met_v6_fmd$`Normalized abundance` <- met_v6_fmd$`Normalized abundance`/met_v0_fmd$`Normalized abundance`

met_all <- rbind(met_v3_md, met_v6_md, met_v3_fmd, met_v6_fmd)

ggplot(filter(met_all, Time == "Month 3"), aes(x = Metabolites, y = `Normalized abundance`, fill = Diet)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  facet_grid(~ Category, scales = "free", space="free") +
  scale_fill_manual(values = list("MD" = "burlywood4", "FMD" = "brown4")) +
  theme_bw() +
  geom_signif(comparisons = list(c("FMD", "MD")), 
              map_signif_level=TRUE)

ggplot(filter(met_all, Time == "Month 6"), aes(x = Metabolites, y = `Normalized abundance`, fill = Diet)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  facet_grid(~ Category, scales = "free", space="free") +
  scale_fill_manual(values = list("MD" = "burlywood4", "FMD" = "brown4")) +
  theme_bw() +
  geom_signif(comparisons = list(c("FMD", "MD")), 
              map_signif_level=TRUE)
```

```{r, fig.width=5, fig.height=2}
#met_v0_norm <- apply(met_v0, 2, function(x) (x - mean(x))/sd(x))
met_v0_md <- met_v0[intersect(names(diet_groups[diet_groups == "M-Diet"]), rownames(met_v0)),
                    intersect(met_annotation$Metabolites[met_annotation$Categories == "Acylcarnitines"],
                              colnames(met_v0))] %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Metabolites", values_to = "Normalized abundance") %>%
  mutate(Time = "Baseline", Diet = "MD", 
         Category = met_annotation$Subcategories[match(.$Metabolites, met_annotation$Metabolites)])

met_v3_md <- met_v3[intersect(names(diet_groups[diet_groups == "M-Diet"]), rownames(met_v3)),
                    intersect(met_annotation$Metabolites[met_annotation$Categories == "Acylcarnitines"],
                              colnames(met_v3))] %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Metabolites", values_to = "Normalized abundance") %>%
  mutate(Time = "Month 3", Diet = "MD", 
         Category = met_annotation$Subcategories[match(.$Metabolites, met_annotation$Metabolites)])
met_v3_md$`Normalized abundance` <- met_v3_md$`Normalized abundance`/met_v0_md$`Normalized abundance`

met_v6_md <- met_v6[intersect(names(diet_groups[diet_groups == "M-Diet"]), rownames(met_v6)),
                    intersect(met_annotation$Metabolites[met_annotation$Categories == "Acylcarnitines"],
                              colnames(met_v6))] %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Metabolites", values_to = "Normalized abundance") %>%
  mutate(Time = "Month 6", Diet = "MD", 
         Category = met_annotation$Subcategories[match(.$Metabolites, met_annotation$Metabolites)])
met_v6_md$`Normalized abundance` <- met_v6_md$`Normalized abundance`/met_v0_md$`Normalized abundance`

met_v0_fmd <- met_v0[intersect(names(diet_groups[diet_groups == "FMD"]), rownames(met_v0)),
                    intersect(met_annotation$Metabolites[met_annotation$Categories == "Acylcarnitines"],
                              colnames(met_v0))] %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Metabolites", values_to = "Normalized abundance") %>%
  mutate(Time = "Baseline", Diet = "FMD", 
         Category = met_annotation$Subcategories[match(.$Metabolites, met_annotation$Metabolites)])

met_v3_fmd <- met_v3[intersect(names(diet_groups[diet_groups == "FMD"]), rownames(met_v3)),
                    intersect(met_annotation$Metabolites[met_annotation$Categories == "Acylcarnitines"],
                              colnames(met_v3))] %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Metabolites", values_to = "Normalized abundance") %>%
  mutate(Time = "Month 3", Diet = "FMD", 
         Category = met_annotation$Subcategories[match(.$Metabolites, met_annotation$Metabolites)])
met_v3_fmd$`Normalized abundance` <- met_v3_fmd$`Normalized abundance`/met_v0_fmd$`Normalized abundance`

met_v6_fmd <- met_v6[intersect(names(diet_groups[diet_groups == "FMD"]), rownames(met_v6)),
                    intersect(met_annotation$Metabolites[met_annotation$Categories == "Acylcarnitines"],
                              colnames(met_v6))] %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Metabolites", values_to = "Normalized abundance") %>%
  mutate(Time = "Month 6", Diet = "FMD", 
         Category = met_annotation$Subcategories[match(.$Metabolites, met_annotation$Metabolites)])
met_v6_fmd$`Normalized abundance` <- met_v6_fmd$`Normalized abundance`/met_v0_fmd$`Normalized abundance`

met_all <- rbind(met_v3_md, met_v6_md, met_v3_fmd, met_v6_fmd)

ggplot(filter(met_all, Time == "Month 3"), aes(x = Metabolites, y = `Normalized abundance`, fill = Diet)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  facet_grid(~ Category, scales = "free", space="free") +
  scale_fill_manual(values = list("MD" = "burlywood4", "FMD" = "brown4")) +
  theme_bw() +
  geom_signif(comparisons = list(c("FMD", "MD")), 
              map_signif_level=TRUE)

ggplot(filter(met_all, Time == "Month 6"), aes(x = Metabolites, y = `Normalized abundance`, fill = Diet)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  facet_grid(~ Category, scales = "free", space="free") +
  scale_fill_manual(values = list("MD" = "burlywood4", "FMD" = "brown4")) +
  theme_bw() +
  geom_signif(comparisons = list(c("FMD", "MD")), 
              map_signif_level=TRUE)
```