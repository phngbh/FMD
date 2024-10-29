#!/usr/bin/env Rscript

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

