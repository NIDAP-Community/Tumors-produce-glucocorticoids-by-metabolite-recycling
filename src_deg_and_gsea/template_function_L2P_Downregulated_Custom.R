# L2P Analysis for Single Comparisons [CCBR] (ff44bcce-10c1-4147-a089-79cadc4ea68a): v103
L2P_Downregulated_Custom <- function(DEG_Analysis,Custom_genes_L2P) {
    
    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    libs <- c("l2p","l2psupp","dplyr","magrittr","ggplot2","stringr","RCurl")
    nidapLoadPackages(libs)
    

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    #Basic Parameters:
    deg_table <- DEG_Analysis 
    gene_name_column <- "Gene"
    species = "Mouse"
    collections_to_include <- c("H")
    select_by_rank <- FALSE
    choose_upregulated_genes = FALSE
    plot_bubble <- TRUE

    #Genelist selected by t-statistic rank parameters:
    column_used_to_rank_genes <- c()[1]
    select_top_percentage_of_genes <- TRUE
    select_top_genes <- 500

    #Genelist selected by fold-change and pval thresholds parameters:
    significance_column <- c("KO-WT_pval")[1]
    significance_threshold <- 0.1
    fold_change_column <- c("KO-WT_FC")[1]
    fold_change_threshold <- 1.2
    minimum_number_of_DEG_genes <- 100
    
    #Visual Parameters:
    number_of_pathways_to_plot <- 12 
    plot_top_pathways <- TRUE
    pathways_to_use <- c()
    pathway_axis_label_font_size <- 10
    pathway_axis_label_max_length <- 50

    #For Bubble Plot:
    sort_bubble_plot_by <- "percent gene hits per pathway"
    plot_bubble_size <- "number of hits"  
    plot_bubble_color <- "Fisher's Exact pval"
    bubble_colors <- "blues"
    

    #For Bar Plot:
    use_fdr_pvals <- TRUE
    color_for_bar <- "GreentoBlue"

    #Advanced Parameters:
    use_built_in_gene_universe <- FALSE
    minimum_pathway_hit_count <- 5
    p_value_threshold_for_output <- 0.05

    pathgenes <- Custom_genes_L2P

    celltype <- colnames(pathgenes)

    cell.list <- list()

    for(cell in celltype) {
        genes <- l2psupp::updategenes(pathgenes %>% select(.data[[cell]]) %>% pull(), trust = 1)
        cell.list[[cell]] <- c(cell,cell,o2o(genes,"mouse","human"))
    }

    ##--------------- ##
    ## Error Messages ##
    ## -------------- ##

    if(plot_top_pathways == FALSE & length(pathways_to_use) == 0){
        stop("ERROR: Enter at least one pathway in 'Pathways to use' when 'Plot top pathways' is set to FALSE")
    }
    if(select_by_rank == FALSE){
      sigcol <- gsub("_pval|p_val_|_adjpval|p_val_adj_","",significance_column)
       fccol <- gsub("_FC|_logFC|avg_logFC_|avg_log2FC_","",fold_change_column)
       if(sigcol != fccol){
           stop("ERROR: when 'select by rank' is FALSE, under Genelist selected by fold-change and pval thresholds, make sure to select significance and fold change columns from the same group comparison")
       }
    }    

    ## --------- ##
    ## Functions ##
    ## --------- ##

    return_org_genes <- function(l2pout){
        l2pgenes <- as.list(l2pout$genesinpathway)
        l2pgenes <- lapply(l2pgenes, function(x) unlist(strsplit(x," ")))
        l2pgenesnew <- lapply(l2pgenes,function(a) o2o(a,"human",species))
        l2pout$org_genes <- l2pgenesnew
        l2pout$org_genes <- sapply(l2pout$org_genes, paste, collapse=' ')
        l2pout %>% arrange(pval) -> l2pout
        return(l2pout)
    }

    return_orig_genes <- function(l2pout){
        l2pgenes <- as.list(l2pout$genesinpathway)
        l2pgenes <- lapply(l2pgenes, function(x) unlist(strsplit(x," ")))
        l2pgenesorig <- lapply(l2pgenes,function(a) names(new_gene_names[a]))
        l2pout$orig_genes <- l2pgenesorig
        l2pout$orig_genes <- sapply(l2pout$orig_genes, paste, collapse=' ')
        l2pout %>% arrange(pval) -> l2pout
        return(l2pout)
    }
    
    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##
    
    if(select_by_rank == TRUE){
        deg_table %>% dplyr::select(.data[[gene_name_column]],.data[[column_used_to_rank_genes]]) -> genesmat
        genesmat %>% dplyr::arrange(desc(.data[[column_used_to_rank_genes]])) -> genesmat
        genesmat %>% dplyr::filter(!is.na(.data[[column_used_to_rank_genes]])) -> genesmat
        if(select_top_percentage_of_genes == TRUE){
            numselect <- ceiling(0.1*dim(deg_table)[1])
        } else {
            numselect <- select_top_genes
        }
    } else {
        numselect <- NULL
        genesmat <- deg_table %>% dplyr::select(.data[[gene_name_column]],.data[[fold_change_column]],.data[[significance_column]]) -> genesmat
    } 

    if (choose_upregulated_genes == TRUE) {
        if(!is.null(numselect)){
            genes_to_include = head(genesmat[[gene_name_column]], numselect)
            } else {
            if(grepl("logFC|log2FC",fold_change_column)){
                logFC_threshold <- log2(fold_change_threshold)
                genes_to_include = deg_table %>% dplyr::arrange(.data[[significance_column]]) %>% dplyr::filter(.data[[significance_column]] <= significance_threshold & .data[[fold_change_column]] >= logFC_threshold) %>% pull(gene_name_column)
                n <- length(genes_to_include)
                cat(sprintf("Number of DEG genes using set cut-offs: %s \n",n))
                if(n < minimum_number_of_DEG_genes){
                    stop(paste("ERROR: DEG gene count:",n,"is not enough to find enriched pathways. Try to loosen criteria to reach",minimum_number_of_DEG_genes))
                } else {
                    lastgene <- tail(genes_to_include,1)
                }
            } else {
                genes_to_include = deg_table %>% dplyr::arrange(.data[[significance_column]]) %>% dplyr::filter(.data[[significance_column]] <= significance_threshold & .data[[fold_change_column]] >= fold_change_threshold) %>% pull(gene_name_column)
                n <- length(genes_to_include)
                cat(sprintf("Number of DEG genes: %s \n",n))
                if(n < minimum_number_of_DEG_genes){
                    stop(paste("ERROR: DEG gene count:",n,"is not enough to find enriched pathways. Try to loosen criteria to reach",minimum_number_of_DEG_genes))
                } else {
                lastgene <- tail(genes_to_include,1)
                }
            }
        }
    } else {
        if(!is.null(numselect)){
            genes_to_include = tail(genesmat[[gene_name_column]], numselect)
            lastgene <- head(genes_to_include,1)
        } else {
            if(grepl("logFC|log2FC",fold_change_column)){
                logFC_threshold <- -1*log2(fold_change_threshold)
                genes_to_include = deg_table %>% dplyr::arrange(.data[[significance_column]]) %>% dplyr::filter(.data[[significance_column]] <= significance_threshold & .data[[fold_change_column]] <= logFC_threshold) %>% pull(gene_name_column)
                n <- length(genes_to_include)
                cat(sprintf("Number of DEG genes: %s \n",n))
                if(n < minimum_number_of_DEG_genes){
                    stop(paste("ERROR: DEG gene count:",n,"is not enough to find enriched pathways. Try to loosen criteria to reach",minimum_number_of_DEG_genes))
                } else { 
                lastgene <- tail(genes_to_include,1)
                }
            } else {
                genes_to_include = deg_table %>% dplyr::arrange(.data[[significance_column]]) %>% dplyr::filter(.data[[significance_column]] <= significance_threshold & .data[[fold_change_column]] <= -1*fold_change_threshold) %>% pull(gene_name_column)
                n <- length(genes_to_include)
                cat(sprintf("Number of DEG genes: %s \n",n))
                if(n < minimum_number_of_DEG_genes){
                    stop(paste("ERROR: DEG gene count:",n,"is not enough to find enriched pathways. Try to loosen criteria to reach",minimum_number_of_DEG_genes))
                } else { 
                lastgene <- tail(genes_to_include,1)
                }
            }
        }
    }

    genes_to_include <- as.vector(unique(unlist(genes_to_include)))
    gene_universe = as.vector(unique(unlist(genesmat[gene_name_column])))
    
    if(species == "Human"){
        new_gene_names <- sapply(genes_to_include, function(x) updategenes(x, trust=1))
        #Print out genes in genelist that are updated:
        updated_genes_idx <- sapply(seq_along(new_gene_names), function(i) names(new_gene_names)[i] != new_gene_names[[i]])
        updated_genes <- new_gene_names[updated_genes_idx]
        updated_genes_num <- sum(updated_genes_idx)
        
        print(paste("Number of updated genes:",updated_genes_num))
        
        cat("Original:Updated\n")
        print(sapply(seq_along(updated_genes), function(i) print(paste0(names(updated_genes)[i],":",updated_genes[i]))))
    
        #Set up genelist and gene universe using updated genes
        genes_to_include <- as.character(new_gene_names)
        gene_universe <- updategenes(gene_universe,trust=1)
        lastgene <- names(tail(new_gene_names,1))

    } else {
        #Get homologs for genelist:
        orth_genes <- sapply(genes_to_include, function(x) o2o(x, species,"human")[1])
        no_orth <- names(orth_genes[unlist(lapply(orth_genes, function(x) is.na(x)))])
        orth <- orth_genes[unlist(lapply(orth_genes, function(x) !is.na(x)))]
        
        #Print out numbers of genes with homologs in genelist:
        num_no_orth <- length(no_orth)
        perc_num_no_orth <- formatC((length(no_orth)/length(genes_to_include))*100,digits=2,format="f")
        num_orth <- length(orth)
        perc_num_orth <- formatC((length(orth)/length(genes_to_include))*100,digits=2,format="f")

        cat(paste("\n\nNumber of genes in the genelist without a homologue:",num_no_orth,", Percentage:",perc_num_no_orth,"%\n"))
        print(no_orth)
        
        cat(paste("\n\nNumber of genes in the genelist with a homologue:",num_orth,",Percentage:",perc_num_orth,"%\n"))
        cat("\nGene:Homolog\n")
        print(sapply(seq_along(orth), function(i) paste0(names(orth)[i],":",orth[i])))
        
        #Set up genelist using gene homologs
        lastgene <- names(tail(orth,1))
        genes_to_include <- as.character(unlist(orth))

        #Get homologs for gene universe:
        orth_gene_universe <- sapply(gene_universe, function(x) o2o(x, species,"human")[1])
        no_orth_gu <- names(orth_gene_universe[unlist(lapply(orth_gene_universe, function(x) is.na(x)))])
        orth_gu <- orth_gene_universe[unlist(lapply(orth_gene_universe, function(x) !is.na(x)))]

        gene_universe <- as.character(unlist(orth_gu))

        #Print out numbers of genes with homologs in genelist:
        num_no_orth <- length(no_orth_gu)
        num_orth <- length(orth_gu)
        cat(paste("\n\nNumber of genes in the gene universe without a homologue:",num_no_orth,"\n"))
        #print(no_orth_gu)
        cat(paste("Number of genes in the gene universe with a homologue:",num_orth,"\n"))
        #print(gene_universe)
    }
    
    sizegenelist <- length(genes_to_include)
    sizeuniv <- length(gene_universe)

    cat("\n\nNumber of genes selected for pathway analysis: ",sizegenelist)
    cat("\nSize of gene universe: ", sizeuniv)
    pctuniv <- (sizegenelist/sizeuniv)*100
    pctuniv <- formatC(pctuniv,digits=2,format="f")
    cat(paste0("\nFinal genelist as percent of gene universe (ideally < 20 %) : ",pctuniv,"%\n"))
    if(select_by_rank == FALSE){
        lastgenedat <- deg_table %>% dplyr::filter(.data[[gene_name_column]] == lastgene) %>% select(.data[[gene_name_column]],fold_change_column,significance_column)
    } else {
        grp <- gsub("_tstat","",column_used_to_rank_genes)
        fold_change_column <- paste0(grp,"_FC")
        significance_column <- paste0(grp,"_pval")
        lastgenedat <- deg_table %>% dplyr::filter(.data[[gene_name_column]] == lastgene) %>% select(.data[[gene_name_column]],column_used_to_rank_genes,fold_change_column,significance_column)
    }
    cat("\n\nCheck p-value for least significant gene in genelist (ideally p <= 0.10) :\n\n")
    print(lastgenedat)

    if (use_built_in_gene_universe == TRUE) {
        x <- l2p(genes_to_include, categories=collections_to_include)
        cat("\n\nUsing built-in gene universe.\n")
        cat(paste0("Total number of pathways tested: ", nrow(x)))
    } else {
        x <- l2p(genes_to_include, custompathways=cell.list, categories=collections_to_include, universe=gene_universe, rawdataframe = 1)
        #x <- l2p(genes_to_include, custompathways=cell.list, categories="", universe=gene_universe, rawdataframe = 1)
        cat("\n\nUsing all genes included in the differential expression analysis as gene universe.\n\n")
        cat(paste0("Total number of pathways tested: ", nrow(x)))
    }

    cat("\n\nOnly pathways that are over-represented are tested by One-Sided Fisher's Exact Test.\nPathways that show under-representation are excluded and not returned")
    
    x <- x %>% select(pathway_name,category,number_hits,percent_gene_hits_per_pathway,enrichment_score,pval,fdr,genesinpathway) %>% 
                mutate(percent_gene_hits_per_pathway = percent_gene_hits_per_pathway) %>% 
                dplyr::filter(number_hits >= minimum_pathway_hit_count) %>% 
                dplyr::filter(pval < p_value_threshold_for_output) %>% 
                dplyr::arrange(pval)
    
    x %>% dplyr::mutate(pathwayname2 = stringr::str_replace_all(pathway_name, "_", " ")) -> goResults
    goResults$pathwayname2 <- str_to_upper(goResults$pathwayname2)
    goResults$pathwayname2 <- trimws(goResults$pathwayname2)
    goResults %>% dplyr::mutate(pathwayname2 = stringr::str_wrap(pathwayname2,pathway_axis_label_max_length)) -> goResults
    goResults <- distinct(goResults, pathwayname2, .keep_all = TRUE)

    if(plot_top_pathways == TRUE){
        goResults %>% top_n(number_of_pathways_to_plot, wt=-log(pval)) -> goResults 
    } else {
        goResults %>% dplyr::filter(.data[["pathway_name"]] %in% pathways_to_use) -> goResults              #User-selected pathways to plot in bubble plot
        if(dim(goResults)[1] < length(pathways_to_use)){
        cat("\nSome selected pathways are not showing in plot, check for spelling errors:\n\n")
        cat(pathways_to_use[!pathways_to_use %in% goResults[["pathway_name"]]] )    
        }
    }

    #Set up plot data and parameters:
    if(plot_bubble == TRUE) {
        leglab = plot_bubble_color
        leglab2 = plot_bubble_size
        x_label <- sort_bubble_plot_by

        plot_bubble_list <- list("Fisher's Exact pval" = "pval", "fdr corrected pval" = "fdr", "number of hits" = "number_hits","percent gene hits per pathway" = "percent_gene_hits_per_pathway","enrichment score" = "enrichment_score") 
        plot_bubble_size <- plot_bubble_list[[plot_bubble_size]]
        plot_bubble_color <- plot_bubble_list[[plot_bubble_color]]
        sort_bubble_plot_by <- plot_bubble_list[[sort_bubble_plot_by]]

        if(plot_bubble_color %in% c("pval","fdr")) {
            goResults$color <- -log10(goResults[[plot_bubble_color]])
            leglab = paste0("-log10(",plot_bubble_color,")")
        } else{
            goResults$color <- goResults[[plot_bubble_color]]
        }

        if(plot_bubble_size %in% c("pval","fdr")) {
            goResults$size <- -log10(goResults[[plot_bubble_size]])
            leglab2 = paste0("-log10(",plot_bubble_size,")")
        } else{
            goResults$size <- goResults[[plot_bubble_size]]
        }

        if(sort_bubble_plot_by %in% c("pval","fdr")) {
            goResults$sort <- -log10(goResults[[sort_bubble_plot_by]])
            x_label <- paste0("-log10(",sort_bubble_plot_by,")")
        } else {
            goResults$sort <- goResults[[sort_bubble_plot_by]]
        }

        goResults$color <- as.numeric(goResults$color)
        minp = floor(min(goResults$color))
        maxp = ceiling(max(goResults$color))
        sizemax = ceiling(max(goResults$size)/10)*10  
        
        goResults %>% dplyr::mutate(percorder = order(goResults$sort)) -> goResults
        goResults$pathwayname2 <- factor(goResults$pathwayname2, levels = goResults$pathwayname2[goResults$percorder])
        xmin = min(goResults$sort) - 0.1*min(goResults$sort)
        xmax = max(goResults$sort) + 0.1*min(goResults$sort)

        
        bubblecols <- list("blues" = c("#56B1F7","#132B43"),
                            "reds" = c("#fddbc7","#b2182b"),
                            "blue to red" = c("dark blue","red"))

        gplot <- goResults %>% 
            ggplot(aes(x=sort, y=pathwayname2, col=color, size=size)) +
                geom_point() +
                theme_classic() +
                labs(col = leglab, size=leglab2, y="Pathway", x=x_label) +
                theme(aspect.ratio = 1, 
                    text = element_text(size=pathway_axis_label_font_size), 
                    legend.position = "right", 
                    legend.key.height = unit(1, "cm"),
                    axis.title.y = element_text(size=20,margin = margin(r = 80)), 
                    axis.title.x = element_text(size=20,margin = margin(t = 20)),
                    axis.text.y = element_text(colour="black"), 
                    axis.text.x = element_text(colour="black", size = 15), 
                    legend.key.size = unit(1, 'cm'), 
                    legend.title = element_text(size=15,margin = margin(l = 20)), 
                    legend.text=element_text(size=15)) +
                xlim(xmin,xmax) +
                scale_colour_gradient(low = bubblecols[[bubble_colors]][1], high = bubblecols[[bubble_colors]][2]) +
                expand_limits(colour = seq(minp, maxp, by = 1), size = seq(0, sizemax,by=10))  +
                guides(colour = guide_colourbar(order = 1), size = guide_legend(order=2))
        } else {
            colpal = list("Blue" = "Blues", "Green" = "Greens", "Grey" = "Greys", 
                            "Red" = "Reds", "Purple" = "Purples", "Orange" = "Oranges",
                            "GreentoBlue" = "GnBu", "OrangetoRed" = "OrRd", "YellowOrangeRed" = "YlOrRd") 
            pal <- colpal[[color_for_bar]]
            if(use_fdr_pvals == TRUE){
                goResults %>% mutate(fdr = -log(fdr)) %>% arrange(desc(fdr))-> df1
                gplot <- ggplot(df1, aes(x=reorder(pathwayname2,fdr), y=fdr)) + 
                    geom_bar(aes(fill = fdr),stat = "identity") +
                    scale_fill_distiller(name = expression(paste("-lo",g[10],"(FDR adj ",italic("p"),"value)")), 
                        palette = pal, direction = 1) +  
                    theme_classic() +
                    labs(y=expression(paste("-lo",g[10],"(FDR adjusted ",italic("p")," value)")), 
                        x="Pathways") + 
                    theme(aspect.ratio = 1, 
                        axis.title.y = element_text(size=20,margin = margin(r = 80)), 
                    axis.title.x = element_text(size=20,margin = margin(t = 20)),
                    axis.text.y = element_text(colour="black"), 
                    axis.text.x = element_text(colour="black", size = pathway_axis_label_font_size), 
                    legend.key.size = unit(1, 'cm'), 
                    legend.title = element_text(size=15,margin = margin(l = 20)), 
                    legend.text=element_text(size=15)) +
                coord_flip()           
            } else {
                goResults %>% mutate(pval = -log(pval)) %>% arrange(desc(pval))-> df1
                gplot <- ggplot(df1, aes(x=reorder(pathwayname2,pval), y=pval)) + 
                    geom_bar(aes(fill = pval), stat = "identity") +
                    scale_fill_distiller(name = expression(paste("-lo",g[10],italic("(p"),"-value)")), 
                        palette = pal, direction = 1, limits = c(min(df1$pval)-1,max(df1$pval))) + 
                    theme_classic() +
                    labs(y=expression(paste("-lo",g[10],italic("(p"),"-value)")), 
                        x="Pathways") + 
                    theme(aspect.ratio = 1, 
                        axis.title.y = element_text(size=20,margin = margin(r = 80)), 
                        axis.title.x = element_text(size=20,margin = margin(t = 20)),
                        axis.text.y = element_text(colour="black"), 
                        axis.text.x = element_text(colour="black", size = pathway_axis_label_font_size), 
                        legend.key.size = unit(1, 'cm'), 
                        legend.title = element_text(size=15,margin = margin(l = 20)), 
                        legend.text=element_text(size=15)) +
                    coord_flip()  
            }
        }

    print(gplot)
    
    if (species != "Human"){
        x <- as.data.frame(return_org_genes(x))
    } else {
        x <- as.data.frame(return_orig_genes(x))
    }

    x$enrichment_score <- as.numeric(formatC(x$enrichment_score,digits=3,format="f"))
    x$percent_gene_hits_per_pathway <- as.numeric(formatC(x$percent_gene_hits_per_pathway,digits=3,format="f"))
   
    return(x)
}

    ## ---------------------------- ##
    ## Global Imports and Functions ##
    ## ---------------------------- ##

    ## Functions defined here will be available to call in
    ## the code for any table.

    ## --------------- ##
    ## End of Template ##
    ## --------------- ##

print("template_function_L2P_Downregulated_Custom.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_DEG_Analysis<-readRDS(paste0(rds_output,"/var_DEG_Analysis.rds"))
Input_is_Seurat_count <- 0
for(item in var_DEG_Analysis){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_DEG_Analysis<-as.data.frame(var_DEG_Analysis)}else{var_DEG_Analysis <- var_DEG_Analysis}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Custom_genes_L2P<-readRDS(paste0(rds_output,"/var_Custom_genes_L2P.rds"))
Input_is_Seurat_count <- 0
for(item in var_Custom_genes_L2P){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Custom_genes_L2P<-as.data.frame(var_Custom_genes_L2P)}else{var_Custom_genes_L2P <- var_Custom_genes_L2P}
invisible(graphics.off())
var_L2P_Downregulated_Custom<-L2P_Downregulated_Custom(var_DEG_Analysis,var_Custom_genes_L2P)
invisible(graphics.off())
saveRDS(var_L2P_Downregulated_Custom, paste0(rds_output,"/var_L2P_Downregulated_Custom.rds"))
