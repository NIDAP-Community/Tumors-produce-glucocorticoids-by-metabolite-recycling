# Expression Heatmap [CCBR] (89a32987-10d9-4233-91c3-e9adf3dcc517): v543
Up_in_Human_TITRs <- function(Batch_Corrected_Counts,ccbr1224_metadata) {
    ## This function uses pheatmap to draw a heatmap, scaling first by rows
    ## (with samples in columns and genes in rows)

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(colorspace)
    library(dendsort)
    library(ComplexHeatmap)
    library(dendextend)
    library(tibble)
    library(stringr)
    library(RColorBrewer)
    library(dplyr)
    library(grid)
    library(gtable)
    library(gridExtra)
    library(gridGraphics)

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    #Basic Parameters:
    counts_matrix <- Batch_Corrected_Counts
    sample_metadata <- ccbr1224_metadata
    gene_column_name <- "Gene"
    group_columns <- c("Group")
    sample_name_column <- "Sample"
    samples_to_include = c("WT_1","WT_2","WT_3","GRfoxp3_1","GRfoxp3_2","GRfoxp3_3")
    include_all_genes <- FALSE
    filter_top_genes_by_variance = FALSE
    top_genes_by_variance_to_include <- 500
    specific_genes_to_include_in_heatmap = "Abhd13 Acaa2 Acot7 Acot9 Acp5 Acsl1 Acsl4 Actr1B Ada Adam12 Adi1 Adprh Agtrap Akirin2 Akr1A1 Ankrd10 Anxa5 Apobec3C Apobec3G Arhgef12 Arhgef9 Arl6Ip6 Aspm Atf3 Atp10D Atp13A3 Atp6V0A1 Azi2 Bach1 Bak1 Batf Bcl2L1 Birc5 Blm Bmpr2 Brca1 Bst2 Bub1 C12Orf75 C15Orf53 C16Orf87 C1Orf112 C3Ar1 Cadm1 Calm3 Calu Capg Capn2 Capzb Casc5 Ccdc88B Ccl20 Ccna2 Ccnb1 Ccnb2 Ccng2 Ccr1 Ccr5 Ccr8 Cd177 Cd276 Cd58 Cd70 Cd79B Cd80 Cd82 Cd83 Cdc20 Cdca4 Cdca7 Cdk1 Cdkn1A Ceacam1 Cenpe Cenpf Ckap2 Cks1B Clec7A Clic5 Coro1C Cotl1 Cpd Cradd Creb3L2 Csf1 Csk Ctnna1 Ctnnb1 Ctsc Ctsh Ctsl1 Cttnbp2Nl Cxcl13 Cxcr6 Dck Dedd2 Dlgap5 Dnph1 Dpysl2 Dtl Dusp1 Dusp16 Dusp4 Ell Elovl5 Entpd1 Epsti1 Eri1 Fah Fam105B Fam129A Fam50A Fam53B Fancg Fanci Fank1 Fbxo5 Fen1 Fkbp4 Fmnl3 Fndc3B Fosl2 Foxp3 Furin Fut8 Gapdh Gbp5 Gcnt1 Gdpd5 Glrx Gng2 Gnly Golga8A Got2 Gpr75-Asb3 Gtpbp1 Gzmb Hapln3 Hexb Hist1H1C Hivep3 Hla-Dma Hla-Dpa1 Hla-Dpb1 Hla-Dqa1 Hla-Dra Hla-Drb1 Hla-Drb5 Hlf Hn1 Hnrpll Hpgd Hsd17B10 Hspa1A Hspa1B Hsph1 Htatip2 Ier2 Ier3 Ifi6 Ikzf4 Il12Rb2 Il17A Il17F Il17Rb Il1R2 Il21R Il21R-As1 Il22 Il23R Il2Ra Il2Rb Il4I1 Irak2 Irf5 Itgae Jakmip1 Kiaa0101 Kif11 Kif14 Kif20B Kif23 Lair2 Lap3 Layn Ldha Lgals1 Lgals3 Lmcd1 Loc100130476 Loc541471 Lrp10 Ltb Lztr1 Macc1 Mapkapk3 Mast4 Mcam Mcm2 Mcm4 Med12 Mettl7A Mgme1 Mical3 Miip Mki67 Mrps6 Mx1 Myo1C Myo7A Nampt Ncaph Ncf4 Ndfip2 Nfkbia Ninj1 Nme1 Nusap1 Oas1 Oas3 Paip1 Parp10 Parp9 Pde4A Pde4D Pdia6 Perp Pex11B Pfkfb3 Phlda1 Pi4Kb Pign Pkm Plagl2 Pld3 Plekhf1 Plekhg1 Plk1 Plscr1 Pmaip1 Polq Pomp Ppa1 Ppp1R16B Prc1 Prdx5 Prkar1A Prnp Psmb6 Pspc1 Ptp4A3 Ptprj Pxk Rab10 Rab11Fip1 Rab33A Ralb Rasgrp1 Rassf5 Rbbp8 Rbpj Rcn2 Rgs1 Rgs19 Rhbdd2 Rhoc Rorc Rpusd3 Rqcd1 Rrm2 Rtkn2 S100A4 Sat1 Sdc4 Sdf4 Sec14L1 Sh2D2A Sh3Tc1 Shcbp1 Slbp Slc35B1 Slc35F2 Slc39A6 Slc6A6 Sox4 Spats2L Spock2 Sppl2A Sqle Ssh1 Stat1 Stat5A Stk17B Stmn1 Suox Surf4 Syngr2 Syt11 Tars Tax1Bp3 Tdp1 Tet3 Tfrc Timeless Tk1 Tmem106C Tmem200A Tmem30A Tmprss3 Tnfrsf18 Tnfrsf1B Tnfrsf4 Tnfrsf8 Tnfrsf9 Tnip1 Tns3 Top2A Tpm4 Tpx2 Traf3 Trpc4Ap Ttc39C Tubb Ube2C Uevld Uhrf1 Uhrf1Bp1L Vdac1 Vdr Vhl Wdr76 Zbed2 Zeb2 Zfp62 Zwint" 
    
    #Visual Parameters:
    heatmap_color_scheme <- "Default"
    autoscale_heatmap_color <- TRUE
    set_min_heatmap_color <- -2
    set_max_heatmap_color <- 2
    group_colors <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
    assign_group_colors <- FALSE
    assign_color_to_sample_groups <- c()
    legend_font_size <- 10 
    display_gene_names <- TRUE
    gene_name_font_size <- 4
    display_sample_names <- TRUE
    sample_name_font_size <- 8
    display_dendrograms <- TRUE
    reorder_dendrogram <- FALSE
    reorder_dendrogram_order <- c()
    manually_rename_samples <- FALSE
    samples_to_rename <- c("")
    display_numbers <- FALSE
    aspect_ratio <- "Auto"

    #Advanced Parameters
    distance_metric <- "correlation"
    clustering_method <- "average"
    center_and_rescale_expression <- TRUE
    cluster_genes <- TRUE
    cluster_samples <- FALSE
    arrange_sample_columns <- TRUE
    order_by_gene_expression <- FALSE
    gene_to_order_columns <- " "
    gene_expression_order <- "low_to_high"
    return_z_scores <- FALSE

    ##--------------- ##
    ## Error Messages ##
    ## -------------- ##

    if(include_all_genes == TRUE && filter_top_genes_by_variance == TRUE){
        stop("ERROR: Choose only one of 'Include all genes' or 'Filter top genes by variance' as TRUE")
    }

    if((cluster_samples == TRUE && arrange_sample_columns == TRUE) | (arrange_sample_columns == TRUE && order_by_gene_expression == TRUE) | 
    (arrange_sample_columns == TRUE && cluster_samples == TRUE) | (cluster_samples == FALSE && arrange_sample_columns == FALSE && order_by_gene_expression == FALSE)) {
     stop("ERROR: Choose only one of 'Cluster Samples', 'Arrange sample columns', or 'order by gene expression' as TRUE")   
    }

    ## --------- ##
    ## Functions ##
    ## --------- ##

    getourrandomcolors<-function(k){
        seed=10
        n <- 2e3
        ourColorSpace <- colorspace::RGB(runif(n), runif(n), runif(n))
        ourColorSpace <- as(ourColorSpace, "LAB")
        currentColorSpace <- ourColorSpace@coords
        # Set iter.max to 20 to avoid convergence warnings.
        set.seed(seed)
        km <- kmeans(currentColorSpace, k, iter.max=20)
        return( unname(hex(LAB(km$centers))))
    }

    ## Begin pal() color palette function∂:
    pal = function (n, h=c(237, 43), c=100, l=c(70, 90), power=1, fixup=TRUE, gamma=NULL, alpha=1, ...) {
        if (n < 1L) {
            return(character(0L))
        }
        h <- rep(h, length.out = 2L)
        c <- c[1L]
        l <- rep(l, length.out = 2L)
        power <- rep(power, length.out = 2L)
        rval <- seq(1, -1, length = n)
        rval <- hex(
            polarLUV(
                L = l[2L] - diff(l) * abs(rval)^power[2L], 
                C = c * abs(rval)^power[1L],
                H = ifelse(rval > 0, h[1L], h[2L])
            ),
            fixup=fixup, ...
        )
        if (!missing(alpha)) {
            alpha <- pmax(pmin(alpha, 1), 0)
            alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
                width = 2L, upper.case = TRUE)
            rval <- paste(rval, alpha, sep = "")
        }
        return(rval)
    } 
    # End pal() color palette function:

    ## Begin doheatmap() function:
    doheatmap <- function(dat, clus, clus2, ht, rn, cn, col, dispnum) {
        #require(pheatmap)
        #require(dendsort)
        col.pal <- np[[col]]
        if (FALSE) {
            col.pal = rev(col.pal)
        }
        # Define metrics for clustering
        drows1 <- distance_metric
        dcols1 <- distance_metric
        minx = min(dat)
        maxx = max(dat)
        if (autoscale_heatmap_color) {
            breaks = seq(minx, maxx, length=100)
            legbreaks = seq(minx, maxx, length=5)
        } else {
            breaks = seq(set_min_heatmap_color, set_max_heatmap_color, length=100)
            legbreaks = seq(set_min_heatmap_color, set_max_heatmap_color, length=5)
        }
        breaks = sapply(breaks, signif, 4)
        legbreaks = sapply(legbreaks, signif, 4)
        # Run cluster method using 
        hc = hclust(dist(t(dat)), method=clustering_method)
        hcrow = hclust(dist(dat), method=clustering_method)
        if (FALSE) {
            sort_hclust <- function(...) as.hclust(rev(dendsort(as.dendrogram(...))))
        } else {
            sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
        }
        if (clus) {
            colclus <- sort_hclust(hc)
        } else {
            colclus = FALSE
        }
        if (clus2) {
            rowclus <- sort_hclust(hcrow)
        } else {
            rowclus = FALSE
        }
        if (display_dendrograms) {
            treeheight <- 25
        } else {
            treeheight <- 0
        }

        hm.parameters <- list(
            dat, 
            color=col.pal,
            legend_breaks=legbreaks,
            legend=TRUE,
            scale="none",
            treeheight_col=treeheight,
            treeheight_row=treeheight,
            kmeans_k=NA,
            breaks=breaks,
            display_numbers=dispnum,
            number_color = "black",
            fontsize_number = 8,
            height=80,
            cellwidth = NA, 
            cellheight = NA, 
            fontsize= legend_font_size,   
            fontsize_row=gene_name_font_size,
            fontsize_col=sample_name_font_size,
            show_rownames=rn, 
            show_colnames=cn,
            cluster_rows=rowclus, 
            cluster_cols=clus,
            clustering_distance_rows=drows1, 
            clustering_distance_cols=dcols1,
            annotation_col = annotation_col,
            annotation_colors = annot_col,
            labels_col = labels_col
        )
        mat = t(dat)
        callback = function(hc, mat) {
            dend = rev(dendsort(as.dendrogram(hc)))
            if(reorder_dendrogram == TRUE) {
                dend %>% dendextend::rotate(reorder_dendrogram_order) -> dend
            } else {
                dend %>% dendextend::rotate(c(1:nobs(dend))) 
            }
            as.hclust(dend)
        }

        ## Make Heatmap
        phm <- do.call("pheatmap", c(hm.parameters, list(clustering_callback=callback)))
        
    }
    # End doheatmap() function.

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

    ## Build different color spectra options for heatmap:
    np0 = pal(100) 
    np1 = diverge_hcl(100, c=100, l=c(30, 80), power=1) # Blue to Red
    np2 = heat_hcl(100, c=c(80, 30), l=c(30, 90), power=c(1/5, 2)) # Red to Vanilla
    np3 = rev(heat_hcl(100, h=c(0, -100), c=c(40, 80), l=c(75, 40), power=1)) # Violet to Pink
    np4 = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100)) #Red to yellow to blue
    np5 = colorRampPalette(c("steelblue","white", "red"))(100)  # Steelblue to White to Red

    ## Gather list of color spectra and give them names for the GUI to show.
    np = list(np0, np1, np2, np3, np4, np5)
    names(np) = c("Default","Blue to Red","Red to Vanilla","Violet to Pink","Bu Yl Rd","Bu Wt Rd")

    ## Parse input counts matrix. Subset by samples.
    df1 <- counts_matrix
    # Swap out Gene Name column name, if it's not 'Gene'.
    if(gene_column_name != "Gene"){
        # Drop original Gene column
        df1 = df1[,!(colnames(df1)%in% c("Gene")) ]
        # Rename column to Gene
        colnames(df1)[which(colnames(df1) == gene_column_name)] <- 'Gene'
    }
    # Get sample columns
    samples_to_include <- samples_to_include[samples_to_include != gene_column_name]
    samples_to_include <- samples_to_include[samples_to_include != "Gene"]
    samples_to_include <- samples_to_include[samples_to_include != "GeneName"]

    # Build new counts matrix containing only sample subset chosen by user.
    df1 <- df1[,append("Gene", samples_to_include)]
    df.orig = df1
    df.orig %>% dplyr::group_by(Gene) %>% summarise_all(funs(mean)) -> df
    df.mat = df[ , (colnames(df) != "Gene" )] %>% as.data.frame
    df %>% dplyr::mutate(Gene = stringr::str_replace_all(Gene, "_", " ")) -> df
    df %>% dplyr::mutate(Gene = stringr::str_wrap(Gene,50)) -> df            ##### MC: added wrapper for looooonnggg names
    row.names(df.mat) <- df$Gene
    rownames(df.mat) <- str_wrap(rownames(df.mat),10)
    df.mat <- as.data.frame(df.mat)

    ## Subset counts matrix by genes.
    # Toggle to include all genes in counts matrix (in addition to any user-submitted gene list).
    if (include_all_genes == FALSE) {
        # Add user-submitted gene list (optional).
        genes_to_include_parsed = c()
        genes_to_include_parsed = strsplit(specific_genes_to_include_in_heatmap, " ")[[1]]
        df.mat[genes_to_include_parsed,] -> df.final.extra.genes
        if(filter_top_genes_by_variance == TRUE) {
            # Want to filter all genes by variance.
            df.final = as.matrix(df.mat)
            var <- matrixStats::rowVars(df.final)
            df <- as.data.frame(df.final)
            rownames(df) <- rownames(df.final)
            df.final <- df
            df.final$var <- var
            df.final %>% rownames_to_column("Gene") -> df.final 
            df.final %>% dplyr::arrange(desc(var)) -> df.final
            df.final.extra.genes = dplyr::filter(df.final, Gene %in% genes_to_include_parsed)
            df.final = df.final[1:top_genes_by_variance_to_include,]
            df.final = df.final[complete.cases(df.final),]
            # Rbind user gene list to variance-filtered gene list and deduplicate.
            df.final <- rbind(df.final, df.final.extra.genes)
            df.final <- df.final[!duplicated(df.final),] 
            rownames(df.final) <- df.final$Gene
            df.final$Gene <- NULL
            df.final$var <- NULL
        } else {
            # Want to use ONLY user-provided gene list.
            df.final <- df.final.extra.genes
            df.final <- df.final[!duplicated(df.final),]
            # Order genes in heatmap by user-submitted order of gene names.
            df.final <- df.final[genes_to_include_parsed,]
            #df.final$Gene <- NULL
        }
    } else {
        df.final <- df.mat
        df.final$Gene <- NULL
    }
    
        ## Optionally apply centering and rescaling (default TRUE).
    if (center_and_rescale_expression == TRUE) {
            tmean.scale = t(scale(t(df.final)))
            tmean.scale = tmean.scale[!is.infinite(rowSums(tmean.scale)),]
            tmean.scale = na.omit(tmean.scale)
    } else {
            tmean.scale = df.final
    }

    if(order_by_gene_expression == TRUE){
        gene_to_order_columns <- gsub(" ","",gene_to_order_columns)
        if(gene_expression_order == "low_to_high"){
        tmean.scale <- tmean.scale[,order(tmean.scale[gene_to_order_columns,])] #order from low to high 
        } else{
        tmean.scale <- tmean.scale[,order(-tmean.scale[gene_to_order_columns,])] #order from high to low  
        }
    }

    df.final <- as.data.frame(tmean.scale)

    ## Parse input sample metadata and add annotation tracks to top of heatmap.
    annot <- sample_metadata
    # Filter to only samples user requests.
    annot %>% dplyr::filter(.data[[sample_name_column]] %in% samples_to_include) -> annot
    # Arrange sample options.
    if(arrange_sample_columns) {
      annot %>% dplyr::arrange_(.dots=group_columns) -> annot
      df.final <- df.final[,match(annot[[sample_name_column]],colnames(df.final))] 
    }
    # Build subsetted sample metadata table to use for figure.

    colorlist <- c("#5954d6","#e1562c","#b80058","#00c6f8","#d163e6","#00a76c","#ff9287","#008cf9","#006e00","#796880","#FFA500","#878500")
    names(colorlist) <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
    group_colors <- colorlist[group_colors]

    annot %>% dplyr::select(group_columns) -> annotation_col    
    annotation_col = as.data.frame(unclass(annotation_col))
    annotation_col[] <- lapply(annotation_col,factor)
    x <- length(unlist(lapply(annotation_col,levels)))
    if(x>length(group_colors)){
        k=x-length(group_colors)
        more_cols<- getourrandomcolors(k) 
        group_colors <- c(group_colors, more_cols)
    }
    rownames(annotation_col) <- annot[[sample_name_column]]
    annot_col = list()
    b=1
    i=1
    while (i <= length(group_columns)){
        nam <- group_columns[i]
        grp <- as.factor(annotation_col[,i])
        c <- b+length(levels(grp))-1
        col = group_colors[b:c]
        names(col) <- levels(grp)
        assign(nam,col)
        annot_col = append(annot_col,mget(nam))
        b = c+1
        i=i+1
    }

    if(assign_group_colors == TRUE){
            colassign <- assign_color_to_sample_groups
            groupname <- c()
            groupcol <- c() 
            for (i in 1:length(colassign)) {
                groupname[i] <- strsplit(colassign[i], ": ?")[[1]][1]
                groupcol[i] <- strsplit(colassign[i], ": ?")[[1]][2]
            }
            annot_col[[1]][groupname] <- groupcol
    }

    ## Setting labels_col for pheatmap column labels.
    if (manually_rename_samples == TRUE) {
        # Use user-provided names to rename samples.
        replacements = samples_to_rename
        old <- c()
        new <- c()
        labels_col <- colnames(df.final)
        for (i in 1:length(replacements)) {
            old <- strsplit(replacements[i], ": ?")[[1]][1]
            new <- strsplit(replacements[i], ": ?")[[1]][2]
            old=gsub("^[[:space:]]+|[[:space:]]+$","",old)
            new=gsub("^[[:space:]]+|[[:space:]]+$","",new)
            labels_col[labels_col==old]=new           
        }
    } else {
        ## Use original column names for samples.
        labels_col <- colnames(df.final)
    }

    ## Print number of genes to log.
    print(paste0("The total number of genes in heatmap: ", nrow(df.final)))

    ## Make the final heatmap.
    p <- doheatmap(dat=df.final, clus=cluster_samples, clus2=cluster_genes, ht=50, rn=display_gene_names, cn=display_sample_names, col=heatmap_color_scheme, dispnum=display_numbers)
    p@matrix_color_mapping@name <- " "
    p@matrix_legend_param$at <- as.numeric(formatC(p@matrix_legend_param$at, 2))
    p@column_title_param$gp$fontsize <- 10
    print(p)

    ## If user sets toggle to TRUE, return Z-scores.
    ## Else return input counts matrix by default (toggle FALSE).
    ## Returned matrix includes only genes & samples used in heatmap.
    if(return_z_scores){
        df.new <- data.frame(tmean.scale) # Convert to Z-scores.
        df.new %>% rownames_to_column("Gene") -> df.new
        return(df.new)
    } else {
        df.final %>% rownames_to_column("Gene") -> df.new
        return(df.new)
    }
}

## ---------------------------- ##
## Global Imports and Functions ##
## ---------------------------- ##

## Functions defined here will be available to call in
## the code for any table.

## --------------- ##
## End of Template ##
## --------------- ##

print("template_function_Up_in_Human_TITRs.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Batch_Corrected_Counts<-readRDS(paste0(rds_output,"/var_Batch_Corrected_Counts.rds"))
Input_is_Seurat_count <- 0
for(item in var_Batch_Corrected_Counts){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Batch_Corrected_Counts<-as.data.frame(var_Batch_Corrected_Counts)}else{var_Batch_Corrected_Counts <- var_Batch_Corrected_Counts}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_ccbr1224_metadata<-readRDS(paste0(rds_output,"/var_ccbr1224_metadata.rds"))
Input_is_Seurat_count <- 0
for(item in var_ccbr1224_metadata){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_ccbr1224_metadata<-as.data.frame(var_ccbr1224_metadata)}else{var_ccbr1224_metadata <- var_ccbr1224_metadata}
invisible(graphics.off())
var_Up_in_Human_TITRs<-Up_in_Human_TITRs(var_Batch_Corrected_Counts,var_ccbr1224_metadata)
invisible(graphics.off())
saveRDS(var_Up_in_Human_TITRs, paste0(rds_output,"/var_Up_in_Human_TITRs.rds"))
