# Filter Low Counts [CCBR] (59e7ee62-1715-4276-b370-e8395168f9d8): v731
Filtered_Counts_1 <- function(Clean_Raw_Counts_1, ccbr1224_rev_response_gsea_meta) {
   
    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(limma)
    library(amap)
    library(colorspace)
    library(dendsort)
    library(dplyr)
    library(edgeR)
    library(ggplot2)
    library(gplots)
    library(gridExtra)
    library(gridGraphics)
    library(lattice)
    library(magrittr)
    library(plotly)
    library(RColorBrewer)
    library(RCurl)
    library(reshape2)
    library(stringr)
    library(tidyverse)
    library(tibble)

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##
    
    #Basic Parameters:
    counts_matrix = Clean_Raw_Counts_1
    sample_metadata <- ccbr1224_rev_response_gsea_meta   
    gene_names_column <- "Gene"
    columns_to_include <- c("Gene","TITR_B16_CLEAN_D16_1","TITR_B16_CLEAN_D16_2","TITR_B16_CLEAN_D16_3","Tu_Tc_B16_CLEAN_D16_1","Tu_Tc_B16_CLEAN_D16_2","Tu_Tc_B16_CLEAN_D16_3")
    sample_names_column <- "Sample"
    groups_column <- "Group"
    labels_column <- "Label"

    #Filtering Parameters:
    outlier_samples_to_remove <- c()
    use_cpm_counts_to_filter <- TRUE
    Minimum_Count_Value_to_be_Considered_Nonzero <- 1
    Minimum_Number_of_Samples_with_Nonzero_Counts_in_Total <- 1
    Use_Group_Based_Filtering <- TRUE
    Minimum_Number_of_Samples_with_Nonzero_Counts_in_a_Group <- 3
    
    #PCA Parameters:
    principal_component_on_x_axis<-1 
    principal_component_on_y_axis<-2 
    legend_position_for_PCA <- "top"
    point_size_for_pca<-1
    add_labels_to_PCA <- TRUE
    label_font_size <- 3
    label_offset_y_ <- 2
    label_offset_x_ <- 2
    samples_to_rename_manually <- c("")

    #Histogram Parameters:
    color_histogram_by_group <- FALSE     
    set_min_max_for_x_axis_for_histogram <- FALSE
    minimum_for_x_axis_for_histogram <- -1
    maximum_for_x_axis_for_histogram <- 1
    legend_position_for_histogram <- 'top'
    legend_font_size_for_histogram <- 10
    number_of_histogram_legend_columns <- 6

    #Visualization Parameters:
    colors_for_plots <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
    number_of_image_rows <- 2
    interactive_plots <- FALSE

    #TCGA:
    plot_correlation_matrix_heatmap <- TRUE
    
    ##--------------- ##
    ## Error Messages ##
    ## -------------- ##

    
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

    make_heatmap <- function(counts_matrix, metadata,colorval) {
        mat <- as.matrix(counts_matrix) 
        tcounts=t(mat)
        tcounts=merge(metadata,tcounts,by.x=sample_names_column,by.y='row.names')
        rownames(tcounts)=tcounts[,labels_column]
        tcounts=tcounts[,!colnames(tcounts)%in%colnames(metadata)]
        d=Dist(tcounts,method="correlation",diag=TRUE)
        dend = rev(dendsort(as.dendrogram(hclust( d,method="average"))))
        m=as.matrix(d)
        sample_metadata <- metadata
        rownames(sample_metadata) = sample_metadata[[labels_column]]
        idx = as.factor(sample_metadata[rownames(m),groups_column])
        col = colorval
        cols <- col[idx]
        new.palette=colorRampPalette(c("blue","green","yellow"),space="rgb")
  
    mk<-function(){
        if(length(colnames(m))>20){
            par(mar=c(0,0,0,0))
            heatmap.2(m,
                        labRow = NA, 
                        labCol = NA,
                        col=new.palette(20),
                        trace="none",
                        colRow = col[idx], 
                        colCol = col[idx],
                        rowDendrogram=dend,
                        colDendrogram=dend,
                        RowSideColors = col[idx],
                        ColSideColors = col[idx],
                        dendrogram = "row",
                        cexRow=3,
                        cexCol=3,
                        margins=c(0,0),   
                        lmat=rbind( c(0,0,2),c(4,1,3) ,c(0,5,6) ), 
                        lhei=c(.2,4,2), 
                        lwid=c(1, .2,4 ), 
                        key.par=list(mgp=c(1.75, .5, 0), 
                        mar=c(7, 2, 3.5, 0), 
                        cex.axis=.1, 
                        cex.lab=3, 
                        cex.main=1, 
                        cex.sub=1),
                        key.xlab = "Correlation",
                        key.ylab="Count",
                        key.title=" ")       
        } else {
            heatmap.2(m,col=new.palette(20),
                        trace="none",
                        colRow = col[idx], 
                        colCol = col[idx],
                        rowDendrogram=dend,
                        colDendrogram=dend,
                        RowSideColors = col[idx],
                        ColSideColors = col[idx],
                        dendrogram = "row",
                        cexRow=3,cexCol=3,margins=c(4,1),  
                        lmat=rbind( c(0,0,2),c(4,1,3) ,c(0,5,6) ), 
                        lhei=c( .2,4,2), 
                        lwid=c(1, .2,4),
                        key.par=list(mgp=c(1.75, .5, 0), mar=c(7, 2, 3.5, 0), cex.axis=.1, cex.lab=3, cex.main=1, cex.sub=1),
                        key.xlab = "Correlation",
                        key.ylab="Count",
                        key.title=" ")
            }
        }
  
            tg<-mk()
            grid.echo(mk)
            gh1<-grid.grab()
            mklegend<-function(){
            plot.new()
            legend(x="top", legend=levels(idx), col=col[as.factor(levels(idx))],pch=15,x.intersp=3,bty ="n",cex=2)
            }
        grid.echo(mklegend )
        gh2<-grid.grab()
        lay <- c(1,3)
        grid.newpage()
        grid.arrange(gh1,gh2,nrow=1,widths=c(unit(1000, "bigpts"),unit(300, "bigpts")))
        gh<-grid.grab()
        return(gh)
    }

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

samples_to_include=columns_to_include[columns_to_include%in%sample_metadata[,sample_names_column,drop=T]]
    anno_col=columns_to_include[columns_to_include%in%sample_metadata[,sample_names_column,drop=T]==F]

    samples_to_include <- samples_to_include[! samples_to_include %in% outlier_samples_to_remove]
    samples_to_include <- samples_to_include[samples_to_include != gene_names_column]
    samples_to_include <- samples_to_include[samples_to_include != "Gene"]
    samples_to_include <- samples_to_include[samples_to_include != "GeneName"]
    samples_to_include <- samples_to_include[samples_to_include %in% sample_metadata[[sample_names_column]]]

##create unique rownames to correctly add back Annocolumns at end of template
counts_matrix[,gene_names_column]=paste0(counts_matrix[,gene_names_column],'_',1:nrow(counts_matrix))

    anno_col=c(anno_col,gene_names_column)%>%unique
    anno_tbl=counts_matrix[,anno_col,drop=F]%>%as.data.frame

    df <- counts_matrix[,c(gene_names_column,samples_to_include)]
    gene_names <- NULL
    gene_names$GeneID <- counts_matrix[,gene_names_column]

#print(colnames(df.final))

    ### This code block does input data validation
    
    sample_metadata <- sample_metadata[match(colnames(df),sample_metadata[[sample_names_column]]),] #First match sample metadata to counts matrix
    sample_metadata <- sample_metadata[rowSums(is.na(sample_metadata)) != ncol(sample_metadata), ] # Remove empty rows
    sample_metadata <- sample_metadata[, colSums(is.na(sample_metadata)) == 0] #Remove empty columns
    rownames(sample_metadata) <- sample_metadata[[sample_names_column]]
    
    
    ### Remove specal characters from Metadata Column. Replace with _
    sample_metadata[,groups_column]=gsub('-| |!|\\*|\\.',"_",sample_metadata[,groups_column])

    
    #### remove low count genes ########
    
    df <- df[complete.cases(df),]
    ## duplicate Rows should be removed in Clean_Raw_Counts template
    #df %>% dplyr::group_by(.data[[gene_names_column]]) %>% summarise_all(sum) %>% as.data.frame() -> df
    print(paste0("Number of features before filtering: ", nrow(df)))

    ## USE CPM Transformation
     if (use_cpm_counts_to_filter == TRUE){
            trans.df=df
            trans.df[, -1]=edgeR::cpm(as.matrix(df[, -1]))
            counts_label="Filtered Counts (CPM)"
     } else {
            trans.df=df
            counts_label="Filtered Counts"

     }

    if (Use_Group_Based_Filtering == TRUE) {
        rownames(trans.df) <- trans.df[,gene_names_column]
        trans.df[,gene_names_column] <- NULL
            
            counts <- trans.df >  Minimum_Count_Value_to_be_Considered_Nonzero # boolean matrix
        
        tcounts <- as.data.frame(t(counts))
        colnum <- dim(counts)[1] # number of genes
        tcounts <- merge(sample_metadata[groups_column], tcounts, by="row.names")
        tcounts$Row.names <- NULL
        melted <- melt(tcounts, id.vars=groups_column)
        tcounts.tot <- dplyr::summarise(dplyr::group_by_at(melted, c(groups_column, "variable")), sum=sum(value))
        tcounts.tot %>% tidyr::spread(variable, sum) -> tcounts.group
        colSums(tcounts.group[(1:colnum+1)]>=Minimum_Number_of_Samples_with_Nonzero_Counts_in_a_Group) >= 1 -> tcounts.keep 
        df.filt <- trans.df[tcounts.keep, ]
        df.filt %>% rownames_to_column(gene_names_column) -> df.filt
    } else {

            trans.df$isexpr1 <- rowSums(as.matrix(trans.df[, -1]) > Minimum_Count_Value_to_be_Considered_Nonzero) >= Minimum_Number_of_Samples_with_Nonzero_Counts_in_Total

            df.filt <- as.data.frame(trans.df[trans.df$isexpr1, ])
    }

    #colnames(df.filt)[colnames(df.filt)==gene_names_column] <- "Gene"
    print(paste0("Number of features after filtering: ", nrow(df.filt)))

    ######## Start PCA ###############

    edf <- log((as.matrix(df.filt[,samples_to_include]+0.5)))
    rownames(edf) <- df.filt[,1]
    tedf <- t(edf)
    tedf <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
    tedf <- tedf[, apply(tedf, 2, var) != 0]
    pca <- prcomp(tedf, scale.=T)
    
    pcx <- paste0("PC",principal_component_on_x_axis)
    pcy <- paste0("PC",principal_component_on_y_axis)
    pca.df <- as.data.frame(pca$x) %>% dplyr::select(.data[[pcx]], .data[[pcy]])
    pca.df$group <- sample_metadata[[groups_column]]
    pca.df$sample <- sample_metadata[[labels_column]]
    perc.var <- (pca$sdev^2/sum(pca$sdev^2))*100
    perc.var <- formatC(perc.var,format = "g",digits=4)
    pc.x.lab <- paste0(pcx," ", perc.var[principal_component_on_x_axis],"%")
    pc.y.lab <- paste0(pcy," ", perc.var[principal_component_on_y_axis],"%")
    labelpos <- pca.df
    labelpos$mean_y <- pca.df[[pcy]]+label_offset_y_
    labelpos$mean_x <- pca.df[[pcx]]+label_offset_x_
    pca.df$xdata <- pca.df[[pcx]]
    pca.df$ydata <- pca.df[[pcy]]

    # Manual changes to sample names
    replacements = samples_to_rename_manually

    if (!is.null(samples_to_rename_manually)) {
        if (replacements != c("")) {
            for (x in replacements) {
                old <- strsplit(x, ": ?")[[1]][1]
                new <- strsplit(x, ": ?")[[1]][2]
                pca.df$sample <- ifelse(pca.df$sample==old, new, pca.df$sample)
            }
        }
    }

    colorlist <- c("#5954d6","#e1562c","#b80058","#00c6f8","#d163e6","#00a76c","#ff9287","#008cf9","#006e00","#796880","#FFA500","#878500")
    names(colorlist) <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
    if(length(colors_for_plots) == 0){
        colors_for_plots <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
    }
    colorval <- colorlist[colors_for_plots]
    colorval <- unname(colorval) #remove names which affect ggplot

    if (length(unique(sample_metadata[[groups_column]])) > length(colorval)) {
        ## Original color-picking code.
        k=length(unique(sample_metadata[[groups_column]]))-length(colorval)
        more_cols<- getourrandomcolors(k) 
        colorval <- c(colorval , more_cols)
    }

    if (add_labels_to_PCA == TRUE){
    g <- ggplot(pca.df, aes(x=xdata, y=ydata)) +
        theme_bw() +
        theme(legend.title=element_blank()) +
        theme(legend.position=legend_position_for_PCA) +
        geom_point(aes(color=group), size=point_size_for_pca) +
        geom_text(data=labelpos, aes(x=labelpos$mean_x, y=labelpos$mean_y, 
            label=sample, color=group, vjust="inward", hjust="inward"), size=label_font_size, show.legend=FALSE) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
        scale_colour_manual(values = colorval) +
        xlab(pc.x.lab) + ylab(pc.y.lab)
    } else {
    g <- ggplot(pca.df, aes(x=xdata, y=ydata)) +
        theme_bw() +
        theme(legend.title=element_blank()) +
        theme(legend.position=legend_position_for_PCA) +
        geom_point(aes(color=group,text=sample), size=point_size_for_pca) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
        scale_colour_manual(values = colorval) +
        xlab(pc.x.lab) + ylab(pc.y.lab)    
    }

    par(mfrow = c(2,1))

    df.m <- melt(edf,id.vars=c(gene_names_column))
    df.m = dplyr::rename(df.m,sample=Var2)

    if(set_min_max_for_x_axis_for_histogram == TRUE){
        xmin = minimum_for_x_axis_for_histogram
        xmax = maximum_for_x_axis_for_histogram
    } else {
        xmin = min(df.m$value)
        xmax = max(df.m$value)
    }

    if(color_histogram_by_group == TRUE){
        df.m %>% mutate(colgroup = sample_metadata[sample,groups_column]) -> df.m
        df.m = df.m[complete.cases(df.m[, "colgroup"]),]
        df.m$colgroup = gsub("\\s","_",df.m$colgroup)
        df.m$colgroup = factor(df.m$colgroup, levels=unique(df.m$colgroup))
        #print(unique(df.m$sample))

        # plot Density 
        g2 = ggplot(df.m, aes(x=value, group=sample)) + 
            geom_density(aes(colour = colgroup)) +
            xlab(counts_label) + ylab("Density") +
            theme_bw() +
            theme(legend.position=legend_position_for_histogram,legend.text = element_text(size = legend_font_size_for_histogram)) + 
            ggtitle("Frequency Histogram") +
            xlim(xmin,xmax) +
            #scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),40)) +
            scale_colour_manual(values=colorval)
            guides(linetype = guide_legend(ncol = number_of_histogram_legend_columns))
    } else {
        
        df.m$sample = sample_metadata[df.m$sample,labels_column]
        n=length(unique(df.m$sample))
        cols<- getourrandomcolors(n) 
        
        g2 = ggplot(df.m, aes(x=value, group=sample)) + 
            geom_density(aes(colour = sample )) +
            xlab(counts_label) + ylab("Density") +
            theme_bw() +
            theme(legend.position=legend_position_for_histogram,legend.text = element_text(size = legend_font_size_for_histogram)) +  
            ggtitle("Frequency Histogram") +
            xlim(xmin,xmax) +
            #scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),n)) +
            scale_colour_manual(values=cols) +
            guides(linetype = guide_legend(ncol = number_of_histogram_legend_columns))
    }
   
    #dev.off()

    imageWidth = 3000
    imageHeight = 1500*2
    dpi = 300

    png(
      filename="Filtered_Counts_1.png",
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")

    
    if(plot_correlation_matrix_heatmap == TRUE){
        if(interactive_plots == TRUE){
            p1=(g)%>%ggplotly(tooltip = c("sample","group"))
            p2=(g2+theme(legend.position = "none")) %>%ggplotly(tooltip = c("sample"))
            fig=subplot(p1,p2,which_layout = 'merge',margin=.05,shareX = F,shareY = F,titleY = T,titleX = T,widths=c(.5,.5),nrows = 1)
            fig=fig %>% layout(title = 'Interactive PCA and Histogram')
            print(fig)
        } else {
            require(gridExtra)
            gh<-make_heatmap(df.filt[,samples_to_include],sample_metadata,colorval)
            grid.arrange(g,g2,gh, nrow=number_of_image_rows)
            #dev.off()
        }  
    } else {
        if(interactive_plots == TRUE){
            p1=(g)%>%ggplotly(tooltip = c("sample","group"))
            p2=(g2+theme(legend.position = "none")) %>%ggplotly(tooltip = "sample" )
            fig=subplot(p1,p2,which_layout = 'merge',margin=.05,shareX = F,shareY = F,titleY = T,titleX = T,widths=c(.5,.5),nrows = 1)
            fig=fig %>% layout(title = 'Interactive PCA and Histogram')
            print(fig)
        } else {
            grid.arrange(g,g2, nrow=number_of_image_rows)
            #dev.off()
        }
        }
    
    df %>% filter(.data[[gene_names_column]] %in% df.filt[,gene_names_column]) -> df.final
   # colnames(df.final)[colnames(df.final)==gene_names_column] <- "Gene"

print('')
print('Sample Columns')
print(colnames(df.final[,!colnames(df.final)%in%gene_names_column]))
print('Annotation Columns')
print(colnames(anno_tbl))

    df.final=merge(anno_tbl,df.final,by=gene_names_column,all.y=T)
   df.final[,gene_names_column]=gsub('_[0-9]+$',"",df.final[,gene_names_column])

    return(df.final)
}

print("template_function_Filtered_Counts_1.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Clean_Raw_Counts_1<-readRDS(paste0(rds_output,"/var_Clean_Raw_Counts_1.rds"))
Input_is_Seurat_count <- 0
for(item in var_Clean_Raw_Counts_1){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Clean_Raw_Counts_1<-as.data.frame(var_Clean_Raw_Counts_1)}else{var_Clean_Raw_Counts_1 <- var_Clean_Raw_Counts_1}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_ccbr1224_rev_response_gsea_meta<-readRDS(paste0(rds_output,"/var_ccbr1224_rev_response_gsea_meta.rds"))
Input_is_Seurat_count <- 0
for(item in var_ccbr1224_rev_response_gsea_meta){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_ccbr1224_rev_response_gsea_meta<-as.data.frame(var_ccbr1224_rev_response_gsea_meta)}else{var_ccbr1224_rev_response_gsea_meta <- var_ccbr1224_rev_response_gsea_meta}
invisible(graphics.off())
var_Filtered_Counts_1<-Filtered_Counts_1(var_Clean_Raw_Counts_1,var_ccbr1224_rev_response_gsea_meta)
invisible(graphics.off())
saveRDS(var_Filtered_Counts_1, paste0(rds_output,"/var_Filtered_Counts_1.rds"))
