# [R] Get Package Info  (bdc9cdf9-c835-42e1-a084-bcf72b756ec1): v3
Session_Info <- function() {
    library(sessioninfo)
    df = as.data.frame(package_info())
    df = df[df$attached == TRUE, c("package", "loadedversion", "date", "source")]
    return(df)
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

print("template_function_Session_Info.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
invisible(graphics.off())
var_Session_Info<-Session_Info()
invisible(graphics.off())
saveRDS(var_Session_Info, paste0(rds_output,"/var_Session_Info.rds"))
