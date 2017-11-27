#' @function survival analysis according to gene expression for primary tumor
#' @param geneSymbol identify gene symbol 
#' @param exp TCGA format gene expression dataset
#' @param  cli Clinical information dataset
#' @param method method use to dive samples into groups. Options are "quantile", "median", "mean". the "quantile" use first and third quartile as threshold
#' @param trans transform the clinical IDs from '-' separate to '.' separate
#' @return include the plot and p value
#' @author Shixiang Wang

# input test data
# exp <- exp.tumor
# cli <- cli_info
# geneSymbol <- "PD1"
# method="quantile"
#
library(survminer)
library(survival)

surv_analysis <- function(geneSymbol=NULL, exp, cli, method="quantile", trans=FALSE, step=10){
    if(any(c(is.null(geneSymbol), is.null(exp), is.null(cli)))){
        stop("Wrong data! Please check your input.")
    }
    
    if(trans){
        cli$sampleID = gsub("-",".",cli$sampleID, fixed = TRUE)
    }
    tumorSample <- cli[cli$sample_type=="Primary Tumor",]$sampleID
    # exp <- exp[, c(1, na.omit(match(tumorSample, colnames(exp) )))]
    
    tumorSampleID <- tumorSample[ na.omit(match(colnames(exp), tumorSample)) ]
    symbol_exp <- exp[exp$sample==geneSymbol, tumorSampleID]
    exp_value <- unlist(symbol_exp)
    
    group_surv <- list()
    
    # compute threshold according to method parameter
    if(method=="quantile"){
        symbol_stat <- summary(exp_value)
        ths1 <- as.numeric(symbol_stat[2])
        ths2 <- as.numeric(symbol_stat[5])
        down_id <- names(exp_value[exp_value<ths1])
        up_id <- names(exp_value[exp_value>ths2])
    }else if(method=="mean"){
        ths <- mean(exp_value)
        down_id <- names(exp_value[exp_value<=ths])
        up_id <- names(exp_value[exp_value>ths])
    }else if(method=="median"){
        ths <- median(exp_value)
        down_id <- names(exp_value[exp_value<=ths])
        up_id <- names(exp_value[exp_value>ths])
    }else{
        stop("The method you input have something wrong! Here only 'quantile', 'mean', 'median' provided.  Please check one of them.")}
    
    # compute the cutoff which p value is minimal
    
    os_mat <- subset(cli, sampleID%in%c(up_id, down_id))   
    os_mat$group <- "high"
    os_mat[os_mat$sampleID%in%down_id,]$group <- "low"
    
    search_cutoff <- function(mat, step=20){
        # mat <- os_mat
        days <- max(na.omit(mat$X_OS))
        
        cutoff_list <- {}
        pval_list <- {}
        cutoff <- list()
        
        while(days > 300){
            mat <- mat[mat$X_OS<=days,]
            sdf <- survdiff(Surv(X_OS, X_EVENT)~group,data = mat)
            p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
            cutoff_list <- c(cutoff_list, days)
            pval_list <- c(pval_list, p.val)
            days <- days - step
        }
        
        min_pval <- min(pval_list)
        min_cutoff <- cutoff_list[which(pval_list==min_pval)]
        
        cutoff$cutoff_list <- cutoff_list
        cutoff$pval_list <- pval_list
        cutoff$min_pval <- min_pval
        cutoff$min_cutoff <- min_cutoff
        
        return(cutoff)
    
    }
    
    cutoff <- search_cutoff(os_mat, step=step)

    surv <- list()
    surv$os_mat <- os_mat
    surv$cutoff <- cutoff
    
    return(surv)
}

plot_surv <- function(os_mat, cutoff, pval=TRUE){
    fit <- survfit(Surv(X_OS, X_EVENT) ~ group, data = os_mat[os_mat$X_OS<=cutoff,])
    ggsurv <- ggsurvplot(
        fit,                     # survfit object with calculated statistics.
        data = os_mat[os_mat$X_OS<=cutoff,],           # data used to fit survival curves.
        risk.table = FALSE,       # show risk table.
        pval = pval,             # show p-value of log-rank test.
        conf.int = FALSE,         # show confidence intervals for
        # point estimates of survival curves.
        palette = c("red", "blue"),
        # palette = c("#E7B800", "#2E9FDF"),
        # xlim = c(0,5000),         # present narrower X axis, but not affect
        # survival estimates.
        xlab = "",   # customize X axis label.
        ylab = "",
        # break.time.by = 100,     # break X axis in time intervals by 500.
        # ggtheme = theme_light(), # customize plot and risk table with a theme.
        risk.table.y.text.col = T,# colour risk table text annotations.
        risk.table.height = 0.25, # the height of the risk table
        risk.table.y.text = FALSE,# show bars instead of names in text annotations
        # in legend of risk table.
        ncensor.plot = FALSE,      # plot the number of censored subjects at time t
        ncensor.plot.height = 0.25,
        # conf.int.style = "step",  # customize style of confidence intervals
        # surv.median.line = "hv",  # add the median survival pointer.
        legend.labs =
            c("APOBEC3B high expression", "APOBEC3B low expression"),    # change legend labels.
        legend.title = "",
        size = 0.5,
        censor.size = 2
     
    )
    
    ggsurv <- ggpar(
        ggsurv,
        font.caption = c(6, "plain", "black"),
        font.tickslab = c(8, "plain", "black"),
        font.legend = c(6, "plain", "black"),
        ggtheme = theme_survminer()+theme(legend.background=element_rect(fill=NA), legend.key.height = unit(3,"mm"),
                                          line = element_line(size=.1), legend.position =c(0.8, 0.9))
    )
    
    ggsurv$plot
}

