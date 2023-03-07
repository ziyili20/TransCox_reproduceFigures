setwd("~/Dropbox/TransCox/TransCox_package/ReproduceSimulation")

alldir <- list.files("~/Dropbox/TransCox/TransCox_package/ReproduceSimulation")
mydir <- alldir[grep(".RData", alldir)]
mydir

ss1 <- sapply(strsplit(mydir, split = "_"), "[[", 3)
ss2 <- gsub(".RData", "", sapply(strsplit(mydir, split = "_"), "[[", 4))

addterm <- c("5Cov_Beta_Figure2", "5Cov_Norm_FigureS3")

## When k = 1, it output boxplots for Figure 2
## When k = 2, it output boxplots for Figure S3
for(k in 1:2) {
    addterm1 <- addterm[k]
    
    load(mydir[k])
    head(fullAbsPersPred[[1]])
    
    fulldat <- c()
    for(i in 1:4) {
        oneres <- fullAbsPersPred[[i]]
        # oneres2 <- oneres[, c(1,3,5,2,4,6)]
        colnames(oneres) <- c("Trans-Cox:APPE1", "Trans-Cox:APPE2", 
                              "Cox-Tonly:APPE1", "Cox-Tonly:APPE2",
                              "Cox-Both:APPE1", "Cox-Both:APPE2",
                              "Cox-Str:APPE1", "Cox-Str:APPE2")
        
        onetmp <- tidyr::gather(as.data.frame(oneres))
        head(onetmp)
        onetmp$setting = paste0("Setting ", i)
        if(i == 1) {
            fulldat = onetmp
        } else {
            fulldat = rbind(fulldat, onetmp)
        }
    }
    
    library(ggplot2)
    mycolor = RColorBrewer::brewer.pal(9, "Paired")[c(1,3,7,5,2,4,8,6)]
    # mycolor2 = rep(rep(rev(mycolor), each = 100), 4)
    fulldat$Method = factor(fulldat$key, levels = c("Cox-Tonly:APPE1",
                                                    "Cox-Both:APPE1",
                                                    "Cox-Str:APPE1",  
                                                    "Trans-Cox:APPE1",
                                                    "Cox-Tonly:APPE2",
                                                    "Cox-Both:APPE2",
                                                    "Cox-Str:APPE2",
                                                    "Trans-Cox:APPE2"))
    fulldat$value <- log(fulldat$value)
    p1 <- ggplot(fulldat, aes(x=Method, y=value, fill=Method)) + 
        geom_boxplot() + xlab("") + ylab("log(Absolute Personalized Prediction Error)")  + 
        facet_wrap(~setting)+
        scale_fill_manual(values=mycolor) + 
        theme_bw()+
        theme(axis.text.x=element_blank(), #remove x axis labels
              axis.ticks.x=element_blank()
        ) + geom_vline(xintercept = 4.5, linetype=1, color = "lightgray")
    p1   #+geom_jitter(color="lightgray", size=0.2, alpha=0.9)
    ggsave(plot = p1, filename = paste0("BoxplotForlogAPPE_", ss1[k], "_", ss2[k], "_", addterm1, ".pdf"), 
           width = 10.5, height = 5.5)
}
