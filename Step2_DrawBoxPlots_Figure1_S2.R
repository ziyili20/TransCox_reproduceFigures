setwd("~/Dropbox/TransCox/TransCox_package/ReproduceSimulation")

alldir <- list.files("~/Dropbox/TransCox/TransCox_package/ReproduceSimulation")
mydir <- alldir[grep(".RData", alldir)]
mydir

ss1 <- sapply(strsplit(mydir, split = "_"), "[[", 3)
ss2 <- gsub(".RData", "", sapply(strsplit(mydir, split = "_"), "[[", 4))

addterm <- c("5Cov_Beta_Figure1", "5Cov_Norm_FigureS2")

## When k = 1, it output boxplots for Figure 1
## When k = 2, it output boxplots for Figure S2
for(k in 1:2) {
    addterm1 <- addterm[k]
    
    load(mydir[k])
    head(fullBetaBias[[1]])
    
    fulldat <- c()
    for(i in 1:4) {
        oneres <- fullBetaBias[[i]]
        
        # oneres2 <- oneres[, c(1,3,5,2,4,6)]
        colnames(oneres) <- c("Trans-Cox:b1", "Trans-Cox:b2", "Trans-Cox:b3",
                              "Cox-Tonly:b1", "Cox-Tonly:b2", "Cox-Tonly:b3",
                              "Cox-Both:b1", "Cox-Both:b2", "Cox-Both:b3",
                              "Cox-Str:b1", "Cox-Str:b2", "Cox-Str:b3")
        
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
    mycolor = c(RColorBrewer::brewer.pal(9, "Paired")[c(1,3,7,5,2,4,8,6)], 
                c("blue","darkolivegreen4","darkgoldenrod","darkred"))
    # mycolor2 = rep(rep(rev(mycolor), each = 100), 4)
    fulldat$Method = factor(fulldat$key, levels = c("Cox-Tonly:b1",
                                                    "Cox-Str:b1",
                                                    "Cox-Both:b1",  
                                                    "Trans-Cox:b1",
                                                    "Cox-Tonly:b2",
                                                    "Cox-Str:b2",
                                                    "Cox-Both:b2",
                                                    "Trans-Cox:b2",
                                                    "Cox-Tonly:b3",
                                                    "Cox-Str:b3",
                                                    "Cox-Both:b3",
                                                    "Trans-Cox:b3"))
    p1 <- ggplot(fulldat, aes(x=Method, y=value, fill=Method)) + 
        geom_boxplot() + xlab("") + ylab("Estimation bias for beta1/beta2/beta3")  + 
        facet_wrap(~setting)+ ylim(c(-1,1))+
        scale_fill_manual(values=mycolor)+ 
        geom_hline(yintercept=0, linetype="dashed", color = "gray") + 
        theme_bw()+
        theme(axis.text.x=element_blank(), #remove x axis labels
              axis.ticks.x=element_blank()
        ) + geom_vline(xintercept = 4.5, linetype=1, color = "lightgray") +
        geom_vline(xintercept = 8.5, linetype=1, color = "lightgray")
    p1   #+geom_jitter(color="lightgray", size=0.2, alpha=0.9)
    ggsave(plot = p1, filename = paste0("BoxplotForBias_", ss1[k], "_", ss2[k], "_b1b2_bias_", addterm1, ".pdf"), 
           width = 11.5, height = 5.5)
    
    
    
    fulldat <- c()
    for(i in 1:4) {
        oneres <- fullHtwoBias[[i]]
        colnames(oneres) <- c("Trans-Cox:H1", "Trans-Cox:H2", 
                              "Cox-Tonly:H1", "Cox-Tonly:H2", 
                              "Cox-Both:H1", "Cox-Both:H2", 
                              "Cox-Str:H1", "Cox-Str:H2")
        
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
    mycolor = c(RColorBrewer::brewer.pal(9, "Paired")[c(1,3,7,5,2,4,8,6)])
    # mycolor2 = rep(rep(rev(mycolor), each = 100), 4)
    fulldat$Method = factor(fulldat$key, levels = c("Cox-Tonly:H1",  
                                                    "Cox-Both:H1", 
                                                    "Cox-Str:H1",
                                                    "Trans-Cox:H1", 
                                                    "Cox-Tonly:H2",
                                                    "Cox-Both:H2",
                                                    "Cox-Str:H2",
                                                    "Trans-Cox:H2"))
    p1 <- ggplot(fulldat, aes(x=Method, y=value, fill=Method)) + 
        geom_boxplot() + xlab("") + ylab("Cumulative baseline hazards bias at t=0.6,1.2")  + 
        facet_wrap(~setting)+
        scale_fill_manual(values=mycolor)+ 
        geom_hline(yintercept=0, linetype="dashed", color = "gray") + 
        theme_bw()+ ylim(c(-1,1)) +
        theme(axis.text.x=element_blank(), #remove x axis labels
              axis.ticks.x=element_blank()
        )+ geom_vline(xintercept = 4.5, linetype=1, color = "lightgray") 
    p1   #+geom_jitter(color="lightgray", size=0.2, alpha=0.9)
    ggsave(plot = p1, filename = paste0("BoxplotForBias_", ss1[k], "_", ss2[k], "_H1H2_bias_", addterm1, ".pdf"), 
           width = 10.5, height = 5.5)
    
}



