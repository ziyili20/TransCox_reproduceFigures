setwd("~/Dropbox/TransCox/TransCox_package/ReproduceSimulation")

alldir <- list.files("~/Dropbox/TransCox/TransCox_package/ReproduceSimulation")
mydir <- alldir[grep(".RData", alldir)]
mydir

ss1 <- sapply(strsplit(mydir, split = "_"), "[[", 3)
ss2 <- gsub(".RData", "", sapply(strsplit(mydir, split = "_"), "[[", 4))

addterm <- c("5Cov_Beta_Figure3", "5Cov_Norm_FigureS4")

## When k = 1, it output boxplots for Figure 3, left and right panel
## When k = 2, it output boxplots for Figure S4, left and right panel
addterm <- c("cov5_Beta", "cov5_Norm")

for(k in 1:2) {
    addterm1 <- addterm[k]
    load(mydir[k])
    
    pdf(paste0("SurvCurv", ss1[k], "_", ss2[k], "_", addterm1, ".pdf"), width = 6, height = 7)
    par(mar = c(2.8, 2.4, 1.3, 0.5), mgp = c(1.4, 0.5, 0), mfrow=c(2,2))
    for(j in 1:4) {
        
        time = seq(0.15, 2, length.out = 50)
        trueH = time^2
        TransM <- colMeans(SurvC_Trans[[j]], na.rm = TRUE)
        qfh <- function(x){
            quantile(x, 0.975, na.rm = TRUE)
        }
        qfl <- function(x){
            quantile(x, 0.025, na.rm = TRUE)
        }
        upper <- apply(SurvC_Trans[[j]], 2, qfh)
        lower <- apply(SurvC_Trans[[j]], 2, qfl)
        
        CB <- colMeans(SurvC_Both[[j]])
        upperB <- apply(SurvC_Both[[j]], 2, qfh)
        lowerB <- apply(SurvC_Both[[j]], 2, qfl)
        
        mycol = RColorBrewer::brewer.pal(5, "Set2")
        
        plot(CB~time, type = "l", col = "blue", main = paste0("Setting ", j),
             xlim = c(0.1, 1.7), ylab = c("Survival Probability"), xlab = "Time")
        points(time, upperB, type = "l", lty = 2, col = mycol[3])
        points(time, lowerB, type = "l", lty = 2, col = mycol[3])
        polygon(x = c(time, rev(time)),
                y = c(upperB, rev(lowerB)),
                col =  adjustcolor(mycol[3], alpha.f = 0.20), border = NA)
        
        lines(TransM~time, type = "l", col = "green")
        points(time, upper, type = "l", lty = 2, col = mycol[5])
        points(time, lower, type = "l", lty = 2, col = mycol[5])
        polygon(x = c(time, rev(time)),
                y = c(upper, rev(lower)),
                col =  adjustcolor(mycol[5], alpha.f = 0.20), border = NA)
        
        TB <- colMeans(SurvC_True[[j]])
        lines(time, TB, col = "red", type = "l")
        
        if(j == 1) {
            legend("topright", bty = "n", col = c("green", "blue", "red"),
                   legend = c("Trans-Cox", "Cox-Both", "True"), lty = 1)
        }
        
    }
    dev.off()
    
}



for(k in 1:2) {
    addterm1 <- addterm[k]
    
    load(mydir[k])
    
    pdf(paste0("CumHaz", ss1[k], "_", ss2[k], "_", addterm[k], ".pdf"), width = 6, height = 7)
    par(mar = c(2.8, 2.4, 1.3, 0.5), mgp = c(1.4, 0.5, 0), mfrow=c(2,2))
    for(j in 1:4) {
        
        time = seq(0.15, 2, length.out = 50)
        trueH = time^2
        TransM <- colMeans(CumHaz_Trans[[j]], na.rm = TRUE)
        qfh <- function(x){
            quantile(x, 0.975, na.rm = TRUE)
        }
        qfl <- function(x){
            quantile(x, 0.025, na.rm = TRUE)
        }
        upper <- apply(CumHaz_Trans[[j]], 2, qfh)
        lower <- apply(CumHaz_Trans[[j]], 2, qfl)
        
        CB <- colMeans(CumHaz_Both[[j]])
        upperB <- apply(CumHaz_Both[[j]], 2, qfh)
        lowerB <- apply(CumHaz_Both[[j]], 2, qfl)
        
        mycol = RColorBrewer::brewer.pal(8, "Set2")
        
        plot(CB~time, type = "l", col = "blue", main = paste0("Setting ", j),
             xlim = c(0.1, 1.7), ylab = c("Cumulative Baseline Hazards"), xlab = "Time")
        points(time, upperB, type = "l", lty = 2, col = mycol[3])
        points(time, lowerB, type = "l", lty = 2, col = mycol[3])
        polygon(x = c(time, rev(time)),
                y = c(upperB, rev(lowerB)),
                col =  adjustcolor(mycol[3], alpha.f = 0.20), border = NA)
        
        lines(TransM~time, type = "l", col = "green")
        points(time, upper, type = "l", lty = 2, col = mycol[5])
        points(time, lower, type = "l", lty = 2, col = mycol[5])
        polygon(x = c(time, rev(time)),
                y = c(upper, rev(lower)),
                col =  adjustcolor(mycol[5], alpha.f = 0.20), border = NA)
        
        lines(time, time^2, col = "red", type = "l")
        
        if(j == 1) {
            legend("topleft", bty = "n", col = c("green", "blue", "red"),
                   legend = c("Trans-Cox", "Cox-Both", "True"), lty = 1)   
        }
    }
    dev.off()
    
}

