setwd("~/Dropbox/TransCox/TransCox_package/ReproduceSimulation")

source("GetData.R")
source("TransCoxRFunctions_auxOnly.R")
library(reticulate)
source_python(system.file("python", "TransCoxFunction.py", package = "TransCox"))
library(TransCox)
library(survival)
set.seed(1345)


t_col <- function(color, percent = 50, name = NULL) {

    rgb.val <- col2rgb(color)

    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100 - percent) * 255 / 100,
                 names = name)
    invisible(t.col)
}

layout.matrix <- matrix(1:8, nrow = 4, ncol = 2, byrow = TRUE)
layout(mat = layout.matrix,
       heights = c(2,2,2,2), # Heights of the two rows
       widths = c(1, 4))
par(mar = c(2.8, 2.3, 0.5, 0.3), mgp = c(1.4, 0.5, 0))

nprim = 250
naux = 6400

learning_rate = 0.0006  
nsteps = 800


pdf(file = "FigureS5A_Simulation_SparsityPlot.pdf", width = 6, height = 4.5)
layout.matrix <- matrix(1:8, nrow = 4, ncol = 2, byrow = TRUE)
layout(mat = layout.matrix,
       heights = c(2,2,2,2), # Heights of the two rows
       widths = c(1, 4))
par(mar = c(2.8, 2.3, 0.5, 0.3), mgp = c(1.4, 0.5, 0))

for(nset in 1:4) {
    
    cat("Setting = ", nset, "\n")
    
    onedata <- GenData_FiveCov(nprim = nprim,
                               naux = naux, 
                               setting = nset,
                               Dist = "Beta")
    pData = onedata$primData
    aData = onedata$auxData

    BestLam1 = BestLam2 = 0.1
    cat(BestLam1, ";", BestLam2, "\n")
    Cout <- GetAuxSurv(aData, cov = c("X1", "X2", "X3", "X4", "X5"))
    Pout <- GetPrimaryParam(pData, q = Cout$q, estR = Cout$estR)
    
    Tres <- runTransCox_one(Pout, l1 = BestLam1, l2 = BestLam2, 
                            learning_rate = learning_rate, nsteps = nsteps, cov = c("X1", "X2", "X3", "X4", "X5"))
    
    fcolor <- ifelse(abs(Tres$eta)<0.001, t_col("blue"), t_col("black"))
    up <- max(c(0.05, Tres$eta))
    lower <- min(c(-0.05, Tres$eta))
    plot(Tres$eta, ylab = "eta", xlim = c(0,6), xaxt = "n", col = fcolor, pch = 19, ylim = c(lower, up))
    axis(1, at = c(1,2,3,4,5), labels = c(1,2,3,4,5))
    abline(h = 0, col = "red", lty = 2)
    scolor <- ifelse(abs(Tres$xi)<0.001, t_col("blue"), t_col("black"))
    plot(Tres$xi, ylab = "xi", col = scolor, pch = 19, ylim = c(-max(abs(Tres$xi)+0.05), max(abs(Tres$xi))+0.05))
    abline(h = 0, col = "red", lty = 2)
    
    if(nset == 1) {
        legend("topleft", legend = c("abs(value)<1e-3", "else"), pch = 19, col = c(t_col("blue"), t_col("black")), bty = "n")
        
    }
}

dev.off()

pdf(file = "FigureS5B_Simulation_SparsityPlot.pdf", width = 3, height = 4.5)
par(mfrow = c(1,1),mar = c(2.8, 2.3, 0.3, 0.3), mgp = c(1.4, 0.5, 0))
hist(Tres$xi[which(Tres$xi>(-0.1))], 100, xlab = "xi", main = "")
dev.off()


