### draw cumulative baseline hazard (FigureS1)

t = seq(0, 2, by = 0.005)

pH <- t^2
sH <- 2/3*t^2


pdf("Target_Source_CumBaseHaz_FigureS1.pdf", width = 5, height = 4)

par(mar = c(2.8, 2.4, 1.3, 0.5), mgp = c(1.4, 0.5, 0), mfrow=c(1,1))

plot(x = t, y = pH, type = "l", lty = 1, col = "darkgreen", lwd = 3, xlab = "Time", 
     ylab = "Cumulative Baseline Hazards", main = "Simulation settings 3 and 4")
points(x = t, y = sH, type = "l", lty = 2, col = "salmon", lwd = 3)
legend("topleft", legend = c("Target cohort", "Source cohort"), col = c("darkgreen", "salmon"), lwd = 3, lty = 1:2,
       bty = "n", text.col = c("darkgreen", "salmon"))

dev.off()
