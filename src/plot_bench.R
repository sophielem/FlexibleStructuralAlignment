Data = read.csv("results/bench_res.csv", sep = ",", header = T)
colorList = c("deepskyblue", "gray31", "darkgoldenrod1", "deeppink", "blue", "pink")
pdf("results/bench_res.pdf")
par(mar = c(8,3,2,2), mgp = c(7,1,0))
plot(Data[,2], type = "b", lty = 1, pch = 1, col = "deepskyblue", ylim = c(0,1), axes = F, xlab = "Names of protein", ylab = "TMscore", main = "TMscore for all algorithms")
lines(Data[,3], type = "b", lty = 1, pch = 1, col = "gray31")
lines(Data[,4], type = "b", lty = 1, pch = 1, col = "darkgoldenrod1")
lines(Data[,5], type = "b", lty = 1, pch = 1, col = "deeppink")
lines(Data[,6], type = "b", lty = 1, pch = 1, col = "blue")
lines(Data[,7], type = "b", lty = 1, pch = 1, col = "pink")
abline(a = 0.17, b = 0, type = "l", lty = 1, col = "red")
legend(x = 0, y = 1, c("X1_2", "X2_1", "gdt1_2", "gdt2_1"),cex = 1,pch = 15, col = colorList)
axis(side = 1,at = 1:36, labels = Data[,1], las = 2)
axis(2)
dev.off()

pdf("results/boxplot.pdf")
boxplot(Data[,-1], col = c("deepskyblue", "gray31", "darkgoldenrod1", "deeppink", "blue", "pink"))
dev.off()
