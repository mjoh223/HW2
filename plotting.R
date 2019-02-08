library(ggplot2)


dat = read.csv("/Users/matt/OneDrive/UCSF/algorithms/hw2/pca_clustering.csv", header = FALSE)
x = dat[1]
y = dat[2]
lab = dat[3]
x2 = data.frame(x,y,lab)
plot(x2[1:2], col = as.numeric(unlist(lab)))

ggplot(dat, aes(x=V2,y=V3,color=factor(V4))) +
  geom_point(size=3) +
  xlab("PCA 1") +
  ylab("PCA 2")
