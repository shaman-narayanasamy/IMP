#!/usr/bin/env Rscript
.libPaths("/mnt/nfs/projects/ecosystem_biology/local_tools/R/library")

args=commandArgs(trailingOnly=TRUE)
data_file = args[1]
out_dir=args[2]
data = read.table(data_file, sep=",", comment.char='', header=T)
colnames(data) = c("BIN_NO", "NumUnique", "NumMultiple")

library(ggplot2)
library(gridExtra)
require(cowplot)

#op=par() #commenting this out prevented spurious Rplots.pdf from being created
out_pdf = paste(out_dir, "analysis_essential_genes.pdf", sep="/")
pdf(file = out_pdf, width=10,height=5.5)

element_text_size = 12

xlab = "Total number of essential genes per Bin"
ylab = "Number of essential genes in\nmultiple copies per Bin"
scatter = ggplot(data = data, aes(x = NumUnique + NumMultiple, y = NumMultiple))
scatter = scatter + 
  geom_point(size=3, alpha=0.5) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.title.y=element_text(vjust=1.4)) + 
  theme(axis.title.x=element_text(vjust=0.1)) + 
  scale_y_continuous(breaks = seq(0,110, by=20)) +
  scale_x_continuous(breaks = seq(0,110, by=20)) +
  xlab(xlab) +
  ylab(ylab)

hist_left = ggplot(data, aes(x = NumMultiple)) +
  stat_density(position="identity", geom="line", trim=T, size=1) + 
  coord_flip() +
  scale_y_reverse() +
  scale_x_continuous(breaks = seq(0,110, by=20)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.title.y=element_text(vjust=1.4)) + 
  theme(axis.title.x=element_text(vjust=0.1)) + 
  xlab(ylab) +
  ylab("Density") +
  theme(legend.position = "none")

hist_top = ggplot(data, aes(x = NumUnique + NumMultiple)) +
  stat_density(position="identity", geom="line", size=1) + 
  scale_x_continuous(breaks = seq(0,110, by=20)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.title.y=element_text(vjust=1.4)) + 
  theme(axis.title.x=element_text(vjust=0.1)) + 
  xlab(xlab) +
  ylab("\nDensity") +
  theme(legend.position = "none")

g = ggdraw() +
  draw_plot(hist_left, 0, 0, 0.25, 0.75) +
  draw_plot(scatter, 0.25, 0, 0.74, 0.75) +
  draw_plot(hist_top, 0.25, 0.74, 0.74, 0.25)
  #draw_plot(hist_left, 0, 0, 0.26, 0.70) +
  #draw_plot(scatter, 0.27, 0, .70, .70) +
  #draw_plot(hist_top, .27, 0.70 , .70 , .30)
g

dev.off()
