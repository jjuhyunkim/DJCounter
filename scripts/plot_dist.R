library(ggplot2)

dat <- read.table("DJ_counts.txt", header = F)
colnames(dat) <- c("Sample", "DJcount", "PeakCount", "Peak_Est")

ggplot(data = dat, aes(x = "DJ", y = Peak_Est)) +
    geom_jitter(shape=1,
     position=position_jitter(seed = 42, width = 0.35, height = 0),
     alpha = 0.3, size = 0.3) +
    geom_violin(trim = F, width = 0.5, alpha = 0, linewidth = 0.2) +
    scale_y_continuous(breaks = seq(0, 14, 1), limits = c(0, 14)) +
    # Set color pallette to set1
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    labs(title="Num. of DJs from k-mers",
	   x="",
	   y="copy number estimate")

ggsave("DJ_counts.png", height = 3.5, width = 3)
ggsave("DJ_counts.pdf", height = 3.5, width = 3)
