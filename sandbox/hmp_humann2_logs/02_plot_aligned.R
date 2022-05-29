# Plot HUMANn2 percent unaligned

setwd("~/github/2020-ibd/")

library(dplyr)
library(ggplot2)
library(scales)
unaligned <- read.table("sandbox/hmp_humann2_logs/OAR.txt", 
                        stringsAsFactors = F, sep = "%")

ggplot(unaligned, aes(x = V1)) +
  geom_density() + 
  scale_y_continuous(labels = percent_format()) +
  theme_minimal() +
  xlim(0, 100) +
  labs(title = "iHMP Metagenomes Percent Reads Aligned by HUMANn2",
       x = "percent of reads mapped",
       y = "percent of samples",
       caption = "mean = 59.6\nsd = 12.7\n1330 metagenome samples from iHMP with HUMANn2 log files") +
  theme(plot.caption = element_text(size=12))

mean(unaligned$V1)
sd(unaligned$V1)

unaligned$fraction <- unaligned$V1 / 100

ggplot(unaligned, aes(x = V1)) +
  geom_density() + 
  scale_y_continuous(labels = percent_format()) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), 
                     labels = c("0%", "25%", "50%", "75%", "100%")) + 
  theme_minimal() +
  labs(x = "percent of reads mapped",
       y = "percent of samples")



ggplot(unaligned, aes(x = fraction)) +
  geom_density() + 
  scale_y_continuous(labels = percent_format()) +
  scale_x_continuous(labels = percent_format(), limits = c(0, 1)) +
  theme_minimal() +
  labs(x = "percent of reads mapped",
       y = "percent of samples")
