library(purrr)
library(dplyr)
library(ggplot2)

files <- list.files("outputs/optimal_rf/", pattern = "acc.csv$", full.names = T)

acc <- files %>%
  map_dfr(read_csv) 

pdf(file = "new_figure_scripts/rf_optimal_accuracy.pdf", width = 7.81, height = 3.14)
ggplot(acc, aes(x = reorder(study, -accuracy), y = accuracy * 100, 
                fill = set, label = round(accuracy*100, digits = 1))) +
  geom_col(position = "dodge") + 
  geom_hline(yintercept = 33.333, linetype = "dashed", color = "grey") +
  theme_minimal() +
  scale_fill_brewer(palette = "Paired", labels = c("training", "validation")) +
  geom_text(size = 4, position = position_dodge(width = 1), vjust = 2, 
            color = "white") +
  labs(y = "percent accuracy", x = "validation study") 
dev.off()
