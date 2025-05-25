
# 清空环境
rm(list=ls())
# 设置Bioconductor镜像
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
library(pROC)
library(ggplot2)

# Nature-style color palette
nature_colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4")

plot_roc_curves <- function(risk_file, method, output_dir = ".") {
  # Read risk matrix
  risk_data <- read.table(risk_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  
  # Extract cohort information
  risk_data$Cohort <- sub("(.*)_.*_.*", "\\1", rownames(risk_data))
  risk_data$Cohort <- sub("(.*)\\.(.*)", "\\1", risk_data$Cohort)
  
  # Function to create ROC plot
  create_roc_plot <- function(data, method, cohort) {
    y <- ifelse(grepl("CS2$", rownames(data)), 0, 1)
    roc_obj <- roc(y, as.numeric(data[[method]]))
    ci <- ci.auc(roc_obj, method = "bootstrap")
    
    roc_data <- data.frame(
      specificity = roc_obj$specificities,
      sensitivity = roc_obj$sensitivities
    )
    
    ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity)) +
      geom_path(color = nature_colors[1], size = 1.5) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      annotate("text", x = 0.65, y = 0.25, 
               label = sprintf("AUC: %.3f\n95%% CI: %.3f-%.3f", 
                               ci[2], ci[1], ci[3]),
               hjust = 0, vjust = 0, size = 4) +
      labs(title = cohort,
           x = "1 - Specificity",
           y = "Sensitivity") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(fill = NA, color = "black", size = 0.5))
  }
  
  # Create ROC curve for each cohort
  for (cohort in unique(risk_data$Cohort)) {
    cohort_data <- risk_data[risk_data$Cohort == cohort, ]
    roc_plot <- create_roc_plot(cohort_data, method, cohort)
    
    ggsave(filename = file.path(output_dir, paste0("ROC_", cohort, ".pdf")),
           plot = roc_plot, width = 5, height = 4.75, units = "in")
  }
}

# Usage
plot_roc_curves("riskMatrix15.txt", "RF")

