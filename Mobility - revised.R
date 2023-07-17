# Loading packages
library(readr)
library(data.table)
library(summarytools)
library(hdm)
library(glmnet)
library(DoubleML)
library(mlr3)
library(ranger)
library(ggplot2)
library(ggpubr)

######################## DATA CLEANING #############################

# Load GDIM at country level
url = "https://datacatalogfiles.worldbank.org/ddh-published/0050771/DR0065670/GDIM_2023_03.csv"
data_raw = read.csv(url)

# Data pre-processing
data_raw = data_raw[,c(1,4,5,8,12:14,26:29,31,40)] #removing unnecessary var
df = data.frame(
  mobility = data_raw$CAT,
  inequality = data_raw$SDc - data_raw$SDp,
  expansion = data_raw$MEANc - data_raw$MEANp,
  dependency = data_raw$COR, # parental dependency
  cohort = as.factor(data_raw$cohort), # generations
  fragile = as.factor(data_raw$fragile),
  developing = as.factor(ifelse(data_raw$incgroup2 == "Developing economies", "Yes", "No")),
  region = as.factor(data_raw$region_noHICgroup),
  mom = as.factor(ifelse(data_raw$parent == "mom", "Yes", "No")),
  daughter = as.factor(ifelse(data_raw$child == "daughter", "Yes", "No"))
)

view(dfSummary(df, 
               plain.ascii  = FALSE, 
               style        = "grid", 
               graph.magnif = 0.75, 
               graph.col = FALSE,
               valid.col    = FALSE))

######################## DESCRIPTIVE ANALYSIS #############################

# Scatter plotting with continuous X
plot_scatter <- function(x_var, y_var) {
  ggscatter(df, x_var, y_var,
            color = "grey",
            add = "reg.line",
            add.params = list(color = "blue", fill = "lightgray"),
            conf.int = TRUE) +
    stat_cor(method = "pearson",label.y = 1.1) + ylim(0, 1.1)
}
ggarrange(plot_scatter("expansion", "mobility"),
          plot_scatter("inequality", "mobility"),
          plot_scatter("dependency", "mobility"),
          ncol = 3, nrow = 1)

# Violin plotting with categorical X
plot_violin <- function(x_var, y_var) {
  ggviolin(df, x = x_var, y = y_var, fill = x_var,
           legend.position = "none",
           add = "boxplot", add.params = list(fill = "white")) +
    stat_compare_means(method = "anova", label.y = 1.1) + guides(fill = FALSE)
}

ggarrange(plot_violin("mom", "mobility"),
          plot_violin("daughter", "mobility"), 
          ncol = 2, nrow = 1)

ggarrange(ggarrange(plot_violin("developing", "mobility"),
                    plot_violin("fragile", "mobility"), 
                    ncol = 2, nrow = 1), 
          plot_violin("region", "mobility") + rotate_x_text(30),
          ncol = 1, nrow = 2, heights=c(1,1.5))

######################## MODELLING #############################

# Convert to the modeled data
data <- df
data$mobility <- with(data_raw, CAT - ave(CAT, country)) # remove fixed effects
offset <- abs(min(data[sapply(data, is.numeric)])) + 1
data[sapply(data, is.numeric)] <- log(data[sapply(data, is.numeric)] + offset)
data <- model.matrix(~ . - 1, data = data) # convert into model data
data <- data[, -which(colnames(data) == "cohort1940")] # transform to dummy-coding
colnames(data) <- gsub("[^[:alnum:]]", "", colnames(data)) # correct var names
view(dfSummary(data, 
               plain.ascii  = FALSE, 
               style        = "grid", 
               graph.magnif = 0.75, 
               graph.col = FALSE,
               valid.col    = FALSE))

# Q1
## Rigorous Lasso without interaction
partiall_out <- function(dt,y, dvar) {
  x <- as.matrix(data)[,-c(which(colnames(dt) == y), which(colnames(dt) == dvar))]
  y <- data[, y]
  d <- data[, dvar]
  effect <- rlassoEffect(x, y, d)
  result <- summary(effect)
  result_df <- as.data.frame(result$coefficients)
  rownames(result_df)[rownames(result_df) == "d1"] = dvar
  return(result_df)
}

Q1_dl = rbind(partiall_out(data,"mobility","inequality"),
              partiall_out(data,"mobility","expansion"),
              partiall_out(data,"mobility","dependency")
)

## Double machine learning - PLR model
fit_dml_plr <- function(dt, y, d) {
  dt = as.data.table(dt)
  dml_data <- DoubleMLData$new(dt, y_col = y, d_cols = d,
                               x_cols = colnames(dt)[!(colnames(dt) %in% c(y,d))])
  set.seed(123) # required to replicate sample split
  learner_l <- lrn("regr.ranger", num.trees = 500, min.node.size = 2, max.depth = 5)
  learner_m <- lrn("regr.ranger", num.trees = 500, min.node.size = 2, max.depth = 5)
  dml_plr <- DoubleMLPLR$new(dml_data,
                             ml_l = learner_l,
                             ml_m = learner_m)
  dml_plr$fit()
  result <- dml_plr$summary()
  return(result)
}
Q1_plr = rbind(fit_dml_plr(data,"mobility","inequality"),
               fit_dml_plr(data,"mobility","expansion"),
               fit_dml_plr(data,"mobility","dependency"))

# Q2 
## Rigorous Lasso with interaction
rlasso_effects <- function(dt, y, x) {
  dt = as.data.frame(dt)
  # Get all column names except x and y
  all_cols <- setdiff(names(dt), c(x, y))
  # Model matrix
  x_names <- paste(all_cols, collapse = "+")
  formula <- as.formula(paste0("~ -1 + ", x, "+", x, ":(", x_names, ")+(", x_names, ")^2"))
  X <- model.matrix(formula, data = dt)
  X <- X[, which(apply(X, 2, var) != 0)] #exclude constant variables 
  index_x <- grep(x, colnames(X)) 
  effects <- rlassoEffects(x = X, y = dt[, y], index = index_x)
  result <- summary(effects)
  result_df <- as.data.frame(result$coefficients)
  return(result_df)
}
Q2_mom_dl = rlasso_effects(data,"mobility","momYes")
Q2_dau_dl = rlasso_effects(data,"mobility","daughterYes")

## Double machine learning - IRM model
fit_dml_irm <- function(dt, y, d) {
  dt = as.data.table(dt)
  dml_data <- DoubleMLData$new(dt, y_col = y, d_cols = d,
                               x_cols = colnames(dt)[!(colnames(dt) %in% c(y,d))])
  set.seed(123) # required to replicate sample split
  learner_g <- lrn("regr.ranger", num.trees = 500, min.node.size = 2, max.depth = 5)
  learner_m <- lrn("classif.ranger", num.trees = 500, min.node.size = 2, max.depth = 5)
  dml_irm <- DoubleMLIRM$new(dml_data,
                             ml_g = learner_g,
                             ml_m = learner_m)
  dml_irm$fit()
  return(dml_irm$summary())
}
Q2_mom_irm = fit_dml_irm(data,"mobility","momYes")
Q2_dau_irm = fit_dml_irm(data,"mobility","daughterYes")

######################## REPORTING #############################
# table formatting
est_tab <- function(dt, model_name) {
  colnames(dt) <- c("estimate", "std_error", "t_value", "p_value")
  dt <- data.frame(model = model_name, var = rownames(dt), dt)
  rownames(dt) <- NULL
  return(dt)
}
Q1 = rbind(est_tab(Q1_dl,"Partialling-out Lasso"), est_tab(Q1_plr,"Double ML"))
Q2_mom = rbind(est_tab(Q2_mom_dl,"Partialling-out Lasso"), est_tab(Q2_mom_irm,"Double ML"))
Q2_dau = rbind(est_tab(Q2_dau_dl,"Partialling-out Lasso"), est_tab(Q2_dau_irm,"Double ML"))

write.csv(Q1, "Q1.csv", row.names = FALSE)
write.csv(Q2_mom, "Q2_mom.csv", row.names = FALSE)
write.csv(Q2_dau, "Q2_dau.csv", row.names = FALSE)

# Plotting results
plot_results <- function(dt) {
  results = dt
  # calculate the lower and upper bounds of the 95% confidence intervals
  results$conf_low <- results$estimate - 1.96*results$std_error
  results$conf_high <- results$estimate + 1.96*results$std_error
  # create a factor variable for the p-value significance levels
  results$p_value_level <- cut(
    results$p_value,
    breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
    labels = c("*** (p < 0.01)", "** (p < 0.05)", "* (p < 0.1)", " (p >= 0.1)")
  )
  # plot the estimates with error bars and size the points according to the p-value levels
  ggplot(results, aes(x = var, y = estimate)) +
    geom_pointrange(aes(ymin = conf_low, ymax = conf_high), color = "darkblue") +
    geom_point(aes(size = p_value_level)) +
    geom_text(aes(label = sprintf("%.3f", estimate)), nudge_x = 0.3, size=3) +
    scale_x_discrete("") + geom_hline(yintercept = 0, color = "red") +
    theme_bw() + ylab(NULL) + 
    scale_size_discrete(name = "p-value",limits = rev(levels(results$p_value_level))) +
    coord_flip() + facet_grid(~model)
}
plot_results(Q1)
plot_results(Q2_mom)
plot_results(Q2_dau)
