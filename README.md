# scaPCoA

# scaPCoA: Supervised and Covariate Adjusted Principal CoOrdinate Analysis

Given that the gut microbiome has been shown to play a critical role in human health and dis-
ease, discerning the relationship between microbial features and disease outcomes is often of interest
in microbiome studies. Most often, this is in terms of visualization of trends in the data and/or
predictive modeling. Dimension reduction techniques such as principal coordinate analysis (PCoA),
an alternative to PCA that relies on non-Euclidean pairwise distances between samples, provides a
lower-dimensional representation for visualization but does not incorporate outcome-related informa-
tion. Therefore, the lower-dimensional principal coordinates may have poor performance visualizing
outcome differences, as well as poor performance in downstream disease prediction tasks. More-
over, covariates, such as demographics, lifestyle factors, and technical variation, can obscure true
biological signals, leading to spurious associations. To address this, we propose a supervised and
covariate-adjusted PCoA algorithm. By minimizing the influence of nuisance variables, our method
produces lower-dimensional representations that better separate disease outcomes in visualizations
and result in improved accuracy in regression and classification models. We demonstrate the effec-
tiveness of the proposed method and highlight its ability to extract meaningful microbiome-disease
associations while mitigating nuisance effects through simulations and application to real 16S rRNA
and metagenomic sequencing studies.


## Installation
You can install the development version of scaPCoA:
```
remotes::install_github("crystalzang/scaPCoA")
```

## Usage

1. Visualization
   `scapcoa` can be used to visualize principal coordinate.
2. Prediction
   `scapcoa` can be used to make prediction on outcome (continuous or binary) given microbiome and covariates.


## Main Function
```
obj = generate_data(seed = 1, n =100, p = 200, l=3, y_type = "bin", zinLDA_min = 100, zinLDA_max = 1000)
lambda_ls = seq(0, 10, 0.5)
results = scapcoa(obj, lambda_ls, nPC = 3, pcoa=TRUE, residualization = TRUE)

# Estimated principal coordinates
pc <- results$PC.scores

# Eigen values
eigen_values <- results$eigen_values
```

## Results
```

pc %>%
  as.data.frame() %>%
  mutate(Y = ifelse(obj$Y == 1, "Case", "Control")) %>%
  ggplot(aes(x = V1, y = V2, color = Y, shape = Y)) +
  geom_point(size = 3.5, alpha = 1) +
  stat_ellipse(aes(group = Y), linetype = 1) +
  labs(
    x = "Principal Coordinate 1",
    y = "Principal Coordinate 2",
    shape = "Outcome (Y)",
    color = "Outcome (Y)",
    title = "scaPCoA: Principal Coordinate Plot"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Case" = "#C84D6A", "Control" = "#A6C36F"))
```

## Authors

**Crystal Zang**, Lu Tang, Rebecca Deek,

Department of Biostatistics and Data Science, University of Pittsburgh School of Public Health
