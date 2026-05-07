library(ggplot2)

L1 <- 1e6
chr_info <- list(
     list(ts_path="dev/testData/msprime_chr0.trees",
     breaks=c(0, L1/2, L1), rates=c(1e-8, 2e-8)))
founderGenomes <- asMapPop(
     chr_info = chr_info,
     inbred = FALSE,
     ploidy = 2L
     )
all_pos <- chrKeptPosBpList
all_pos <- unlist(all_pos)
breaks <- seq(min(all_pos), max(all_pos), length.out = 11)
bg_counts <- cut(all_pos, breaks = breaks, include.lowest = TRUE) %>%
     table() %>%
     as.numeric()

n_iterations <- 50
n_bins <- 10
segSites <- 77

chr_info <- list(
  list(ts_path="dev/testData/msprime_chr0.trees",
       breaks=c(0, L1/2, L1), rates=c(1e-8, 2e-8), segSites=segSites))

all_bp_positions <- list()
all_gen_map_positions <- list()

#set.seed(42)
random_seeds <- sample(1:1000000, size = n_iterations)

for (i in 1:n_iterations) {
  current_seed <- random_seeds[i]

  founderGenomes <- asMapPop(
    chr_info = chr_info,
    inbred = FALSE,
    ploidy = 2L,
    site_sampling_seed = current_seed
  )

  bp_pos <- unlist(chrKeptPosBpList)
  all_bp_positions[[i]] <- bp_pos

  gen_map_pos <- unlist(founderGenomes@genMap)
  all_gen_map_positions[[i]] <- gen_map_pos
}


df_bp <- data.frame(pos = unlist(all_bp_positions))
df_gen <- data.frame(pos = unlist(all_gen_map_positions))

sampled_counts <- cut(df_bp$pos, breaks = breaks, include.lowest = TRUE) %>%
     table() %>%
     as.numeric()

plot_data <- data.frame(
   bin_mid = (breaks[-1] + breaks[-length(breaks)]) / 2,
   sampling_ratio = (sampled_counts) / bg_counts
)

p1 <- ggplot(plot_data, aes(x = bin_mid/1000, y = sampling_ratio)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "white") +
  geom_hline(yintercept = mean(plot_data$sampling_ratio, na.rm = TRUE),
             linetype = "dashed", color = "red") +
  labs(title = paste(segSites, "segSites (", length(all_pos), "in total) over", n_iterations, "runs"),
       x = "Position (kbp)",
       y = "Frequency (Count) / Background") +
  theme_minimal()
print(p1)
