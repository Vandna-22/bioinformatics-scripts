# =====================================
# Correlation Network Plot (Community Style)
# =====================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(scales)
})

# =====================================
# Input / Output Settings
# =====================================
input_file <- "your_correlation_file.csv"
out_edges <- "network_edges.csv"
out_plot <- "network_plot.png"

# Minimum absolute correlation to keep an edge
threshold <- 0.6

# Community detection method: "louvain" or "walktrap"
community_method <- "louvain"

# =====================================
# Load and Filter Correlation Table
# =====================================
# Expected columns: Source | Target | Correlation
data <- read.csv(input_file, stringsAsFactors = FALSE)
colnames(data) <- c("Source", "Target", "Correlation")

edges <- data %>%
  filter(!is.na(Correlation)) %>%
  mutate(Correlation = as.numeric(Correlation)) %>%
  filter(abs(Correlation) >= threshold) %>%
  transmute(from = Source, to = Target, corr = Correlation, abs_corr = abs(Correlation))

if (nrow(edges) == 0) {
  stop("No edges left after thresholding. Lower 'threshold' and rerun.")
}

# =====================================
# Build Network + Compute Metrics
# =====================================
network <- graph_from_data_frame(edges, directed = FALSE)
network <- simplify(network, remove.multiple = TRUE, remove.loops = TRUE,
                    edge.attr.comb = list(corr = "mean", abs_corr = "mean", "ignore"))

# Degree for node size
V(network)$degree <- degree(network)

# Community detection
if (community_method == "walktrap") {
  comm <- cluster_walktrap(network, weights = E(network)$abs_corr)
} else {
  comm <- cluster_louvain(network, weights = E(network)$abs_corr)
}

V(network)$community <- factor(membership(comm))

# Label logic: only genus names are emphasized (bold/italic).
# Heuristic for genus label: one capitalized word (e.g., "Selenomonas").
V(network)$is_genus <- str_detect(V(network)$name, "^[A-Z][a-z]+$")

# =====================================
# Plot
# =====================================
# Style target: similar to user reference
# - Node color by community
# - Edge color by correlation direction (blue negative, red positive)
# - Edge width by magnitude
# - Degree-scaled node size
# - Dense labels with only genus names emphasized
set.seed(123)

p <- ggraph(network, layout = "fr") +
  geom_edge_link(
    aes(width = abs_corr, color = corr),
    alpha = 0.35,
    show.legend = TRUE
  ) +
  scale_edge_color_gradient2(
    low = "#4575B4",
    mid = "#F7F7F7",
    high = "#D73027",
    midpoint = 0,
    name = "Correlation\nDirection"
  ) +
  scale_edge_width(
    range = c(0.2, 1.8),
    name = "Correlation\nMagnitude"
  ) +
  geom_node_point(
    aes(size = degree, fill = community),
    shape = 21,
    color = "white",
    stroke = 0.35,
    alpha = 0.95
  ) +
  # Genus nodes are overlaid with genus-specific colors so different genera
  # (e.g., Moraxella vs Selenomonas) are always visually distinct.
  geom_node_point(
    data = function(x) x %>% filter(is_genus),
    aes(color = name),
    size = 4.2,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  scale_color_discrete() +
  geom_node_text(
    aes(label = name),
    size = 2.9,
    color = "#111111",
    repel = TRUE,
    point.padding = unit(0.14, "lines"),
    box.padding = unit(0.20, "lines"),
    max.overlaps = Inf,
    segment.size = 0.2,
    segment.alpha = 0.5
  ) +
  geom_node_text(
    data = function(x) x %>% filter(is_genus),
    aes(label = name),
    fontface = "bold.italic",
    size = 5.0,
    color = "#1F1F1F",
    repel = TRUE,
    point.padding = unit(0.25, "lines"),
    box.padding = unit(0.35, "lines"),
    max.overlaps = Inf
  ) +
  scale_fill_brewer(palette = "Set1", name = "Functional\nCommunity") +
  scale_size_continuous(range = c(2.5, 14), name = "Connectivity\n(Degree)") +
  guides(
    edge_color = guide_colorbar(order = 2),
    edge_width = guide_legend(order = 1),
    fill = guide_legend(order = 3, override.aes = list(size = 6)),
    size = guide_legend(order = 4)
  ) +
  labs(
    title = "Community Analysis with Correlation Dynamics",
    subtitle = "Node colors: Functional Clusters | Edge colors: Correlation Direction (Red+, Blue-)",
    caption = "Spearman Correlation | All interaction names displayed | Node size = connectivity"
  ) +
  theme_void(base_size = 13) +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 12, color = "#444444", hjust = 0.5, margin = margin(b = 12)),
    plot.caption = element_text(size = 11, color = "#555555", hjust = 0.5, face = "italic", margin = margin(t = 8)),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.background = element_rect(fill = "#F4F4F4", color = NA),
    panel.background = element_rect(fill = "#F4F4F4", color = NA)
  )

print(p)

# =====================================
# Save Outputs
# =====================================
write.csv(edges %>% rename(Source = from, Target = to, Correlation = corr), out_edges, row.names = FALSE)
ggsave(out_plot, plot = p, width = 18, height = 14, dpi = 320, bg = "#F4F4F4")
