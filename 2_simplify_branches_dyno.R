library(dynwrap)
library(dynutils)
library(dynplot)
library(tidyverse)
library(dyno)


setwd("/home/cang/Dropbox/Projects_UCI/ZebraFishNeuralCrest/dev/publication")
zebraNC <- readRDS(file = './data/WT/Arch1_new_aggregate_filtered_noNT.RDS')

object_counts <- Matrix::t(as(as.matrix(zebraNC@assays$RNA@counts), 'sparseMatrix'))
object_expression <- Matrix::t(as(as.matrix(zebraNC@assays$RNA@data), 'sparseMatrix'))
dataset <- wrap_expression(
  counts = object_counts, 
  expression = object_expression
)
branches_fname = "./data/branches.txt"
branches <-read.csv(branches_fname, stringsAsFactors = FALSE, colClasses=c("branch_id"="character", "directed"="logical"))
branch_network_fname = "./data/branch_network.txt"
branch_network <-read.csv(branch_network_fname, stringsAsFactors = FALSE, colClasses=c("from"="character", "to"="character"))
branch_progressions_fname = "./data/branch_progressions.txt"
branch_progressions <-read.csv(branch_progressions_fname, stringsAsFactors = FALSE, colClasses=c("branch_id"="character"))
grouping_fname = "./data/grouping.txt"
grouping <-read.csv(grouping_fname, stringsAsFactors = FALSE, colClasses=c("group_id"="character"))

model_paga_tree <- add_branch_trajectory(dataset,
  grouping = grouping,
  branch_progressions = branch_progressions,
  branches = branches,
  branch_network = branch_network
)

celltime <- as.vector( read.csv("./data/WT/cell_time.csv")[,"time"] )
celltype <- as.vector( read.csv("./data/cellassign/celltypes.csv")[,"x"])

model_paga_tree <- model_paga_tree %>% add_root(root_milestone_id = "2")
model_paga_tree_simplify <- simplify_trajectory(model_paga_tree)

plot_dimred(model_paga_tree, color_density = "grouping", grouping = dynwrap::group_onto_nearest_milestones(model_paga_tree))
plot_dendro(model_paga_tree, color_cells="grouping", grouping = dynwrap::group_onto_nearest_milestones(model_paga_tree))

pseudotime <- calculate_pseudotime(model_paga_tree_simplify)

cellassign_probs <- read.csv("./data/cellassign/cellprobs_bin.csv")
cellassign_probs <- as(as.matrix(cellassign_probs), 'sparseMatrix')
rownames(cellassign_probs) <- dataset$cell_ids
dataset_cellassign <- wrap_expression(expression=cellassign_probs, counts=cellassign_probs)
plot_dendro(model_paga_tree_simplify, color_cells="feature", feature_oi="earlyNC", expression=dataset_cellassign, size_cells=1, alpha_cells=1.0)
plot_dendro(model_paga_tree_simplify, color_cells="feature", feature_oi="neural_glial", expression=dataset_cellassign, size_cells=1, alpha_cells=1.0)
plot_dendro(model_paga_tree_simplify, color_cells="feature", feature_oi="pigment", expression=dataset_cellassign, size_cells=1, alpha_cells=1.0)
plot_dendro(model_paga_tree_simplify, color_cells="feature", feature_oi="skeletal", expression=dataset_cellassign, size_cells=1, alpha_cells=1.0)
plot_dendro(model_paga_tree_simplify, color_cells="feature", feature_oi="unassigned", expression=dataset_cellassign, size_cells=1, alpha_cells=1.0)
plot_dendro(model_paga_tree_simplify, color_cells="grouping", grouping=celltime, color_milestones="auto")
plot_dendro(model_paga_tree_simplify, color_cells="grouping", grouping=celltype)
plot_dendro(model_paga_tree_simplify, color_cells="grouping", grouping=dynwrap::group_onto_nearest_milestones(model_paga_tree_simplify))

plot_dimred(model_paga_tree_simplify, color_cells="grouping", grouping=celltype)
plot_dimred(model_paga_tree_simplify, color_cells="grouping", grouping=celltime)


write.csv(model_paga_tree_simplify$milestone_network, file="./data/simplified_milestone_network.csv")
write.csv(model_paga_tree_simplify$milestone_percentages, file="./data/simplified_milestone_percentages.csv")
write.csv(model_paga_tree_simplify$progressions, file="./data/simplified_progressions.csv")

# Branching de gene analysis
branching_milestone <- model_paga_tree_simplify$milestone_network %>% group_by(from) %>% filter(n() > 1) %>% pull(from) %>% first()
branch_feature_importance <- calculate_branching_point_feature_importance(model_paga_tree_simplify, expression_source=dataset$expression, milestones_oi = branching_milestone)
branching_point_features <- branch_feature_importance %>% top_n(50, importance) %>% pull(feature_id)
plot_heatmap(
  model_paga_tree_simplify,
  expression_source = dataset$expression,
  features_oi = branching_point_features
)

# Overall de gene analysis
overall_feature_importances <- dynfeature::calculate_overall_feature_importance(model_paga_tree_simplify, expression_source=dataset)
overall_features <- overall_feature_importances %>% top_n(50, importance) %>% pull(feature_id)
plot_heatmap(
  model_paga_tree_simplify,
  expression_source = dataset$expression,
  features_oi = overall_features
)

write.csv(pseudotime, file = "./pseudotime.csv")
write.csv(model_paga_tree_simplify$cell_ids, file="./cellids.csv")
