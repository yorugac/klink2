# three semantic relations according to Klink-2 paper
# similarityLink is a temporal structure for computation of relatedEquivalent
semantic <- c("relatedEquivalent", "broaderGeneric", "contributesTo", "similarityLink")

# input relations taken into consideration
relations <- c("publication", "author", "venue", "area")
quantified <- c(FALSE, FALSE, FALSE, FALSE)
rn <- length(relations)

# verbosity level
# 0 - no messages
# 1 - main statistics per iteration
# 2 - notifications for the start of key procedures
# 3 - word-by-word messaging
verbosity <- 2

# what is the minimum connection strength for keywords to be considered related?
relkeyT <- 1

# weights for linear combination of n measure (string similarity)
nweights <- c(1, 1, 1, 1)

## metric params
# threshold for hierarchical metrics
tR <- c(1, 1, 1, 1)
# threshold for hierarchical indicators, i.e. how many should point in the same direction
th <- 2
# threshold for relatedEquivalent metric
tS <- 0.5
# threshold for relatedEquivalent indicators, i.e. how many should be positive
tre <- 2
# coefficient for T metric
gamma <- 2 # must be > 0

## clustering params
# clustering threshold for mergeSimilarWords
merge_t <- 2
# clustering threshold for intersectBasedClustering
intersect_t <- 1
# clustering threshold for quickHierarchicalClustering
quick_t <- 0.8

## filter params
# number of main keywords
nmain <- 20
# co-occurrence coverage by main keywords
maincover <- 0.15
