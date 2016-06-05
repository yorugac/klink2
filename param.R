# three semantic relations according to Klink-2 paper
# similarityLink is a temporal structure for computation of relatedEquivalent
semantic <- c("relatedEquivalent", "broaderGeneric", "contributesTo", "similarityLink")

# Input relations taken into consideration:
# relation 1: if 2 keywords are used in the same publication
# relation 2: if 2 keywords are used by the same author in the same year
# relation 3: if 2 keywords are used in the same venue (name of journal and so on) in the same year
# relation 4: if 2 keywords are classified as belonging to the same research area in the same year
# then there is a co-occurrence in regards of corresponding relation.
relations <- c("publication", "author", "venue", "area")
quantified <- c(FALSE, FALSE, FALSE, FALSE)
rn <- length(relations)

# verbosity level
# 0 - no messages
# 1 - main statistics per iteration
# 2 - notifications for the start and end of key procedures
# 3 - word-by-word messaging
verbosity <- 2

# what is the minimum connection strength for keywords to be considered related?
relkeyT <- 1

## metric params
# weights for linear combination of n measure (string similarity)
nweights <- c(1, 1, 1, 1)
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
