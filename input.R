# Prepare dataset as input to Klink-2.

# solar dataset
d <- read.csv('ac2012.tsv', sep="\t", header=TRUE, stringsAsFactors=FALSE)

empty <- function(elem) {
    is.null(elem) || elem==""
}

# need only items with defined keywords fields (DE, ID)
w <- c()
for(i in 1:dim(d)[1])
    if(!(empty(d$DE[i]) && empty(d$ID[i]))) w=c(i,w)
d = d[w,]

# fields used: 2 keywords fields, document title, authors, publication name, research areas, year
# NOTE: these fields are already checked for empty in case of solar dataset,
# but must be tailored for another dataset
d = d[c("DE", "ID", "TI", "AU", "SO", "SC", "PY")]

# input relations taken into consideration
relations <- c("publication", "author", "venue", "area")
# number of relations
rn <- length(relations)

keywordsdb <- new.env(parent=globalenv(), hash=TRUE)
reldb_df <- list()
reldb_l <- list()
k <- 1

# for testing purposes
d = d[1:5000,]
# d = d[ceiling(runif(5000, min=0, max=dim(d)[1]-1)),]

# length of saved relation values
m <- 100

process_item <- function(item) {
    keywords <- unique(unlist(sapply(strsplit(tolower(item[c("DE", "ID")]), ";"), trimws)))
    newkeywords = setdiff(keywords, ls(keywordsdb))
    oldkeywords = setdiff(keywords, newkeywords)

    authors <- unlist(sapply(strsplit(item["AU"], ";"), trimws))
    areas <- unique(unlist(sapply(strsplit(tolower(item["SC"]), ";"), trimws)))
    venues <- tolower(item["SO"])
    relation <- factor(c("publication", rep("author", length(authors)), "venue", rep("area", length(areas))))
    entity <- c(item["TI"], authors, venues, areas)
    quantity <- rep(NA_integer_, length(entity))
    year <- as.numeric(rep(item["PY"], length(entity)))

    for(i in seq_along(newkeywords)) {
        keywordsdb[[newkeywords[i]]] <<- k
        reldb_l[[k]] <<- list()
        reldb_l[[k]]$publication <<- item["TI"]
        reldb_l[[k]]$author <<- authors
        reldb_l[[k]]$venue <<- venues
        reldb_l[[k]]$area <<- areas
        names(reldb_l)[k] <<- newkeywords[i]
        reldb_df[[k]] <<- data.frame(relation, entity, quantity, year, stringsAsFactors=FALSE)
        names(reldb_df)[k] <<- newkeywords[i]
        k <<- k + 1
    }
    for(i in seq_along(oldkeywords)) {
        index <- keywordsdb[[oldkeywords[i]]]
        reldb_l[[index]]$publication <<- c(reldb_l[[index]]$publication, item["TI"])
        reldb_l[[index]]$author <<- unique(c(reldb_l[[index]]$author, authors))
        reldb_l[[index]]$venue <<- unique(c(reldb_l[[index]]$venue, venues))
        reldb_l[[index]]$area <<- unique(c(reldb_l[[index]]$area, areas))
        reldb_df[[index]] <<- rbind(reldb_df[[index]], data.frame(relation, entity, quantity, year, stringsAsFactors=FALSE))
    }
}

apply(d, 1, process_item)

# number of keywords
n = length(ls(keywordsdb))

inputm <- matrix(0, nrow=m, ncol=n*2*rn)
source('relations.R')
# O(n^2):
for(i in 1:length(reldb_l)) {
    cat("process ", i, "\n")
    irel = reldb_l[[i]]
    for(r in 1:rn) {
        rname <- relations[r]
        co_m <- calc.cooccurrence(rname, i, irel)
        inputm[, 2*r-1 + 2*rn*(i-1)] = co_m[,1]
        inputm[, 2*r + 2*rn*(i-1)] = co_m[,2]
    }
}
# sorting in decreasing order of relation values
for(i in seq(1,ncol(inputm)-1,by=2)) {
    ind <- order(inputm[, i+1], decreasing=TRUE)
    inputm[, i] = inputm[, i][ind]
    inputm[, i+1] = inputm[, i+1][ind]
}
save('reldb_df', 'reldb_l', 'keywordsdb', 'inputm', file="input5000.Rdata")
