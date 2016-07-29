# Prepare dataset as input to Klink-2.

source('relations.R')
library(Rcpp)
sourceCpp('utils.cpp')

# input relations taken into consideration
relations <- c("publication", "author", "venue", "area")
# number of relations
rn <- length(relations)
# length of saved cooccurrence values
m <- 100

# 4 input variables to Klink-2:
keywordsdb <- new.env(parent=globalenv(), hash=TRUE)
reldb_df <- list()
reldb_l <- list()
inputm <- NULL

# reads dataset into global raw input variables: keywordsdb, reldb_df, reldb_l
# limit - read only limited number of articles
# NOTE: on full data set can take up to 2 hours
read_dataset <- function(limit=-1) {
    # solar data set
    d <- read.csv('ac2012.tsv', sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=limit*2)

    empty <- function(elem) is.null(elem) || elem==""

    # need only items with defined keywords fields (DE, ID)
    w <- c()
    for(i in 1:dim(d)[1])
        if(!(empty(d$DE[i]) && empty(d$ID[i]))) w=c(i,w)
    d = d[w,]

    d = d[1:limit,]

    # fields used: 2 keywords fields, document title, authors, publication name, research areas, year
    # NOTE: these fields are already checked for empty in case of solar dataset,
    # but must be tailored for another dataset
    d = d[c("DE", "ID", "TI", "AU", "SO", "SC", "PY")]

    k <- 1
    a <- 1
    n <- nrow(d)
    process_item <- function(item) {
        cat('process article ', a, ' out of ', n, '\n'); a <<- a + 1
        keywords <- unique(unlist(sapply(strsplit(tolower(item[c("DE", "ID")]), ";"), trimws)))
        newkeywords = setdiff(keywords, ls(keywordsdb))
        oldkeywords = setdiff(keywords, newkeywords)

        authors <- unlist(sapply(strsplit(item["AU"], ";"), trimws))
        areas <- unique(unlist(sapply(strsplit(tolower(item["SC"]), ";"), trimws)))
        venues <- tolower(item["SO"])
        relation <- c(1, rep(2, length(authors)), 3, rep(4, length(areas)))
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

    invisible(apply(d, 1, process_item))
}

# populates inputm
# NOTE: costly function O(n^2); currently on full data set can take several days
cache_cooccurrence <- function() {
    # number of keywords
    n = length(ls(keywordsdb))
    inputm <<- matrix(0, nrow=m, ncol=n*2*rn)

    for(i in 1:n) {
        cat('process keyword ', i, ' out of ', n, '\n')
        irel = reldb_l[[i]]
        for(r in 1:rn) {
            co_m <- calc.cooccurrence(relations[r], i, irel)
            inputm[, cached.keys(i, r)] <<- co_m[,1]
            inputm[, cached.values(i, r)] <<- co_m[,2]
        }
    }

    # sorting in decreasing order of relation values
    for(i in seq(1,ncol(inputm)-1,by=2)) {
        ind <- order(inputm[, i+1], decreasing=TRUE)
        inputm[, i] <<- inputm[, i][ind]
        inputm[, i+1] <<- inputm[, i+1][ind]
    }
}

# How to prepare all input variables for Klink-2 in one go:
run_all <- function(limit=-1) {
    read_dataset(limit)
    cache_cooccurrence()
    save('reldb_df', 'reldb_l', 'keywordsdb', 'inputm', file=paste('input', limit, '.Rdata', sep=''))
}
# To preserve time, functions can be run separately so long as input variables are saved and loaded.

# prints main information about input variables
# filename must contain input variables
inspect_dataset <- function(filename) {
    load(filename)
    n <- length(reldb_df)
    cat('Number of keywords: ', n, '\n')
    for(i in 1:rn) {
        cat('Value of co-occurrence,', relations[i], 'relation: [', 
            min(inputm[, cached.values(1:n, i)]), ', ', 
            max(inputm[, cached.values(1:n, i)]), ']\n')
    }
}
