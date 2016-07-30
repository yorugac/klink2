# Prepare dataset as input to Klink-2.

source('relations.R')
Rcpp::sourceCpp('utils.cpp')

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
inputm <- matrix(, nrow=m, ncol=0)

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

    if(limit > 0) d = d[1:limit,]

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

        authors <- as.vector(unlist(sapply(strsplit(item["AU"], ";"), trimws)))
        areas <- unique(unlist(sapply(strsplit(tolower(item["SC"]), ";"), trimws)))
        venues <- tolower(item["SO"])
        relation <- c(1, rep(2, length(authors)), 3, rep(4, length(areas)))
        entity <- c(item["TI"], authors, venues, areas)
        quantity <- rep(NA_integer_, length(entity))
        year <- as.numeric(rep(item["PY"], length(entity)))

        for(i in seq_along(newkeywords)) {
            keywordsdb[[newkeywords[i]]] <<- k
            reldb_l[[k]] <<- list()
            reldb_l[[k]]$publication <<- paste(item["TI"], year[1], sep="_")
            reldb_l[[k]]$author <<- vapply(authors, paste, "", year[1], sep="_")
            reldb_l[[k]]$venue <<- vapply(venues, paste, "", year[1], sep="_")
            reldb_l[[k]]$area <<- vapply(areas, paste, "", year[1], sep="_")
            names(reldb_l)[k] <<- newkeywords[i]
            reldb_df[[k]] <<- data.frame(relation, entity, quantity, year, stringsAsFactors=FALSE)
            names(reldb_df)[k] <<- newkeywords[i]
            k <<- k + 1
        }
        for(i in seq_along(oldkeywords)) {
            index <- keywordsdb[[oldkeywords[i]]]
            reldb_l[[index]]$publication <<- c(reldb_l[[index]]$publication, paste(item["TI"], year[1], sep="_"))
            reldb_l[[index]]$author <<- unique(c(reldb_l[[index]]$author, vapply(authors, paste, "", year[1], sep="_")))
            reldb_l[[index]]$venue <<- unique(c(reldb_l[[index]]$venue, vapply(venues, paste, "", year[1], sep="_")))
            reldb_l[[index]]$area <<- unique(c(reldb_l[[index]]$area, vapply(areas, paste, "", year[1], sep="_")))
            reldb_df[[index]] <<- rbind(reldb_df[[index]], data.frame(relation, entity, quantity, year, stringsAsFactors=FALSE))
        }
    }

    apply(d, 1, process_item)

    for(i in 1:length(reldb_l)) {
        if(!is.null(reldb_l[[i]]$publication)) reldb_l[[i]]$publication <<- sort(reldb_l[[i]]$publication)
        if(!is.null(reldb_l[[i]]$author)) reldb_l[[i]]$author <<- sort(reldb_l[[i]]$author)
        if(!is.null(reldb_l[[i]]$venue)) reldb_l[[i]]$venue <<- sort(reldb_l[[i]]$venue)
        if(!is.null(reldb_l[[i]]$area)) reldb_l[[i]]$area <<- sort(reldb_l[[i]]$area)
    }
}

# range of number of entities one keyword can be associated with,
# regardless of relation
entities_range <- function(reldb_l) {
    a <- Inf; b <- 0
    for(i in 1:length(names(reldb_l)))
        for(r in 1:rn) {
            t <- length(reldb_l[[i]][[r]])
            if(t < a) a = t
            if(t > b) b = t
        }
    c(a, b)
}

# returns populated inputm
# NOTE: costly function O(n^2); currently on full data set can take up to day
cache_cooccurrence <- function() {
    # number of keywords
    n = length(ls(keywordsdb))
    inputm <- matrix(0, nrow=m, ncol=n*2*rn)
    maxsize <- entities_range(reldb_l)[2]

    for(i in 1:n) {
        cat('process keyword ', i, ' out of ', n, '\n')
        irel = reldb_l[[i]]
        for(r in 1:rn) {
#             co_m <- calc.cooccurrence(relations[r], i, irel)
            co_m <- calc_cooccurrence_C(n, m, i, r, reldb_l, maxsize)
            inputm[, cached.keys(i, r)] = co_m[, 1]
            inputm[, cached.values(i, r)] = co_m[, 2]
        }
    }
    inputm
}

# How to prepare all input variables for Klink-2 in one go:
# with negative limit all articles will be read,
# with posiive only this number of articles will be read.
run_all <- function(limit=-1) {
    # just in case
    keywordsdb <<- new.env(parent=globalenv(), hash=TRUE)
    reldb_df <<- list()
    reldb_l <<- list()

    read_dataset(limit)
    inputm <<- cache_cooccurrence()
    if (limit) fname <- paste('input', limit, '.Rdata', sep='')
    else fname <- 'input.Rdata'
    save('reldb_df', 'reldb_l', 'keywordsdb', 'inputm', file=fname)
    cat('Input variables saved to', fname, '\n')
}
# To preserve time, functions can be run separately so long as input variables are saved and loaded.

# prints main information about input variables
# filename must contain input variables
inspect_dataset <- function(filename) {
    load(filename)
    n <- length(reldb_df)
    entrange <- entities_range(reldb_l)
    cat('Number of keywords: ', n, '\n')
    cat('Range of entities associated with a keyword: [', entrange[1], ',', entrange[2], ']\n')
    for(i in 1:rn) {
        cat('Value of co-occurrence,', relations[i], 'relation: [',
            min(inputm[, cached.values(1:n, i)]), ',',
            max(inputm[, cached.values(1:n, i)]), ']\n')
    }
}
