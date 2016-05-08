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

# length of saved relation values
m <- 100

process_item <- function(item) {
    keywords <- unique(unlist(sapply(strsplit(tolower(item[c("DE", "ID")]), ";"), trimws)))
    newkeywords = setdiff(keywords, ls(keywordsdb))
    oldkeywords = setdiff(keywords, newkeywords)

    authors <- unlist(sapply(strsplit(item["AU"], ";"), trimws))
    areas <- unique(unlist(sapply(strsplit(tolower(item["SC"]), ";"), trimws)))
    venues <- tolower(item["SO"])
    relation <- c("publication", rep("author", length(authors)), "venue", rep("area", length(areas)))
    entity <- c(item["TI"], authors, venues, areas)
    quantity <- rep(NA_integer_, length(entity))
    years <- rep(item["PY"], length(entity))

    for(i in seq_along(newkeywords)) {
        keywordsdb[[newkeywords[i]]] <<- k
        reldb_l[[k]] <<- list()
        reldb_l[[k]]$publication <<- item["TI"]
        reldb_l[[k]]$author <<- authors
        reldb_l[[k]]$area <<- areas
        reldb_l[[k]]$venue <<- venues
        names(reldb_l)[k] <<- newkeywords[i]
        reldb_df[[k]] <<- data.frame(relation, entity, quantity, years, stringsAsFactors=FALSE)
        names(reldb_df)[k] <<- newkeywords[i]
        k <<- k + 1
    }
    for(i in seq_along(oldkeywords)) {
        index <- keywordsdb[[oldkeywords[i]]]
        reldb_l[[index]]$publication <<- c(reldb_l[[index]]$publication, item["TI"])
        reldb_l[[index]]$author <<- unique(c(reldb_l[[index]]$author, authors))
        reldb_l[[index]]$area <<- unique(c(reldb_l[[index]]$area, areas))
        reldb_l[[index]]$venue <<- unique(c(reldb_l[[index]]$venue, venues))
        reldb_df[[index]] <<- rbind(reldb_df[[index]], data.frame(relation, entity, quantity, years, stringsAsFactors=FALSE))
    }
    # for(i in seq_along(keywords)) {
    #     index <- keywordsdb[[keywords[i]]]
    #     for(j in seq_along(keywords)){
    #         if(i != j) {
    #             jindex <- keywordsdb[[keywords[j]]]
    #             # publication relation
    #             r <- 1
    #
    #             # is there j word in the stored list?
    #             jm = match(j, inputm[, 2*r-1 + 2*rn*(i-1)])
    #             if(!is.na(jm)) {
    #                 inputm[jm, 2*r + 2*rn*(i-1)] <<- inputm[jm, 2*r + 2*rn*(i-1)] + 1
    #             } else {
    #                 # is there any 0 values to store publication with j?
    #                 jm = match(0, inputm[, 2*r + 2*rn*(i-1)])
    #                 if(!is.na(jm)) {
    #                      inputm[jm, 2*r-1 + 2*rn*(i-1)] <<- j
    #                      inputm[jm, 2*r + 2*rn*(i-1)] <<- 1
    #                 }
    #             }
    #         }
    #     }
    # }
    # if(k == n) {
    #     inputm <<- cbind(inputm, matrix(0, nrow=m, ncol=n*2*rn))
    # }
}

apply(d, 1, process_item)

# number of keywords
n = length(ls(keywordsdb))

inputm <- matrix(0, nrow=m, ncol=n*2*rn)

# O(n^2):
for(i in seq_along(reldb)) {
    cat("process ", i, "\n")
    irel = reldb_l[[i]]
    for(j in seq_along(reldb)) {
        if(i != j) {
            jrel = reldb_l[[j]]
            for(r in 1:rn) {
                rel = relations[r]
                cooccur = irel[[rel]] %in% jrel[[rel]]
                v = length(which(cooccur))
                # is there j word in the stored list?
                im = match(j, inputm[, 2*r-1 + 2*rn*(i-1)])
                if(!is.na(im)) {
                    inputm[im, 2*r + 2*rn*(i-1)] = inputm[im, 2*r + 2*rn*(i-1)] + v
                } else {
                    # is there place in the store list for the value with j?
                    im = which.min(inputm[, 2*r + 2*rn*(i-1)])
                    if(v > inputm[im, 2*r + 2*rn*(i-1)]) {
                         inputm[im, 2*r-1 + 2*rn*(i-1)] = j
                         inputm[im, 2*r + 2*rn*(i-1)] = v
                    }
                }
            }
        }
    }
}
# sorting in ascending order of keywords
for(i in seq(1,ncol(inputm)-1,by=2)) {
    ind <- order(inputm[, i])
    inputm[, i] = inputm[, i][ind]
    inputm[, i+1] = inputm[, i+1][ind]
}
save('reldb_df', 'keywordsdb', 'inputm', file="input5000.Rdata")
