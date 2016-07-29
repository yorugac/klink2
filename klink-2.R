source('param.R')
source("relations.R")
source("graph.R")
library(stringi)
library(hash)
library(fastcluster)
library(Rcpp)
sourceCpp('utils.cpp')

# Var naming here:
# r - relation, string
# k, x, y - keywords, string

# output semantic relations
triples <- data.frame(k1=character(0), k2=character(0), relation=numeric(0), stringsAsFactors=FALSE)
# working semantic relations
semrel <- list()
# iteration regulator
continue <- TRUE

# for S measure scaling
largestS <- rep(0, rn)

# filled during inference
cachedS <- list()
cachedCV <- list() # conn.vec structure

# I_r(x,y) conditional probability that
# an element associated with x will be associated with y
I.prob <- function(r, x, y, diachronic=FALSE) {
    x = keyword.index(x); y = keyword.index(y)
    v <- rel.value(x, y, r)
    if(v==0) return(0)
    if(!quantified[r] && diachronic) {
        ce <- common.entities(x, y, r)
        v = v * ((ent_year_C(get.reldf(y), ce) - min(ent_year_C(get.reldf(x), ce)) + 1)^(-gamma))
    }
    # TODO quantified relation
    sum(v)
}

# longest common substring
lcs <- function(x, y) {
    trafos <- as.vector(attr(adist(x, y, counts=T), "trafos"))
    subs <- stri_sub(y, stri_locate_all_regex(trafos, "M+")[[1]])
    subs[which.max(nchar(subs))]
}

# percentage of identical words
identical.words <- function(x, y) {
    wordsx <- stri_split_fixed(x, " ")[[1]]
    wordsy <- stri_split_fixed(y, " ")[[1]]
    length(intersect(wordsx, wordsy)) / max(length(wordsx), length(wordsy))
}

# search for acronyms, separated by dots
# NOTE: does not check if acronyms are matching
# All keywords are lowered during input handling,
# so only presence of dots can serve as acronym detector
have.acronym <- function(x, y) {
    # acronym must be at least 2 words, so find at least two dots
    length(stri_split_fixed(x, ".")[[1]]) > 1 &&
      length(stri_split_fixed(y, ".")[[1]]) > 1
}

# n(x,y) measure
string.similarity <- function(x, y) {
    sum(nweights * c(length(lcs(x,y)), identical.words(x,y), common_chars_C(x,y), have.acronym(x,y)))
}

# c_r(x,y) measure
# vx, vy as arguments for optimization purpose, matrices returned by conn.vector
semantic.similarity <- function(vx, vy) {
    if(length(vx) == 0 || length(vy) == 0) return(0)
    s <- semantic_similarity_C(vx, vy) / (sqrt(sum(vx[,2]^2)) * sqrt(sum(vy[,2]^2)))
    if(is.nan(s)) s = 0
    s
}

# vx, vy - matrices
H.metric <- function(r, x, y, vx, vy, diachronic=FALSE, stringsim=NULL) {
    m <- (I.prob(r, x, y, diachronic) / I.prob(r, x, x, diachronic) - I.prob(r, y, x, diachronic) / I.prob(r, y, y, diachronic)) *
        semantic.similarity(vx, vy) * stringsim
    if(is.nan(m)) m = 0
    m
}

# vx, vy - lists with 3 matrices
S.metric <- function(r, x, y, vx=NULL, vy=NULL) {
    if(is.null(vx)) vx = conn.vector(r, x)
    if(is.null(vy)) vy = conn.vector(r, y)
    semantic.similarity(vx[[1]], vy[[1]]) /
        (max(semantic.similarity(vx[[2]], vy[[2]]), semantic.similarity(vx[[3]], vy[[3]])) + 1)
}

# infer relations
infer <- function(x, y) {
    # returns whether metric1 prevail over metric2
    prevalence <- function(metric1, metric2) {
        # is this correct understanding of prevalence?
        cmp <- mapply(function(a, b){
                if(a > 0 && b > 0 && a > b) return(TRUE)
                if(a < 0 && b < 0 && a < b) return(TRUE)
                FALSE
            }, metric1, metric2)
        # metric vectors are expected to be of the same size
        length(which(cmp)) > floor(length(metric1) / 2)
    }

    # for each realtion holds {-1, 0, 1} for
    # {x super of y, no hierarchy, x sub-area of y} respectively
    hierarchy <- c()
    # for each relation holds metric values
    hmetrics <- c()
    tmetrics <- c()
    smetrics <- c()
    stringsim <- string.similarity(x, y)
    for(r in 1:rn) {
        vx <- conn.vector(r, x)
        vy <- conn.vector(r, y)
        h <- H.metric(r, x, y, vx[[1]], vy[[1]], stringsim=stringsim)
        hmetrics = c(hmetrics, h)
        if (h >= tR[r]) hierarchy = c(hierarchy, 1)
        else if (h <= -tR[r]) hierarchy = c(hierarchy, -1)
        else hierarchy = c(hierarchy, 0)
        tmetrics = c(tmetrics, H.metric(r, x, y, vx[[1]], vy[[1]], diachronic=TRUE, stringsim=stringsim))
        s <- S.metric(r, x, y, vx, vy)
        if(s > largestS[r]) largestS[r] <<- s
        smetrics = c(smetrics, s)
    }

    if(is.null(cachedS[[x]])) cachedS[[x]] <<- hash()
    cachedS[[x]][y] <<- smetrics

    # is there enough to infer hierarchical relation?
    if (length(which(hierarchy == 1)) >= th)
        hierarchy = 1
    else if (length(which(hierarchy == -1)) >= th)
        hierarchy = -1
    else
        hierarchy = 0
    # determine type of hierarchical relation
    if (hierarchy == 1 || hierarchy == -1) {
        if (age(x) > age(y) && nentities(x) > nentities(y) && prevalence(tmetrics, hmetrics)) {
            # broaderGeneric
            if(hierarchy == 1)
                set.semantic(x, semantic[2], y)
            else
                set.semantic(y, semantic[2], x)
        } else {
            # contributesTo
            if(hierarchy == 1)
                set.semantic(y, semantic[3], x)
            else
                set.semantic(x, semantic[3], y)
        }
    }
    # "similarityLink"
    if(length(which(smetrics > tS)) >= tre) {
        set.semantic(x, semantic[4], y)
    }
}

# break loops for broaderGeneric
fix.loops <- function() {
    cleanup.semrel()
    semmatrix <- get.semantic(2)
    if(nrow(semmatrix) > 1) {
        semmatrix = delete.cycles(semmatrix)
        semrel[[2]][1:nrow(semmatrix), ] <<- semmatrix
        semrel$sizes[2] <<- nrow(semmatrix)
    }
}

fast.expand <- function(v1, v2)
    cbind(rep.int(v1, length(v2)),
          rep.int(v2, rep.int(length(v1), length(v2))))

# keywords: vector of keyword ids;
# based on cachedS value
distance.matrix.cached <- function(keywords) {
    keynames <- vapply(keywords, keyword.name, "")
    n <- length(keywords)
    distances <- matrix(0, nrow=n, ncol=n)
    for(i1 in 1:n) {
        for(i2 in i1:n) {
            if(i1 != i2) {
                if(!is.null(keynames) && !is.null(cachedS[[keynames[i1]]]) && has.key(keynames[i2], cachedS[[keynames[i1]]])) {
                    s <- sum(values(cachedS[[keynames[i1]]], keys=keynames[i2]) / largestS)
                } else {
                    s <- sum(vapply(1:rn, function(r) {
                            S.metric(r, keywords[i1], keywords[i2]) / largestS[r]
                        }, 0))
                }
                s = 1 - s / rn
                distances[i1,i2] = s
                distances[i2,i1] = s
            }
        }
    }
    diag(distances) = Inf
    distances
}

# keywords: list of keyword objects; does not rely on cachedS: for use with pseudo keywords
distance.matrix <- function(keywords) {
    n <- length(keywords)
    distances <- matrix(0, nrow=n, ncol=n)
    for(i1 in 1:n) {
        for(i2 in i1:n) {
            if(i1 != i2) {
                s <- sum(vapply(1:rn, function(r) {
                        S.metric(r, keywords[[i1]], keywords[[i2]]) / largestS[r]
                    }, 0))
                s = 1 - s / rn
                distances[i1,i2] = s
                distances[i2,i1] = s
            }
        }
    }
    diag(distances) = Inf
    distances
}

# merge i and j in clusters
# auxiliary method for hand-written clustering
merge.cluster <- function(clusters, i, j) {
    ncl <- c(clusters[[i]], clusters[[j]])
    clusters[[i]] = ncl
    clusters = clusters[-j]
    clusters
}

# mergeSimilarKeywords
similar <- function() {
    # similarityLink - candidates for relatedEquivalent
    links <- get.semantic(4)
    keywords <- unique(as.vector(links))
    if(verbosity>=2) cat("mergeSimilarKeywords for", nrow(links), "links or", length(keywords), "keywords.\n")
    if(!length(keywords)) return()

    cluster_v <- cutree(hclust(as.dist(distance.matrix.cached(keywords)), method="single"), h=merge_t)
    nclusters <- max(cluster_v)
    for(k in 1:nclusters) {
        cl = keywords[which(cluster_v == k)]
        pairs <- combn(cl, 2)
        # relatedEquivalent relation is a symmetric one
        set.semantic(pairs[1, ], semantic[1], pairs[2, ])
        merge.keywords(cl)
    }
    if(verbosity>=2) cat("Merging resulted in", nclusters, "clusters.\n")

    # reset similarityLink
    semrel[[4]][1:nrow(links), ] <<- 0
    semrel$sizes[4] <<- 0
    cleanup.semrel()

    continue <<- TRUE
}

harm.mean <- function(x) 1 / mean(1/x)

# chooses higher-level keyword from cluster with respect to k
# that will be used to name pseudo-keywords
high.in.cluster <- function(k, other) {
    hmeans <- mapply(function(x,y) harm.mean(c(x,y)),
        vapply(other, cooccur, 0, k),
        vapply(other, npapers, 0))
    other[which.max(hmeans)]
}

# returns list of pseudo-keyword objects intersected with k
# number of pseudo-keywords is equal to number of clusters
gen.pseudos <- function(k, clusters) {
    pseudos <- list()
    for(i in seq_along(clusters)) {
        pseudos[[i]] = create.pseudo(k, clusters[[i]])
    }
    pseudos
}

quick.clustering <- function(keywords) {
    if(length(keywords) <= 1) return(list())
    distances <- distance.matrix.cached(keywords)
    weights <- vapply(keywords, npapers, 0)
    clusters <- as.list(1:length(keywords))
    d <- distances
    while(TRUE) {
        i <- which(d==min(d, na.rm=TRUE), arr.ind=T)
        if(length(i) > 2) i = i[1,]
        if(length(clusters) > 1 && d[i[1],i[2]] <= quick_t) {
            clusters = merge.cluster(clusters, i[1], i[2])
            d = update_dist_C(distances, clusters, weights)
        } else break
    }
    for(i in seq_along(clusters)) {
        clusters[[i]] = keywords[clusters[[i]]]
    }
    clusters
}

# ambk - potentially ambiguous keyword
# keywords - set of related to k keywords to clusterize
# returns whether split happened or not
intersect.clustering <- function(ambk, keywords) {
    k = keyword.index(ambk)
    clusters <- as.list(keywords)
    pseudos <- gen.pseudos(k, clusters)
    d <- distance.matrix(pseudos)
    while(TRUE) {
        i <- which(d==min(d, na.rm=TRUE), arr.ind=T)
        if(length(i) > 2) i = i[1,]
        if(length(clusters) > 1 && d[i[1],i[2]] < intersect_t) {
            clusters = merge.cluster(clusters, i[1], i[2])
            # pseudos = gen.pseudos(k, clusters)
            pseudos = update.pseudos(pseudos, i[1], i[2])
            d = distance.matrix(pseudos)
        } else break
    }
    if(length(clusters) > 1) {
        if(verbosity>=3) cat("splitting", ambk, "into ", length(clusters), "keywords\n")
        # add pseudos to global vars
        add.pseudos(k, pseudos, clusters)
        # delete ambiguous keyword
        delete.keyword(k)
        # split was done
        continue <<- TRUE
        return(TRUE)
    }
    FALSE
}

# splitAmbiguousKeywords
ambiguous <- function() {
    if(verbosity>=2) cat("Seeking ambiguous keywords.\n")
    totalsplit <- 0
    for(k in all.keywords()) {
        if(verbosity>=3) cat("splitAmbiguousKeywords for", k, "\n")
        rk <- related.keywords(k, threshold=relkeyAmbig)
        clusters <- quick.clustering(rk)
        if(length(clusters) > 1) {
            if(verbosity>=3) cat("intersect clustering: #clusters =", length(clusters), "\n")
            s <- intersect.clustering(k, rk)
            if(s) totalsplit = totalsplit + 1
        }
    }
    if(verbosity>=2) cat("Splitting (intersect clustering) was done", totalsplit, "times.\n")
}

# filterNotAcademicKeywords
academic <- function() {
    if(verbosity>=2) cat("Filtering keywords.\n")
    # 1: keywords without relations
    # Semantic relations are stored in a separate data structure,
    # so no need in this check.

    # 2: distribution check
    keywords <- unique(triples$k1, triples$k2)
    for(k in keywords) {
        mo <- main.cooccur(k, nmain, index=FALSE)
        p <- apply(mo, 2, sum) / total.cooccur(k)
        p[is.nan(p)] = 1
        # or all?
        if(any(p < maincover)) {
            # delete keyword from output semantic relations
            triples <<- triples[!triples$k1==k | triples$k2==k,]
        }
    }

    # 3: only one source at the moment, so no need
}

# empties global output variables and loads input
prepare.globals <- function(inputfile) {
    load(inputfile)
    reldb_df <<- reldb_df
    reldb_l <<- reldb_l
    keywordsdb <<- keywordsdb
    inputm <<- inputm

    triples <<- data.frame(k1=character(0), k2=character(0), relation=numeric(0), stringsAsFactors=FALSE)
    semrel <<- list()
    prepare.semrel()
    continue <<- TRUE

    # defined in relations.R
    maxindex <<- NULL
}

# must be called once per iteration
update.caches <- function() {
    cachedS <<- list()
    cachedCV <<- list()
    cached.names <<- NULL
}

klink2 <- function(inputfile) {
    # ensure correctness of global variables
    prepare.globals(inputfile)

    split_merge <- TRUE
    iter <- 1
    while(continue) {
        update.caches()
        if(verbosity>=1) cat("Iteration", iter, "\nNumber of keywords =", nkeywords(), "\n")
        if(verbosity>=2) cat("Keyword inference.\n")
        # set to true only if there was splitting / merging done
        continue <<- FALSE
        for(k in all.keywords()) {
            if(verbosity>=3) cat("Infering keyword:", k, "\n")
            rk <- related.keywords(k)
            for(k2 in rk) {
                infer(k, keyword.name(k2))
            }
        }
        fix.loops()
        if(verbosity>=1)
            cat("Number of working semantic relations after inference:\n\trelatedEquivalent: ", semrel$sizes[1],
                    "\n\tbroaderGeneric: ", semrel$sizes[2],
                    "\n\tcontributesTo:", semrel$sizes[3],
                    "\n\tsimilarityLink:", semrel$sizes[4], "\n")
        if(split_merge)
            ambiguous()
        else
            similar()
        split_merge = !split_merge
        iter = iter + 1
    }
    academic()

    rm('reldb_df', 'reldb_l', 'keywordsdb', 'inputm', envir=globalenv())
    output.stats()
}

output.keywords <- function() {
    unique(c(triples$k1, triples$k2))
}

output.stats <- function() {
    cat("Output\n\tnumber of keywords: ", length(output.keywords()),
        "\n\tnumber of relations: ", dim(triples)[1], "\n")
}
