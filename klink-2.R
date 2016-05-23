source('param.R')
source("relations.R")
source("graph.R")
library(stringi)

# Var naming here:
# r - relation, string
# k, x, y - keywords, string

# I_r(x,y) conditional probability that
# an element associated with x will be associated with y
I.prob <- function(r, x, y, diachronic=FALSE) {
    v <- rel.value(x, y, r)
    if(v==0) return(0)
    w <- ifelse(diachronic,
            (rel.year(y, r) - min(rel.year(x, r)) + 1)^(-gamma),
            1)
    w * v
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

# number of common characters
common.chars <- function(x, y) {
    charsx <- trimws(strsplit(x, "")[[1]])
    charsx = unique(charsx[which(nchar(charsx) > 0)])
    charsy <- trimws(strsplit(y, "")[[1]])
    charsy = unique(charsy[which(nchar(charsy) > 0)])
    length(intersect(charsx, charsy))
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

# c_r(x,y) measure
# vx, vy as arguments for optimization purpose
semantic.similarity <- function(r, x, y, is_super=FALSE, is_sib=FALSE, vx=NULL, vy=NULL) {
    if(is.null(vx)) vx <- conn.vector(r, x, is_super, is_sib)
    if(is.null(vy)) vy <- conn.vector(r, y, is_super, is_sib)
    s <- as.double(vx %*% vy) / (sqrt(sum(vx^2)) * sqrt(sum(vy^2)))
    ifelse(is.nan(s), 0, s)
}

# n(x,y) measure
string.similarity <- function(x, y) {
    sum(nweights * c(length(lcs(x,y)), identical.words(x,y), common.chars(x,y), have.acronym(x,y)))
}

H.metric <- function(r, x, y, diachronic=FALSE) {
    m <- (I.prob(r, x, y, diachronic) / I.prob(r, x, x, diachronic) - I.prob(r, y, x, diachronic) / I.prob(r, y, y, diachronic)) *
        semantic.similarity(r, x, y) * string.similarity(x, y)
    ifelse(is.nan(m), 0, m)
}

T.metric <- function(r, x, y) {
    H.metric(r, x, y, diachronic=TRUE)
}

S.metric <- function(r, x, y) {
    vx <- conn.vector(r, x, is_super=TRUE, is_sib=TRUE)
    vy <- conn.vector(r, y, is_super=TRUE, is_sib=TRUE)
    semantic.similarity(r, x, y, vx=vx[,1], vy=vy[,1]) /
        (max(semantic.similarity(r, x, y, vx=vx[,2], vy=vy[,2]), semantic.similarity(r, x, y, vx=vx[,3], vy=vy[,3])) + 1)
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
    for(i in seq_along(relations)) {
        r <- relations[i]
        h <- H.metric(r, x, y)
        hmetrics = c(hmetrics, h)
        if (h >= tR[i]) hierarchy = c(hierarchy, 1)
        else if (h <= -tR[i]) hierarchy = c(hierarchy, -1)
        else hierarchy = c(hierarchy, 0)
        tmetrics = c(tmetrics, T.metric(r, x, y))
        smetrics = c(smetrics, S.metric(r, x, y))
    }
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
            ifelse(hierarchy == 1,
                set.semantic(x, semantic[2], y),
                set.semantic(y, semantic[2], x))
        } else {
            # contributesTo
            ifelse(hierarchy == 1,
                set.semantic(y, semantic[3], x),
                set.semantic(x, semantic[3], y))
        }
    }
    # "similarityLink"
    if(length(which(smetrics > tS)) >= tre) {
        set.semantic(x, semantic[4], y)
    }
}

# break loops for broaderGeneric
fix.loops <- function() {
    if(nrow(semrel[[2]]) > 1)
        semrel[[2]] <<- delete.cycles(semrel[[2]])
}

# keywords: list of keyword objects or vector of keyword ids
distance.matrix <- function(keywords) {
    n <- length(keywords)
    distances <- matrix(0, nrow=n, ncol=n)
    for(i1 in 1:n) {
        for(i2 in i1:n) {
            if(i1 != i2) {
                distances[i1,i2] = sum(sapply(relations, function(r) {
                    ifelse(is.list(keywords),
                                S.metric(r, keywords[[i1]], keywords[[i2]]),
                                S.metric(r, keywords[i1], keywords[i2]))
                    }))
                distances[i2,i1] = distances[i1,i2]
            }
        }
    }
    diag(distances) = Inf
    distances
}

# merge i and j in clusters
merge.cluster <- function(clusters, i, j) {
    ncl <- c(clusters[[i]], clusters[[j]])
    clusters[[i]] = ncl
    clusters = clusters[-j]
    clusters
}

# mergeSimilarKeywords
similar <- function() {
    links <- semrel[[4]]
    keywords <- unique(as.vector(links))
    if(verbosity>=2) cat("mergeSimilarKeywords for", nrow(links), "links or", length(keywords), "keywords.\n")
    if(!length(keywords)) return()
    distances <- distance.matrix(keywords)
    clusters <- as.list(keywords)
    update.dist <- function() {
        d <- matrix(0, nrow=length(clusters), ncol=length(clusters))
        for(i in seq_along(clusters)) {
            for(j in seq_along(clusters)) {
                if(i != j) {
                    p1 <- which(keywords %in% clusters[[i]])
                    p2 <- which(keywords %in% clusters[[j]])
                    d[i,j] = sum(distances[as.matrix(expand.grid(p1,p2))])
                }
            }
        }
        diag(d) = Inf
        d
    }
    d <- distances
    while(TRUE) {
        i <- which(d==min(d, na.rm=TRUE), arr.ind=T)
        if(length(i) > 2) i = i[1,]
        if(length(clusters) > 1 && d[i[1],i[2]] < merge_t) {
            clusters = merge.cluster(clusters, i[1], i[2])
            d = update.dist()
        } else break
    }
    for(cl in clusters) {
        for(i in cl)
            for(j in cl)
                if(i!=j) set.semantic(i, semantic[1], j)
        merge.keywords(cl)
    }
    semrel[[4]] <<- matrix(, nrow=0, ncol=2)
    continue <<- TRUE
}

harm.mean <- function(x) 1 / mean(1/x)

# chooses higher-level keyword from cluster with respect to k
# that will be used to name pseudo-keywords
high.in.cluster <- function(k, other) {
    hmeans <- mapply(function(x,y) harm.mean(c(x,y)),
        sapply(other, cooccur, k),
        sapply(other, npapers))
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
    distances <- distance.matrix(keywords)
    clusters <- as.list(keywords)
    update.dist <- function(keywords) {
        d <- matrix(0, nrow=length(clusters), ncol=length(clusters))
        for(i in seq_along(clusters)) {
            for(j in seq_along(clusters)) {
                if(i != j) {
                    p1 <- which(keywords %in% clusters[[i]])
                    p2 <- which(keywords %in% clusters[[j]])
                    w <- rep(sapply(p1, npapers), length(p2))
                    d[i,j] = sum((w * distances[as.matrix(expand.grid(p1,p2))]) / sum(w))
                }
            }
        }
        diag(d) = Inf
        d
    }
    d <- distances
    while(TRUE) {
        i <- which(d==min(d, na.rm=TRUE), arr.ind=T)
        if(length(i) > 2) i = i[1,]
        if(length(clusters) > 1 && d[i[1],i[2]] < quick_t) {
            clusters = merge.cluster(clusters, i[1], i[2])
            d = update.dist(keywords)
        } else break
    }
    clusters
}

# ambk - potentially ambiguous keyword
# keywords - set of related to k keywords to clusterize
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
            pseudos = gen.pseudos(k, clusters)
            d = distance.matrix(pseudos)
        } else break
    }
    if(length(clusters) > 1) {
        if(verbosity>=3) cat("splitting", ambk, "into ", length(clusters), "keywords\n")
        # add pseudo-keywords
        for(i in seq_along(clusters)) {
            add.pseudo(pseudos[[i]],
                paste(keyword.name(k), " (", keyword.name(high.in.cluster(k, clusters[[i]])), ")", sep=""),
                clusters[[i]])
        }
        # delete ambiguous keyword
        delete.keyword(k)
        # split was done
        continue <<- TRUE
    }
}

# splitAmbiguousKeywords
ambiguous <- function() {
    if(verbosity>=2) cat("Seeking ambiguous keywords.\n")
    for(k in all.keywords()) {
        if(verbosity>=3) cat("splitAmbiguousKeywords for", k, "\n")
        rk <- related.keywords(k, threshold=relkeyT * 2)
        clusters <- quick.clustering(rk)
        if(length(clusters) > 1) {
            if(verbosity>=3) cat("intersect clustering: #clusters =", length(clusters), "\n")
            intersect.clustering(k, rk)
        }
    }
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

# output semantic relations
triples <- data.frame(k1=character(0), k2=character(0), relation=numeric(0), stringsAsFactors=FALSE)
# working semantic relations
semrel <- list()
# iteration regulator
continue <- TRUE

prepare.semrel <- function() {
    for(i in seq_along(semantic)) {
        semrel[[i]] <<- matrix(, nrow=0, ncol=2)
    }
    names(semrel) <<- semantic
}

klink2 <- function() {
    prepare.semrel()

    split_merge <- TRUE
    iter <- 1
    while(continue) {
        if(verbosity>=1) cat("Iteration", iter, "\nNumber of keywords =", nkeywords(), "\n")
        if(verbosity>=2) cat("Keyword inferrence.\n")
        # set to true only if there was splitting / merging done
        continue <<- FALSE
        for(k in all.keywords()) {
            if(verbosity>=3) cat("Inferring keyword:", k, "\n")
            rk <- related.keywords(k)
            for(k2 in rk) {
                infer(k, keyword.name(k2))
            }
        }
        fix.loops()
        if(verbosity>=1)
            cat("Number of working semantic relations after inference:\n\trelatedEquivalent: ", nrow(semrel[[1]]),
                "\n\tbroaderGeneric: ", nrow(semrel[[2]]), "\n\tcontributesTo:", nrow(semrel[[3]]),
                "\n\tsimilarityLink:", nrow(semrel[[4]]), "\n")
        if(split_merge)
            ambiguous()
        else
            similar()
        split_merge = !split_merge
        iter = iter + 1
    }
    academic()
}
