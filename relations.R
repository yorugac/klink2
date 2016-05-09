# Handling keywords and their relations.

# It is assumed that there are global objects keywordsdb, reldb_df and inputm
# that hold all the keywords, relations and co-occurrences calculations;
# as well as semrel for output

# returns dataframe with all relations for the given keyword
get.reldf <- function(keyword) {
    index <- keyword.index(keyword)
    reldb_df[[index]]
}

keyword.index <- function(keyword) {
    ifelse(is.character(keyword),
        keywordsdb[[keyword]],
        as.integer(keyword))
}

keyword.name <- function(index) {
    names(reldb_df)[[index]]
}

relation.index <- function(rel) {
    if(is.character(rel)) rel = which(relations == rel)
    as.integer(rel)
}

nkeywords <- function() {
    length(ls(keywordsdb))
}

# strength of connection between k1 and k2 in regards to relation rel
# returns number of co-occurrences in case of unquantified relation and
# sum of minimum values in case of quantified
# rel can be an index or a string in relations
rel.value <- function(k1, k2, r) {
    if(k1==k2) return(length(rel.entity(k1, r)))
    rel = relation.index(r)
    ik1 <- keyword.index(k1)
    ik2 <- keyword.index(k2)
    i <- match(ik2, inputm[, 2*rel-1 + 2*rn*(ik1-1)])
    if(is.na(i)) {
         0
    } else {
        ifelse(quantified[rel], {
            xdf <- get.reldf(k1)
            ydf <- get.reldf(k1)
            xv <- rel.entity(k1, r, xdf)
            xq <- rel.quantity(k1, r, xdf)
            yv <- rel.entity(k2, r, ydf)
            yq <- rel.quantity(k2, r, ydf)
            cooccur <- intersect(xv, yv)
            sum(pmin(xq[which(xv==cooccur)], yq[which(yv==cooccur)]))
        }, inputm[i, 2*rel + 2*rn*(ik1-1)])
    }
}

set.semantic <- function(k1, semtype, k2) {
    if(!semtype %in% semantic)
        warning(paste("used unknown semantic relation ", semtype, sep=""))
    ik1 <- keyword.index(k1)
    ik2 <- keyword.index(k2)
    semrel[[semtype]] <<- unique(rbind(semrel[[semtype]], c(ik1, ik2)))
}

# returns keywords that are super according to hierarchical relations
super <- function(keyword) {
    ik <- keyword.index(keyword)
    ibr <- which(ik == semrel[[2]][,1])
    ict <- which(ik == semrel[[3]][,2])
    res <- c()
    if(length(ibr)) res = c(res, semrel[[2]][ibr, 2])
    if(length(ict)) res = c(res, semrel[[3]][ict, 1])
    res
}

# returns keywords that are siblings according to similarity relation
sib <- function(keyword) {
    ik <- keyword.index(keyword)
    ire1 <- which(ik == semrel[[1]][,1])
    ire2 <- which(ik == semrel[[1]][,2])
    res <- c()
    if(length(ire1)) res = c(res, semrel[[1]][ire1, 2])
    if(length(ire2)) res = c(res, semrel[[1]][ire2, 1])
    res
}

# total co-occurrence between k1 and k2
# without inputm
cooccur <- function(k1, k2) {
    ik1 <- keyword.index(k1)
    ik2 <- keyword.index(k2)
    rel1 <- reldb_df[[k1]]
    rel2 <- reldb_df[[k2]]
    v <- 0
    for(r in 1:rn) {
        rel = relations[r]
        cooccur = rel.entity(k1, rel, rel1) %in% rel.entity(k2, rel, rel2)
        v = v + length(which(cooccur))
    }
    v
}

# values of connections for keyword and rel;
# with length equal to number of keywords
# With super / sib only super or sib keywords are taken into account.
conn.vector <- function(rel, keyword, is_super=FALSE, is_sib=FALSE) {
    ik <- keyword.index(keyword)
    rel = relation.index(rel)
    v <- inputm[, 2*rel + 2*rn*(ik-1)]
    iv <- inputm[, 2*rel-1 + 2*rn*(ik-1)]
    ind <- which(v > 0)
    v = v[ind]
    iv = iv[ind]
    if(is_super) {
        leave = iv %in% super(ik)
        iv = iv[leave]
        v = v[leave]
    }
    if(is_sib) {
        leave = iv %in% sib(ik)
        iv = iv[leave]
        v = v[leave]
    }
    bv = rep(0, nkeywords())
    bv[iv] = v
    bv
}

total.cooccur <- function(k) {
    ik <- keyword.index(k)
    rel <- 1:rn
    v <- inputm[, 2*rel + 2*rn*(ik-1)]
    apply(v, 2, sum)
}

# returns subset of inputm matrix for keyword k and n first keywords
main.cooccur <- function(k, n, index=FALSE) {
    ik <- keyword.index(k)
    if(index)
        inputm[1:n, sort(c(2*1:rn-1 + 2*rn*(ik-1), 2*1:rn + 2*rn*(ik-1)))]
     else
        inputm[1:n, 2*1:rn + 2*rn*(ik-1)]
}

# number of entities associated with keyword
nentities <- function(keyword) {
    length(unique(get.reldf(keyword)$entity))
}

# number of papers associated with keyword
npapers <- function(keyword) {
    reldf <- get.reldf(keyword)
    length(reldf[reldf$relation == "publication",]$entity)
}

# returns entity(s) for the given keyword and relation name
rel.entity <- function(keyword, relation, reldf=NULL) {
    if(is.null(reldf)) reldf = get.reldf(keyword)
    reldf[reldf$relation %in% relation,]$entity
}

# returns quantity(s) for the given keyword and relation name
rel.quantity <- function(keyword, relation, reldf=NULL) {
    if(is.null(reldf)) reldf = get.reldf(keyword)
    reldf[reldf$relation %in% relation,]$quantity
    # get rid of NAs?
}

# returns year(s) for the given keyword and relation name
rel.year <- function(keyword, relation, reldf=NULL) {
    if(is.null(reldf)) reldf = get.reldf(keyword)
    as.numeric(reldf[reldf$relation %in% relation,]$year)
}

# when keyword was first used
age <- function(keyword) {
    min(get.reldf(keyword)$year, na.rm=TRUE)
}

# returns vector of related keywords to the given keyword
# threshold - min number of entities to be shared by related keywords
# (NOTE: among ones that are saved in input matrix)
related.keywords <- function(keyword, threshold=relkeyT) {
    ik <- keyword.index(keyword)
    v <- inputm[, 2*1:rn + 2*rn*(ik-1)]
    iv <- inputm[, 2*1:rn-1 + 2*rn*(ik-1)]
    unique(iv[which(v > threshold)])
}

merge.keywords <- function(cluster) {
    newk <- paste(keyword.name(cluster[1]), " merged", sep="")
    index <- nkeywords() + 1
    rel <- get.reldf(cluster[1])
    for(k in 2:length(cluster)) {
        rel = rbind(rel, get.reldf(cluster[k]))
    }
    m <- nrow(inputm)
    for(r in 1:rn) {
        v <- as.vector(inputm[, 2*r + 2*rn*(cluster-1)])
        iv <- as.vector(inputm[, 2*r-1 + 2*rn*(cluster-1)])
        ind <- order(v)
        v = v[ind][1:m]
        iv = iv[ind][1:m]
        if(ncol(inputm) < 2*r-1 + 2*rn*(index-1))
            inputm <<- cbind(inputm, matrix(0, nrow=m, ncol=2*rn))
        inputm[, 2*r-1 + 2*rn*(index-1)] <<- iv
        inputm[, 2*r + 2*rn*(index-1)] <<- v
    }
    reldb_df[[index]] <<- rel
    names(reldb_df)[index] <<- newk
    keywordsdb[[newk]] <<- index
    for(k in cluster) delete.keyword(k)
}

add.pseudo <- function(keyword, pseudo, cluster) {
    newk <- paste(keyword.name(keyword), " (", keyword.name(pseudo), ")", sep="")
    index <- nkeywords() + 1
    relk <- get.reldf(keyword)

    m <- nrow(inputm)
    for(r in 1:rn) {
        v <- as.vector(inputm[, 2*r + 2*rn*(cluster-1)])
        iv <- as.vector(inputm[, 2*r-1 + 2*rn*(cluster-1)])
        ind <- order(v)
        v = v[ind][1:m]
        iv = iv[ind][1:m]
        if(ncol(inputm) < 2*r-1 + 2*rn*(index-1))
            inputm <<- cbind(inputm, matrix(0, nrow=m, ncol=2*rn))
        inputm[, 2*r-1 + 2*rn*(index-1)] <<- iv
        inputm[, 2*r + 2*rn*(index-1)] <<- v
    }
    entities <- c()
    for(k in cluster) {
        rel <- get.reldf(k)
        entities = unique(c(entities, intersect(relk$entity, rel$entity)))
    }
    reldb_df[[index]] <<- relk[relk$entity %in% entities,]
    names(reldb_df)[index] <<- newk
    keywordsdb[[newk]] <<- index
    index
}

# keyword must be character
delete.keyword <- function(keyword) {
    if(is.character(keyword)) {
        k = keyword
        ik <- keyword.index(keyword)
    } else {
        k = keyword.name(keyword)
        ik <- keyword.index(keyword)
    }
    rm(list=k, envir=keywordsdb)
    for(i in seq_along(semrel)) {
        if(nrow(semrel[[i]]) > 0) for(r in 1:nrow(semrel[[i]])) {
            if(any(semrel[[i]][r,]==ik)) semrel[[i]] <<- semrel[[i]][-r,]
        }
    }
    # delete cols in inputm?
}
