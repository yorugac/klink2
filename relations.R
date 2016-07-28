# Handling keywords and their relations; internal functions behind Klink-2.

# It is assumed that there are global objects keywordsdb, reldb_df, reldb_l and inputm
# that hold all the keywords, relations and co-occurrences calculations;
# as well as semrel and triples for output.
#
# No keywords are truly deteled. Whether keyword should be included in analysis
# can be checked with keyword.exists -- whether there is a name for it in reldb_df.

maxindex <- NULL
zeromatrix <- NULL # used for semrel structure growth

semantic.index <- function(semtype) {
    if(is.character(semtype)) semtype = which(semantic == semtype)
    as.integer(semtype)
}

prepare.semrel <- function() {
    zeromatrix <<- matrix(0, nrow=nkeywords()*2, ncol=2)
    for(i in seq_along(semantic)) {
        semrel[[i]] <<- zeromatrix
    }
    names(semrel) <<- semantic
    semrel$sizes <<- rep(0, length(semantic))
}

cleanup.semrel <- function() {
    for(i in seq_along(semantic)) {
        t <- unique(get.semantic(i))
        semrel[[i]] <<- rbind(t, matrix(0, nrow=nkeywords()*2, ncol=2))
        semrel$sizes[i] <<- nrow(t)
    }
}

# returns dataframe with all relations for the given keyword
get.reldf <- function(keyword) {
    index <- keyword.index(keyword)
    reldb_df[[index]]
}

keyword.index <- function(keyword) {
    if(is.character(keyword))
        keywordsdb[[keyword]]
    else as.integer(keyword)
}

keyword.name <- function(index) {
    names(reldb_df)[[index]]
}

keyword.exists <- function(index) {
    !is.na(keyword.name(index))
}

relation.index <- function(rel) {
    if(is.character(rel)) rel = which(relations == rel)
    as.integer(rel)
}

nkeywords <- function() {
    length(all.keywords())
}

all.keywords <- function() {
    names(reldb_df)[!is.na(names(reldb_df))]
}

new.index <- function() {
    if(is.null(maxindex)) maxindex <<- nkeywords()
    maxindex <<- maxindex + 1
    maxindex
}

cached.keys <- function(ik, r) 2*r-1 + 2*rn*(ik-1)
cached.values <- function(ik, r) 2*r + 2*rn*(ik-1)

# representation of keyword as a single object, without relation dataframe
keyword.object <- function(k) {
    ik <- keyword.index(k)
    list(irel=inputm[, cached.keys(ik, 1)],
        rel=inputm[, cached.values(ik, 1)],
        super=super(ik),
        sib=sib(ik),
        df=NULL,
        rellist=NULL)
}

# calculation of m largest co-occurrences for keyword i and relation rname;
# compares i against keywords in words.
# Returns matrix with column of indexes and column of values.
# rname - relation name
# i - keyword index
# irel_l - object of reldb_l for i
# irel_df - object reldb_df for i
calc.cooccurrence <- function(rname, i, irel_l, irel_df=NULL, words=all.keywords()) {
    if(is.null(irel_df)) irel_df <- get.reldf(i)
    m <- nrow(inputm)
    res <- matrix(0, nrow=m, ncol=2)
    for(jk in words) {
        j <- keyword.index(jk)
        if(i != j) {
            jrel_l <- reldb_l[[j]]
            cooccur <- intersect_C(irel_l[[rname]], jrel_l[[rname]])
            value <- cooccurrence_C(irel_df, get.reldf(j), cooccur)
            # is there j word in the stored list?
            im <- match_C(j, res[,1])
            if(!is.na(im)) {
                res[im,2] = res[im,2] + value
            } else {
                # is there place in the stored list for the value with j?
                im = which.min(res[,2])
                if(value > res[im,2]) {
                     res[im,1] = j
                     res[im,2] = value
                }
            }
        }
    }
    # sorting in decreasing order of relation values
    ind <- order(res[,2], decreasing=TRUE)
    res[,1] = res[ind,1]
    res[,2] = res[ind,2]
    res
}

# retrieval of m largest co-occurrences for given keywords and relation rname
# Returns matrix with column of indexes and column of values.
# rname - relation name
# words - keywords
calc.cooccurrence.cached <- function(rname, words) {
    m <- nrow(inputm)
    r <- relation.index(rname)
    iks <- vapply(words, keyword.index, 0)
    iv <- as.vector(inputm[, cached.keys(iks, r)])
    v <- as.vector(inputm[, cached.values(iks, r)])

    ind <- order(iv)
    ordiv <- iv[ind]
    ordv <- v[ind]
    outiv <- unique(ordiv)
    outv <- unlist(lapply(split(ordv, ordiv), sum))

    # sorting in decreasing order of relation values
    ind = order(outv, decreasing=TRUE)
    # forming final result
    res <- matrix(0, nrow=m, ncol=2)
    res[, 1] = outiv[ind][1:m]
    res[, 2] = outv[ind][1:m]
    if(length(outiv) < m) res[(length(outiv)+1):m, ] = 0
    res
}

# strength of connection between k1 and k2 in regards to relation rel
# returns number of co-occurrences in case of unquantified relation and
# vector of minimum values in case of quantified
# rel can be an index or a string in relations
rel.value <- function(k1, k2, r) {
    if(k1==k2) return(length(rel.entity(k1, r)))
    rel = relation.index(r)
    ik1 <- keyword.index(k1)
    ik2 <- keyword.index(k2)
    i <- match(ik2, inputm[, cached.keys(ik1, rel)])
    if(is.na(i)) {
         0
    } else {
        if(quantified[rel]) {
            xdf <- get.reldf(k1)
            ydf <- get.reldf(k1)
            xv <- rel.entity(k1, r, xdf)
            xq <- rel.quantity(k1, r, xdf)
            yv <- rel.entity(k2, r, ydf)
            yq <- rel.quantity(k2, r, ydf)
            cooccur <- intersect(xv, yv)
            pmin(xq[which(xv==cooccur)], yq[which(yv==cooccur)])
        } else  {
            inputm[i, cached.values(ik1, rel)]
        }
    }
}

# returns matrix of all working relations for semtype
get.semantic <- function(semtype) {
    i <- semantic.index(semtype)
    size <- semrel$sizes[i]
    if(size == 0) matrix(, nrow=0, ncol=2)
    else matrix(semrel[[i]][1:size, ], ncol=2)
    # additional matrix call to ensure that size=1 is returned as matrix, not vector
}

# k1 and k2 can be vectors
set.semantic <- function(k1, semtype, k2) {
    if(length(k1) != length(k2)) stop("set.semantic got vectors of different size")

    i <- semantic.index(semtype)
    ik1 <- k1; ik2 <- k2
    if(!is.numeric(k1)) ik1 = vapply(k1, keyword.index, 0)
    if(!is.numeric(k2)) ik2 = vapply(k2, keyword.index, 0)
    prevsize <- semrel$sizes[i]
    s <- prevsize + length(ik1)
    # check if there is place to put new relation
    if(s >= nrow(semrel[[i]])) {
        semrel[[i]] <<- rbind(semrel[[i]][1:prevsize, ], matrix(c(ik1, ik2), ncol=2), zeromatrix)
    } else {
        semrel[[i]][(prevsize+1):s, 1] <<- ik1
        semrel[[i]][(prevsize+1):s, 2] <<- ik2
    }
    semrel$sizes[i] <<- s
}

# returns keywords that are super according to hierarchical relations
super <- function(keyword) {
    ik <- keyword.index(keyword)
    semmatrix2 <- get.semantic(2)
    semmatrix3 <- get.semantic(3)
    ibr <- which(ik == semmatrix2[,1])
    ict <- which(ik == semmatrix3[,2])
    res <- c()
    if(length(ibr)) res = c(res, semmatrix2[ibr, 2])
    if(length(ict)) res = c(res, semmatrix3[ict, 1])
    res
}

# returns keywords that are siblings according to similarity relation
sib <- function(keyword) {
    ik <- keyword.index(keyword)
    semmatrix <- get.semantic(1)
    ire1 <- which(ik == semmatrix[,1])
    ire2 <- which(ik == semmatrix[,2])
    res <- c()
    if(length(ire1)) res = c(res, semmatrix[ire1, 2])
    if(length(ire2)) res = c(res, semmatrix[ire2, 1])
    res
}

# rel - string
# returns vector of entities shared by k1 and k2 via relation rel
common.entities <- function(k1, k2, rel) {
    rel1 <- reldb_df[[k1]]
    rel2 <- reldb_df[[k2]]
    intersect_C(rel.entity(k1, rel, rel1), rel.entity(k2, rel, rel2))
}

# total co-occurrence between k1 and k2
# without inputm
cooccur <- function(k1, k2) {
    v <- 0
    for(r in 1:rn) {
        rel = relations[r]
        cooccur = common.entities(k1, k2, rel)
        v = v + length(common.entities(k1, k2, rel))
    }
    v
}

# values of connections for keyword and rel in form of matrix;
# with length equal to number of keywords
# keyword: id, string or keyword object
# With is_super / is_sib only super or sib keywords are taken into account.
# If both super and sib flags are set, list with 3 matrices is returned.
conn.vector <- function(rel, keyword, is_super=FALSE, is_sib=FALSE) {
    rel <- relation.index(rel)
    if(!is.list(keyword)) {
        ik <- keyword.index(keyword)
        v <- inputm[, cached.values(ik, rel)]
        iv <- inputm[, cached.keys(ik, rel)]
        if(is_super) superv <- super(ik)
        if(is_sib) sibv <- sib(ik)
    } else {
        v <- keyword$rel[, rel]
        iv <- keyword$irel[, rel]
        superv <- keyword$super
        sibv <- keyword$sib
    }
    ind <- which(v > 0)
    # if keywords were already deleted
    ind = which(iv[ind] <= nkeywords())
    v = v[ind]
    iv = iv[ind]
    if(is_super && is_sib) {
        bv = list()
        bv[[1]] = matrix(c(iv, v), ncol=2)
        leave = iv %in% superv
        bv[[2]] = matrix(c(iv[leave], v[leave]), ncol=2)
        leave = iv %in% sibv
        bv[[3]] = matrix(c(iv[leave], v[leave]), ncol=2)
        return(bv)
    }
    if(is_super) {
        leave = iv %in% superv
        iv = iv[leave]
        v = v[leave]
    } else if(is_sib) {
        leave = iv %in% sibv
        iv = iv[leave]
        v = v[leave]
    }
    matrix(c(iv, v), ncol=2)
}

total.cooccur <- function(k) {
    ik <- keyword.index(k)
    rel <- 1:rn
    v <- inputm[, cached.values(ik, rel)]
    apply(v, 2, sum)
}

# returns subset of inputm matrix for keyword k and n first keywords
main.cooccur <- function(k, n, index=FALSE) {
    ik <- keyword.index(k)
    if(index)
        inputm[1:n, sort(c(cached.keys(ik, 1:rn), cached.values(ik, 1:rn)))]
     else
        inputm[1:n, cached.values(ik, 1:rn)]
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
    reldf[reldf$relation == relation,]$entity
}

# returns quantity(s) for the given keyword and relation name
rel.quantity <- function(keyword, relation, reldf=NULL) {
    if(is.null(reldf)) reldf = get.reldf(keyword)
    reldf[reldf$relation %in% relation,]$quantity
    # get rid of NAs?
}

# returns year(s) for the given keyword and entities
ent.year <- function(keyword, entities, reldf=NULL) {
    if(is.null(reldf)) reldf = get.reldf(keyword)
    sort(reldf[reldf$entity %in% entities,]$year)
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
    v <- inputm[, cached.values(ik, 1:rn)]
    iv <- inputm[, cached.keys(ik, 1:rn)]
    rk <- unique(iv[which(v > threshold)])
    # some words might have been already deleted
    validrk <- vapply(rk, keyword.exists, TRUE)
    if(length(rk) == 0) c()
    else rk[validrk]
}

merge.keywords <- function(cluster) {
    newk <- paste(keyword.name(cluster[1]), " merged", sep="")
    index <- new.index()
    rel <- get.reldf(cluster[1])
    rel_l <- list()
    for(k in 2:length(cluster)) {
        rel = rbind(rel, get.reldf(cluster[k]))
    }
    m <- nrow(inputm)
    for(r in 1:rn) {
        v <- as.vector(inputm[, cached.values(cluster, r)])
        iv <- as.vector(inputm[, cached.keys(cluster, r)])
        ind <- order(v)
        v = v[ind][1:m]
        iv = iv[ind][1:m]
        if(ncol(inputm) < 2*r-1 + 2*rn*(index-1))
            inputm <<- cbind(inputm, matrix(0, nrow=m, ncol=2*rn))
        inputm[, cached.keys(index, r)] <<- iv
        inputm[, cached.values(index, r)] <<- v
        rel_l[[relations[r]]] = rel.entity(newk, relations[r], rel)
    }
    reldb_df[[index]] <<- rel
    names(reldb_df)[index] <<- newk
    reldb_l[[index]] <<- rel_l
    names(reldb_l)[index] <<- newk
    keywordsdb[[newk]] <<- index
    for(k in cluster) delete.keyword(k)
}

# keyword - ambiguous keyword
# creates an intersection of relations between keyword and cluster,
# as used in intersect clustering
intersect.reldf <- function(keyword, cluster) {
    relk <- get.reldf(keyword)
    entities <- c()
    for(k in cluster) {
        rel <- get.reldf(k)
        entities = unique(c(entities, intersect_C(relk$entity, rel$entity)))
    }
    relk[relk$entity %in% entities,]
}


# creates pseudo-keyword from cluster (set of ids) of keywords
# keyword - ambiguous keyword
# returns keyword object
create.pseudo <- function(keyword, cluster) {
    rel_df <- intersect.reldf(keyword, cluster)
    m <- nrow(inputm)
    rel <- matrix(0, nrow=m, ncol=rn)
    irel <- matrix(0, nrow=m, ncol=rn)
    rel_l <- list()
    words <- unique(unlist(lapply(c(keyword.index(keyword), cluster), related.keywords, threshold=0)))

    # NOTE: a full re-calculation of co-occurrence values
    for(r in 1:rn) {
        rname <- relations[r]
        # index does not matter here since it is a pseudo keyword and data frame is provided
        rel_l[[rname]] = rel.entity(1, rname, rel_df)
        # NOTE: not full cooccurrence calculation
        co_m <- calc.cooccurrence.cached(rname, words)
        irel[,r] = co_m[,1]
        rel[,r] = co_m[,2]
    }
    cluster = c(keyword, cluster)
    list(irel=irel,
        rel=rel,
        super=unique(unlist(lapply(cluster, super))),
        sib=unique(unlist(lapply(cluster, sib))),
        df=rel_df,
        rellist=rel_l)
}

# ko - keyword object
# name - string for new keyword
# cluster - all keywords used in intersect clustering
add.pseudo <- function(ko, name, cluster) {
    index <- new.index()
    m <- nrow(inputm)
    for(r in 1:rn) {
        if(ncol(inputm) < cached.keys(index, r))
            inputm <<- cbind(inputm, matrix(0, nrow=m, ncol=2*rn))
        inputm[, cached.keys(index, r)] <<- ko$irel[,r]
        inputm[, cached.values(index, r)] <<- ko$rel[,r]
    }
    reldb_df[[index]] <<- ko$df
    reldb_l[[index]] <<- ko$rellist
    names(reldb_df)[index] <<- name
    names(reldb_l)[index] <<- name
    keywordsdb[[name]] <<- index
    index
}

# keyword must be character
# marks keyword as deleted or as not to be analysed, preserves its semantic
# relations and cleans up working semantic relations
delete.keyword <- function(keyword) {
    if(is.character(keyword)) {
        k = keyword
        ik <- keyword.index(keyword)
    } else {
        k = keyword.name(keyword)
        ik <- keyword.index(keyword)
    }
    for(i in seq_along(semantic)) {
        if(semrel$sizes[i] > 0) {
            todel <- c()
            semmatrix <- get.semantic(i)
            for(r in 1:nrow(semmatrix))
                if(any(semmatrix[r, ] == ik)) todel = c(todel, r)
            if(i < 4) {
                # preserve three semantic relations
                triples <<- rbind(triples, data.frame(list(
                    k1=unlist(vapply(semmatrix[todel, 1], keyword.name, "")),
                    k2=unlist(vapply(semmatrix[todel, 2], keyword.name, "")),
                    relation=rep(i, length(todel))), stringsAsFactors=FALSE))
            }
            if(length(todel) > 0) {
                semmatrix = semmatrix[-todel,]
                semrel[[i]][] <<- 0
                if(length(semmatrix) == 2) {
                    semrel[[i]][1, ] <<- semmatrix
                    semrel$sizes[i] <<- 1
                } else if(length(semmatrix) == 0) {
                    semrel$sizes[i] <<- 0
                } else { # nrow(semmatrix) > 1
                    semrel[[i]][1:nrow(semmatrix), ] <<- semmatrix
                    semrel$sizes[i] <<- nrow(semmatrix)
                }
            }
        }
    }

    # rm(list=k, envir=keywordsdb)
    # cannot delete element of reldb_df, instead:
    names(reldb_df)[ik] <<- NA
    names(reldb_l)[ik] <<- NA
    reldb_df[[ik]] <<- NA
    reldb_l[[ik]] <<- NA
}
