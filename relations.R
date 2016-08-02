# Handling keywords and their relations; internal functions behind Klink-2.

# It is assumed that there are global objects keywordsdb, reldb_df, reldb_l and inputm
# that hold all the keywords, relations and co-occurrences calculations;
# as well as semrel and triples for output.
#
# No keywords are truly deteled. Whether keyword should be included in analysis
# can be checked with keyword.exists -- whether there is a name for it in reldb_df.

maxindex <- NULL

# keywords are added/deleted only by interserct.clustering and similar -- once per klink-2 iteration;
# therefore they can be cached for all other times
cached.names <- NULL

zeromatrix <- NULL # used for semrel structure growth

# a short-hand for relation index; TODO
publicationi <- 1

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
    if(is.null(cached.names)) {
        n <- names(reldb_df)
        cached.names <<- n[!is.na(n)]
    }
    cached.names
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
calc.cooccurrence <- function(rname, i, irel_l, words=all.keywords()) {
    i_entities <- irel_l[[rname]]
    m <- nrow(inputm)
    res <- matrix(0, nrow=m, ncol=2)
    for(jk in words) {
        j <- keyword.index(jk)
        if(i != j) {
            v <- match(i_entities, reldb_l[[j]][[rname]], NA)
            value <- sum(!is.na(v))
            # is there place in the stored list for the value with j?
            im = which.min(res[,2])
            if(value > res[im,2]) {
                res[im,1] = j
                res[im,2] = value
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
# limiter - total number of unique entities for given words and rname; cooccurrence cannot be larger than limiter
calc.cooccurrence.cached <- function(rname, words, limiter=Inf) {
    m <- nrow(inputm)
    if(limiter == 0) return(matrix(0, nrow=m, ncol=2))
    iks <- vapply(words, keyword.index, 0)
    iv <- as.vector(inputm[, cached.keys(iks, r)])
    v <- as.vector(inputm[, cached.values(iks, r)])

    ind <- order(iv)
    ordiv <- iv[ind]
    ordv <- v[ind]
    outiv <- unique(ordiv)
    outv <- unlist(lapply(split(ordv, ordiv), sum), use.names=FALSE)

    # sorting in decreasing order of relation values
    ind = order(outv, decreasing=TRUE)
    # forming final result
    res <- matrix(0, nrow=m, ncol=2)
    res[, 1] = outiv[ind][1:m]
    res[, 2] = outv[ind][1:m]
    res[, 2] = pmin(res[, 2], rep(limiter, m))
    if(length(outiv) < m) res[(length(outiv)+1):m, ] = 0
    res
}

# strength of connection between k1 and k2 in regards to relation rel
# returns number of co-occurrences in case of unquantified relation and
# vector of minimum values in case of quantified
# k1, k2, rel - indexes
rel.value <- function(k1, k2, r) {
    if(k1==k2) return(length(rel.entity(k1, r)))
    i <- match(k2, inputm[, cached.keys(k1, r)])
    if(is.na(i)) {
         0
    } else {
        if(quantified[r]) {
            xdf <- get.reldf(k1)
            ydf <- get.reldf(k1)
            xv <- rel.entity(k1, r)
            xq <- rel.quantity(k1, r, xdf)
            yv <- rel.entity(k2, r)
            yq <- rel.quantity(k2, r, ydf)
            cooccur <- intersect(xv, yv)
            pmin(xq[which(xv==cooccur)], yq[which(yv==cooccur)])
        } else  {
            inputm[i, cached.values(k1, r)]
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
    i <- which(ik == c(semmatrix2[, 1], semmatrix3[,2]))
    if(length(i)) c(semmatrix2[, 2], semmatrix3[, 1])[i]
    else c()
}

# returns keywords that are siblings according to similarity relation
sib <- function(keyword) {
    ik <- keyword.index(keyword)
    semmatrix <- get.semantic(1)
    i <- which(ik == c(semmatrix[, 1], semmatrix[, 2]))
    if(length(i)) c(semmatrix[, 2], semmatrix[, 1])[i]
    else c()
}

# rel - string
# returns vector of entities shared by k1 and k2 via relation rel
common.entities <- function(k1, k2, rel) {
    intersect(rel.entity(k1, rel), rel.entity(k2, rel))
}

# total co-occurrence between k1 and k2
# without inputm
cooccur <- function(k1, k2) {
    v <- 0
    for(r in 1:rn) {
        rel = relations[r]
        v = v + length(common.entities(k1, k2, rel))
    }
    v
}

# values of connections for keyword and rel in form of matrix;
# with length equal to number of keywords
# keyword: id, string or keyword object
# List with 3 matrices is returned: 1st is filtered inputm,
# 2nd filtered for super keywords, 3rd filtered for sib keywords.
conn.vector <- function(rel, k) {
    n <- nkeywords()
    if(!is.list(k)) {
        if(!is.character(k)) kname = keyword.name(k)
        else kname = k
        if(!is.null(cachedCV[[kname]])) return(cachedCV[[kname]])
        k = keyword.index(k)
        super = super(k); sib = sib(k)
    } else {
        super = k$super; sib = k$sib
    }
    if(is.null(super)) super = numeric(0)
    if(is.null(sib)) sib = numeric(0)

    if(is.list(k)) {
        conn_vector_C(k$irel[, rel], k$rel[, rel], n, super, sib)
    } else {
        cachedCV[[kname]] <<-
            conn_vector_C(inputm[, cached.keys(k, rel)], inputm[, cached.values(k, rel)], n, super, sib)
        cachedCV[[kname]]
    }
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
npapers <- function(k) {
    length(reldb_l[[keyword.index(k)]][[publicationi]])
}

# returns entities from reldb_l for the given keyword and relation name
rel.entity <- function(k, r) {
    reldb_l[[keyword.index(k)]][[r]]
}

# construct entities (as they are stored in reldb_l) from given data frame
entity.from.df <- function(df, r) {
    sort(do.call(paste, c(df[df$relation==r, ][c('entity', 'year')], sep="_")))
}

# returns quantity(s) for the given keyword and relation name
rel.quantity <- function(keyword, relation, reldf=NULL) {
    if(is.null(reldf)) reldf = get.reldf(keyword)
    reldf[reldf$relation %in% relation,]$quantity
    # get rid of NAs?
}

# extract year(s) from given entities
ent.year <- function(entities) {
    as.integer(vapply(entities, stri_sub, "", -4))
}

# when keyword was first used
debut <- function(keyword) {
    min(get.reldf(keyword)$year, na.rm=TRUE)
}

# returns vector of related keywords to the given keyword
# threshold - min number of entities to be shared by related keywords
# (NOTE: among ones that are saved in input matrix)
related.keywords <- function(keyword, threshold=relkeyT) {
    ik <- keyword.index(keyword)
    v <- inputm[, cached.values(ik, 1:rn)]
    iv <- inputm[, cached.keys(ik, 1:rn)]
    rk <- list()
    for(i in 1:rn) rk[[i]] = iv[which(v[, i] > threshold[i]), i]
    rk = unique(unlist(rk))
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
    if(ncol(inputm) < cached.keys(index, 1))
        inputm <<- cbind(inputm, matrix(0, nrow=m, ncol=2*rn))
    for(r in 1:rn) {
        v <- as.vector(inputm[, cached.values(cluster, r)])
        iv <- as.vector(inputm[, cached.keys(cluster, r)])
        ind <- order(v)
        v = v[ind][1:m]
        iv = iv[ind][1:m]
        inputm[, cached.keys(index, r)] <<- iv
        inputm[, cached.values(index, r)] <<- v
        rel_l[[relations[r]]] = entity.from.df(rel, r)
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
        entities = unique(c(entities, intersect(relk$entity, rel$entity)))
    }
    df = relk[relk$entity %in% entities,]
    df
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
    words <- unique(unlist(lapply(c(keyword.index(keyword), cluster), related.keywords, threshold=0), use.names=FALSE))
    # words <- c(keyword, cluster)

    for(r in 1:rn) {
        rname <- relations[r]
        rel_l[[rname]] = entity.from.df(rel_df, relation.index(rname))
        # NOTE: not full cooccurrence calculation
        co_m <- calc.cooccurrence(rname, 1, rel_l, words)
        irel[,r] = co_m[,1]
        rel[,r] = co_m[,2]
    }
    cluster = c(keyword, cluster)
    list(irel=irel,
        rel=rel,
        super=unique(unlist(lapply(cluster, super), use.names=FALSE)),
        sib=unique(unlist(lapply(cluster, sib), use.names=FALSE)),
        df=rel_df,
        rellist=rel_l)
}

# pseudos - list of keyword objects
# i1 and i2 are to be "merged"
# cluster - keywords of "merged" cluster
# list with updated pseudos is to be returned
update.pseudos <- function(pseudos, i1, i2, cluster) {
    pseudos[[i1]] = create.pseudo(cluster[1], tail(cluster, -1))
    pseudos[-i2]
}

# add pseudo keywords in a batch
# k - ambiguous keyword (index)
# clusters - list of clusters
# size of pseudos and clusters must be equal
add.pseudos <- function(k, clusters) {
    n <- length(clusters)
    m <- nrow(inputm)
    keysm <- matrix(0, nrow=m, ncol=rn * n)
    valuesm <- matrix(0, nrow=m, ncol=rn * n)
    stored_index <- numeric(n)

    # if(ncol(inputm) < cached.keys(nkeywords() + n, 1))
        inputm <<- cbind(inputm, matrix(0, nrow=m, ncol=(2 * rn * n)))

    for(i in seq_along(clusters)) {
        pseudoname = paste(keyword.name(k), " (", keyword.name(high.in.cluster(k, clusters[[i]])), ")", sep="")
        index = new.index()
        pseudo = create.pseudo(k, clusters[[i]])

        where <- seq((i-1)*4+1, (i-1)*4+4)
        keysm[, where] = pseudo$irel
        valuesm[, where] = pseudo$rel

        keywordsdb[[pseudoname]] <<- index
        reldb_df[[index]] <<- pseudo$df

        reldb_l[[index]] <<- pseudo$rellist
        names(reldb_df)[index] <<- pseudoname
        names(reldb_l)[index] <<- pseudoname
        stored_index[i] = index
    }
    ikeys <- as.vector(vapply(stored_index, cached.keys, numeric(4), 1:rn))
    ivalues <- as.vector(vapply(stored_index, cached.values, numeric(4), 1:rn))
    inputm[, ikeys] <<- keysm
    inputm[, ivalues] <<- valuesm
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
                    k1=unlist(vapply(semmatrix[todel, 1], keyword.name, ""), use.names=FALSE),
                    k2=unlist(vapply(semmatrix[todel, 2], keyword.name, ""), use.names=FALSE),
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
