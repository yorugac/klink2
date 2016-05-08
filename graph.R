library(igraph)
# edges - matrix n x 2 with n > 2
delete.cycles <- function(edges) {
    find.edge <- function(edges, e) {
        which(apply(t(edges) == e, 2, all))
    }
    # find cycle starting with given vertex
    detect.cycle <- function(g, v) {
        starts <- tail_of(g, incident(g, v, mode="out"))
        ends <- head_of(g, incident(g, v, mode="in"))
        for(i in starts) {
            for(j in ends) {
                p <- vertex.disjoint.paths(g, i, j)
                if(p) g[v,i] <<- 0
            }
        }
    }

    g <- graph(edges=t(edges), directed=T)
    if(is.dag(g)) return(edges)
    i <- 1
    while(TRUE) {
        detect.cycle(g, i)
        i = i + 1
        if(i == vcount(g)) break
    }
    as_edgelist(g)
}

tograph <- function(semrel) {
    # different types of edges?
    e <- rbind(semrel[[2]], semrel[[3]])
    g <- graph(edges=t(e), directed=T)
    g
}
