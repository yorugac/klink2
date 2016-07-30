#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>

// TODO what is going on with character cases in R vs std sorting?

// [[Rcpp::export]]
// v1 and v2 are previously sorted
std::vector<std::string> intersect_C(const std::vector<std::string> & v1, const std::vector<std::string> & v2) {
    std::vector<std::string> v(v1.size() + v2.size());
    std::vector<std::string>::iterator it = std::set_intersection(
        v1.begin(), v1.end(),
        v2.begin(), v2.end(), v.begin());
    v.resize(std::distance(v.begin(), it));
    return v;
}

// [[Rcpp::export]]
int cooccurrence_C(List df1, List df2, CharacterVector entities) {
    if(!entities.size()) return 0;

    CharacterVector ent1 = as<CharacterVector>(df1["entity"]);
    NumericVector y1 = as<NumericVector>(df1["year"]);
    CharacterVector ent2 = as<CharacterVector>(df2["entity"]);
    NumericVector y2 = as<NumericVector>(df2["year"]);
    int n1 = ent1.size(), n2 = ent2.size(), n = entities.size();

    std::vector<int> matched_y1, matched_y2;
    for(int j = 0; j < n; ++j) {
        size_t i = std::distance(ent1.begin(), std::find(ent1.begin(), ent1.end(), entities[j]));
        if(i != n1) matched_y1.push_back(y1[i]);
        i = std::distance(ent2.begin(), std::find(ent2.begin(), ent2.end(), entities[j]));
        if(i != n2) matched_y2.push_back(y2[i]);
    }
    std::vector<int> res(matched_y1.size() + matched_y2.size());
    std::sort(matched_y1.begin(), matched_y1.end());
    std::sort(matched_y2.begin(), matched_y2.end());
    std::vector<int>::iterator it = std::set_intersection(
        matched_y1.begin(), matched_y1.end(),
        matched_y2.begin(), matched_y2.end(), res.begin());
    return std::distance(res.begin(), it);
}

bool keyvalue_cmp(const std::pair<int, int>& a, const std::pair<int, int>& b) {
    return a.second > b.second;
}

// [[Rcpp::export]]
// for use in data pre-processing only, since it considers only keywords with indexes 1..n
// n - number of keywords
// m - number of cached values to save
// i - index of keyword which is analyzed
// r - number of relation
// maxsize - size of largest possible entities vector
NumericMatrix calc_cooccurrence_C(const int n, const int m, const int i, int r, const List& reldb_l, int maxsize=1000) {
    CharacterVector ientities = as<List>(reldb_l[i-1])[r-1], jentities;
    std::vector<std::pair<int, int>> v(m);
    std::vector<std::string> interv(maxsize + ientities.size());
    auto zeros = std::make_pair(0, 0);
    int value, j;
    for(j = 0; j < n; ++j) {
        if(i-1 != j) {
            jentities = as<List>(reldb_l[j])[r-1];
            auto it = std::set_intersection(
                ientities.begin(), ientities.end(),
                jentities.begin(), jentities.end(), interv.begin());
            value = std::distance(interv.begin(), it);
            if(value) {
                auto p = std::make_pair(j+1, value);
                auto it = std::lower_bound(v.begin(), v.end(), zeros, keyvalue_cmp);
                if(it == v.end()) // no more free space
                    *(std::lower_bound(v.begin(), v.end(), p, keyvalue_cmp)) = p;
                else if(it != v.end()) *it = p;
            }
        }
    }
    std::sort(v.begin(), v.end(), keyvalue_cmp);
    NumericMatrix res(m, 2);
    for(j = 0; j < m; ++j) { res(j, 0) = v[j].first; res(j, 1) = v[j].second; }
    return res;
}

// [[Rcpp::export]]
// number of common characters
int common_chars_C(CharacterVector w1, CharacterVector w2) {
    // w1 and w2 consist only of one word
    std::string s1(w1[0]);
    std::string s2(w2[0]);
    std::sort(s1.begin(), s1.end());
    std::sort(s2.begin(), s2.end());
    s1.erase(unique(s1.begin(), s1.end()), s1.end());
    s2.erase(unique(s2.begin(), s2.end()), s2.end());
    std::string v(s1.size() + s2.size(), ' ');
    std::string::iterator it = std::set_intersection(
        s1.begin(), s1.end(),
        s2.begin(), s2.end(), v.begin());
    return std::distance(v.begin(), it);
}

// [[Rcpp::export]]
List conn_vector_C(NumericVector k, NumericVector v, int n, NumericVector super, NumericVector sib) {
    // std::sort(super.begin(), super.end());
    // std::sort(sib.begin(), sib.end());
    std::vector<int> ind1, ind2, ind3;
    int i;
    // k and v are of the same size
    for(i = 0; i < k.size(); ++i) {
        if(v[i] > 0 && k[i] <= n) {
            ind1.push_back(i);
            if(std::find(super.begin(), super.end(), k[i]) != super.end()) {
                ind2.push_back(i);
            }
            if(std::find(sib.begin(), sib.end(), k[i]) != sib.end()) {
                ind3.push_back(i);
            }
        }
    }
    List l(3);
    int is = ind1.size();
    NumericMatrix m(is, 2);
    for(i = 0; i < is; ++i) {
        m(i, 0) = k[ind1[i]];
        m(i, 1) = v[ind1[i]];
    }
    l[0] = m;
    is = ind2.size();
    m = NumericMatrix(is, 2);
    for(i = 0; i < is; ++i) {
        m(i, 0) = k[ind2[i]];
        m(i, 1) = v[ind2[i]];
    }
    l[1] = m;
    is = ind3.size();
    m = NumericMatrix(is, 2);
    for(i = 0; i < is; ++i) {
        m(i, 0) = k[ind3[i]];
        m(i, 1) = v[ind3[i]];
    }
    l[2] = m;
    // faster way? more idiomatic?
    return l;
}

// [[Rcpp::export]]
NumericMatrix update_dist_C(NumericMatrix d, List clusters, NumericVector weights) {
    int s = clusters.size();
    double t, w;
    NumericMatrix d2(s, s);
    NumericVector cl1, cl2;
    for(int i = 0; i < s; ++i)
        for(int j = i; j < s; ++j) {
            if(i != j) {
                cl1 = as<NumericVector>(clusters[i]);
                cl2 = as<NumericVector>(clusters[j]);
                t = 0; w = 0;
                for(int k1 = 0; k1 < cl1.size(); ++k1) {
                    for(int k2 = 0; k2 < cl2.size(); ++k2) {
                        t += d(cl1[k1]-1, cl2[k2]-1) * weights[cl1[k1]-1];
                        w += weights[cl1[k1]-1];
                    }
                }
                // weights cannot be 0
                d2(i, j) = t / w;
                d2(j, i) = d2(i, j);
            } else d2(i, j) = R_PosInf;
        }
    return d2;
}

double ss_accumf(const int& a, const int& b) {
    return a + std::pow(b, 2);
}
// [[Rcpp::export]]
double semantic_similarity_C(NumericMatrix x, NumericMatrix y, bool run_complete=false) {
    int xn = x.nrow(), j;
    NumericMatrix::Column kx = x.column(0);
    NumericMatrix::Column vx = x.column(1);
    NumericMatrix::Column ky = y.column(0);
    NumericMatrix::Column vy = y.column(1);
    double t = 0;
    // smth faster?
    for(int i = 0; i < xn; ++i) {
        j = std::distance(ky.begin(), std::find(ky.begin(), ky.end(), kx[i]));
        if(j != ky.size()) t += vx[i] * vy[j];
    }
    if(run_complete) {
        double d = std::sqrt(std::accumulate(vx.begin(), vx.end(), 0, ss_accumf)) *
            std::sqrt(std::accumulate(vy.begin(), vy.end(), 0, ss_accumf));
        if(d == 0) return 0;
        return t / d;
    }
    return t;
}

// [[Rcpp::export]]
double S_metric_C(const List& vx, const List& vy) {
    return semantic_similarity_C(vx[0], vy[0], true) /
        (std::max(semantic_similarity_C(vx[1], vy[1], true),
                    semantic_similarity_C(vx[2], vy[2], true)) + 1);
}

// // [[Rcpp::export]]
// NumericMatrix distance_matrix_C(int n, const List& connvectors, const NumericVector& largestS, int rn) {
//     double t, s;
//     NumericMatrix d(n, n);
//     for(int i = 0; i < n; ++i)
//         for(int j = i; j < n; ++j) {
//             if(i != j) {
//                 t = 0;
//                 for(int r = 1; i <= rn; ++r) {
//                     s = S_metric_C(as<List>(connvectors[i*rn + r]), as<List>(connvectors[j*rn + r]));
//                     t += s[0] / largestS[r-1];
//                 }
//                 t = 1 - t / rn;
//                 d(i, j) = t;
//                 d(j, i) = t;
//             } else d(i, j) = R_PosInf;
//         }
//     return d;
// }
