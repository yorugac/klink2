#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>

// [[Rcpp::export]]
int match_C(int j, NumericVector v) {
    for(int i = 0; i < v.size(); ++i)
        if(v[i] == j) return i + 1;
    return NA_INTEGER;
}

// [[Rcpp::export]]
NumericVector ent_year_C(List df, CharacterVector entities) {
    CharacterVector ent = as<CharacterVector>(df["entity"]);
    NumericVector y = as<NumericVector>(df["year"]);
    int n1 = ent.size(), n = entities.size();
    std::vector<int> matched;
    for(int j = 0; j < n; ++j) {
        size_t i = std::distance(ent.begin(), std::find(ent.begin(), ent.end(), entities[j]));
        if(i != n1) matched.push_back(y[i]);
    }
    return wrap(matched);
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
std::vector<std::string> intersect_C(CharacterVector v1, CharacterVector v2) {
    std::vector<std::string> v(v1.size() + v2.size());
    std::vector<std::string> s1 = as<std::vector<std::string> >(v1);
    std::vector<std::string> s2 = as<std::vector<std::string> >(v2);
    std::sort(s1.begin(), s1.end());
    std::sort(s2.begin(), s2.end());

    std::vector<std::string>::iterator it = std::set_intersection(
        s1.begin(), s1.end(),
        s2.begin(), s2.end(), v.begin());
    v.resize(it - v.begin());
    return v;
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

// bool keys_compare(const NumericMatrix::Row& a, const NumericMatrix::Row& b){
//     return a[0] < b[0];
// }

// [[Rcpp::export]]
double semantic_similarity_C(NumericMatrix x, NumericMatrix y) {
    NumericMatrix::Column kx = x.column(0);
    NumericMatrix::Column vx = x.column(1);
    NumericMatrix::Column ky = y.column(0);
    NumericMatrix::Column vy = y.column(1);
    int xn = x.nrow(), j;
    double t = 0;
    // smth faster?
    for(int i = 0; i < xn; ++i) {
        j = std::distance(ky.begin(), std::find(ky.begin(), ky.end(), kx[i]));
        if(j != ky.size()) t += vx[i] * vy[j];
    }
    return t;
}
