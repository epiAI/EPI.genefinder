#include <Rcpp.h>
#pragma once
// [[Rcpp::plugins(cpp11)]]



template<typename T>
class emissionFromIndex_t {
public:
  // pointer to begin of a matrix
  T* const ptr;
  // # of columns of this matrix
  const int ncol;
  // # of rows of this matrix
  const int nrow;
  // a map that links "positions" with the columns of the matrix
  int* const map;

  emissionFromIndex_t(T* D, int nc, int nr, int* m) :
    ptr(D), ncol(nc), nrow(nr), map(m) {};
  inline T* vecPtr(int pos) { return ptr + map[pos] * nrow;};


};

Rcpp::IntegerVector toInt(const Rcpp::IntegerMatrix& mat){
  Rcpp::IntegerVector res(mat.nrow());
  int tmp = 1;

  for (int col = 0; col < mat.ncol(); col++){


    for (int row = 0; row < mat.nrow(); row++){
      if (mat(row, col) == 1){
        res[row] += tmp;
      }
    }
    tmp *= 2;
  }
  return res;
}

// index data
// [[Rcpp::export]]
Rcpp::List index(const Rcpp::IntegerMatrix& mat){

  Rcpp::IntegerVector states = toInt(mat);
  Rcpp::IntegerVector values = Rcpp::sort_unique(states);
  Rcpp::IntegerVector map = Rcpp::match(states, values);
  Rcpp::IntegerMatrix observation(mat.ncol(), values.size());
  for (int col = 0; col < observation.ncol(); col++){
    int value = values[col];
    for (int row = 0; row < observation.nrow(); row++){
      observation(row, col) = value % 2;
      value /= 2;

    }

  }
  return Rcpp::List::create(
    Rcpp::Named("observation")= observation,
    Rcpp::Named("map") = map - 1 // 0-based index
  );

}

// calculate emissions
// [[Rcpp::export]]
Rcpp::NumericVector emission(const Rcpp::List& index, const Rcpp::NumericVector& emissionProb){
  Rcpp::NumericVector logp = Rcpp::log(emissionProb);
  Rcpp::NumericVector logq = Rcpp::log(1 - emissionProb);
  Rcpp::IntegerMatrix observation = Rcpp::as<Rcpp::IntegerMatrix>(index["observation"]);
  Rcpp::NumericVector res(observation.ncol());
  for (int i = 0; i < res.size(); i++){
    for (int row = 0; row < observation.nrow(); row++){
      if (observation(row, i) == 1)
        res[i] += logp(row);
      else
        res[i] += logq(row);
    }


  }

  return res;


}


