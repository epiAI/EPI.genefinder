// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// index
Rcpp::List index(const Rcpp::IntegerMatrix& mat);
RcppExport SEXP EPI_genefinder_index(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(index(mat));
    return rcpp_result_gen;
END_RCPP
}
// emission
Rcpp::NumericVector emission(const Rcpp::List& index, const Rcpp::NumericVector& emissionProb);
RcppExport SEXP EPI_genefinder_emission(SEXP indexSEXP, SEXP emissionProbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type index(indexSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type emissionProb(emissionProbSEXP);
    rcpp_result_gen = Rcpp::wrap(emission(index, emissionProb));
    return rcpp_result_gen;
END_RCPP
}
// baumWelch
Rcpp::List baumWelch(const Rcpp::List& index, size_t nFeatures, size_t nObservations, const Rcpp::IntegerVector& starts, const Rcpp::IntegerVector& ends, size_t nContigs, Rcpp::NumericMatrix& emissionProb, const Rcpp::LogicalVector& updateEmission, Rcpp::NumericMatrix& transitionMat, const Rcpp::LogicalMatrix& updateTransition, Rcpp::NumericVector& initial, size_t nHiddenStates, size_t nPos, size_t maxIter);
RcppExport SEXP EPI_genefinder_baumWelch(SEXP indexSEXP, SEXP nFeaturesSEXP, SEXP nObservationsSEXP, SEXP startsSEXP, SEXP endsSEXP, SEXP nContigsSEXP, SEXP emissionProbSEXP, SEXP updateEmissionSEXP, SEXP transitionMatSEXP, SEXP updateTransitionSEXP, SEXP initialSEXP, SEXP nHiddenStatesSEXP, SEXP nPosSEXP, SEXP maxIterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type index(indexSEXP);
    Rcpp::traits::input_parameter< size_t >::type nFeatures(nFeaturesSEXP);
    Rcpp::traits::input_parameter< size_t >::type nObservations(nObservationsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type starts(startsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type ends(endsSEXP);
    Rcpp::traits::input_parameter< size_t >::type nContigs(nContigsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type emissionProb(emissionProbSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector& >::type updateEmission(updateEmissionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type transitionMat(transitionMatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalMatrix& >::type updateTransition(updateTransitionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type initial(initialSEXP);
    Rcpp::traits::input_parameter< size_t >::type nHiddenStates(nHiddenStatesSEXP);
    Rcpp::traits::input_parameter< size_t >::type nPos(nPosSEXP);
    Rcpp::traits::input_parameter< size_t >::type maxIter(maxIterSEXP);
    rcpp_result_gen = Rcpp::wrap(baumWelch(index, nFeatures, nObservations, starts, ends, nContigs, emissionProb, updateEmission, transitionMat, updateTransition, initial, nHiddenStates, nPos, maxIter));
    return rcpp_result_gen;
END_RCPP
}
// viterbi
Rcpp::IntegerVector viterbi(const Rcpp::List& index, size_t nFeatures, size_t nObservations, const Rcpp::IntegerVector& starts, const Rcpp::IntegerVector& ends, size_t nContigs, Rcpp::NumericMatrix& emissionProb, Rcpp::NumericMatrix& transitionMat, Rcpp::NumericVector& initial, size_t nHiddenStates, size_t nPos);
RcppExport SEXP EPI_genefinder_viterbi(SEXP indexSEXP, SEXP nFeaturesSEXP, SEXP nObservationsSEXP, SEXP startsSEXP, SEXP endsSEXP, SEXP nContigsSEXP, SEXP emissionProbSEXP, SEXP transitionMatSEXP, SEXP initialSEXP, SEXP nHiddenStatesSEXP, SEXP nPosSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type index(indexSEXP);
    Rcpp::traits::input_parameter< size_t >::type nFeatures(nFeaturesSEXP);
    Rcpp::traits::input_parameter< size_t >::type nObservations(nObservationsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type starts(startsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type ends(endsSEXP);
    Rcpp::traits::input_parameter< size_t >::type nContigs(nContigsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type emissionProb(emissionProbSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type transitionMat(transitionMatSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type initial(initialSEXP);
    Rcpp::traits::input_parameter< size_t >::type nHiddenStates(nHiddenStatesSEXP);
    Rcpp::traits::input_parameter< size_t >::type nPos(nPosSEXP);
    rcpp_result_gen = Rcpp::wrap(viterbi(index, nFeatures, nObservations, starts, ends, nContigs, emissionProb, transitionMat, initial, nHiddenStates, nPos));
    return rcpp_result_gen;
END_RCPP
}
