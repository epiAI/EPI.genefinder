# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

index <- function(mat) {
    .Call('EPI_genefinder_index', PACKAGE = 'EPI.genefinder', mat)
}

emission <- function(index, emissionProb) {
    .Call('EPI_genefinder_emission', PACKAGE = 'EPI.genefinder', index, emissionProb)
}

baumWelch <- function(index, nFeatures, nObservations, starts, ends, nContigs, emissionProb, updateEmission, transitionMat, updateTransition, initial, nHiddenStates, nPos, maxIter = 20L) {
    .Call('EPI_genefinder_baumWelch', PACKAGE = 'EPI.genefinder', index, nFeatures, nObservations, starts, ends, nContigs, emissionProb, updateEmission, transitionMat, updateTransition, initial, nHiddenStates, nPos, maxIter)
}

viterbi <- function(index, nFeatures, nObservations, starts, ends, nContigs, emissionProb, transitionMat, initial, nHiddenStates, nPos) {
    .Call('EPI_genefinder_viterbi', PACKAGE = 'EPI.genefinder', index, nFeatures, nObservations, starts, ends, nContigs, emissionProb, transitionMat, initial, nHiddenStates, nPos)
}

