#include <Rcpp.h>
#include "emission_t.h"

inline
double logSumExp(const Rcpp::NumericVector& x){
  double maxVal = Rcpp::max(x);
  return maxVal + log(Rcpp::sum(Rcpp::exp(x - maxVal)));
}

Rcpp::NumericMatrix toLogTransitionMat(const Rcpp::NumericMatrix& transitionMat, bool transpose){
  Rcpp::NumericMatrix res(transitionMat.nrow(), transitionMat.ncol());
  if (transpose){

    for (size_t row = 0; row < transitionMat.nrow(); row++){
      for (size_t col = 0; col < transitionMat.ncol(); col++){
        res(row, col) = log(transitionMat(col, row));
      }
    }

  }
  else{
    for (size_t row = 0; row < transitionMat.nrow(); row++){
      for (size_t col = 0; col < transitionMat.ncol(); col++){
        res(col, row) = log(transitionMat(col, row));
      }
    }

  }
  return res;

}


double forward(
  // emission related
  Rcpp::NumericMatrix& emissionTab,
  Rcpp::IntegerVector& map,
  // transition related
  const Rcpp::NumericMatrix& transitionMat,
  const Rcpp::NumericVector& initial,
  // dimensions
  // the contig
  const size_t start,
  const size_t end,
  // the model
  const size_t nHiddenStates,

  // return
  Rcpp::NumericMatrix& f
){
  emissionFromIndex_t<double> emissionMat(emissionTab.begin(),
    emissionTab.ncol(), emissionTab.nrow(), map.begin());

  Rcpp::NumericMatrix logTransitionMat = toLogTransitionMat(transitionMat, false);
  Rcpp::NumericVector tmp(nHiddenStates);
  size_t pos = start;
  size_t hiddenState;

  // first position transitions only from start
  // with transition probability initial
  double* emissionPtr = emissionMat.vecPtr(pos);

  for (hiddenState = 0; hiddenState < nHiddenStates; hiddenState++){

    f(hiddenState, pos) = log(initial(hiddenState)) + *(emissionPtr + hiddenState);
  }

  pos++;
  for (; pos <= end; pos++){
    emissionPtr = emissionMat.vecPtr(pos);

    for (hiddenState = 0; hiddenState < nHiddenStates; hiddenState++){
      tmp = logTransitionMat(hiddenState, Rcpp::_) + f(Rcpp::_, pos - 1);

      f(hiddenState, pos) = logSumExp(tmp) + *(emissionPtr + hiddenState);

    }
    //std::cout << pos << '\t' << logSumExp(f(Rcpp::_, pos)) << '\n';

  }
  // return loglikelihood
  return logSumExp(f(Rcpp::_, end));
}

void backward(
    // emission related
    Rcpp::NumericMatrix& emissionTab,
    Rcpp::IntegerVector& map,
    // transition related
    const Rcpp::NumericMatrix& transitionMat,
    // dimensions
    // the contig
    const size_t start,
    const size_t end,
    // the model
    const size_t nHiddenStates,

    // return
    Rcpp::NumericMatrix& b
){
  emissionFromIndex_t<double> emissionMat(emissionTab.begin(),
                                          emissionTab.ncol(), emissionTab.nrow(), map.begin());

  Rcpp::NumericMatrix logTransitionMat = toLogTransitionMat(transitionMat, true);
  Rcpp::NumericVector tmp(nHiddenStates);
  size_t pos;
  size_t hiddenState;
  double* emissionPtr;

  for (pos = end - 1; pos >= start && pos < end; pos--){

    emissionPtr = emissionMat.vecPtr(pos + 1);
    for (hiddenState = 0; hiddenState < nHiddenStates; hiddenState++){
      tmp[hiddenState] = *(emissionPtr + hiddenState);
    }
    for (hiddenState = 0; hiddenState < nHiddenStates; hiddenState++){
      b(hiddenState, pos) = logSumExp(tmp + logTransitionMat(hiddenState, Rcpp::_) + b(Rcpp::_, pos + 1));
    }
    //std::cout << pos << '\t' << logSumExp(b(Rcpp::_, pos)) << '\n';
  }
}

// baum-welch

void baumWelch_core(
    const Rcpp::List& index,
    size_t nFeatures,
    size_t nObservations,
    const Rcpp::IntegerVector& starts,
    const Rcpp::IntegerVector& ends,
    size_t nContigs,
    Rcpp::NumericMatrix& emissionProb,
    const Rcpp::LogicalVector& updateEmission,
    Rcpp::NumericMatrix& transitionMat,
    const Rcpp::LogicalMatrix& updateTransition,

    Rcpp::NumericVector& initial,
    size_t nHiddenStates,
    Rcpp::NumericMatrix& f,
    Rcpp::NumericMatrix& b,
    size_t maxIter


){

  size_t contig;
  size_t feature;
  size_t hiddenState;
  size_t fromHiddenState;
  size_t toHiddenState;
  Rcpp::NumericVector likelihood(nContigs);
  Rcpp::IntegerVector map = Rcpp::as<Rcpp::IntegerVector>(index["map"]);
  Rcpp::IntegerMatrix observation = Rcpp::as<Rcpp::IntegerMatrix>(index["observation"]);
  emissionFromIndex_t<int> observationMat(observation.begin(),
    observation.ncol(), observation.nrow(), map.begin());
  int* observationPtr;
  double currentPost;
  double lastLlik = -1e-300;
  for (size_t iter = 0; iter < maxIter; iter++){
    // calculate emissions
    // iterate over hiddenStates
    Rcpp::NumericMatrix emissionTab(nHiddenStates, nObservations);
    std::cout << "creating emission Tab\n";
    for (hiddenState = 0; hiddenState < nHiddenStates; hiddenState++){
      emissionTab(hiddenState, Rcpp::_)= emission(index, emissionProb(Rcpp::_, hiddenState));
    }
    std::cout << emissionTab.nrow() << '\t' << emissionTab.ncol() << '\n';
    //Rcpp::Rcout << emissionTab << std::endl;
    std::cout << "forward/backward\n";
    // forward/baclward
    for (contig = 0; contig < nContigs; contig++){

      likelihood(contig) = forward(emissionTab, map, transitionMat, initial,
        starts[contig], ends[contig], nHiddenStates, f);

      backward(emissionTab, map, transitionMat,
        starts[contig], ends[contig], nHiddenStates, b);
      std::cout << contig << '\t' << likelihood[contig] << '\n';

    }
    std::cout << "update emissions\n";
    // update emissions
    contig = 0;
    double currentLikelihood = likelihood[0];
    Rcpp::NumericMatrix newEmissionProb(nFeatures, nHiddenStates);
    Rcpp::NumericVector postHiddenStates(nHiddenStates);

    for (size_t pos = 0; pos < f.ncol(); pos++){
      if (pos > ends[contig]){
        contig++;
        currentLikelihood = likelihood[contig];
      }
      observationPtr = observationMat.vecPtr(pos);

      for (hiddenState = 0; hiddenState < nHiddenStates; hiddenState++){
        currentPost = exp(f(hiddenState, pos) + b(hiddenState, pos) - currentLikelihood);
        if (currentPost > 1){
          std::cout << pos << '\t' << currentPost << '\n';
        }
        postHiddenStates[hiddenState] += currentPost;
        if (updateEmission[hiddenState]){
          //currentPost = exp(f(hiddenState, pos) + b(hiddenState, pos) - currentLikelihood);
          //postHiddenStates[hiddenState] += currentPost;
          for (feature = 0; feature < nFeatures; feature++){
            if (*(observationPtr + feature) == 1){
              newEmissionProb(feature, hiddenState) += currentPost;

            }
          }

        }
      }

    }
    Rcpp::Rcout << "hiddenStates: " << postHiddenStates << std::endl;

    for (hiddenState = 0; hiddenState < nHiddenStates; hiddenState++){
      if (updateEmission[hiddenState]){
        newEmissionProb(Rcpp::_, hiddenState) = newEmissionProb(Rcpp::_, hiddenState) / postHiddenStates[hiddenState];
      }
      else{
        newEmissionProb(Rcpp::_, hiddenState) = emissionProb(Rcpp::_, hiddenState);
      }
    }

    double llik = Rcpp::sum(likelihood);
    if (llik - lastLlik < 0.01){

    }
    std::cout << iter << '\t' << llik << std::endl;
    std::cout << "update transitions\n";
    //update transitions

    Rcpp::NumericVector newInitial(initial.size(), 0);
    Rcpp::NumericMatrix newTransitionMat(transitionMat.nrow(), transitionMat.ncol());

    contig = 0;
    currentLikelihood = likelihood[0];
    emissionFromIndex_t<double> emissionMat(emissionTab.begin(),
                                            emissionTab.ncol(), emissionTab.nrow(), map.begin());
    double* emissionPtr;
    Rcpp::NumericMatrix logTransitionMat = toLogTransitionMat(transitionMat, false);
    for (size_t pos = 0; pos < f.ncol(); pos++){

      if (pos > ends[contig]){
        contig++;
        currentLikelihood = likelihood[contig];
      }
      if (pos == starts[contig]){
        newInitial = newInitial + Rcpp::exp(f(Rcpp::_, pos) + b(Rcpp::_, pos) - currentLikelihood);
      }
      else{
        emissionPtr = emissionMat.vecPtr(pos);

        for (fromHiddenState = 0; fromHiddenState < nHiddenStates; fromHiddenState++){
          for (toHiddenState = 0; toHiddenState < nHiddenStates; toHiddenState++){
            if (updateTransition(toHiddenState, fromHiddenState)){
              newTransitionMat(toHiddenState, fromHiddenState)
                += exp(
                    f(fromHiddenState, pos - 1)
                    + logTransitionMat(toHiddenState, fromHiddenState)
                    + *(emissionPtr + toHiddenState)
                    + b(toHiddenState, pos)
                    - currentLikelihood
                   );

            }
          }
        }
      }

    } // for (size_t pos = 0; pos < f.ncol(); pos++)
    Rcpp::NumericVector transitionSum(nHiddenStates);
    for (hiddenState = 0; hiddenState < nHiddenStates; hiddenState++){
      transitionSum[hiddenState] = Rcpp::sum(newTransitionMat(Rcpp::_, hiddenState));
    }

    for (fromHiddenState = 0; fromHiddenState < nHiddenStates; fromHiddenState++){
      for (toHiddenState = 0; toHiddenState < nHiddenStates; toHiddenState++){
        if (updateTransition(toHiddenState, fromHiddenState)){
          newTransitionMat(toHiddenState, fromHiddenState) /= transitionSum(fromHiddenState);
        }
        else{
          newTransitionMat(toHiddenState, fromHiddenState) = transitionMat(toHiddenState, fromHiddenState);
        }
      }

    }
    Rcpp::Rcout << newInitial << std::endl;
    // copy
    emissionProb = newEmissionProb;
    initial = newInitial / Rcpp::sum(newInitial);
    transitionMat = newTransitionMat;
    lastLlik = llik;

  } // for (size_t iter = 0; iter < maxIter; iter++)

}
// [[Rcpp::export]]
Rcpp::List baumWelch(
    const Rcpp::List& index,
    size_t nFeatures,
    size_t nObservations,
    const Rcpp::IntegerVector& starts,
    const Rcpp::IntegerVector& ends,
    size_t nContigs,

    Rcpp::NumericMatrix& emissionProb,
    const Rcpp::LogicalVector& updateEmission,
    Rcpp::NumericMatrix& transitionMat,
    const Rcpp::LogicalMatrix& updateTransition,

    Rcpp::NumericVector& initial,
    size_t nHiddenStates,
    size_t nPos,
    size_t maxIter = 20
){
  Rcpp::NumericMatrix f(nHiddenStates, nPos);
  Rcpp::NumericMatrix b(nHiddenStates, nPos);
  baumWelch_core(index, nFeatures, nObservations, starts, ends, nContigs, emissionProb, updateEmission,
                 transitionMat, updateTransition, initial, nHiddenStates, f, b, maxIter);
  return Rcpp::List::create(
    Rcpp::Named("emissionProb")= emissionProb,
    Rcpp::Named("transitionMat") = transitionMat,
    Rcpp::Named("initial") = initial
  );

}
// [[Rcpp::export]]
Rcpp::IntegerVector viterbi(
  const Rcpp::List& index,
  size_t nFeatures,
  size_t nObservations,
  const Rcpp::IntegerVector& starts,
  const Rcpp::IntegerVector& ends,
  size_t nContigs,
  Rcpp::NumericMatrix& emissionProb,
  Rcpp::NumericMatrix& transitionMat,

  Rcpp::NumericVector& initial,
  size_t nHiddenStates,
  size_t nPos
){
  Rcpp::IntegerVector map = Rcpp::as<Rcpp::IntegerVector>(index["map"]);
  Rcpp::NumericMatrix v(nHiddenStates, nPos);
  Rcpp::IntegerMatrix ptr(nHiddenStates, nPos);
  size_t pos;
  size_t hiddenState;
  size_t fromHiddenState;
  size_t maxHiddenState;
  Rcpp::NumericMatrix logTransitionMat = toLogTransitionMat(transitionMat, false);
  Rcpp::NumericMatrix emissionTab(nHiddenStates, nObservations);
  std::cout << "creating emission Tab\n";
  for (hiddenState = 0; hiddenState < nHiddenStates; hiddenState++){
    emissionTab(hiddenState, Rcpp::_)= emission(index, emissionProb(Rcpp::_, hiddenState));
  }

  emissionFromIndex_t<double> emissionMat(emissionTab.begin(),
                                          emissionTab.ncol(), emissionTab.nrow(), map.begin());
  size_t contig;
  double* emissionPtr;
  Rcpp::NumericVector tmp(nHiddenStates);
  Rcpp::IntegerVector stateSequence(nPos);
  for (contig = 0; contig < nContigs; contig++){
    // first iteration
    pos = starts[contig];
    emissionPtr = emissionMat.vecPtr(pos);
    for (hiddenState = 0; hiddenState < nHiddenStates; hiddenState++){
      v(hiddenState, pos) = log(initial(hiddenState)) + *(emissionPtr + hiddenState);
      ptr(hiddenState,pos) = -1;
    }
    pos++;

    for ( ; pos <= ends[contig]; pos++){
      emissionPtr = emissionMat.vecPtr(pos);
      for (hiddenState = 0; hiddenState < nHiddenStates; hiddenState++){
        tmp = logTransitionMat(hiddenState, Rcpp::_) + v(Rcpp::_, pos - 1);
        maxHiddenState = 0;

        for (fromHiddenState = 1; fromHiddenState < nHiddenStates; fromHiddenState++){
          if (tmp[fromHiddenState] > tmp[maxHiddenState]){
            maxHiddenState = fromHiddenState;
          }

        }
        v(hiddenState, pos) = *(emissionPtr + hiddenState) + tmp[maxHiddenState];
        ptr(hiddenState,pos) = maxHiddenState;
      }
    }
    // backtrace
    pos = ends[contig];
    maxHiddenState = 0;
    for (hiddenState = 0; hiddenState < nHiddenStates; hiddenState++){
      if (v(hiddenState, pos) > v(maxHiddenState, pos)){
        maxHiddenState = hiddenState;
      }
    }
    stateSequence[pos] = maxHiddenState;

    for (; pos > starts[contig]; pos--){
      maxHiddenState = ptr(maxHiddenState, pos);
      stateSequence[pos - 1] = maxHiddenState;
    }

  }
  return stateSequence + 1;
}
