#define VCL_NAMESPACE vcl
#include "vectorclass.h"
using namespace vcl;

#include "vectormath_exp.h"
#include <math.h>


extern "C" double vcl_logsum(int N, double alpha, double *x){

  int i=0, j=0;
  double total = 0.0;
  Vec4d xvec, yvec;
  if (N>4){
    for(i=0; i<N-4+1; i+=4){
      xvec.load(x+i);
      yvec = vcl::log(xvec * alpha);
      total += vcl::horizontal_add(yvec);
    }
  }
  for(j=i; j<N; j++){
    total += log(x[j]*alpha);
  }
  return total;
}


extern "C" double vcl_sum(int N, double *x){

  int i=0, j=0;
  double total = 0.0;
  Vec4d xvec;
  if (N>4){
    for(i=0; i<N-4+1; i+=4){
      xvec.load(x+i);
      total += vcl::horizontal_add(xvec);
    }
  }
  for(j=i; j<N; j++){
    total += x[j];
  }
  return total;
}
