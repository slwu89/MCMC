#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


// 2D random walker
// [[Rcpp::export]]
NumericMatrix random_walk2d(int n){
  NumericMatrix out(Dimension(n,2));
  out(0,_) = NumericVector(2,0.0);
  for(int i=1; i<n; i++){
    double r_num = 1 + floor(R::runif(0,1) * 4);
    if(r_num == 1){ //move +1 in x
      out(i,_) = out(i-1,_);
      out(i,0) = out(i-1,0) + 1;
    }
    if(r_num == 2){ //move -1 in x
      out(i,_) = out(i-1,_);
      out(i,0) = out(i-1,0) - 1;
    }
    if(r_num == 3){ //move +1 in y
      out(i,_) = out(i-1,_);
      out(i,1) = out(i-1,1) + 1;
    }
    if(r_num == 4){ //move -1 in y
      out(i,_) = out(i-1,_);
      out(i,1) = out(i-1,1) - 1;
    }
    if(i % 10 == 0){ //print diagnostics
      Rcout << "Currently at step: " << i << std::endl;
    }
  }
  return(out);
}


// 3D random walker
// [[Rcpp::export]]
NumericMatrix random_walk3d(int n){
  NumericMatrix out(Dimension(n,3));
  out(0,_) = NumericVector(3,0.0);
  for(int i=1; i < n; i++){
    //6 movement possibilities in 3 axes at each step
    double r_num = 1 + floor(R::runif(0,1) * 6);
    if(r_num == 1){ //move +1 in x
      out(i,_) = out(i-1,_);
      out(i,0) = out(i-1,0) + 1;
    }
    if(r_num == 2){ //move -1 in x
      out(i,_) = out(i-1,_);
      out(i,0) = out(i-1,0) - 1;
    }
    if(r_num == 3){ //move +1 in y
      out(i,_) = out(i-1,_);
      out(i,1) = out(i-1,1) + 1;
    }
    if(r_num == 4){ //move -1 in y
      out(i,_) = out(i-1,_);
      out(i,1) = out(i-1,1) - 1;
    }
    if(r_num == 5){ //move +1 in z
      out(i,_) = out(i-1,_);
      out(i,2) = out(i-1,2) + 1;
    }
    if(r_num == 6){ //move -1 in z
      out(i,_) = out(i-1,_);
      out(i,2) = out(i-1,2) - 1;
    }
    if(i % 10 == 0){ //print diagnostics
      Rcout << "Currently at step: " << i << std::endl;
    }
  }
  return(out);
}


// Random walker on graph
// [[Rcpp::export]]
NumericVector random_walkGraph(NumericMatrix graph, int init_pos, int n){
  NumericVector out(n);
  out(0) = init_pos;
  int pos_begin;
  int pos_end;
  int num_nodes = graph.ncol();
  IntegerVector node_indices = seq_len(num_nodes) - 1;
  for(int i=1; i<n; i++){
    pos_begin = out(i-1);
    NumericVector trans_i = graph(pos_begin,_);
    //sample a node to transition to
    IntegerVector samp = RcppArmadillo::sample(node_indices,1,false,trans_i);
    pos_end = samp(0);
    out(i) = pos_end;
    if(i % 10 == 0){ //print diagnostics
      Rcout << "Currently at step: " << i << std::endl;
    }
    // Rcout << "Currently at step: " << i << std::endl;
  }
  return(out);
}
