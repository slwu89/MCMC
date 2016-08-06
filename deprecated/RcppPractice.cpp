#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

//Various examples, practice, and functions in C++ and Armadillo from coding the RW and Adaptive Metropolis-Hastings
//and realizing that I am just using RStudio as my C++ IDE at this point...


// test creating and subsetting an arma vec
// [[Rcpp::export]]
arma::vec armaVec_create(){
  arma::vec mu(2);
  mu(0) = 0;
  mu(1) = 0;
  return(mu);
}

/*** R
message("Running armaVec_create!")
armaVec_create()
  */


// test feeding C++ function an R function
// [[Rcpp::export]]
double targetEval(NumericVector x, Function target){
  return(as<double>(target(x)));
}

/*** R
message("Running targetEval!")
target <- function(vec){
  x <- vec[1]
  y <- vec[2]
  z <- vec[3]
  return(x^2 + y^2 / z^2)
}
targetEval(x=c(1,2,3),target)
*/


// test feeding C++ function an R function with arma rowvec as input
// [[Rcpp::export]]
double targetEvalArma(arma::rowvec parm, Function target){
  return(as<double>(target(parm)));
}
/*** R
message("Running targetEvalArma!")
target <- function(vec){
  x <- vec[1]
  y <- vec[2]
  z <- vec[3]
  return(x^2 + y^2 / z^2)
}
targetEvalArma(parm=c(1,2,3),target)
*/


// test outputting a R list from C++ function
// [[Rcpp::export]]
List output_test(arma::rowvec theta_init, arma::mat sigma, int iterations){
  return(List::create(Named("theta")=theta_init,Named("sigma")=sigma,Named("iter")=iterations));
}
/*** R
message("Running output_test!")
output_test(c(1,2,3),diag(c(1,2,3)),5)
*/


// test subsetting an arma matrix by row
// [[Rcpp::export]]
arma::rowvec subsetRow_armaMat(){
  arma::mat A = arma::randu(10,10);
  return(A.row(0));
}
/*** R
message("Running subsetRow_armaMat!")
subsetRow_armaMat()
*/


// test setting up empty objects and output R list from C++
// [[Rcpp::export]]
List output_testArma(){
  //da matrix
  arma::mat A;
  A.zeros(10,3);
  //dat vector
  arma::rowvec B;
  B.zeros(20);
  //return that shit
  return(List::create(Named("da_matrix")=A,Named("yo_vector")=B));
}
/*** R
message("Running output_testArma!")
output_testArma()
*/


//test converting arma randu to double
// [[Rcpp::export]]
double randTestArma(){
  arma::vec armaRand = arma::randu(1);
  return(as<double>(wrap(armaRand(0))));
}
/*** R
message("Running randTestArma!")
randTestArma()
*/


//function to give dimensions of an arma matrix
// [[Rcpp::export]]
NumericVector dimensions(arma::mat x) {
  NumericVector out(2);
  out[0] = x.n_rows;
  out[1] = x.n_cols;
  return(out);
}
/*** R
message("Running dimensions!")
dimensions(diag(rep(5,3)))
dimensions(diag(rep(5,10)))
*/


//function to return empty Arma square matrix
// [[Rcpp::export]]
arma::mat zeroMat(int dim){
  arma::mat sigma_empirical = arma::zeros(dim,dim);
  return(sigma_empirical);
}
/*** R
message("Running zeroMat!")
zeroMat(10)
*/


//function to return 3d array; each matrix is drawn from runif
// [[Rcpp::export]]
arma::cube arrayC(int number, int dim){
  arma::cube out = arma::zeros(dim,dim,number);
  for(int i=0; i < number; i++){
    out.slice(i) = arma::randu(dim,dim);
  }
  return(out);
}
/*** R
message("Running arrayC!")
arrayC(5,4)
*/


// testing arithmetic with C++ bool class
// [[Rcpp::export]]
int boolTest(bool input){
  int out;
  out = input + 1;
  return(out);
}
/*** R
message("Running boolTest!")
boolTest(FALSE)
boolTest(TRUE)
*/


//testing "any" function in Rcpp sugar
// [[Rcpp::export]]
RObject anyTest(NumericVector input){
  LogicalVector out = Rcpp::any(input < 0.5);
  return(out);
}
/*** R
message("Running anyTest!")
anyTest(c(0.76,.8,1,5,3))
anyTest(runif(10))
*/


//testing modifying the diagonal of an arma matrix
// [[Rcpp::export]]
arma::mat diagTest(arma::mat input){
  int dim = input.n_rows;
  arma::mat out = input;
  out.diag() = arma::randu(dim);
  return(out);
}
/*** R
message("Running diagTest!")
diagTest(matrix(data=c(1,1,1,1,1,1,1,1,1),ncol=3))
*/


//testing finding length of arma vectors
// [[Rcpp::export]]
List lengthTest(arma::vec input){
  int elem = input.n_elem;
  int row = input.n_rows;
  int col = input.n_cols;
  return(List::create(Named("elem")=elem,Named("row")=row,Named("col")=col));
}
/*** R
message("Running lengthTest!")
lengthTest(rep(1,10))
*/


// testing matrix multiplication in arma
// [[Rcpp::export]]
arma::mat multiplyTest(arma::mat a, arma::mat b){
  arma::mat c;
  c = a*b;
  return(c);
}
/*** R
message("Running multiplyTest!")
a = matrix(data=c(4,5,3,6),2,2)
b = matrix(data=c(2,1,7,4),2,2)
a %*% b
multiplyTest(a,b)
*/

//testing armadillo vector printing to cout
// [[Rcpp::export]]
void printTest(){
  for(int i=0; i < 10; i++){
    arma::vec someVector = arma::randu(10);
    //need to fix, doesnt look good.
    Rcout << "at iteration i: " << i << ", someVector: " << someVector.t() << std::endl;
  }
}
/*** R
message("Running printTest!")
printTest()
*/


//testing matrix multiplication by a scalar in arma
// [[Rcpp::export]]
arma::mat scalarTest(arma::mat a, int b){
  arma::mat c;
  c = a * b;
  return(c);
}
/*** R
message("Running scalarTest!")
a = matrix(c(1,2,3,4),2,2)
b = 5
scalarTest(a,b)
*/


//testing vector to vector transpose to get matrix in arma
// [[Rcpp::export]]
arma::mat vecXvecTest(arma::rowvec a){
  arma::mat c;
  c = a.t() * a;
  return(c);
}
/*** R
message("Running vecXvecTest!")
a = c(1,2,3,4,5)
vecXvecTest(a)
d <- 1:5
d %*% t(d)
*/

//testing vector to vector transpose to get matrix in arma
// [[Rcpp::export]]
arma::mat colXcolcTest(arma::vec a){
  arma::mat c;
  c = a * a.t();
  return(c);
}
/*** R
message("Running colXcolcTest!")
a = c(1,2,3,4,5)
colXcolcTest(a)
d <- 1:5
d %*% t(d)
*/


//testing logical subsetting over a vector
// [[Rcpp::export]]
LogicalVector boolVec(int n, NumericVector low, NumericVector high){
  arma::vec vec = arma::randu(n);
  LogicalVector out(n);
  for(int i=0; i < n; i++){
    if((vec(i) <= high[i]) && (vec(i) >= low[i])){
      out[i] = true;
    } else {
      out[i] = false;
    }
  }
  return(out);
}
/*** R
message("Running boolVec!")
n <- 10
low <- c(.1,.2,.3,.4,.5,.1,.2,.3,.4,.5)
high <- c(.9,.8,.7,.6,.55,.9,.8,.7,.6,.55)
boolVec(n,low,high)
*/


//testing using LogicalVector inside of C++ function because std::vector<bool> is bad
// [[Rcpp::export]]
void boolVecTest(LogicalVector vector){
  int n = vector.size();
  for(int i=0; i < n; i++){
    if(vector[i]==TRUE){
      Rcout << "at position " << i+1 << "it is TRUE" << std::endl;
    } else {
      Rcout << "at position " << i+1 << "it is FALSE" << std::endl;
    }
  }
}
/*** R
message("Running boolVecTest!")
boolVec <- c(T,F,T,T,T,F,F,T,T,T,T,F,T,F,F,F,F,T,T,F)
boolVecTest(boolVec)
*/


//what the hell does arma::sum(x%x) do?
// [[Rcpp::export]]
List armaSumTest(arma::vec x){
  arma::vec armaMultiply;
  armaMultiply = x%x;
  double armaSum;
  armaSum = arma::sum(x%x);
  return(List::create(Named("mult")=armaMultiply,Named("sum")=armaSum));
}
/*** R
message("Running armaSumTest!")
armaSumTest(c(1,2,3))
*/


// //trying to call an R function from a package within C++
// // [[Rcpp::export]]
// NumericVector packageFuncTest(int n, NumericVector mu, NumericMatrix Sigma){
//   NumericVector out(n);
//   Environment MASS("package:MASS"); 
//   Function mvrnorm = MASS["mvrnorm"];
//   out = mvrnorm(n,mu,Sigma);
//   return(out);
// }
// /*** R
// message("Running packageFuncTest!")
// packageFuncTest(5,c(1,2,3,4,5),diag(c(1,1,1,1,1)))
// */
// 
// 
// //trying to call an R function with confusing output from a package within C++
// // [[Rcpp::export]]
// double packageFuncWeirdShit(NumericVector low, NumericVector high, NumericVector mean, NumericMatrix sigma){
//   SEXP out;
//   Environment mvtnorm("package:mvtnorm");
//   Function pmvnorm = mvtnorm["pmvnorm"];
//   out = pmvnorm(low,high,mean,R_NilValue,sigma);
//   double real_out = as<double>(out);
//   return(real_out);
// }
// /*** R
// message("Running packageFuncWeirdShit!")
// require(mvtnorm)
// packageFuncWeirdShit(-Inf,c(2,2),c(1,1),diag(2)*2)
//   */


//infinity in C++
// [[Rcpp::export]]
double inf_test(){
  double out = log(0.0);
  return(out);
}
/***R
message("Running inf_test")
inf_test()
*/


//roll a 6-sided dice in C++
// [[Rcpp::export]]
double d6_roll(){
  double out = 1 + (rand() % 6);
  return(out);
}
/*** R
message("Running d6_roll")
replicate(10,d6_roll())
*/


//roll a 6-sided dice in C++ using Rcpp sugar
// [[Rcpp::export]]
double d6_rollSugar(){
  double roll = R::runif(0,1);
  double out = 1 + floor(roll * 6);
  return(out);
}
/***R
message("Running d6_rollSugar")
replicate(10,d6_rollSugar())
*/


//testing RcppArmadillo sample function
// [[Rcpp::export]]
double cpp_sample(NumericVector x, NumericVector probs){
  double out;
  NumericVector samp = RcppArmadillo::sample(x,1,false,probs);
  out = samp(0);
  return(out);
}
/***R
message("Running cpp_sample")
vector <- 1:5
rand <- runif(5)
rand <- rand/sum(rand)
cpp_sample(vector,rand)
*/


