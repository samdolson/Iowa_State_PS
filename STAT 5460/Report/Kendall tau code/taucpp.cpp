#include <RcppArmadillo.h>
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]


//[[Rcpp::export]]
mat bubbleSortX(vec x,vec y) {
  int n=x.size();
  mat  A = ones<mat>(n,2);
  int i,j;
  double tmpx,tmpy;
  for(i = n - 1; i > 0; i--) {
    for(j = 0; j < i; j++) {
      if(x(j+1)<x(j)){
        tmpx = x(j+1);
        x(j+1) = x(j);
        x(j) = tmpx;
        tmpy = y(j+1);
        y(j+1) = y(j);
        y(j) = tmpy;
      }
    }
    A(i,0)=x(i);
    A(i,1)=y(i);
    
  }
  
  A(0,0)=x(0);
  A(0,1)=y(0);
  return(A);
} 

//[[Rcpp::export]]
double mtau(vec x,vec y) {
  int n=x.size();
  mat  A = bubbleSortX(x,y);
  double c=0;
  for(int i=0;i<n;i++){
    double Ri=0;
    for(int j=0;j<i;j++){
      if(A(j,1)<A(i,1))
        Ri=Ri+1;
    }
    c+=Ri;
  }
  return( 4*c/n/(n-1)-1.0);
} 


// [[Rcpp::export]]
vec taucpp(mat X){
  vec TN = zeros<vec>(2);
  int p=X.n_cols, n=X.n_rows;
  mat Sigma2=zeros(p,p);
  for(int i=0;i<(p-1);i++){
    for(int j=(i+1);j<p;j++){
      Sigma2(i,j)= pow(mtau(X.col(i),X.col(j)),2);
    }
  }
  double Tn=0;
  Tn=sum(sum( Sigma2 ));
  double mu=p*(p-1)*(2*n+5)/9.0/n/(n-1);
  double sigmanp=4.0*p*(p-1)*(n-2)*(100*pow(n,3)+492*pow(n,2)+731*n+279)/2025.0/(pow(n,3)*pow(n-1,3));
  TN(0)=(Tn-mu)/sqrt(sigmanp);
  sigmanp=2.0*(2*n+5)/9.0/n/(n-1);
  TN(1)=max(max(Sigma2))/(sigmanp);
  return(TN);
}






