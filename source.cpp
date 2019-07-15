#include <R.h>
#include <Rmath.h>
#include <math.h>
// sudo apt-get install r-cran-rcpp r-cran-rcpparmadillo r-cran-rcppeigen
// [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
using namespace Rcpp;
// using namespace RcppArmadillo;

/* function to test the independence of two snps */
/* include standard C math library and R internal function declarations     */

//[[Rcpp::export]]
double mychisq(NumericVector x, NumericVector y)
{
  NumericVector p (x.size());
  p=rep(1/float(x.size()), x.size());
//   if (size(x) != size(y)){ 
   if (x.size() != y.size()){ 
    stop("'x' and 'y' must have the same length");
  }
  bool NAx, NAy;
  NAx=is_true(any(is_na(x)));
  NAy=is_true(any(is_na(y)));
  if (NAx || NAy){
	  stop("'x' and 'y' should not contain NAs");
  }
  if (x.size() != y.size()){ 
    stop("'x' and 'y' must have the same length");
  }
  NumericVector uniX,uniY;
  int LuniX,LuniY;
  uniX=unique(x);
  LuniX=uniX.size();
  uniY=unique(y);	
  LuniY=uniY.size();
  if(LuniX < 2 || LuniY < 2){
    stop("'x' and 'y' must have at least 2 levels");
  }
  NumericMatrix tabla(LuniX,LuniY);
  for (int i=0;i<LuniX;i++){
	  for(int j=0;j<LuniY;j++){
		int cnt =0;
		for(int k=0;k<x.size();k++){
		if (x[k]==uniX[i] && y[k]==uniY[j]){
			cnt++;
		}
		tabla(i,j)=cnt;
		}
		
	}
}
	if (is_true(any(is_na(tabla)))){
		stop("D");
	}
	unsigned int nr,nc,n;
	
	nr = tabla.nrow();
	nc = tabla.ncol();
	IntegerVector sr(nr);
	IntegerVector sc(nc);
	
	for (int i=0;i < nr;i++){
		sr[i] = sum(tabla(i,_));
	}
	for (int j=0;j < nc;j++){
		sc[j] = sum(tabla(_,j));
	}
	NumericMatrix E(nr,nc);
	n=sum(sr);
	for (unsigned int i=0; i<nr; i++) {
    for (unsigned int j=0; j<nc; j++)  {
      E(i,j)=sr[i]*sc[j];
    }
  }
	E=E/n;
	NumericVector times,nr2;
	NumericVector  Y(nc*nr);
	NumericVector  X(nc*nr);
	times=rep(float(nr), int(nc));
	IntegerVector times2, idx;
	times2=as<IntegerVector>(times);
	int l=0;
	int i,k;
  	for(i =0;i<nc;i++){
	for(k=l;k<l+times2[i];k++){
		
	Y[k]=sc[i];
 	}
 	l=l+times2[i];
	}
	float ax;
	ax=float(nc*nr) / float(nr);
	int aux;
 	aux=ceil(ax);
	X=rep(sr,aux);
	NumericVector V;
	float n3, YATES;
	n3=pow(n,3);
	V=Y * X * (n - X) * (n - Y)/n3;	
	NumericVector minimo(2);
	minimo[0]=0.5;
	minimo[1]=0.5;
	NumericMatrix prt(tabla.nrow(),tabla.ncol());
	
		for (int i =0;i<tabla.nrow();i++){
			for (int j =0;j<tabla.ncol();j++){
				prt(i,j)=fabs(tabla(i,j)-E(i,j));
				if (prt(i,j) < minimo[1]){
					minimo[1]=prt(i,j);
				}
		}
	}
	if(nr ==2 & nc == 2){
		YATES=min(minimo);
	}
else{
 	YATES=0;
	}
	double stat;
	double PVAL;
	NumericMatrix auxilio(tabla.nrow(),tabla.ncol());
	for (int i =0;i<tabla.nrow();i++){
		for (int j =0;j<tabla.ncol();j++){
			auxilio(i,j)=pow((prt(i,j)-YATES),2)/E(i,j);
		}
	}
	
	stat=sum(auxilio);
   double PARAMETER = (nr - 1L) * (nc - 1L);
   PVAL=Rf_pchisq(stat, PARAMETER,0,0);

 return PVAL;
	
}

//[[Rcpp::export]]
NumericVector callmychi(  NumericVector x, NumericMatrix y){
	NumericVector pval(y.ncol());
    NumericVector uniX,uniY;
    int LuniX,LuniY;
    
    for (int i =0;i<y.ncol();i++){
        uniX=unique(x);
        LuniX=uniX.size();
        uniY=unique(y(_,i));	
        LuniY=uniY.size();
        if(LuniX < 2 || LuniY < 2){
            pval[i]=1;
        }else{
		    pval[i]=mychisq(x,y(_,i));
	    }
    }
    return pval;
}

//[[Rcpp::export]]
Rcpp::DataFrame numVecEx2(Rcpp::NumericVector xs){
    Rcpp::NumericVector x1(xs);
    Rcpp::NumericVector x2(Rcpp::clone(xs));
    x1[0] = 22;
    x2[1] = 44;
    return(Rcpp::DataFrame::create(Named("orig", xs),
                                   Named("x1", x1),
                                   Named("x2", x2)));
}
