#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double Likelihood(List model, arma::cube dat,Function compute_states, Function gammas,Function ddens1, Function ddens0, arma::mat P_S,double ts) {
  
  NumericMatrix m = model["Q"];
  NumericMatrix theta =model["prior.theta"];
  int nstates = m.ncol();
  int times = dat.n_slices;
  
  arma::cube mem_ex;
  arma::rowvec cub_subset;
  NumericVector dat_vec1;
  NumericVector dat_vec0;
  
  double a= 0;
  double small_a= 0;
  double  ll=0 ;
  
  NumericVector temp(nstates) ;
  NumericMatrix gamma_val;
  NumericMatrix R_vec;
  
  
  //iterate over all perturbation experiments
  for(int q = 0 ; q < (m.nrow()); q = q + 1){
    
    R_vec = compute_states(model["W"], m(q ,_), model["r0"], (times),"ODE",P_S,ts);
    
    for( int y = 0 ; y<(times); y++){
      for( int x = 0 ; x<(nstates); x++)
      {if(R_vec(y,x) > 700)
        R_vec(y,x)=700; //* sign(R_vec[i]) in case of - states R_vec = R_vec % sign(R_vec);
      if(R_vec(y,x) < (-700))
        R_vec(y,x)= -700; //* sign(R_vec[i]) in case of - states R_vec = R_vec % sign(R_vec);
      }
    }
    gamma_val=gammas(R_vec,model["r0"]);
    //compute Eq. (3), ensuring you never take log(0)  plug into formula in supplements:
    //iterate over all E-nodes
    //add the second formula in the Supplements to ll
    for(int i=0; i < (dat.n_rows); i++){
      
      mem_ex=dat.tube(i,q);
      cub_subset= arma::rowvec(mem_ex.memptr(),mem_ex.n_elem,1,false );
      
      dat_vec1= ddens1(as<NumericVector>(wrap(cub_subset)));
      dat_vec0= ddens0(as<NumericVector>(wrap(cub_subset)));
      
      //iterate over all S-nodes
      for(int s=0; s<nstates;s++)
      {
        NumericVector temp_t(times);    
        for (int t=0; t<times;t++)
        {
          //temp[s]=sum(log((gamma_val(_ ,s)*dat_vec1) + ((1-gamma_val(_ ,s))*dat_vec0)+1));   //ddens_0 normal round 0, but we worked with absolute log change
          
          temp_t[t]= log(((gamma_val(t,s)*dat_vec1[t]) + ((1-gamma_val(t,s))*dat_vec0[t]))*(theta(i,s))+1);
          //Rcout<<"-------"<<log(((gamma_val(t,s)*dat_vec1[t]) + ((1-gamma_val(t,s))*dat_vec0[t]))*(theta(i,s))+1)<<"\n";
        }
        temp[s] = sum(temp_t);
        //       
      }
      a= max(temp);
      small_a=log(sum(exp(temp-a)));
      ll=ll+a+small_a;
    }
  }
  
  return (ll);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
