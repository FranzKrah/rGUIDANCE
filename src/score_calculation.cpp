// ## This code is part of the polenta package
// ## © F.-S. Krah (last update 2017-08-18)

// This Rcpp code implements the basic steps to calculate
// the residue pair score between a REF (reference) MSA and
// one or more ALT (alternative) MSA(s).

#include <Rcpp.h>
#include <cstdlib>
using namespace std;
using namespace Rcpp;


//[[Rcpp::export]]
int nChoosek( int n, int k )
{
  if (k > n) return 0;
  if (k * 2 > n) k = n-k;
  if (k == 0) return 1;

  int result = n;
  for( int i = 2; i <= k; ++i ) {
    result *= (n-i+1);
    result /= i;
  }
  return result;
}

/////////////////////////////////////////////////////////////////
// Cmatrix //
// For more info about the Cmatrix recoding see page 3 on the supplemantary material of:
// Satija R, Novak A., Mikls I., Lyngs R., and Hein J. (2009) BigFoot:
// Bayesian alignment and phylogenetic footprinting with MCMC,
// BMC Evolutionary Biology, 9, 217.
// http://www.biomedcentral.com/content/supplementary/1471-2148-9-217-s1.pdf
// Here implemented as in GUIDANCE program msa_set_score
// Penn et al. (2010). An alignment confidence score capturing robustness to guide tree
// uncertainty. Molecular Biology and Evolution. 27:1759--1767
// [[Rcpp::export]]
List msa_recode(NumericMatrix msa){

  int n_row = msa.nrow();
  int n_col = msa.ncol();

  NumericMatrix col2res(n_row, n_col);
  NumericMatrix res2col(n_row, n_col);

  for(int row=0; row<n_row; row++){
    int lastnongap = -1;
    int res = 0;

    for(int col=0; col < n_col; col++){
      if (msa(row,col) == 0){
        col2res(row,col) = lastnongap+1;
      } else {
        // 'res2col': indicies of the characters in the alignment
        res2col(row,res) = col;
        // 'col2res': characters are represented by odd numbers
        //            and gaps by even numbers
        col2res(row,col) = (2*res)+1;
        lastnongap = (2*res)+1;
        res++;
        // Rcout << count << endl;
      }
    }
  }
  return List::create(Rcpp::Named("col2res") = col2res,
    Rcpp::Named("res2col") = res2col);
}

/////////////////////////////////////////////////////////////////
// Produce res_pair_hit object as basis for compare residue pairs
// this matrix contains all residue pairs in the reference MSA and the information
// if they are bases (hit = 0) or gabs (-1); the hits are later updated by
// comparisons with the alternative MSAs (see add_msa)
// The function 'init_counters' in 'set_msa_score' uses 3D matrices,
// which is not straight forewardly implemented in Rcpp, which is why I
// work with a nxk matrix, with k many columns as the REF MSA and n rows for all
// residue pairs.
// This is almost equally fast, however, the subsequent calculation of derived
// scores is not as handy. These are done in the function 'daughter_scores'.
// [[Rcpp::export]]
NumericMatrix res_pair_hit(NumericMatrix col2res){

  // double n_row = col2res.nrow();
  // double n_col = col2res.ncol();

  double nrows = nChoosek(col2res.nrow(), 2);
  NumericMatrix respairhit(nrows, col2res.ncol());

  for(int col=0; col < col2res.ncol()-1; col++){

    double count = 0;
    for(int row1 = 0; row1 <= col2res.nrow()-2; row1++){
      int resup = row1+1;
      for(int row2 = resup; row2 <= col2res.nrow()-1; row2++){

        int res1 = col2res(row1, col);
        int res2 = col2res(row2, col);

        int oddornot1 = res1 % 2;
        int oddornot2 = res2 % 2;

        // respairhit(count,0) = col+1;
        // respairhit(count,1) = row1+1;
        // respairhit(count,2) = row2+1;

        // if both residues are bases => 0
        // if one of the basis is a gab => -1
        if ( oddornot1 == 1 && oddornot2 == 1 ) {
          respairhit(count,col) = 0;
        }else{
          respairhit(count,col) = -1;
        }
        count++;
      }
    }
  }
  return(respairhit);
}

/////////////////////////////////////////////////////////////////
// Based on the output of msa_recode and res_pair_hit the MSAs are compared
// and scores are added to the scores matrix (res_pair_hit matrix)
// [[Rcpp::export]]
IntegerMatrix add_msa(IntegerMatrix ref_col2res, IntegerMatrix alt_col2res,
  IntegerMatrix alt_res2col, IntegerMatrix hit_mat){

  // Rcout << "initializing done"<< endl;
  // int ncol = ref_col2res.ncol();
  // int nrow = ref_col2res.nrow();

  for(int col=0; col <= ref_col2res.ncol()-1; col++){

    double count = 0;
    for(int row1 = 0; row1 <= ref_col2res.nrow()-2; row1++){

      int res;
      // in the Cmatrix (col2res) the value for each non gap residue is 2*res+1
      // => this trasformes it back to the residue number
      // We need this number to look up the position of the residue in the
      // alternative MSA (by using res2col)
      if(ref_col2res(row1,col)%2 == 1){
        // coution: modulo works only with type int
        res=0.5*(ref_col2res(row1,col)+1)-1;
      }
      else{ res = -1;}

      double raise = row1+1;
      for(int row2 = raise; row2 <= ref_col2res.nrow()-1; row2++){

        // if the residue pair in the REF MSA is not a gab
        // and the second residue associated with the first residue is equal
        // => hit and rise the hit matrix (res_pair_hit) by one for this residue pair
        if(res > -1){
          if(hit_mat(count,col)>-1 && (ref_col2res(row2,col) == alt_col2res(row2, alt_res2col(row1, res)))){
            hit_mat(count,col)++;
          }
        }
        count++;
      }
    }
  }
  return(hit_mat);
}
