#include <Rcpp.h>
using namespace Rcpp;
//' @title Calculate C-index (consistency index)
//' @description Calculate C-index (consistency index)
//' @param data_1 data in matrix form
//' @param beta_6 the estimator that has been obtained from AFT model
//' @return Returns a number between 0-1
//' @examples
//' \dontrun{
//'     result <- computeCindex(as.numeric(beta_0), as.matrix(data_r))
//' }
//' @export
 // [[Rcpp::export]]
 double computeCindex(NumericVector beta_6, NumericMatrix data_1) {
   int n = data_1.nrow();
   int p = data_1.ncol() - 2; // Exclude "Y" and "status" columns
   
   NumericVector Y_time = log(data_1(_, 0));
   NumericMatrix data_x = data_1(_, Range(1, p));
   NumericVector status_in = data_1(_, p + 1);
   
   NumericMatrix W(n, n);
   NumericMatrix V(n, n);
   NumericMatrix U(n, n);
   NumericMatrix S(n, n);
   
   for (int i = 0; i < n; ++i) {
     for (int j = i + 1; j < n; ++j) {
       W(i, j) = (Y_time[i] - Y_time[j] > 0) ? 1 : -1;
       NumericVector diff_x = data_x(i, _) - data_x(j, _);
       V(i, j) = sum(diff_x * beta_6);
       V(i, j) = (V(i, j) > 0) ? 1 : -1;
       U(i, j) = V(i, j) - W(i, j);
       S(i, j) = status_in[i] - status_in[j];
     }
   }
   
   double C1 = 0;
   double C2 = 0;
   
   for (int i = 0; i < n; ++i) {
     for (int j = i + 1; j < n; ++j) {
       if (i == n - 1) {
         break;
       }
       if (S(i, j) == 0 && status_in[i] == 1 && U(i, j) == 0) {
         C1++;
       } else if (S(i, j) == 0 && status_in[i] == 1 && U(i, j) != 0) {
         C2++;
       } else if (S(i, j) != 0 && status_in[i] == 1 && W(i, j) < 0 && U(i, j) == 0) {
         C1++;
       } else if (S(i, j) != 0 && status_in[i] == 1 && W(i, j) < 0 && U(i, j) != 0) {
         C2++;
       } else if (S(i, j) != 0 && status_in[i] == 0 && W(i, j) > 0 && U(i, j) == 0) {
         C1++;
       } else if (S(i, j) != 0 && status_in[i] == 0 && W(i, j) > 0 && U(i, j) != 0) {
         C2++;
       }
     }
   }
   
   double res = C1 / (C1 + C2);
   return res;
 }



