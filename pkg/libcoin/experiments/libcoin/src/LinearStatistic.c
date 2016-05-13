
#include "libcoin_internal.h"
#include "Sums.h"
#include "Utils.h"
#include "Tables.h"


/* Variables
  x:		a double N x P matrix
  y:		a double N x Q matrix / a double M x Q matrix
  ix:		integer vector of length N with elements 0...(Lx - 1)
  iy:		integer vector of length N with elements 0...(Ly - 1)
  weights:	an integer vector of length N
  subset:       an integer Nsubset vector with elements 0...(N - 1)
  Lb:		number of levels of block
  PQ_ans:	return value, a double P x Q matrix
  P_ans:	return value, a double P vector
  PP_sym_ans:	return value, a symmetric double P x P matrix in lower packed format
  Q_ans:	return value, a double Q vector
  QQ_sym_ans:	return value, a symmetric double Q x Q matrix in lower packed format
  LbQQ_sym_ans:	return values, Lb symmetric double Q x Q matrices in lower packed format
  PQPQ_sym_ans:	return value, a symmetric double PQ x PQ matrix in lower packed format
  ExpInf:	expectation of influence function y, a double Q vector
  CovInf:	covariance of influence function y, a double Q x (Q + 1) / 2 matrix
  ExpX:		"expectation" of x, a double P vector
  CovX:		"covariance" of x, a double P * (P + 1) / 2 matrix
  VarX:		"variance" of x, a double P vector
  PP_tmp:	temp variable, a symmetric double P x P matrix in lower packed format
  P_tmp:	temp variable, a double P vector
  add:		integer; 0 means init return value with 0 and 1 means add to existing values
*/   
   

void C_LinearStatistic_(double *x, int N, int P, double *y, int Q, 
                        double *PQ_ans) 
{
     C_KronSums_(x, N, P, y, Q, PQ_ans);
}

void C_LinearStatistic_weights(double *x, int N, int P, double *y, int Q, 
                               int *weights, double *PQ_ans)
{                               
     C_KronSums_weights(x, N, P, y, Q, weights, PQ_ans);
}
     
void C_LinearStatistic_subset(double *x, int N, int P, double *y, int Q, 
                              int *subset, int Nsubset, double *PQ_ans) 
{
     C_KronSums_subset(x, N, P, y, Q, subset, subset, Nsubset, PQ_ans);
}

void C_LinearStatistic_weights_subset(double *x, int N, int P, double *y, int Q, 
                                      int *weights, int *subset, int Nsubset, 
                                      double *PQ_ans) 
{
     C_KronSums_weights_subset(x, N, P, y, Q, weights, subset, Nsubset, PQ_ans);
}


void C_PermutedLinearStatistic(double *x, int N, int P, double *y, int Q, 
                               int *perm, int *original, int Nperm, 
                               double *PQ_ans) 
{
     C_KronSums_subset(x, N, P, y, Q, perm, original, Nperm, PQ_ans);
}

void C_LinearStatistic_2d(double *x, int N, int P, double *y, int M, int Q, 
                          int *weights, double *PQ_ans) 
{
    C_KronSums_2dweights(x, N, P, y, M, Q, weights, PQ_ans);
}

void C_LinearStatisticXfactor_(int *ix, int N, int P, double *y, int Q, 
                               double *PQ_ans) 
{
    C_tapplySum_(y, N, Q, ix, P, PQ_ans);
}

void C_LinearStatisticXfactor_weights(int *ix, int N, int P, double *y, int Q, 
                                      int *weights, double *PQ_ans)
{                               
    C_tapplySum_weights(y, N, Q, ix, P, weights, PQ_ans);
}
     
void C_LinearStatisticXfactor_subset(int *ix, int N, int P, double *y, int Q, 
                              int *subset, int Nsubset, double *PQ_ans) 
{
    C_tapplySum_subset(y, N, Q, ix, P, subset, Nsubset, PQ_ans);
}

void C_LinearStatisticXfactor_weights_subset(int *ix, int N, int P, double *y, 
                                             int Q, int *weights, int *subset, 
                                             int Nsubset, double *PQ_ans) 
{
    C_tapplySum_weights_subset(y, N, Q, ix, P, weights, subset, Nsubset, 
                               PQ_ans);
}
     
void C_LinearStatisticXfactor_2d(int N, double *y, int M, int Q, 
                                 int *weights, double *Nm1Q_ans) 
{
    C_tapplySum_2d(y, M, Q, N, weights, Nm1Q_ans);
}

void C_LinearStatisticXfactor(int *x, int N, int P, double* y, int Q, 
                              int *weights, int *sumweights,
                              int *subset, int *Nsubset, int Lb, 
                              double *PQ_ans) 
{

    int sw = 0, ns = 0;
    
    for (int b = 0; b < Lb; b++) {
        sw = sw + sumweights[b];
        ns = ns + Nsubset[b];
    }

    if (ns == 0) {
        if (sw == 0) {
              C_LinearStatisticXfactor_(x, N, P, y, Q, PQ_ans);
        } else {
              C_LinearStatisticXfactor_weights(x, N, P, y, Q, weights, PQ_ans);
        }
    } else {
        if (sw == 0) {
            C_LinearStatisticXfactor_subset(x, N, P, y, Q, subset, ns, PQ_ans);
        } else {
            C_LinearStatisticXfactor_weights_subset(x, N, P, y, Q, weights, 
                                                    subset, ns, PQ_ans);
        }
    }
}

void C_LinearStatistic(SEXP x, int N, int P, double* y, int Q, 
                       int *weights, int *sumweights,
                       int *subset, int *Nsubset, int Lb, 
                       double *PQ_ans) 
{

    int sw = 0, ns = 0;

    if (isInteger(x)) {
        C_LinearStatisticXfactor(INTEGER(x), N, P, y, Q, 
                                 weights, sumweights, subset, Nsubset,
                                 Lb, PQ_ans);
    } else {
        for (int b = 0; b < Lb; b++) {
            sw = sw + sumweights[b];
            ns = ns + Nsubset[b];
        }

        if (ns == 0) {
            if (sw == 0) {
                  C_LinearStatistic_(REAL(x), N, P, y, Q, PQ_ans);
            } else {
                  C_LinearStatistic_weights(REAL(x), N, P, y, Q, weights, PQ_ans);
        }
        } else {
            if (sw == 0) {
                C_LinearStatistic_subset(REAL(x), N, P, y, Q, subset, ns, PQ_ans);
            } else {
                C_LinearStatistic_weights_subset(REAL(x), N, P, y, Q, weights, 
                                                 subset, ns, PQ_ans);
            }
        }
    }
}


void C_ExpectationInfluence_(double* y, int N, int Q, double *Q_ans) 
{
     C_colSums_(y, N, Q, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / N;
}

void C_ExpectationInfluence_weights(double* y, int N, int Q, int *weights, 
                                    int sumweights, double *Q_ans) 
{
     C_colSums_weights(y, N, Q, weights, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / sumweights;
}

void C_ExpectationInfluence_subset(double* y, int N, int Q, int *subset, 
                                   int Nsubset, double *Q_ans) 
{
     C_colSums_subset(y, N, Q, subset, Nsubset, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / Nsubset;
}

void C_ExpectationInfluence_weights_subset(double* y, int N, int Q, 
                                           int *weights, int sumweights, 
                                           int *subset, int Nsubset, 
                                           double *Q_ans) 
{
     C_colSums_weights_subset(y, N, Q, weights, subset, Nsubset, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / sumweights;
}

     
void C_CovarianceInfluence_(double* y, int N, int Q, double *ExpInf, 
                            double *QQ_sym_ans) 
{
     C_KronSums_sym_center_(y, N, Q, ExpInf, QQ_sym_ans);
     for (int q = 0; q < Q * (Q + 1) / 2; q++) 
         QQ_sym_ans[q] = QQ_sym_ans[q] / N;
}

void C_CovarianceInfluence_weights(double* y, int N, int Q, int *weights, 
                                   int sumweights, double *ExpInf, 
                                   double *QQ_sym_ans) 
{
     C_KronSums_sym_center_weights(y, N, Q, weights, ExpInf, QQ_sym_ans);
     for (int q = 0; q < Q * (Q + 1) / 2; q++) 
         QQ_sym_ans[q] = QQ_sym_ans[q] / sumweights;
}

void C_CovarianceInfluence_subset(double* y, int N, int Q, int *subset, 
                                  int Nsubset, double *ExpInf, 
                                  double *QQ_sym_ans) 
{
     C_KronSums_sym_center_subset(y, N, Q, subset, Nsubset, 
                                  ExpInf, QQ_sym_ans);
     for (int q = 0; q < Q * (Q + 1) / 2; q++) 
         QQ_sym_ans[q] = QQ_sym_ans[q] / Nsubset;
}

void C_CovarianceInfluence_weights_subset(double* y, int N, int Q, int *weights, 
                                          int sumweights, int *subset, 
                                          int Nsubset, double *ExpInf, 
                                          double *QQ_sym_ans) 
{
     C_KronSums_sym_center_weights_subset(y, N, Q, weights, subset, Nsubset, 
                                          ExpInf, QQ_sym_ans);
     for (int q = 0; q < Q * (Q + 1) / 2; q++) 
         QQ_sym_ans[q] = QQ_sym_ans[q] / sumweights;
}

void C_VarianceInfluence_(double* y, int N, int Q, double *ExpInf, 
                          double *Q_ans) 
{
     C_colSums2_center_(y, N, Q, ExpInf, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / N;
}

void C_VarianceInfluence_weights(double* y, int N, int Q, int *weights, 
                                 int sumweights, double *ExpInf, 
                                 double *Q_ans)
{
     C_colSums2_center_weights(y, N, Q, weights, ExpInf, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / sumweights;
}

void C_VarianceInfluence_subset(double* y, int N, int Q, int *subset, 
                                int Nsubset, double *ExpInf, double *Q_ans) 
{
     C_colSums2_center_subset(y, N, Q, subset, Nsubset, 
                              ExpInf, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / Nsubset;
}

void C_VarianceInfluence_weights_subset(double* y, int N, int Q, int *weights, 
                                        int sumweights, int *subset, 
                                        int Nsubset, double *ExpInf, 
                                        double *Q_ans) 
{
     C_colSums2_center_weights_subset(y, N, Q, weights, subset, Nsubset, 
                                      ExpInf, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / sumweights;
}

void C_ExpectationCovarianceInfluence(double* y, int N, int Q, 
                                      int *weights, int *sumweights, 
                                      int *subset, int *Nsubset, int Lb, 
                                      int varonly, double *LbQ_ans, 
                                      double *LbQQ_sym_ans) 
{
     int ns = 0;
     double *ExpInf, *CovInf, *VarInf;

     for (int b = 0; b < Lb; b++) {
         ExpInf = LbQ_ans + b * Q;
         VarInf = LbQQ_sym_ans + b * Q;
         CovInf = LbQQ_sym_ans + b * Q * (Q + 1) / 2;
         if (Nsubset[b] == 0) {
             if (sumweights[b] == 0) {
                 C_ExpectationInfluence_(y, N, Q, ExpInf);
                 if (varonly) {
                     C_VarianceInfluence_(y, N, Q, ExpInf, VarInf);
                 } else {
                     C_CovarianceInfluence_(y, N, Q, ExpInf, CovInf);
                 }
             } else {
                 C_ExpectationInfluence_weights(y, N, Q, weights, sumweights[b], 
                                                ExpInf);
                 if (varonly) {
                     C_VarianceInfluence_weights(y, N, Q, weights, 
                                                 sumweights[b], ExpInf, 
                                                 VarInf);
                 } else {
                     C_CovarianceInfluence_weights(y, N, Q, weights, 
                                                   sumweights[b], ExpInf, 
                                                   CovInf);
                 }
             }
         } else {
             if (sumweights[b] == 0) {
                 C_ExpectationInfluence_subset(y, N, Q, subset + ns, Nsubset[b], 
                                               ExpInf);
                 if (varonly) {
                     C_VarianceInfluence_subset(y, N, Q, subset + ns, 
                                                Nsubset[b], ExpInf, VarInf);
                 } else {
                     C_CovarianceInfluence_subset(y, N, Q, subset + ns, 
                                                  Nsubset[b], ExpInf, CovInf);
                 }
             } else {
                 C_ExpectationInfluence_weights_subset(y, N, Q, weights, 
                                                       sumweights[b], subset + ns, 
                                                       Nsubset[b], ExpInf);
                 if (varonly) {
                     C_VarianceInfluence_weights_subset(y, N, Q, weights, 
                                                        sumweights[b], subset + ns, 
                                                        Nsubset[b], ExpInf, VarInf);
                 } else {
                     C_CovarianceInfluence_weights_subset(y, N, Q, weights, 
                                                          sumweights[b], 
                                                          subset + ns, 
                                                          Nsubset[b], 
                                                          ExpInf, CovInf);
                 }
             }
             ns = ns + Nsubset[b];
         }
     }
}


void C_ExpectationX_(double* x, int N, int P, double *P_ans) 
{
     C_colSums_(x, N, P, P_ans);
}

void C_ExpectationX_weights(double* x, int N, int P, int *weights, 
                            double *P_ans) 
{
     C_colSums_weights(x, N, P, weights, P_ans);
}

void C_ExpectationX_subset(double* x, int N, int P, int *subset, int Nsubset, 
                           double *P_ans) 
{
     C_colSums_subset(x, N, P, subset, Nsubset, P_ans);
}

void C_ExpectationX_weights_subset(double* x, int N, int P, int *weights, 
                                   int *subset, int Nsubset, 
                                   double *P_ans) 
{
     C_colSums_weights_subset(x, N, P, weights, subset, Nsubset, P_ans);
}

void C_CovarianceX_(double *x, int N, int P, double *PP_sym_ans) 
{
     C_KronSums_sym_(x, N, P, PP_sym_ans);
}
     
void C_CovarianceX_weights(double *x, int N, int P, 
                           int *weights, double *PP_sym_ans) 
{
     C_KronSums_sym_weights(x, N, P, weights, PP_sym_ans);
}     
     
void C_CovarianceX_subset(double *x, int N, int P, int *subset, 
                          int Nsubset, double *PP_sym_ans) 
{
     C_KronSums_sym_subset(x, N, P, subset, Nsubset, PP_sym_ans);
}
     
void C_CovarianceX_weights_subset(double *x, int N, int P, int *weights,
                                  int *subset, int Nsubset, double *PP_sym_ans) 
{
     C_KronSums_sym_weights_subset(x, N, P, weights, subset, Nsubset, 
                                   PP_sym_ans);
}

void C_VarianceX_(double *x, int N, int P, double *P_ans) 
{
     C_colSums2_(x, N, P, P_ans);
}
     
void C_VarianceX_weights(double *x, int N, int P, int *weights, double *P_ans) 
{
     C_colSums2_weights(x, N, P, weights, P_ans);
}     
     
void C_VarianceX_subset(double *x, int N, int P, int *subset, 
                        int Nsubset, double *P_ans) 
{
     C_colSums2_subset(x, N, P, subset, Nsubset, P_ans);
}
     
void C_VarianceX_weights_subset(double *x, int N, int P, int *weights,
                                int *subset, int Nsubset, double *P_ans) 
{
     C_colSums2_weights_subset(x, N, P, weights, subset, Nsubset, P_ans);
}

void C_ExpectationLinearStatistic(int P, int Q, double *ExpInf, double *ExpX, 
                                  int add, double *PQ_ans)
{

    if (!add)
        for (int p = 0; p < P * Q; p++) PQ_ans[p] = 0.0;
        
    for (int p = 0; p < P; p++) {
        for (int q = 0; q < Q; q++)
            PQ_ans[q * P + p] += ExpX[p] * ExpInf[q];
    }
}          

void C_CovarianceLinearStatistic(int P, int Q, double *CovInf, double *ExpX, 
                                 double *CovX, int sumweights, 
                                 double *PP_sym_tmp, int add,
                                 double *PQPQ_sym_ans) 
{
    double f1 = (double) sumweights / (sumweights - 1);
    double f2 = 1.0 / (sumweights - 1);
        
    if (P * Q == 1) {
        PQPQ_sym_ans[0] = f1 * CovInf[0] * CovX[0];
        PQPQ_sym_ans[0] -= f2 * CovInf[0] * ExpX[0] * ExpX[0];
    } else {
        C_KronSums_sym_(ExpX, 1, P, PP_sym_tmp);
        for (int p = 0; p < P * (P + 1) / 2; p++)
            PP_sym_tmp[p] = f1 * CovX[p] - f2 * PP_sym_tmp[p];
        C_kronecker_sym(CovInf, Q, PP_sym_tmp, P, 1 - add, PQPQ_sym_ans);
    }
}

void C_VarianceLinearStatistic(int P, int Q, double *VarInf, double *ExpX, 
                               double *VarX, int sumweights, double *P_tmp, 
                               int add, double *PQ_ans) 
{
    if (P * Q == 1) {
        C_CovarianceLinearStatistic(P, Q, VarInf, ExpX, VarX, 
                                    sumweights, P_tmp, add, PQ_ans);
    } else {

        double f1 = (double) sumweights / (sumweights - 1);
        double f2 = 1.0 / (sumweights - 1);
        for (int p = 0; p < P; p++)
            P_tmp[p] = f1 * VarX[p] - f2 * ExpX[p] * ExpX[p];
        C_kronecker(VarInf, 1, Q, P_tmp, 1, P, 1 - add, PQ_ans);
    }
}

void C_ExpectationCovarianceLinearStatistic(SEXP x, int N, int P, int Q,
                                            int *weights, int *sumweights, 
                                            int *subset, int *Nsubset, int Lb, 
                                            double *ExpInf, double *CovInf, 
                                            double *work, double *PQ_ans,  
                                            double *PQPQ_sym_ans) 
{
     int bQ, ns = 0, PQ = P * Q, sw = 0, *ix, *subtmp;
     double *ExpX, *CovX, *PPtmp;


     for (int i = 0; i < P + 2 * P * (P + 1) / 2 + 1; i++) work[i] = 0.0;

     /* work[0] counts NAs in ix (ix[i] == 0)  */
     ExpX = work + 1;
     CovX = ExpX + P;
     PPtmp = CovX + P * (P + 1) / 2;

     for (int b = 0; b < Lb; b++) {
         bQ = b * PQ * (PQ + 1) / 2;
         if (Nsubset[b] == 0) {
             if (sumweights[b] == 0) {
                 if (isInteger(x)) {
                     ix = INTEGER(x);
                     /* work[0] counts NAs */
                     for (int i = 0; i < N; i++) work[ix[i]]++; 
                     /* CovX = diag(ExpX) */
                     for (int p = 0; p < P; p++)
                         CovX[S(p, p, P)] = work[p];
                 } else {
                     C_ExpectationX_(REAL(x), N, P, ExpX);
                     C_CovarianceX_(REAL(x), N, P, CovX);
                 }
                 sw = N;
             } else {
                 if (isInteger(x)) {
                     ix = INTEGER(x);
                     for (int i = 0; i < N; i++) 
                         work[ix[i]] += (double) weights[ix[i]];
                     for (int p = 0; p < P; p++)
                         CovX[S(p, p, P)] = work[p];
                 } else {
                     C_ExpectationX_weights(REAL(x), N, P, weights, ExpX);
                     C_CovarianceX_weights(REAL(x), N, P, weights, CovX);
                 }
                 sw = sumweights[b];
             }
         } else {
             if (sumweights[b] == 0) {
                 if (isInteger(x)) {
                 ix = INTEGER(x);
                 subtmp = subset + ns;
                 for (int i = 0; i < Nsubset[b]; i++) 
                     work[ix[subtmp[i]]]++; 
                 for (int p = 0; p < P; p++)
                      CovX[S(p, p, P)] = work[p];
                 } else {
                     C_ExpectationX_subset(REAL(x), N, P, subset + ns, 
                                           Nsubset[b], ExpX);
                     C_CovarianceX_subset(REAL(x), N, P, subset + ns, 
                                          Nsubset[b], CovX);
                 }
                 sw = Nsubset[b];
             } else {
                 if (isInteger(x)) {
                     ix = INTEGER(x);
                     subtmp = subset + ns;
                     for (int i = 0; i < Nsubset[b]; i++) 
                         work[ix[subtmp[i]]] += (double) weights[subtmp[i]]; 
                     for (int p = 0; p < P; p++)
                          CovX[S(p, p, P)] = work[p];
                 } else {
                     C_ExpectationX_weights_subset(REAL(x), N, P, weights, 
                                                   subset + ns, Nsubset[b], 
                                                   ExpX);
                     C_CovarianceX_weights_subset(REAL(x), N, P, weights, 
                                                  subset + ns, Nsubset[b], 
                                                  CovX);
                 }
                 sw = sumweights[b];
             }
         }
         C_ExpectationLinearStatistic(P, Q, ExpInf + b * Q, ExpX, b, PQ_ans);
         C_CovarianceLinearStatistic(P, Q, CovInf + b * Q * (Q + 1) / 2, 
                                     ExpX, CovX, sw, PPtmp, b, PQPQ_sym_ans);
         ns = ns + Nsubset[b];
     }
}

void C_ExpectationVarianceLinearStatistic(SEXP x, int N, int P, int Q, 
                                          int *weights, int *sumweights, 
                                          int *subset, int *Nsubset, int Lb, 
                                          double *ExpInf, double *VarInf, 
                                          double *work, double *PQ_ans_Exp, 
                                          double *PQ_ans_Var) 
{
     int bQ, ns = 0, PQ = P * Q, sw = 0, *ix, *subtmp;
     double *ExpX, *VarX, *PPtmp;

     for (int i = 0; i < 3 * P + 1; i++) work[i] = 0.0;

     /* work[0] counts NAs in ix (ix[i] = 0)*/
     ExpX = work + 1;
     VarX = ExpX + P;
     PPtmp = VarX + P;
     
     for (int b = 0; b < Lb; b++) {
         bQ = b * PQ;
         if (Nsubset[b] == 0) {
             if (sumweights[b] == 0) {
                 if (isInteger(x)) {
                     ix = INTEGER(x);
                     /* work[0] counts NAs */
                     for (int i = 0; i < N; i++) work[ix[i]]++;
                     VarX = ExpX;
                 } else {
                     C_ExpectationX_(REAL(x), N, P, ExpX);
                     C_VarianceX_(REAL(x), N, P, VarX);
                 }
                 sw = N;
             } else {
                 if (isInteger(x)) {
                     ix = INTEGER(x);
                     for (int i = 0; i < N; i++) 
                         work[ix[i]] += (double) weights[ix[i]]; 
                     VarX = ExpX;
                 } else {
                     C_ExpectationX_weights(REAL(x), N, P, weights, ExpX);
                     C_VarianceX_weights(REAL(x), N, P, weights, VarX);
                 }
                 sw = sumweights[b];
             }
         } else {
             if (sumweights[b] == 0) {
                 if (isInteger(x)) {
                     ix = INTEGER(x);
                     subtmp = subset + ns;
                     for (int i = 0; i < Nsubset[b]; i++) 
                         work[ix[subtmp[i]]]++; 
                     VarX = ExpX;
                 } else {
                     C_ExpectationX_subset(REAL(x), N, P, subset + ns, 
                                           Nsubset[b], ExpX);
                     C_VarianceX_subset(REAL(x), N, P, subset + ns, 
                                        Nsubset[b], VarX);
                 }
                 sw = Nsubset[b];
             } else {
                 if (isInteger(x)) {
                     ix = INTEGER(x);
                     subtmp = subset + ns;
                     for (int i = 0; i < Nsubset[b]; i++) 
                         work[ix[subtmp[i]]] += (double) weights[subtmp[i]]; 
                     VarX = ExpX;
                 } else {
                     C_ExpectationX_weights_subset(REAL(x), N, P, weights, 
                                                   subset + ns, Nsubset[b], 
                                                   ExpX);
                     C_VarianceX_weights_subset(REAL(x), N, P, weights, 
                                                subset + ns, Nsubset[b], 
                                                VarX);
                 }
                 sw = sumweights[b];
             }
         }
         C_ExpectationLinearStatistic(P, Q, ExpInf + b * Q, ExpX, b, 
                                      PQ_ans_Exp);
         C_VarianceLinearStatistic(P, Q, VarInf + b * Q, ExpX, VarX, 
                                   sw, PPtmp, b, PQ_ans_Var);
         ns = ns + Nsubset[b];
     }
}
