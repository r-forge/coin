
/**
    Exact Distribution of Two-Sample Permutation Tests
    van de Wiel split-up Algorithm

    Author: Mark van de Wiel (2001-2005) <m.a.v.d.wiel@TUE.nl>
            with modifications for R by Torsten Hothorn <Torsten.Hothorn@R-project.org>
      
    *\file $RCSfile$
    *\author $Author$
    *\date $Date$
*/
                    
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
                    

/**
    The probability distribution for the independent two sample problem.
    
    REFERENCE
    
    M.A. van de Wiel (2001). The split-up algorithm: a fast symbolic method 
    for computing p-values of rank statistics, Computational Statistics, 16, 519-538.
    
*/
                                        
typedef struct {
    long length;
    double *c;
    double *x;
} celW;

double binomi(int m, int n) { 

    double bin = 1;
    double bin1 = 1;
    double bin2 = 1;
    int i, j;
        
    for (i = 1; i <= n; i++) bin1 = bin1 * (m + 1 -i);
    for (j = 1; j <= n; j++) bin2 = bin2 * j;
    bin = bin1/bin2;
    
    return(bin);
}

celW** reserveerW(int a, int b)
{
    long res = 0;
    int i, j;
    celW** W;
    
    W = Calloc(a + 1, celW*);

    for (i = 0; i <= a; i++)
        W[i] = Calloc(b + 1, celW);
        
    for (i = 0; i <= a; i++) {
        for (j = i; j <= b; j++) {
            res = (long) binomi(j,i);
            W[i][j].c = Calloc(res, double);
            W[i][j].x = Calloc(res, double);
        }
    }
    return(W);
}

void initW(int a, int b, celW **W) {

    int i, j;

    for (i = 1; i <= a; i++)
    for (j = 0; j <= b; j++) {
        W[i][j].length = 0;
    }
    for (j = 0; j <= b; j++) {
        W[0][j].length = 1;
        W[0][j].c[0] = 1;  
        W[0][j].x[0] = 0;  
    }
}

void mult(celW *tem, int a, int b, int rank, double *rs) {

    int j;
    for (j = 0; j < tem[0].length; j++)
        tem[0].x[j] += rs[rank];
}

void plus(celW **W, celW *tempie, int a, int b) {

    int elep = 0;
    int k = 0;
    int test = 1;
    int i, j;
    
    for (i = 0; i < W[a][b-1].length; i++) {

        test = 1;
        
        for (j = elep; j < tempie[0].length && test==1; j++) {
        
            if (tempie[0].x[j] - 0.000001 <= W[a][b-1].x[i]
                && W[a][b-1].x[i] <= tempie[0].x[j] + 0.000001) {

                tempie[0].c[j] += W[a][b-1].c[i];
                test = 0;
                elep = j;             
            }
        }
         
        if (test == 1) {
            tempie[0].c[tempie[0].length + k] = W[a][b-1].c[i];
            tempie[0].x[tempie[0].length + k] = W[a][b-1].x[i];
            k++;
        }
    }
    tempie[0].length += k;
}

void mergesort(celW temptw, long tijd)
{
    celW copiep;
    int t1 = 0; 
    int t2 = 0;
    int i, j;

    copiep.c = Calloc(temptw.length, double);
    copiep.x = Calloc(temptw.length, double);
        
    for (i = 0; i < temptw.length; i++) {
        copiep.c[i] = temptw.c[i];
        copiep.x[i] = temptw.x[i];
    }
    
    for (j = 0; j < temptw.length; j++) {
        if (t1 <= tijd-1 && t2 <= temptw.length - tijd - 1) {
            if (copiep.x[t1] < copiep.x[tijd + t2]) {
                temptw.x[j] = copiep.x[t1];
                temptw.c[j] = copiep.c[t1];
                t1++;
            } else {
                temptw.x[j] = copiep.x[tijd + t2];
                temptw.c[j] = copiep.c[tijd + t2];
                t2++;
            }
        } else {
            if (t1 > tijd - 1) {
                temptw.x[j] = copiep.x[tijd + t2];
                temptw.c[j] = copiep.c[tijd + t2];
                t2++; 
            } else {   
                temptw.x[j] = copiep.x[t1];
                temptw.c[j] = copiep.c[t1];
                t1++;
            }
        }          
    } 
    Free(copiep.c);
    Free(copiep.x);
}

void vulcel(celW **W, int i1, int j1, int r, double *rs) {
    
    long tijd;
    celW temp2;
    int j, j2;

    temp2.c = Calloc(W[i1 - 1][j1 - 1].length + W[i1][j1 - 1].length, double);
    temp2.x = Calloc(W[i1 - 1][j1 - 1].length + W[i1][j1 - 1].length, double);
    temp2.length = W[i1 - 1][j1 - 1].length;

    for (j = 0; j < temp2.length; j++) {
       temp2.c[j] = W[i1 - 1][j1 - 1].c[j];
       temp2.x[j] = W[i1 - 1][j1 - 1].x[j];
    }

    if (i1 == j1) {       
        mult(&temp2, i1 - 1, j1 - 1, r, rs); 
    } else {           
        mult(&temp2, i1 - 1, j1 - 1, r, rs);                        
        tijd = temp2.length;                                
        plus(W, &temp2, i1, j1);                            
        mergesort(temp2, tijd);                              
    }

    W[i1][j1].length = temp2.length;

    for (j2 = 0; j2 < temp2.length; j2++) {
        W[i1][j1].c[j2] = temp2.c[j2];
        W[i1][j1].x[j2] = temp2.x[j2];
    }          

    Free(temp2.c);
    Free(temp2.x);
}

void spiegelW(celW **W,int ce, int bep, int start, double *rs) {   

    double totsum = 0;
    long len;
    int r, h;
    
    for (r = 0; r < bep; r++) totsum += rs[start + r];
    
    len = W[bep-ce][bep].length;
        
    for (h = 0; h < len; h++) {
        W[ce][bep].length = W[bep-ce][bep].length;
        W[ce][bep].c[len-1-h] = W[bep-ce][bep].c[h];
        W[ce][bep].x[len-1-h] = totsum - W[bep-ce][bep].x[h];
    }
}

void maakW(celW **W, int a, int b, int start, double *rs) {

    long i,j;
    int rank;
    int hulp;

    for (j = 1; j <= b; j++) {  /* verander naar 0!! */

        if (j < a) {
            hulp = j; 
        } else {
            hulp = a;
        }

        for (i=1; i <= hulp; i++) {   
            if (i <= j/2 || j == 1) {
                rank = start+j;
                vulcel(W, i, j, rank - 1, rs);
            } else {
                spiegelW(W, i, j, start, rs);
            }                               
        }
    }
}

void cumulcoef(celW **W, int i1, int j1) {

     double coef = 0;
     int i;
     
     for(i = 0; i < W[i1][j1].length; i++) {
         W[i1][j1].c[i] += coef;
         coef = W[i1][j1].c[i];
     }
}

double aantalklein(int c, int b, double ob, celW **W1, celW **W2) {

    double tot = 0;
    long le;
    int test = 1;
    int be = b/2;
    int bp = (b+1)/2;
    int tempel = 0;
    int i, j, h;
    
    for (h = 0; h <= c; h++) {

        tempel = 0;
        le = W2[c-h][bp].length;
        
        for (i = 0; i < W1[h][be].length; i++) {
            test = 1;
            for (j = tempel; j < le && test == 1; j++) {
                if (W1[h][be].x[i] + W2[c-h][bp].x[le-j-1] <= ob) {
                    tot += W1[h][be].c[i] * W2[c - h][bp].c[le - j -1];
                    tempel = j;
                    test = 0;
                }
            }
        }
    }
    return(tot);
}

       
SEXP R_split_up_2sample(SEXP scores, SEXP m, SEXP obs) {

    int b, c, d, u;
    double tot, bino, prob;
    double ob;  
    SEXP ans;

    celW **W1;
    celW **W2;
    double *rs;

    b = LENGTH(scores);
    rs = REAL(scores);
    c = INTEGER(m)[0];
    d = b - INTEGER(m)[0];
    ob = REAL(obs)[0];

    /* allocate and initialise memory */
    W1 = reserveerW(c, (b+1)/2);
    initW(c, (b+1)/2, W1);
    W2 = reserveerW(c, (b+1)/2);
    initW(c, (b+1)/2, W2);
    
    maakW(W1, c, b/2, 0, rs);  
    maakW(W2, c, (b+1)/2, b/2, rs);

    for (u = 0; u <= c; u++) cumulcoef(W2, u, (b+1)/2);

    /* number of permutations <= ob */
    tot = aantalklein(c, b, ob, W1, W2);
    
    /* total number of possible permutations */
    bino = binomi(b, c);
    
    /* probability */
    prob = tot/bino; 
    
    Free(W1);
    Free(W2);

    /* return to R */
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0] = prob;
    UNPROTECT(1);
    return(ans);
}
