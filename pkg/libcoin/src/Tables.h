
/* Variables:
  ix:		integer vector of length N with elements 0...(Lx - 1)
  iy:		integer vector of length N with elements 0...(Ly - 1)
  weights:	integer vector of length N
  subset:       an integer Nsubset vector with elements 0...(N - 1)
  block:	an integer N vector with elements 1...Lb
  LxLyLb_ans:	return value, integer array Lx x Ly x Lb 
*/

extern void RC_2dtable
(
    const SEXP ix, 
    const SEXP iy, 
    const SEXP weights, 
    const SEXP subset, 
    const SEXP block, 
    int *LxLyLb_ans
); 

/* table(ix) */
extern void C_1dtable_
(
    const int *ix, 
    const int Lx, 
    const int N, 
    int *Lx_ans
);
    
/* table(ix[subset]) */
extern void C_1dtable_subset
(
    const int *ix, 
    const int Lx, 
    const int *subset, 
    const int Nsubset, 
    int *Lx_ans
);

/* xtabs(weights ~ ix) */
extern void C_1dtable_weights
(
    const int *ix, 
    const int Lx, 
    const int *weights, 
    const int N, 
    int *Lx_ans
);

/* xtabs(weights ~ ix, subset = subset) */
extern void C_1dtable_weights_subset
(
    const int *ix, 
    const int Lx, 
    const int *weights, 
    const int *subset, 
    const int Nsubset, 
    int *Lx_ans
);
