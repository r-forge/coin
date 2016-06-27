
/* Variables:
  ix:		integer vector of length N with elements 0...(Lx - 1)
  iy:		integer vector of length N with elements 0...(Ly - 1)
  weights:	integer vector of length N
  subset:       an integer Nsubset vector with elements 0...(N - 1)
  block:	an integer N vector with elements 1...Lb
  LxLyLb_ans:	return value, integer array Lx x Ly x Lb 
*/

void RC_2dtable(SEXP ix, SEXP iy, SEXP weights, SEXP subset, SEXP block, 
                int *LxLyLb_ans); 

/* table(ix) */
void C_1dtable_(int *ix, int Lx, int N, int *Lx_ans);
/* table(ix[subset]) */
void C_1dtable_subset(int *ix, int Lx, int *subset, int Nsubset, int *Lx_ans);
/* xtabs(weights ~ ix) */
void C_1dtable_weights(int *ix, int Lx, int *weights, int N, int *Lx_ans);
/* xtabs(weights ~ ix, subset = subset) */
void C_1dtable_weights_subset(int *ix, int Lx, int *weights, int *subset, 
                              int Nsubset, int *Lx_ans);
