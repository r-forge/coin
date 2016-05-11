
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
