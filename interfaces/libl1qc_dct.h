
typedef int l1c_int;

typedef struct LBResult{
  double l1;
  int    total_newton_iter;
  int    total_cg_iter;
  int    status;

}LBResult;

typedef struct L1qcDctOpts{
  double epsilon;
  double mu;
  double lbtol;
  double tau;
  double lbiter;
  double newton_tol;
  int newton_max_iter;
  int verbose;
  double l1_tol;
  double cgtol;
  int cgmaxiter;
  int warm_start_cg;
}L1qcDctOpts;

int l1qc_dct(int N, double *x_out, int M, double *b, l1c_int *pix_idx,
             L1qcDctOpts opts, LBResult *lb_res);

int l1qc_dct2(int Nrow, int Ncol, double *x_out, int M, double *b, l1c_int *pix_idx,
              L1qcDctOpts opts, LBResult *lb_res);
