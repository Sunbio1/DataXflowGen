#ifndef _MY_neu_model_900BC41043F09DAAEAD3703013836320
#define _MY_neu_model_900BC41043F09DAAEAD3703013836320

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_sparse.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <cvodes/cvodes_klu.h>
#include <sundials/sundials_sparse.h>
#include <udata.h>
#include <math.h>
#include <mex.h>
#include <arInputFunctionsC.h>



 void fu_neu_model_900BC41043F09DAAEAD3703013836320(void *user_data, double t);
 void fsu_neu_model_900BC41043F09DAAEAD3703013836320(void *user_data, double t);
 void fv_neu_model_900BC41043F09DAAEAD3703013836320(realtype t, N_Vector x, void *user_data);
 void dvdx_neu_model_900BC41043F09DAAEAD3703013836320(realtype t, N_Vector x, void *user_data);
 void dvdu_neu_model_900BC41043F09DAAEAD3703013836320(realtype t, N_Vector x, void *user_data);
 void dvdp_neu_model_900BC41043F09DAAEAD3703013836320(realtype t, N_Vector x, void *user_data);
 int fx_neu_model_900BC41043F09DAAEAD3703013836320(realtype t, N_Vector x, N_Vector xdot, void *user_data);
 void fxdouble_neu_model_900BC41043F09DAAEAD3703013836320(realtype t, N_Vector x, double *xdot_tmp, void *user_data);
 void fx0_neu_model_900BC41043F09DAAEAD3703013836320(N_Vector x0, void *user_data);
 int dfxdx_neu_model_900BC41043F09DAAEAD3703013836320(long int N, realtype t, N_Vector x,N_Vector fx, DlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
 int dfxdx_out_neu_model_900BC41043F09DAAEAD3703013836320(realtype t, N_Vector x, realtype* J, void *user_data); int dfxdx_sparse_neu_model_900BC41043F09DAAEAD3703013836320(realtype t, N_Vector x,N_Vector fx, SlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
 int fsx_neu_model_900BC41043F09DAAEAD3703013836320(int Ns, realtype t, N_Vector x, N_Vector xdot,int ip, N_Vector sx, N_Vector sxdot, void *user_data,N_Vector tmp1, N_Vector tmp2);
 int subfsx_neu_model_900BC41043F09DAAEAD3703013836320(int Ns, realtype t, N_Vector x, N_Vector xdot,int ip, N_Vector sx, N_Vector sxdot, void *user_data,N_Vector tmp1, N_Vector tmp2);
 void fsx0_neu_model_900BC41043F09DAAEAD3703013836320(int ip, N_Vector sx0, void *user_data);
 void subfsx0_neu_model_900BC41043F09DAAEAD3703013836320(int ip, N_Vector sx0, void *user_data);
 void csv_neu_model_900BC41043F09DAAEAD3703013836320(realtype t, N_Vector x, int ip, N_Vector sx, void *user_data);
 void dfxdp0_neu_model_900BC41043F09DAAEAD3703013836320(realtype t, N_Vector x, double *dfxdp0, void *user_data);

 void dfxdp_neu_model_900BC41043F09DAAEAD3703013836320(realtype t, N_Vector x, double *dfxdp, void *user_data);

 void fz_neu_model_900BC41043F09DAAEAD3703013836320(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *z, double *p, double *u, double *x);
 void fsz_neu_model_900BC41043F09DAAEAD3703013836320(double t, int nt, int it, int np, double *sz, double *p, double *u, double *x, double *z, double *su, double *sx);

 void dfzdx_neu_model_900BC41043F09DAAEAD3703013836320(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *dfzdxs, double *z, double *p, double *u, double *x);
#endif /* _MY_neu_model_900BC41043F09DAAEAD3703013836320 */



