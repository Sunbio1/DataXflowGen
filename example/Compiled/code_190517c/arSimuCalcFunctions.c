#include "neu_model_900BC41043F09DAAEAD3703013836320.h"
#include "neu_data_15A1F055C1DB84993186F02596F89DDC.h"
#include "neu_data_246828910E5EBCA41A80196ADDB1F18A.h"
#include "neu_data_BE42CE855E5AD5AEABA6DBC01469C272.h"
#include "neu_data_1BD4D844BEBA067E729747DD9F896963.h"
#include "neu_data_2C6C200814D6F45EB0E9F79EC5041D83.h"
#include "neu_data_FAF3A73978E04AC9C0B11BF8EA78E748.h"
#include "neu_data_765854CC8C58746184D7F2BA78817CA4.h"

 int AR_CVodeInit(void *cvode_mem, N_Vector x, double t, int im, int ic){
  if((im==0) & (ic==0)) return CVodeInit(cvode_mem, fx_neu_model_900BC41043F09DAAEAD3703013836320, RCONST(t), x);
  return(-1);
}

 void fx(realtype t, N_Vector x, double *xdot, void *user_data, int im, int ic){
  if((im==0) & (ic==0)) fxdouble_neu_model_900BC41043F09DAAEAD3703013836320(t, x, xdot, user_data);
}

 void fx0(N_Vector x0, void *user_data, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fx0_neu_model_900BC41043F09DAAEAD3703013836320(x0, data);
}

 int AR_CVDlsSetDenseJacFn(void *cvode_mem, int im, int ic, int setSparse){
  if((im==0) & (ic==0) & (setSparse==0)){ 
 return CVDlsSetDenseJacFn(cvode_mem, dfxdx_neu_model_900BC41043F09DAAEAD3703013836320);
 
 }else if((im==0) & (ic==0) & (setSparse==1)){ 
 return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_neu_model_900BC41043F09DAAEAD3703013836320);

}
  return(-1);
}

 void getdfxdx(int im, int ic, realtype t, N_Vector x, realtype *J, void *user_data){
  if((im==0) & (ic==0)) { dfxdx_out_neu_model_900BC41043F09DAAEAD3703013836320(t, x, J, user_data); return; }

}

 void fsx0(int is, N_Vector sx_is, void *user_data, int im, int ic, int sensitivitySubset) {
  UserData data = (UserData) user_data;
  if ( sensitivitySubset == 0 ) {
    if((im==0) & (ic==0)) fsx0_neu_model_900BC41043F09DAAEAD3703013836320(is, sx_is, data);
  } else {
    if((im==0) & (ic==0)) subfsx0_neu_model_900BC41043F09DAAEAD3703013836320(is, sx_is, data);
  }
}

 void csv(realtype t, N_Vector x, int ip, N_Vector sx, void *user_data, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) csv_neu_model_900BC41043F09DAAEAD3703013836320(t, x, ip, sx, data);
}

 int AR_CVodeSensInit1(void *cvode_mem, int nps, int sensi_meth, int sensirhs, N_Vector *sx, int im, int ic, int sensitivitySubset){
  if (sensirhs == 1) {
    if (sensitivitySubset == 0) {
      if((im==0) & (ic==0)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_neu_model_900BC41043F09DAAEAD3703013836320, sx);
    } else {
      if((im==0) & (ic==0)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_neu_model_900BC41043F09DAAEAD3703013836320, sx);
  }
  } else {
    if((im==0) & (ic==0)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, NULL, sx);
  }
  return(-1);
}

 void fu(void *user_data, double t, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fu_neu_model_900BC41043F09DAAEAD3703013836320(data, t);
}

 void fsu(void *user_data, double t, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fsu_neu_model_900BC41043F09DAAEAD3703013836320(data, t);
}

 void fv(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) fv_neu_model_900BC41043F09DAAEAD3703013836320(t, x, data);
}

 void fsv(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) {
	dvdp_neu_model_900BC41043F09DAAEAD3703013836320(t, x, data);
	dvdu_neu_model_900BC41043F09DAAEAD3703013836320(t, x, data);
	dvdx_neu_model_900BC41043F09DAAEAD3703013836320(t, x, data);
}
}

 void dfxdp0(void *user_data, double t, N_Vector x, double *dfxdp0, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) dfxdp0_neu_model_900BC41043F09DAAEAD3703013836320(t, x, dfxdp0, data);
}

 void dfxdp(void *user_data, double t, N_Vector x, double *dfxdp, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) dfxdp_neu_model_900BC41043F09DAAEAD3703013836320(t, x, dfxdp, data);
}

void fz(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *z, double *p, double *u, double *x, int im, int ic){
  if((im==0) & (ic==0)) fz_neu_model_900BC41043F09DAAEAD3703013836320(t, nt, it, nz, nx, nu, iruns, z, p, u, x);
}

void dfzdx(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *dfzdx, double *z, double *p, double *u, double *x, int im, int ic){
  if((im==0) & (ic==0)) dfzdx_neu_model_900BC41043F09DAAEAD3703013836320(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);
}

void fsz(double t, int nt, int it, int np, double *sz, double *p, double *u, double *x, double *z, double *su, double *sx, int im, int ic){
  if((im==0) & (ic==0)) fsz_neu_model_900BC41043F09DAAEAD3703013836320(t, nt, it, np, sz, p, u, x, z, su, sx);
}

 void fy(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y, double *p, double *u, double *x, double *z, int im, int id){
  if((im==0) & (id==0)) fy_neu_data_15A1F055C1DB84993186F02596F89DDC(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==1)) fy_neu_data_246828910E5EBCA41A80196ADDB1F18A(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==2)) fy_neu_data_BE42CE855E5AD5AEABA6DBC01469C272(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==3)) fy_neu_data_1BD4D844BEBA067E729747DD9F896963(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==4)) fy_neu_data_2C6C200814D6F45EB0E9F79EC5041D83(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==5)) fy_neu_data_FAF3A73978E04AC9C0B11BF8EA78E748(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
  if((im==0) & (id==6)) fy_neu_data_765854CC8C58746184D7F2BA78817CA4(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);
}

 void fy_scale(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y_scale, double *p, double *u, double *x, double *z, double *dfzdx, int im, int id){
  if((im==0) & (id==0)) fy_scale_neu_data_15A1F055C1DB84993186F02596F89DDC(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==1)) fy_scale_neu_data_246828910E5EBCA41A80196ADDB1F18A(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==2)) fy_scale_neu_data_BE42CE855E5AD5AEABA6DBC01469C272(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==3)) fy_scale_neu_data_1BD4D844BEBA067E729747DD9F896963(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==4)) fy_scale_neu_data_2C6C200814D6F45EB0E9F79EC5041D83(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==5)) fy_scale_neu_data_FAF3A73978E04AC9C0B11BF8EA78E748(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
  if((im==0) & (id==6)) fy_scale_neu_data_765854CC8C58746184D7F2BA78817CA4(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);
}

 void fystd(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x, double *z, int im, int id){
  if((im==0) & (id==0)) fystd_neu_data_15A1F055C1DB84993186F02596F89DDC(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==1)) fystd_neu_data_246828910E5EBCA41A80196ADDB1F18A(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==2)) fystd_neu_data_BE42CE855E5AD5AEABA6DBC01469C272(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==3)) fystd_neu_data_1BD4D844BEBA067E729747DD9F896963(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==4)) fystd_neu_data_2C6C200814D6F45EB0E9F79EC5041D83(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==5)) fystd_neu_data_FAF3A73978E04AC9C0B11BF8EA78E748(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
  if((im==0) & (id==6)) fystd_neu_data_765854CC8C58746184D7F2BA78817CA4(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);
}

 void fsy(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *z, double *su, double *sx, double *sz, int im, int id){
  if((im==0) & (id==0)) fsy_neu_data_15A1F055C1DB84993186F02596F89DDC(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==1)) fsy_neu_data_246828910E5EBCA41A80196ADDB1F18A(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==2)) fsy_neu_data_BE42CE855E5AD5AEABA6DBC01469C272(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==3)) fsy_neu_data_1BD4D844BEBA067E729747DD9F896963(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==4)) fsy_neu_data_2C6C200814D6F45EB0E9F79EC5041D83(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==5)) fsy_neu_data_FAF3A73978E04AC9C0B11BF8EA78E748(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
  if((im==0) & (id==6)) fsy_neu_data_765854CC8C58746184D7F2BA78817CA4(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);
}

 void fsystd(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *z, double *sy, double *su, double *sx, double *sz, int im, int id){
  if((im==0) & (id==0)) fsystd_neu_data_15A1F055C1DB84993186F02596F89DDC(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==1)) fsystd_neu_data_246828910E5EBCA41A80196ADDB1F18A(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==2)) fsystd_neu_data_BE42CE855E5AD5AEABA6DBC01469C272(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==3)) fsystd_neu_data_1BD4D844BEBA067E729747DD9F896963(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==4)) fsystd_neu_data_2C6C200814D6F45EB0E9F79EC5041D83(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==5)) fsystd_neu_data_FAF3A73978E04AC9C0B11BF8EA78E748(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
  if((im==0) & (id==6)) fsystd_neu_data_765854CC8C58746184D7F2BA78817CA4(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);
}

/* for arSSACalc.c */

 void fvSSA(void *user_data, double t, N_Vector x, int im, int ic){
  UserData data = (UserData) user_data;
  if((im==0) & (ic==0)) {
    fu_neu_model_900BC41043F09DAAEAD3703013836320(data, t);
    fv_neu_model_900BC41043F09DAAEAD3703013836320(t, x, data);
  }
}

