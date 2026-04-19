#include "neu_data_FAF3A73978E04AC9C0B11BF8EA78E748.h"
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_sparse.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_sparse.h>
#include <cvodes/cvodes_klu.h>
#include <udata.h>
#include <math.h>
#include <mex.h>
#include <arInputFunctionsC.h>





 void fy_neu_data_FAF3A73978E04AC9C0B11BF8EA78E748(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y, double *p, double *u, double *x, double *z){
  y[ny*nt*iruns+it+nt*0] = x[nx*ntlink*iruns+itlink+ntlink*0];
  y[ny*nt*iruns+it+nt*1] = x[nx*ntlink*iruns+itlink+ntlink*1];
  y[ny*nt*iruns+it+nt*2] = x[nx*ntlink*iruns+itlink+ntlink*2];
  y[ny*nt*iruns+it+nt*3] = x[nx*ntlink*iruns+itlink+ntlink*3];
  y[ny*nt*iruns+it+nt*4] = x[nx*ntlink*iruns+itlink+ntlink*4];
  y[ny*nt*iruns+it+nt*5] = x[nx*ntlink*iruns+itlink+ntlink*5];

  return;
}


 void fystd_neu_data_FAF3A73978E04AC9C0B11BF8EA78E748(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x, double *z){
  ystd[it+nt*0] = 3.41304516704898E-2;
  ystd[it+nt*1] = 4.418516079357334E-2;
  ystd[it+nt*2] = 5.403351879598717E-2;
  ystd[it+nt*3] = 4.984687254590423E-2;
  ystd[it+nt*4] = 3.919868612840364E-2;
  ystd[it+nt*5] = 3.427852417761523E-2;

  return;
}


 void fsy_neu_data_FAF3A73978E04AC9C0B11BF8EA78E748(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *z, double *su, double *sx, double *sz){
  sy[it+nt*0] = sx[itlink+ntlink*0];
  sy[it+nt*1] = sx[itlink+ntlink*1];
  sy[it+nt*2] = sx[itlink+ntlink*2];
  sy[it+nt*3] = sx[itlink+ntlink*3];
  sy[it+nt*4] = sx[itlink+ntlink*4];
  sy[it+nt*5] = sx[itlink+ntlink*5];
  sy[it+nt*6] = sx[itlink+ntlink*6];
  sy[it+nt*7] = sx[itlink+ntlink*7];
  sy[it+nt*8] = sx[itlink+ntlink*8];
  sy[it+nt*9] = sx[itlink+ntlink*9];
  sy[it+nt*10] = sx[itlink+ntlink*10];
  sy[it+nt*11] = sx[itlink+ntlink*11];
  sy[it+nt*12] = sx[itlink+ntlink*12];
  sy[it+nt*13] = sx[itlink+ntlink*13];
  sy[it+nt*14] = sx[itlink+ntlink*14];
  sy[it+nt*15] = sx[itlink+ntlink*15];
  sy[it+nt*16] = sx[itlink+ntlink*16];
  sy[it+nt*17] = sx[itlink+ntlink*17];
  sy[it+nt*18] = sx[itlink+ntlink*18];
  sy[it+nt*19] = sx[itlink+ntlink*19];
  sy[it+nt*20] = sx[itlink+ntlink*20];
  sy[it+nt*21] = sx[itlink+ntlink*21];
  sy[it+nt*22] = sx[itlink+ntlink*22];
  sy[it+nt*23] = sx[itlink+ntlink*23];
  sy[it+nt*24] = sx[itlink+ntlink*24];
  sy[it+nt*25] = sx[itlink+ntlink*25];
  sy[it+nt*26] = sx[itlink+ntlink*26];
  sy[it+nt*27] = sx[itlink+ntlink*27];
  sy[it+nt*28] = sx[itlink+ntlink*28];
  sy[it+nt*29] = sx[itlink+ntlink*29];
  sy[it+nt*30] = sx[itlink+ntlink*30];
  sy[it+nt*31] = sx[itlink+ntlink*31];
  sy[it+nt*32] = sx[itlink+ntlink*32];
  sy[it+nt*33] = sx[itlink+ntlink*33];
  sy[it+nt*34] = sx[itlink+ntlink*34];
  sy[it+nt*35] = sx[itlink+ntlink*35];
  sy[it+nt*36] = sx[itlink+ntlink*36];
  sy[it+nt*37] = sx[itlink+ntlink*37];
  sy[it+nt*38] = sx[itlink+ntlink*38];
  sy[it+nt*39] = sx[itlink+ntlink*39];
  sy[it+nt*40] = sx[itlink+ntlink*40];
  sy[it+nt*41] = sx[itlink+ntlink*41];
  sy[it+nt*42] = sx[itlink+ntlink*42];
  sy[it+nt*43] = sx[itlink+ntlink*43];
  sy[it+nt*44] = sx[itlink+ntlink*44];
  sy[it+nt*45] = sx[itlink+ntlink*45];
  sy[it+nt*46] = sx[itlink+ntlink*46];
  sy[it+nt*47] = sx[itlink+ntlink*47];
  sy[it+nt*48] = sx[itlink+ntlink*48];
  sy[it+nt*49] = sx[itlink+ntlink*49];
  sy[it+nt*50] = sx[itlink+ntlink*50];
  sy[it+nt*51] = sx[itlink+ntlink*51];
  sy[it+nt*52] = sx[itlink+ntlink*52];
  sy[it+nt*53] = sx[itlink+ntlink*53];
  sy[it+nt*54] = sx[itlink+ntlink*54];
  sy[it+nt*55] = sx[itlink+ntlink*55];
  sy[it+nt*56] = sx[itlink+ntlink*56];
  sy[it+nt*57] = sx[itlink+ntlink*57];
  sy[it+nt*58] = sx[itlink+ntlink*58];
  sy[it+nt*59] = sx[itlink+ntlink*59];
  sy[it+nt*60] = sx[itlink+ntlink*60];
  sy[it+nt*61] = sx[itlink+ntlink*61];
  sy[it+nt*62] = sx[itlink+ntlink*62];
  sy[it+nt*63] = sx[itlink+ntlink*63];
  sy[it+nt*64] = sx[itlink+ntlink*64];
  sy[it+nt*65] = sx[itlink+ntlink*65];
  sy[it+nt*66] = sx[itlink+ntlink*66];
  sy[it+nt*67] = sx[itlink+ntlink*67];
  sy[it+nt*68] = sx[itlink+ntlink*68];
  sy[it+nt*69] = sx[itlink+ntlink*69];
  sy[it+nt*70] = sx[itlink+ntlink*70];
  sy[it+nt*71] = sx[itlink+ntlink*71];
  sy[it+nt*72] = sx[itlink+ntlink*72];
  sy[it+nt*73] = sx[itlink+ntlink*73];
  sy[it+nt*74] = sx[itlink+ntlink*74];
  sy[it+nt*75] = sx[itlink+ntlink*75];
  sy[it+nt*76] = sx[itlink+ntlink*76];
  sy[it+nt*77] = sx[itlink+ntlink*77];
  sy[it+nt*78] = sx[itlink+ntlink*78];
  sy[it+nt*79] = sx[itlink+ntlink*79];
  sy[it+nt*80] = sx[itlink+ntlink*80];
  sy[it+nt*81] = sx[itlink+ntlink*81];
  sy[it+nt*82] = sx[itlink+ntlink*82];
  sy[it+nt*83] = sx[itlink+ntlink*83];
  sy[it+nt*84] = sx[itlink+ntlink*84];
  sy[it+nt*85] = sx[itlink+ntlink*85];
  sy[it+nt*86] = sx[itlink+ntlink*86];
  sy[it+nt*87] = sx[itlink+ntlink*87];
  sy[it+nt*88] = sx[itlink+ntlink*88];
  sy[it+nt*89] = sx[itlink+ntlink*89];
  sy[it+nt*90] = sx[itlink+ntlink*90];
  sy[it+nt*91] = sx[itlink+ntlink*91];
  sy[it+nt*92] = sx[itlink+ntlink*92];
  sy[it+nt*93] = sx[itlink+ntlink*93];
  sy[it+nt*94] = sx[itlink+ntlink*94];
  sy[it+nt*95] = sx[itlink+ntlink*95];
  sy[it+nt*96] = sx[itlink+ntlink*96];
  sy[it+nt*97] = sx[itlink+ntlink*97];
  sy[it+nt*98] = sx[itlink+ntlink*98];
  sy[it+nt*99] = sx[itlink+ntlink*99];
  sy[it+nt*100] = sx[itlink+ntlink*100];
  sy[it+nt*101] = sx[itlink+ntlink*101];
  sy[it+nt*102] = sx[itlink+ntlink*102];
  sy[it+nt*103] = sx[itlink+ntlink*103];
  sy[it+nt*104] = sx[itlink+ntlink*104];
  sy[it+nt*105] = sx[itlink+ntlink*105];
  sy[it+nt*106] = sx[itlink+ntlink*106];
  sy[it+nt*107] = sx[itlink+ntlink*107];
  sy[it+nt*108] = sx[itlink+ntlink*108];
  sy[it+nt*109] = sx[itlink+ntlink*109];
  sy[it+nt*110] = sx[itlink+ntlink*110];
  sy[it+nt*111] = sx[itlink+ntlink*111];
  sy[it+nt*112] = sx[itlink+ntlink*112];
  sy[it+nt*113] = sx[itlink+ntlink*113];
  sy[it+nt*114] = sx[itlink+ntlink*114];
  sy[it+nt*115] = sx[itlink+ntlink*115];
  sy[it+nt*116] = sx[itlink+ntlink*116];
  sy[it+nt*117] = sx[itlink+ntlink*117];
  sy[it+nt*118] = sx[itlink+ntlink*118];
  sy[it+nt*119] = sx[itlink+ntlink*119];
  sy[it+nt*120] = sx[itlink+ntlink*120];
  sy[it+nt*121] = sx[itlink+ntlink*121];
  sy[it+nt*122] = sx[itlink+ntlink*122];
  sy[it+nt*123] = sx[itlink+ntlink*123];
  sy[it+nt*124] = sx[itlink+ntlink*124];
  sy[it+nt*125] = sx[itlink+ntlink*125];
  sy[it+nt*126] = sx[itlink+ntlink*126];
  sy[it+nt*127] = sx[itlink+ntlink*127];
  sy[it+nt*128] = sx[itlink+ntlink*128];
  sy[it+nt*129] = sx[itlink+ntlink*129];
  sy[it+nt*130] = sx[itlink+ntlink*130];
  sy[it+nt*131] = sx[itlink+ntlink*131];
  sy[it+nt*132] = sx[itlink+ntlink*132];
  sy[it+nt*133] = sx[itlink+ntlink*133];
  sy[it+nt*134] = sx[itlink+ntlink*134];
  sy[it+nt*135] = sx[itlink+ntlink*135];
  sy[it+nt*136] = sx[itlink+ntlink*136];
  sy[it+nt*137] = sx[itlink+ntlink*137];
  sy[it+nt*138] = sx[itlink+ntlink*138];
  sy[it+nt*139] = sx[itlink+ntlink*139];
  sy[it+nt*140] = sx[itlink+ntlink*140];
  sy[it+nt*141] = sx[itlink+ntlink*141];
  sy[it+nt*142] = sx[itlink+ntlink*142];
  sy[it+nt*143] = sx[itlink+ntlink*143];

  return;
}


 void fsystd_neu_data_FAF3A73978E04AC9C0B11BF8EA78E748(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *z, double *sy, double *su, double *sx, double *sz){


  return;
}


 void fy_scale_neu_data_FAF3A73978E04AC9C0B11BF8EA78E748(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y_scale, double *p, double *u, double *x, double *z, double *dfzdx){
  y_scale[ny*nt*iruns+it+nt*0] = 1.0;
  y_scale[ny*nt*iruns+it+nt*1] = 0.0;
  y_scale[ny*nt*iruns+it+nt*2] = 0.0;
  y_scale[ny*nt*iruns+it+nt*3] = 0.0;
  y_scale[ny*nt*iruns+it+nt*4] = 0.0;
  y_scale[ny*nt*iruns+it+nt*5] = 0.0;
  y_scale[ny*nt*iruns+it+nt*6] = 0.0;
  y_scale[ny*nt*iruns+it+nt*7] = 1.0;
  y_scale[ny*nt*iruns+it+nt*8] = 0.0;
  y_scale[ny*nt*iruns+it+nt*9] = 0.0;
  y_scale[ny*nt*iruns+it+nt*10] = 0.0;
  y_scale[ny*nt*iruns+it+nt*11] = 0.0;
  y_scale[ny*nt*iruns+it+nt*12] = 0.0;
  y_scale[ny*nt*iruns+it+nt*13] = 0.0;
  y_scale[ny*nt*iruns+it+nt*14] = 1.0;
  y_scale[ny*nt*iruns+it+nt*15] = 0.0;
  y_scale[ny*nt*iruns+it+nt*16] = 0.0;
  y_scale[ny*nt*iruns+it+nt*17] = 0.0;
  y_scale[ny*nt*iruns+it+nt*18] = 0.0;
  y_scale[ny*nt*iruns+it+nt*19] = 0.0;
  y_scale[ny*nt*iruns+it+nt*20] = 0.0;
  y_scale[ny*nt*iruns+it+nt*21] = 1.0;
  y_scale[ny*nt*iruns+it+nt*22] = 0.0;
  y_scale[ny*nt*iruns+it+nt*23] = 0.0;
  y_scale[ny*nt*iruns+it+nt*24] = 0.0;
  y_scale[ny*nt*iruns+it+nt*25] = 0.0;
  y_scale[ny*nt*iruns+it+nt*26] = 0.0;
  y_scale[ny*nt*iruns+it+nt*27] = 0.0;
  y_scale[ny*nt*iruns+it+nt*28] = 1.0;
  y_scale[ny*nt*iruns+it+nt*29] = 0.0;
  y_scale[ny*nt*iruns+it+nt*30] = 0.0;
  y_scale[ny*nt*iruns+it+nt*31] = 0.0;
  y_scale[ny*nt*iruns+it+nt*32] = 0.0;
  y_scale[ny*nt*iruns+it+nt*33] = 0.0;
  y_scale[ny*nt*iruns+it+nt*34] = 0.0;
  y_scale[ny*nt*iruns+it+nt*35] = 1.0;

  return;
}


