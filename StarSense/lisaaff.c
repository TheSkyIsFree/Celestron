/*******************************************************************************************************************************
** FILE: lisaaff.c
** DATE: 2010.02.23 DJM
** DESC: Implements the affine transformation routines used with the LISA algorithm
**
** NOTE: Copywrite (c) 2009 by Celestron.
*******************************************************************************************************************************/
#define __SKYAFF_C__
#include "lisaaff.h"                                           // interface to LISA affine transformation routines



//==============================================================================================================================
// global variable declarations and subsystem function prototypes
//==============================================================================================================================



//==============================================================================================================================
// static variable declarations
//==============================================================================================================================



/*******************************************************************************************************************************
********                                       static function implementations                                          ********
*******************************************************************************************************************************/

//==============================================================================================================================
// FUNC: (static) GJ()
//==============================================================================================================================
// DESC:
//==============================================================================================================================
static void GJ(
  double
    *Am,
  double
    *BM,
  int
    N
)
{
  int
    indxc[4],
    indxr[4],
    ipiv[4];
  int
    i, icol=0, irow=0,
    j, k, L, ll;
  double
    big, dum, pivinv, Temp;

  for(j = 1; j <= N; j++) {
    ipiv[j] = 0;
  }
  for(i = 1; i <= N; i++) {
    big = 0;
    for(j = 1; j <= N; j++) {
      if(ipiv[j] != 1) {
        for(k = 1; k <= N; k++) {
          if(ipiv[k] == 0) {
            if(fabs(Am[j + k * 4]) >= big) {
              big = (double)fabs((double)Am[j + k * 4]);
              irow = j;
              icol = k;
            }
          }
        }
      }
    }
    ipiv[icol]++;
    if(irow != icol) {
      for(L = 1; L <= N; L++) {
        Temp = Am[irow + L * 4];
        Am[irow + L * 4] = Am[icol + L * 4];
        Am[icol + L * 4] = Temp;
      }
      Temp = BM[irow];
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if(Am[icol + icol * 4] == 0) {
      //return FALSE;
    return;
    }
    pivinv = 1 / Am[icol + icol * 4];
    Am[icol + icol * 4] = 1;
    for(L = 1; L <= N; L++) {
      Am[icol + L * 4] *= pivinv;
    }
    BM[icol] *= pivinv;
    for(ll = 1; ll <= N; ll++) {
      if(ll != icol) {
        dum = Am[ll + icol * 4];
        Am[ll + icol * 4] = 0;
        for(L = 1; L <= N; L++) {
          Am[ll + L * 4] -= Am[icol + L *  4] * dum;
        }
        BM[ll] -= BM[icol] * dum;
      }
    }
  }
  for(L = N; L >= 1; L--) {
    if(indxr[L] != indxc[L]) {
      for(k = 1; k <= N; k++) {
        Temp = Am[k + indxr[L] * 4];
        Am[k + indxr[L] * 4] = Am[k + indxc[L] * 4];
        Am[k + indxc[L] * 4] = Temp;
      }
    }
  }
}



/*******************************************************************************************************************************
********                                      global function implementations                                           ********
*******************************************************************************************************************************/

//==============================================================================================================================
// FUNC: LisaAff_AffToGeo()
//==============================================================================================================================
// DESC:
//==============================================================================================================================
LisaGeoTransformTy LisaAff_AffToGeo(
  LisaAffTransformTy
    Aff
)
{
  LisaGeoTransformTy
    Geo;

  double
    Alfa,
    Beta,
    Gamma,
    Theta,
    Sn,
    Cs;


  Alfa = sqrt(Aff.A * Aff.A + Aff.B * Aff.B);
  Beta = (Aff.A * Aff.D - Aff.B * Aff.C) / (Alfa * Alfa);
  Theta = atan2(-Aff.B, Aff.A);
  Sn = sin(Theta);
  Cs = cos(Theta);
  if(fabs(Sn) > fabs(Cs)) {
    Gamma = (Cs - Aff.D / (Alfa * Beta)) / Sn;
  }
  else {
    Gamma = (-Sn + Aff.C / (Alfa * Beta)) / Cs;
  }
  Geo.Alfa = Alfa;
  Geo.Beta = Beta;
  Geo.Gamma = Gamma;
  Geo.Theta = Theta;
  Geo.U0 = Aff.U0;
  Geo.V0 = Aff.V0;
  return Geo;
}



//==============================================================================================================================
// FUNC: LisaAff_Fit()
//==============================================================================================================================
// DESC: using the first NumMatch star matches in MatchList(), this routine calculates the LSF to the Affine transformation that
//       transforms StarList1 coordinates into StarList2 coordinates. Star MatchList(j) is used only if MatchList(j).Used is
//       true the error vector DistError() contains the distance matching error in StarList1 coords, for each star the routine
//       returns the RMS fit error in StarList2 coords
//==============================================================================================================================
double LisaAff_Fit(
  LisaStarTy
    *StarList1,
  LisaStarTy
    *StarList2,
  LisaMatchListTy
    *MatchList,
  int
    NumMatch,
  LisaAffTransformTy
    *Aff,
  double
    *DistError
)
{
  double
    LisaAffFit;

  int
    i, j, I1, I2;
  double
    X, Y, U, V,
    S1, Sx, Sy, Sxx, Syy, SXY,
    Su, SV, Sux, Suy, Svx, Svy,
    Err;
  double
    dU[4],
    dV[4],
    Mu[6 * 4],
    Mv[6 * 4];

  S1 = 0;
  Sx = 0;
  Sy = 0;
  Sxx = 0;
  Syy = 0;
  SXY = 0;
  Su = 0;
  SV = 0;
  Sux = 0;
  Suy = 0;
  Svx = 0;
  Svy = 0;

  for(i = 0; i < NumMatch; i++) {
    if(MatchList[i].Used) {
      I1 = MatchList[i].StarList1Index;
      I2 = MatchList[i].StarList2Index;
      X = StarList1[I1].Xcen;
      Y = StarList1[I1].Ycen;
      U = StarList2[I2].Xcen;
      V = StarList2[I2].Ycen;
      S1 = S1 + 1;
      Sx = Sx + X;
      Sy = Sy + Y;
      Sxx = Sxx + X * X;
      Syy = Syy + Y * Y;
      SXY = SXY + X * Y;
      Su = Su + U;
      SV = SV + V;
      Sux = Sux + U * X;
      Suy = Suy + U * Y;
      Svx = Svx + V * X;
      Svy = Svy + V * Y;
    }
  }

  for(i = 0; i <= 4; i++) {
    for(j = 0; j <= 4; j++) {
      Mv[i + j * 4] = 0.0;
    }
  }

  Mv[1 + 1 * 4] = Sxx;
  Mv[1 + 2 * 4] = SXY;
  Mv[1 + 3 * 4] = Sx;
  Mv[2 + 1 * 4] = SXY;
  Mv[2 + 2 * 4] = Syy;
  Mv[2 + 3 * 4] = Sy;
  Mv[3 + 1 * 4] = Sx;
  Mv[3 + 2 * 4] = Sy;
  Mv[3 + 3 * 4] = S1;

  for(i = 0; i <= 4; i++) {
    for(j = 0; j <= 4; j++) {
      Mu[i + j * 4] = Mv[i + j * 4];
    }
  }

  dU[1] = Sux;
  dU[2] = Suy;
  dU[3] = Su;

  dV[1] = Svx;
  dV[2] = Svy;
  dV[3] = SV;

  GJ(Mu, dU, 3);
  GJ(Mv, dV, 3);

  Aff->A = dU[1];
  Aff->B = dU[2];
  Aff->U0 = dU[3];
  Aff->C = dV[1];
  Aff->D = dV[2];
  Aff->V0 = dV[3];

  Err = 0;
  S1 = 0;

  for (i = 0; i < NumMatch; i++) {
    I1 = MatchList[i].StarList1Index;
    I2 = MatchList[i].StarList2Index;
    X = StarList1[I1].Xcen;
    Y = StarList1[I1].Ycen;
    U = StarList2[I2].Xcen;
    V = StarList2[I2].Ycen;
    DistError[i] =
      pow((U - (Aff->A * X + Aff->B * Y + Aff->U0)), 2);
    DistError[i] +=
      pow((V - (Aff->C * X + Aff->D * Y + Aff->V0)), 2);
    if (MatchList[i].Used) {
      Err = Err + DistError[i];
      S1 = S1 + 1;
    }
    DistError[i] = sqrt(DistError[i]);
  }
  LisaAffFit = sqrt(Err / S1);
  return LisaAffFit;
}



//==============================================================================================================================
// FUNC: LisaAff_Transform()
//==============================================================================================================================
// DESC:
//==============================================================================================================================
void LisaAff_Transform(
  LisaAffTransformTy
    Aff,
  double
    X,
  double
    Y,
  double
    *U,
  double
    *V
)
{
  *U = Aff.A * X + Aff.B * Y + Aff.U0;
  *V = Aff.C * X + Aff.D * Y + Aff.V0;
}



//==============================================================================================================================
// FUNC: LisaAff_TransformInv()
//==============================================================================================================================
// DESC:
//==============================================================================================================================
void LisaAff_TransformInv(
  LisaAffTransformTy
    Aff,
  double
    U,
  double
    V,
  double
    *X,
  double
    *Y
)
{
  double
    Det;
  Det = Aff.A * Aff.D - Aff.C * Aff.B;
  if(Det != 0) {
    *X = (Aff.D * (U - Aff.U0) - Aff.B * (V - Aff.V0)) / Det;
    *Y = (-Aff.C * (U - Aff.U0) + Aff.A * (V - Aff.V0)) / Det;
  }
  else {
    *X = 0;
    *Y = 0;
  }
}
