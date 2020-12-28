/*******************************************************************************************************************************
** FILE: Nxp.c -- embedded implementation of PointXP (originally by Dave Rowe)
** DATE: 2009.03.23 DJM
** DESC:
**
** NOTE: Copywrite (c) 2009 by Celestron. Licensce for portions of this code that are Copyright by Keil Software is assumed via
**       purchase of Keil MDK and RL RTL packages.
*******************************************************************************************************************************/
#define __NXP_C__
#include "DmTypes.hpp"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "nxp.hpp"                                                // NexPoint model interface
#include "target.h"

//#ifndef ANDROID_NDK
//extern char sneakypath[1024];
//#define DAN_DBG
//#define DAN_DBG2
//#define DAN_DBGPRINT
//#endif


//==============================================================================================================================
// global variable declarations and subsystem function prototypes
//==============================================================================================================================
AlignStarTy
  AlignStarASt[NXP_MAXCALSTARS];
CelestialRefTy
  SkyAlignStarASt[3];
ModelTy
  NxpModelSt;
BOOL_D
  PolarAlignGotoB;
uint32_t 
  NxpNumSkyAlignStarsUT,
  PolarAlignStepUT;
char
  AlignedStars[ALIGNED_STARS_BUFSIZE] = { 0 };
AlignStarTy
  NamedDbStar1[NUM_NAMED_STARS],
  NamedDbStar2[NUM_NAMED_STARS],
  NamedDbStar3[NUM_NAMED_STARS];
int
  Ts1, Ts2, Ts3;


/*******************************************************************************************************************************
********                                       static function implementations                                          ********
*******************************************************************************************************************************/

//==============================================================================================================================
// FUNC: (static) Nxp_MdlCalStarAddHelper() -- do the add to reference
//==============================================================================================================================
// DESC: does the actual adding of the data to the model, and then requests the update to the model
//==============================================================================================================================
static void Nxp_MdlCalStarAddHelper(                            // RETR: none
  const uint32_t                                                     //
    IndexUT,                                                    // PASS: index at which to add star data
  CelestialRefTy                                                //
    StarSt                                                      // PASS: cal star data
)                                                               //
{

  //if(NxpModelSt.LatRadD < 0) {
  //  StarSt.alt += PI;                                           //
  //  if(StarSt.alt > PI_MUL_2) StarSt.alt -= PI_MUL_2;           //
  //  StarSt.azm += PI;                                           //
  //  if(StarSt.azm > PI_MUL_2) StarSt.azm -= PI_MUL_2;           //
  //}

  AlignStarASt[IndexUT].RaRadD = StarSt.ra;                       // update the reference data array
  AlignStarASt[IndexUT].DecRadD = StarSt.dec;                     //
  AlignStarASt[IndexUT].AzmRadD = StarSt.azm;                     //
  AlignStarASt[IndexUT].AltRadD = StarSt.alt;                     //
  AlignStarASt[IndexUT].LstMarkD = StarSt.lst;                    //
  AlignStarASt[IndexUT].JdMarkD = StarSt.jd;                      //
  AlignStarASt[IndexUT].DistErrRadD = 0.0;                        //
  AlignStarASt[IndexUT].RaErrRadD = 0.0;                          //
  AlignStarASt[IndexUT].DecErrRadD = 0.0;                         //
  AlignStarASt[IndexUT].IsUsedB = TRUE_D;                         //

  Nxp_MdlUpdate();                                              // update the pointing model after adding new star

#ifdef TAKI
  Align_AddTakiRef(StarSt, IndexUT);                            //
#endif

  //Align_NxpCalStarSync(&StarSt);                                //
  //os_sem_send(NxpMdlUpdateSEM);                                 // signal to NXP task to update the model
  //os_sem_wait(NxpMdlUpdateDoneSEM, 0xffff);                     // wait for Nxp model to be udpated
}



/*******************************************************************************************************************************
********                                      global function implementations                                           ********
*******************************************************************************************************************************/



//==============================================================================================================================
// FUNC: Nxp_FindErrors()
//==============================================================================================================================
// DESC: solves for the mount error parameters by LSF over all used cal stars, and sets the model instability factor
//==============================================================================================================================
BOOL_D Nxp_FindErrors(
  void
)
{
  int
    i,
    j,
    k,
    L;
  double
    D[NXP_MAXTERMS];
  int
    ErrorIndex[NXP_MAXTERMS];
  double
    Xs[3];
  double
    Xm[3];
  double
    P[NXP_MAXTERMS * NXP_MAXTERMS];
  double
    M[NXP_MAXTERMS * NXP_MAXTERMS];
  Vect3Ty
    R;
  double
    AltEnc,
    AzEnc,
    LstMarkD,
    JDmark;
  double
    RA,
    Dec,
    Theta,
    Phi;
  int
    NumErrorsUsed;                                              // the number of mount errors computed
  int
    J1,
    k1;
  double
    MaxElement = 0,
    CalStarsUsed = 0;

  for (i = 0; i < NxpModelSt.NumCalStarsUT; i++) {
    if (AlignStarASt[i].IsUsedB == TRUE_D) CalStarsUsed++;
  }
  for (j = 0; j < NxpModelSt.MaxTermsUT; j++) {
    NxpModelSt.MountErrorsD[j] = 0;
    ErrorIndex[j] = j;
    D[j] = 0;
    for (k = 0; k < NxpModelSt.MaxTermsUT; k++) {
      M[k + j * NXP_MAXTERMS] = 0;
    }
  }
  if(CalStarsUsed < 2) {
    return FALSE_D;
  }
  if(CalStarsUsed < 3) {
    NumErrorsUsed = 4;
  }
  else if(CalStarsUsed < 4) {
    NumErrorsUsed = 5;
  }
  else {
    NumErrorsUsed = 6;
  }
  for(J1 = 0; J1 < 6; J1++) {
    for(k1 = 0; k1< 6; k1++) {
      P[J1 + k1 * NXP_MAXTERMS] = 0;
    }
  }

  for (i = 0; i < NxpModelSt.NumCalStarsUT; i++) {
    if(AlignStarASt[i].IsUsedB) {
      AltEnc = AlignStarASt[i].AltRadD;
      AzEnc = AlignStarASt[i].AzmRadD;
      LstMarkD = AlignStarASt[i].LstMarkD;
      JDmark = AlignStarASt[i].JdMarkD;
      RA = AlignStarASt[i].RaRadD;
      Dec = AlignStarASt[i].DecRadD;

      R = Nxp_FnCtoX(RA, Dec, LstMarkD, JDmark);                // the true location of star in apparent 3-vect coords

      R = Nxp_FnXtoT(R);                                        // the location of the star in T-coords:
      Xs[0] = R.X;
      Xs[1] = R.Y;
      Xs[2] = R.Z;

      Nxp_FnEnctoTAA(AltEnc, AzEnc, &Theta, &Phi);              // find measured location of star in the 3-vector telescope
      R = Nxp_FnAzmAlttoX(Theta, Phi);                          // coordinate system (T-coords)

      Xm[0] = R.X;                                              // 3-vector T-coords, measured
      Xm[1] = R.Y;
      Xm[2] = R.Z;

      Nxp_PartialDerivatives(Theta, Phi, P);                    // the partial derivatives in T-coords:

      for (L = 0; L < 3; L++) {
        for (J1 = 0; J1 < NumErrorsUsed; J1++) {
          j = ErrorIndex[J1];
          //D[J1] = D[J1] + (Xs[L] - Xm[L]) * P[L + j * NXP_MAXTERMS];
          D[J1] = D[J1] + (Xs[L] - Xm[L]) * P[j + L * NXP_MAXTERMS];
          for (k1 = 0; k1 < NumErrorsUsed; k1++) {
            k = ErrorIndex[k1];
            //M[J1 + k1 * NXP_MAXTERMS] += P[L + k * NXP_MAXTERMS] * P[L + j * NXP_MAXTERMS];
            M[k1 + J1 * NXP_MAXTERMS] += P[k + L * NXP_MAXTERMS] * P[j + L * NXP_MAXTERMS];
          }
        }
      }
    }
  }

  Nxp_GJ(M, D, NumErrorsUsed);

  MaxElement = 0;                                               // find the stability parameter, ModelStability
  for (i = 0; i < NumErrorsUsed; i++) {
    for (j = 0; j < NumErrorsUsed; j++) {
      if (fabs(M[i + j * NXP_MAXTERMS]) > MaxElement) {
        //MaxElement = fabs(M[i + j * NXP_MAXTERMS]);
        MaxElement = fabs(M[j + i * NXP_MAXTERMS]);
      }
    }
  }
  NxpModelSt.InstabilityD = MaxElement;

  for (k1 = 0; k1 < NumErrorsUsed; k1++) {                      // correctly place the error parameters
    k = ErrorIndex[k1];
    NxpModelSt.MountErrorsD[k] = D[k1];
  }

  return TRUE_D;
}



//==============================================================================================================================
// FUNC: Nxp_FnAzmAlttoX()
//==============================================================================================================================
// DESC: converts from theta and phi to 3-vector coords; theta is positive downward from Z-axis toward X-axis, phi increases
//       from X-axis to Y-axis; theta, phi in radians; this is simple coordinate transform, no error terms are involved
//==============================================================================================================================
Vect3Ty Nxp_FnAzmAlttoX(
  double
    Theta,
  double
    Phi
)
{
  Vect3Ty
    P;

  P.X = cos(Phi) * sin(Theta);
  P.Y = sin(Phi) * sin(Theta);
  P.Z = cos(Theta);
  P = Nxp_Normalize(P);

  return P;
}



//==============================================================================================================================
// FUNC: Nxp_FnCtoEnc()
//==============================================================================================================================
// DESC: this is the main telescope pointing routine RA,Dec are star coords in radians JdMarkD is the Julian Day at time of
//       measurement LSTMark is the LST at time of measurement AltEnc and AzEnc are the corrected and scaled encoder outputs
//       full error model is used generates negative theta solution if Tsign = -1
//==============================================================================================================================
void Nxp_FnCtoEnc(
  double
    LstMarkD,
  double
    JdMarkD,
  double
    RA,
  double
    Dec,
  int
    Tsign,
  double
    *AltEnc,
  double
    *AzEnc
)
{
  Vect3Ty
    R;

  R = Nxp_FnCtoX(RA, Dec, LstMarkD, JdMarkD);                   // find topocentric coords of star
  R = Nxp_FnXtoT(R);                                            // transform to telescope coords
  Nxp_FnTtoEnc(R, Tsign, AltEnc, AzEnc);                        // find encoder readings that will point to R
}



//==============================================================================================================================
// FUNC: Nxp_FnCtoX()
//==============================================================================================================================
// DESC: given the J2000 celestial coordinates of the object, observer's Latitude, local sidereal time of measurement, Julian
//       Day of the measurement, this routine routine finds the X,Y,Z coordinates of object in Topocentric coordinate system
//       (X-coords) which include the effects of precession and refraction
//==============================================================================================================================
Vect3Ty Nxp_FnCtoX(
  double
    RaRad0D,
  double
    DecRad0,
  double
    LstMarkD,
  double
    JdMarkD
)
{
  double
    LHAc,
    Lc,
    Dc,
    LHA,
    JD2000 = 2451545,
    RaRad1D = 0.0,
    DecRad1D = 0.0,
    Theta,
    Phi;
  Vect3Ty
    P;

  Nxp_Precess(JD2000, JdMarkD, RaRad0D, DecRad0, &RaRad1D, &DecRad1D);

  Lc = PI_DIV_2 - NxpModelSt.LatRadD;
  Dc = PI_DIV_2 - DecRad1D;

  LHA = LstMarkD - RaRad1D;
  LHAc = -(PI_DIV_2 + LHA);

  P.Y = -cos(LHAc) * sin(Dc);
  P.X = sin(LHAc) * sin(Dc) * cos(Lc) + cos(Dc) * sin(Lc);
  P.Z = cos(Dc) * cos(Lc) - sin(LHAc) * sin(Dc) * sin(Lc);
  P = Nxp_Normalize(P);

  Nxp_FnXtoAzmAlt(P, 1, &Theta, &Phi);                          // convert to alt-az coordinates using positive theta
  Theta = Nxp_Refract(Theta, TRUE_D);                             // refract
  return Nxp_FnAzmAlttoX(Theta, Phi);                           // convert back to X-coords
}



//==============================================================================================================================
// FUNC: Nxp_EnctoC()
//==============================================================================================================================
// DESC: this the main Telescope to Celestial pointing routine given the encoder readings, it calculates where telescope is
//       pointing in RA and Dec it uses the mount error model LST is the current Local Sidereal Time JD is the current Julian
//       Day RA and Dec are in radians
//==============================================================================================================================
void Nxp_FnEnctoC(
  double
    LST,
  double
    JD,
  double
    AltEnc,
  double
    AzEnc,
  double
    *RARad,
  double
    *DecRad
)
{
  Vect3Ty
    R;

  R = Nxp_FnEnctoT(AltEnc, AzEnc);
  R = Nxp_FnTtoX(R);

  Nxp_FnXtoC(R, LST, JD, RARad, DecRad);
}



//==============================================================================================================================
// FUNC: Nxp_FnEnctoT()
//==============================================================================================================================
// DESC: given the encoder readings and the mount errors, this routine calculates the 3-vector pointing direction of the
//       telescope in the telescopes coordinate system (T-coords)
//==============================================================================================================================
Vect3Ty Nxp_FnEnctoT(
  double
    AltEnc,
  double
    AzEnc
)
{
  double
    Theta,
    Phi;
  int
    k;
  double
    P[NXP_MAXTERMS * NXP_MAXTERMS];
  Vect3Ty
    R;

  Nxp_FnEnctoTAA(AltEnc, AzEnc, &Theta, &Phi);                  // from encoders, find Theta and Phi
  R = Nxp_FnAzmAlttoX(Theta, Phi);                              // find uncorrected pointing direction 3-vector in T-coords
  Nxp_PartialDerivatives(Theta, Phi, P);                        // calculate mount error parial derivatives

  for (k = 0; k < NxpModelSt.MaxTermsUT; k++) {
    R.X = R.X + P[k + 0 * NXP_MAXTERMS] * NxpModelSt.MountErrorsD[k];
    R.Y = R.Y + P[k + 1 * NXP_MAXTERMS] * NxpModelSt.MountErrorsD[k];
    R.Z = R.Z + P[k + 2 * NXP_MAXTERMS] * NxpModelSt.MountErrorsD[k];
  }
  R = Nxp_Normalize(R);
  return R;
}



//==============================================================================================================================
// FUNC: Nxp_FnEnctoTAA()
//==============================================================================================================================
// DESC: scales encoder readings to theta and phi in radians removes encoder offsets theta and phi are T-coords
//==============================================================================================================================
void Nxp_FnEnctoTAA(
  double
    AltEnc,
  double
    AzEnc,
  double
    *Theta,
  double
    *Phi
)
{
  int
    AzSign;

  *Theta = (AltEnc - NxpModelSt.AltEncZeroD);
  if (NxpModelSt.AltIncUpB) {
    *Theta = PI_D / 2 - *Theta;
  }

  if (NxpModelSt.AzIncCwB) {
    AzSign = -1;
  }
  else {
    AzSign = 1;
  }
  *Phi = AzSign * (AzEnc - NxpModelSt.AzEncZeroD);
  if (*Phi < 0) {
    *Phi = *Phi + (2 * PI_D);
  }
  if (*Phi > (2 * PI_D)) {
    *Phi = *Phi - (2 * PI_D);
  }
}



//==============================================================================================================================
// FUNC: Nxp_FnTAAtoEnc()
//==============================================================================================================================
// DESC: scales Theta and Phi in radians to encoder ticks adds encoder offsets theta and phi are T-coords
//==============================================================================================================================
void Nxp_FnTAAtoEnc(
  double
    Theta,
  double
    Phi,
  double
    *AltEnc,
  double
    *AzEnc
)
{
  if (NxpModelSt.AltIncUpB)
  {
    //Model.AltEnc0 = AltEnc - (Math.PI / 2 - Theta);
    //AltEnc = ((Math.PI / 2) - Theta) + Model.AltEnc0;
    *AltEnc = (PI_DIV_2 - Theta) + NxpModelSt.AltEncZeroD;

  }
  else
  {
    //Model.AltEnc0 = AltEnc - Theta;
    //AltEnc = Theta + Model.AltEnc0;
    *AltEnc = Theta + NxpModelSt.AltEncZeroD;
  }
  //if (AltEnc < 0) AltEnc = AltEnc + (2 * Math.PI);
  //if (AltEnc > (2 * Math.PI)) AltEnc = AltEnc - (2 * Math.PI);
  if (*AltEnc < 0) *AltEnc = *AltEnc + PI_MUL_2;
  if (*AltEnc > PI_MUL_2) *AltEnc = *AltEnc - PI_MUL_2;


  //AzEnc = -1.0 * Phi + Model.AzEnc0;
  if (NxpModelSt.AzIncCwB)
  {
    //Model.AzEnc0 = AzEnc + Phi; //
    //AzEnc = Model.AzEnc0 - Phi;
    *AzEnc = NxpModelSt.AzEncZeroD - Phi;
  }
  else
  {
    //Model.AzEnc0 = AzEnc - Phi; //
    //AzEnc = Model.AzEnc0 + Phi;
    *AzEnc = NxpModelSt.AzEncZeroD + Phi;
  }
  //if (AzEnc < 0) AzEnc = AzEnc + (2 * Math.PI);
  //if (AzEnc > (2 * Math.PI)) AzEnc = AzEnc - (2 * Math.PI);
  if (*AzEnc < 0) *AzEnc = *AzEnc + PI_MUL_2;
  if (*AzEnc > PI_MUL_2) *AzEnc = *AzEnc - PI_MUL_2;
}



//==============================================================================================================================
// FUNC: Nxp_FnTtoEnc()
//==============================================================================================================================
// DESC: routine calculates the encoder readings needed to point to the 3-vector R in T-coords assuming mount has errors finds
//       negative theta solution if Tsign = -1
//==============================================================================================================================
void Nxp_FnTtoEnc(
  Vect3Ty
    R,
  int
    Tsign,
  double
    *AltEnc,
  double
    *AzEnc
)
{
  double
    Theta,
    Phi;

  Nxp_FnTtoTAA(R, Tsign, &Theta, &Phi);                         // using the error model, compute telescope pointing direction
  // (Theta,Phi) from telescope 3-vector coords, R
  // if Tsign = -1 use negative theta solution
  Nxp_FnTAAtoEnc(Theta, Phi, AltEnc, AzEnc);                          // from (Theta,Phi) scale to encoder readings
}



//==============================================================================================================================
// FUNC: Nxp_FnTtoTAA
//==============================================================================================================================
// DESC: given the mount errors, Point.MountErorrs(), this routine calculates the telescope pointing direction, Theta and Phi,
//       in the telescope coord system given the uncorrected, desired telescope coordinate three-vector, R Ne is the number of
//       errors to use in the solution generates negative theta solution if Tsign = -1
//==============================================================================================================================
void Nxp_FnTtoTAA(
  Vect3Ty
    R,
  int
    Tsign,
  double
    *Theta,
  double
    *Phi
)
{

  double
    X, Y, Z, Den;
  double
    X1, Y1, Z1;
  double
    Rho, Dist2;
  int
    k, Iter;
  double
    P[NXP_MAXTERMS * NXP_MAXTERMS];

  X = R.X;
  Y = R.Y;
  Z = R.Z;

  // find solution assuming all errors are zero
  *Phi = atan2(Y, X);
  Rho = sqrt(X * X + Y * Y);
  *Theta = atan2(Rho, Z);
  if (Tsign < 0){
    *Theta = -*Theta;
    *Phi = *Phi + PI_D;
    if (*Phi > (2 * PI_D)) *Phi = *Phi - (2 * PI_D);
  }

  // find telescope pointing direction assuming errors
  // by iterating until converged
  // there is no guarantee that a solution exists when
  // trying to point near the polar (azimuth) axis

  Dist2 = 1.0;
  Iter = 0;
  while ((Iter < 7) && (Dist2 > 0.000000000001)) {
    Iter++;
    // add mount errors to desired T-coord vector
    Nxp_PartialDerivatives(*Theta, *Phi, P);
    X1 = X; Y1 = Y; Z1 = Z; X = R.X; Y = R.Y; Z = R.Z;

    for (k = 0; k < NxpModelSt.MaxTermsUT; k++) {
      //if (NxpModelSt.UseError[k]) {
        X = X - P[k + 0 * NXP_MAXTERMS] * NxpModelSt.MountErrorsD[k];
        Y = Y - P[k + 1 * NXP_MAXTERMS] * NxpModelSt.MountErrorsD[k];
        Z = Z - P[k + 2 * NXP_MAXTERMS] * NxpModelSt.MountErrorsD[k];
      //}
    }

    //normalize
    Den = sqrt(X * X + Y * Y + Z * Z);
    X = X / Den;
    Y = Y / Den;
    Z = Z / Den;

    // calculate new telescope pointing direction, phi, theta
    *Phi = atan2(Y, X);
    Rho = sqrt(X * X + Y * Y);
    *Theta = atan2(Rho, Z);

    // flip 180 degrees if Tsign is negative
    if (Tsign < 0) {
      *Theta = -*Theta;
      *Phi = *Phi + PI_D;
      if (*Phi > (2 * PI_D)) *Phi = *Phi - (2 * PI_D);
    }

    // find error distance
    Dist2 = (X1 - X) * (X1 - X) + (Y1 - Y) * (Y1 - Y) + (Z1 - Z) * (Z1 - Z);
  }
}



//==============================================================================================================================
// FUNC: Nxp_FnTtoX()
//==============================================================================================================================
// DESC: transforms telescope coords, T, to X-coords telescope axis is assumed to have a pointing error Model.(ErrorWest,
//       ErrorNorth) if equatorial mount then coords are rotated by the colatitude of User.Latitude
//==============================================================================================================================
Vect3Ty Nxp_FnTtoX(
  Vect3Ty
    T
)
{
  double
    Lc;
  Vect3Ty
    P;

  memcpy(&P, &T, sizeof(Vect3Ty));
  if(!PolarAlignGotoB) {
    P = Nxp_RotAboutY(P, NxpModelSt.ErrorNorthD);
    P = Nxp_RotAboutX(P, NxpModelSt.ErrorWestD);
  }
  //else if(PolarAlignStepUT == 0) {
  //  P = Nxp_RotAboutY(P, NxpModelSt.ErrorNorthD);
  //}
  //else {
  //  P = Nxp_RotAboutX(P, NxpModelSt.ErrorWestD);
  //}
  if (NxpModelSt.IsEqAlignB) {
    Lc = (PI_DIV_2) - NxpModelSt.LatRadD;                       // get co-latitude of observer in radians
    P = Nxp_RotAboutY(P, Lc);                                   // rotate by co-latitude
  }
  return P;
}



//==============================================================================================================================
// FUNC: Nxp_FnXtoAzmAlt()
//==============================================================================================================================
// DESC: finds the AltAz coords in radians from the X,Y,Z coords given in P phi is measured from the X-axis increasing toward
//       Y-axis theta is measured from Z-axis downward if Tsign is -1 then negative theta solution is calculated this is simple
//       coordinate transform, no error terms are involved
//==============================================================================================================================
void Nxp_FnXtoAzmAlt(
  Vect3Ty
    P,
  int
    Tsign,
  double
    *Theta,
  double
    *Phi
)
{
  double Rho = sqrt(P.X * P.X + P.Y * P.Y);
  *Phi = atan2(P.Y, P.X);
  if (*Phi < 0) *Phi = *Phi + 2 * PI_D;
  *Theta = atan2(Rho, P.Z);
  if (Tsign < 0) {
    *Theta = -*Theta;
    *Phi = *Phi + PI_D;
    if (*Phi > (2 * PI_D)) *Phi = *Phi - (2 * PI_D);
  }
}



//==============================================================================================================================
// FUNC: Nxp_FnXtoC()
//==============================================================================================================================
// DESC: converts Topocentric X-coords at LSTMark LSTMark is local sidereal time celestial coords are in radians star position
//       is de-refracted and de-precessed to JD2000
//==============================================================================================================================
void Nxp_FnXtoC(
  Vect3Ty
    R,
  double
    LstMarkD,
  double
    JdMarkD,
  double
    *RARad,
  double
    *DecRad
)
{
  Vect3Ty
    P;
  double
    Den;
  double
    Xp, Yp, Zp;
  double
    Rho, Lc, Q;
  double
    Theta, Phi;
  double
    RaRad0D, DecRad0;
  double
    JD2000 = 2451545;

  *RARad = 0.0;
  *DecRad = 0.0;

  Lc = PI_DIV_2 - NxpModelSt.LatRadD;                               // co-latitude in radians

  Nxp_FnXtoAzmAlt(R, 1, &Theta, &Phi);                          // convert to alt-az coordinates using positive theta
  Theta = Nxp_Refract(Theta, FALSE_D);                            // refraction
  P = Nxp_FnAzmAlttoX(Theta, Phi);                              // convert back to X-coords

  Xp = P.X * cos(Lc) - P.Z * sin(Lc);                           // convert to x,y,z coords
  Yp = P.Y;
  Zp = P.X * sin(Lc) + P.Z * cos(Lc);

  Den = sqrt(Xp * Xp + Yp * Yp + Zp * Zp);                      // normalize
  Xp = Xp / Den;
  Yp = Yp / Den;
  Zp = Zp / Den;

  // convert to alt-ax coords
  Rho = sqrt(Xp * Xp + Yp * Yp);
  DecRad0 = atan2(Zp, Rho);
  Q = atan2(-Yp, -Xp);
  RaRad0D = LstMarkD + Q;

  Nxp_Precess(JdMarkD, JD2000, RaRad0D, DecRad0, RARad, DecRad);

  if (*RARad > (2 * PI_D)) *RARad = *RARad - (2 * PI_D);
  if (*RARad < 0.0) *RARad = *RARad + (2 * PI_D);
}



//==============================================================================================================================
// FUNC: Nxp_FnXtoT()
//==============================================================================================================================
// DESC:
//==============================================================================================================================
Vect3Ty Nxp_FnXtoT(
  Vect3Ty
    P
)
{
  //   transform x-coords to telescope coords
  //   telescope axis is assumed to have a pointing error NxpModelSt.(ErrorWestD, ErrorNorthD)
  //   if equatorial mount then coords are rotated by the colatitude of User.Latitude
  double
    Lc;
  Vect3Ty
    T;

  memcpy(&T, &P, sizeof(Vect3Ty));
  if (NxpModelSt.IsEqAlignB) {
    //Lc = (PI_DIV_2) - NxpModelSt.LatRadD * PI / 180;          // get co-latitude of observer in radians
    Lc = (PI_DIV_2) - NxpModelSt.LatRadD;                       // get co-latitude of observer in radians
    T = Nxp_RotAboutY(T, -Lc);                                  // rotate by co-latitude
  }
  if(!PolarAlignGotoB) {
    T = Nxp_RotAboutX(T, -NxpModelSt.ErrorWestD);               // in T-coord system, take out first order axis misalignment
    T = Nxp_RotAboutY(T, -NxpModelSt.ErrorNorthD);
  }
  //else if(PolarAlignStepUT == 0) {
  //  P = Nxp_RotAboutY(P, NxpModelSt.ErrorNorthD);
  //}
  //else {
  //  P = Nxp_RotAboutX(P, NxpModelSt.ErrorWestD);
  //}

  return T;
}



//==============================================================================================================================
// FUNC: Nxp_GetTSignFromRa()
//==============================================================================================================================
// DESC: if using GEM mount, determines which side of merdian object is on based on RA position, and LST
//==============================================================================================================================
int32_t  Nxp_GetTSignFromRa(                                        // RETR: TSign value (1, or -1)
  double                                                        //
    LstRadD,                                                    // PASS: Local Sidereal Time (Radians)
  double                                                        //
    RaRadD                                                      // PASS: Ra Position (Radians)
)
{
  //   this routine calculates which side of the meridian the object is on, based on RA
  //   and sets Tsign appropriately
  //   if the mount is not a GEM, Tsign is positive
  //   LSTRad is the local sidereal time, in radians
  double
    Meridian;

  if (!NxpModelSt.IsGemB) return 1;
  Meridian = LstRadD - RaRadD;
  if ((Meridian < 0 && Meridian > -PI_D) || (Meridian > PI_D && Meridian < 2 * PI_D)) return 1;
  return -1;
}



//==============================================================================================================================
// FUNC: Nxp_GetTSignFromAltEnc()
//==============================================================================================================================
// DESC: if using GEM mount, determines which side of merdian object is on based on encoder position
//==============================================================================================================================
int32_t  Nxp_GetTSignFromAltEnc(
  double
    AltEnc
)
{
  if(!NxpModelSt.IsGemB) {
    return 1;
  }
  if(AltEnc > PI_DIV_2 && AltEnc < (PI_DIV_2 * 3)) {
    return -1;
  }
  return 1;
}



//==============================================================================================================================
// FUNC: Nxp_GJ()
//==============================================================================================================================
// DESC: Adds a star to the calibration list
//==============================================================================================================================
BOOL_D Nxp_GJ(
  double
    *Am,
  double
    *Bm,
  int
    N
)
{
  // GJ stands for Gauss-Jordan, the two dudes who made the matrix inversion algorithm
  //   solves the equation b=ax
  //   N is matrix order
  //   on entry, Bm(1 to N) is the left-hand-side column vector
  //   on exit, Bm(1 to N) is the solution vector, x
  //   on entry, Am(1 to N, 1 to N) is the matrix coefficients
  //   on exit, Am is the inverse matrix
  //   returns false if matrix is singular

  int
    indxc[NXP_MAXTERMS],
    indxr[NXP_MAXTERMS],
    ipiv[NXP_MAXTERMS];
  int
    i, icol = 0, irow = 0,
    j, k, L, ll;
  double
    big, dum, pivinv, temp;
  BOOL_D Singular = FALSE_D;

  for (j = 0; j < N; j++) {
    ipiv[j] = 0;
  }
  for (i = 0; i < N; i++) {
    big = 0;
    for (j = 0; j < N; j++) {
      if (ipiv[j] != 1) {
        for (k = 0; k < N; k++) {
          if (ipiv[k] == 0) {
            //if (fabs(Am[j + k * N]) >= big) {
            //  big = fabs(Am[j + k * N]);
            if (fabs(Am[k + j * NXP_MAXTERMS]) >= big) {
              big = fabs(Am[k + j * NXP_MAXTERMS]);
              irow = j;
              icol = k;
            }
          }
        }
      }
    }
    ipiv[icol] = ipiv[icol] + 1;
    if (irow != icol) {
      for (L = 0; L < N; L++) {
        temp = Am[irow + L * NXP_MAXTERMS];
        //Am[irow + L * N] = Am[icol + L * N];
        //Am[icol + L * N] = temp;
        Am[L + irow * NXP_MAXTERMS] = Am[L + icol * NXP_MAXTERMS];
        Am[L + icol * NXP_MAXTERMS] = temp;
      }
      temp = Bm[irow];
      Bm[irow] = Bm[icol];
      Bm[icol] = temp;
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if(Am[icol + icol * NXP_MAXTERMS] != 0) {
      //pivinv = 1.0 / Am[icol + icol * 4];                     // ModV2.4
      pivinv = 1.0 / Am[icol + icol * NXP_MAXTERMS];            // ModV2.4
    }
    else {
      Singular = TRUE_D;                                          // ModV2.4
      pivinv = 1000000000000.0;                                 // ModV2.4
    }
    Am[icol + icol * NXP_MAXTERMS] = 1.0;
    for (L = 0; L < N; L++) {
      //Am[icol + L * N] = Am[icol + L * N] * pivinv;
      Am[L + icol * NXP_MAXTERMS] = Am[L + icol * NXP_MAXTERMS] * pivinv;
    }
    Bm[icol] = Bm[icol] * pivinv;
    for (ll = 0; ll < N; ll++) {
      if (ll != icol) {
        //dum = Am[ll + icol * N];
        dum = Am[icol + ll * NXP_MAXTERMS];
        //Am[ll + icol * N] = 0;
        Am[icol + ll * NXP_MAXTERMS] = 0;
        for (L = 0; L < N; L++) {
          //Am[ll + L * N] = Am[ll + L * N] - Am[icol + L * N] * dum;
          Am[L + ll * NXP_MAXTERMS] = Am[L + ll * NXP_MAXTERMS] - Am[L + icol * NXP_MAXTERMS] * dum;
        }
        Bm[ll] = Bm[ll] - Bm[icol] * dum;
      }
    }
  }

  for (L = N - 1; L >= 0; L--) {
    if (indxr[L] != indxc[L]) {
      for (k = 0; k < N; k++) {
        //temp = Am[k + indxr[L] * N];
        //Am[k + indxr[L] * N] = Am[k + indxc[L] * N];
        //Am[k + indxc[L] * N] = temp;
        temp = Am[indxr[L] + k * NXP_MAXTERMS];
        Am[indxr[L] + k * NXP_MAXTERMS] = Am[indxc[L] + k * NXP_MAXTERMS];
        Am[indxc[L] + k * NXP_MAXTERMS] = temp;
      }
    }
  }
  return !Singular;
}



#ifdef DAN_DBG2
void GenSkyAlignFileName2(char *fname, int fnamesz)
{
  time_t		timer;
  struct tm	t;

  time ( &timer );
  t = *localtime ( &timer );
  sprintf(
    fname,
    "/%04d%02d%02d_%02d%02d%02d_%02d_Stars.dat",
    t.tm_year + 1900,
    t.tm_mon + 1,
    t.tm_mday,
    t.tm_hour,
    t.tm_min,
    t.tm_sec,
    NxpModelSt.NumCalStarsUT
  );
}
#endif


//==============================================================================================================================
// FUNC: Nxp_MdlCalStarAdd() -- add's celestial reference data to NXP model
//==============================================================================================================================
// DESC: Adds a reference point to the calibration list, does not include reference for pointing if already have three refs
//==============================================================================================================================
BOOL_D Nxp_MdlCalStarAdd(                                         // RETR: result of add reference point operation
  CelestialRefTy                                                //
    StarSt                                                      // PASS: celestial refernce struct with ref data
)
{
  uint32_t                                                           //
    IndexUT;                                                    //
  if(NxpModelSt.NumCalStarsUT >= NXP_MAXCALSTARS) {             // if: maxed out the number of reference slots
    return FALSE_D;                                             //    return -- did not add reference
  }                                                             //
  IndexUT = NxpModelSt.NumCalStarsUT++;                         //
  Nxp_MdlCalStarAddHelper(IndexUT, StarSt);                     //
#ifdef DAN_DBG2
  char path[255], fname[32];
  getdocsdir ( path, sizeof ( path ) );
  GenSkyAlignFileName2(fname, 255);
  strlcat ( path, fname, sizeof ( path ) );
  FILE *fp = fopen ( path, "wb" );
  int sizeofa = sizeof(NxpModelSt);
  int sizeofb = sizeof(AlignStarTy) * NxpModelSt.NumCalStarsUT;
  if ( fp != NULL ) {
    fwrite (&NxpModelSt, sizeofa, 1, fp );
    fwrite (&AlignStarASt, sizeofb, 1, fp );
  }
  fclose ( fp );
#endif
  return TRUE_D;                                                // return -- reference point was added
}



//==============================================================================================================================
// FUNC: Nxp_MdlCalStarAddAt() -- add reference at specified index
//==============================================================================================================================
// DESC: Adds a reference point to the calibration list; replaces if the slot is occupied, if index > num refs, it appends
//==============================================================================================================================
BOOL_D Nxp_MdlCalStarAddAt(                                       // RETR: result of add reference point operation
  uint32_t                                                           //
    IndexUT,                                                    // PASS: 'at' index
  CelestialRefTy                                                //
    StarSt                                                      // PASS: celestial refernce struct with ref data
)
{
  if(IndexUT > NxpModelSt.NumCalStarsUT) {                      // if: index is out of range of list
    if(NxpModelSt.NumCalStarsUT >= NXP_MAXCALSTARS) {           //   if: maxed out the number of reference slots
      return FALSE_D;                                           //      return -- did not add reference
    }                                                           //
    IndexUT = NxpModelSt.NumCalStarsUT;                         //   (else) append to list
  }                                                             //
  Nxp_MdlCalStarAddHelper(IndexUT, StarSt);                     //
  return TRUE_D;                                                // return -- reference point was added
}



//==============================================================================================================================
// FUNC: Nxp_MdlCalStarClear() -- clears the reference array at index
//==============================================================================================================================
// DESC: clears the fields associated with the reference index (does not re-calc. model)
//==============================================================================================================================
void Nxp_MdlCalStarClear(                                       //
  uint32_t                                                           //
    RefIdxUT                                                    //
)
{
  AlignStarASt[RefIdxUT].AzmRadD = 0.0;                         //
  AlignStarASt[RefIdxUT].AltRadD = 0.0;                         //
  AlignStarASt[RefIdxUT].RaRadD = 0.0;                          //
  AlignStarASt[RefIdxUT].DecRadD = 0.0;                         //
  AlignStarASt[RefIdxUT].DistErrRadD = 0.0;                     //
  AlignStarASt[RefIdxUT].LstMarkD = 0.0;                        //
  AlignStarASt[RefIdxUT].JdMarkD = 0.0;                         //
  AlignStarASt[RefIdxUT].RaErrRadD = 0.0;                       //
  AlignStarASt[RefIdxUT].DecErrRadD = 0.0;                      //
  AlignStarASt[RefIdxUT].IsUsedB = FALSE_D;                     //
}



//==============================================================================================================================
// FUNC: Nxp_MdlCalStarDel()
//==============================================================================================================================
// DESC: Adds a star to the calibration list
//==============================================================================================================================
BOOL_D Nxp_MdlCalStarDel(                                         // RETR: result of delete star from model operation
  uint32_t                                                           //
    RefIdxUT                                                    // PASS: index of reference data in array
)
{
  if(NxpModelSt.NumCalStarsUT <= RefIdxUT) {                    // if: range check the RefIdx fails
    return FALSE_D;                                             //   return -- did not remove index
  }                                                             //
  memcpy(                                                       // remove the data at ref index location, and shift any
    (uint8_t *)AlignStarASt + RefIdxUT * sizeof(AlignStarTy),     //   data above down in the array
    (uint8_t *)AlignStarASt + (RefIdxUT + 1) * sizeof(AlignStarTy),    //
    (NxpModelSt.NumCalStarsUT - 1 - RefIdxUT) * sizeof(AlignStarTy)   //
  );                                                            //
  Nxp_MdlCalStarClear(--NxpModelSt.NumCalStarsUT);              // clear last array spot after data was shifted down
  return TRUE_D;                                                // return -- reference at index was removed
}



//==============================================================================================================================
// FUNC: Nxp_MdlCalibrate()
//==============================================================================================================================
// DESC: Performs complete mount calibration
//==============================================================================================================================
BOOL_D Nxp_MdlCalibrate(
  void
)
{
  Vect3Ty
    R;
  double
    Theta,
    Phi,
    AltEnc,
    AzEnc,
    RA,
    Dec,
    LstMarkD,
    JdMarkD,
    CalStarsUsed;
  uint32_t 
    IdxUT;

  CalStarsUsed = 0;
  for(IdxUT = 0; IdxUT < NxpModelSt.NumCalStarsUT; IdxUT++) {
    if(AlignStarASt[IdxUT].IsUsedB == TRUE_D) CalStarsUsed++;
  }

  // disable error use based on number of cal stars available
  if (CalStarsUsed < 3) NxpModelSt.MaxTermsUT = 4;
  else if (CalStarsUsed < 4) NxpModelSt.MaxTermsUT = 5;
  else NxpModelSt.MaxTermsUT = 6;

  NxpModelSt.ErrorWestD = NxpModelSt.ErrorNorthD = 0.0;
  NxpModelSt.AzEncZeroD = 0.0;
  if(NxpModelSt.IsEqAlignB == true) {
    NxpModelSt.AltEncZeroD =
      NxpModelSt.LatRadD < 0 ? 270 * PI_D / 180 : 0;
  }
  else {
    NxpModelSt.AltEncZeroD = 0.0;
  }


  //starting estimate of encoder zeros from first selected cal star
  AzEnc = 0.0;
  AltEnc = 0.0;
  RA = 0.0;
  Dec = 0.0;
  LstMarkD = 0.0;
  JdMarkD = 0.0;
  for(IdxUT = 0; IdxUT < NxpModelSt.NumCalStarsUT; IdxUT++) {   // find the first star that's used
    if (AlignStarASt[IdxUT].IsUsedB == TRUE_D) {
      AltEnc = AlignStarASt[IdxUT].AltRadD;
      AzEnc = AlignStarASt[IdxUT].AzmRadD;
      RA = AlignStarASt[IdxUT].RaRadD;
      Dec = AlignStarASt[IdxUT].DecRadD;
      LstMarkD = AlignStarASt[IdxUT].LstMarkD;
      JdMarkD = AlignStarASt[IdxUT].JdMarkD;
      break;
    }
  }

  R = Nxp_FnCtoX(RA, Dec, LstMarkD, JdMarkD);                   // calculates theta and phi in telescope coords
  R = Nxp_FnXtoT(R);                                            // transform to T-coords
  Nxp_FnXtoAzmAlt(R, Nxp_GetTSignFromAltEnc(AltEnc), &Theta, &Phi);    // calculate pointint direction in T-coords

  //one-star cal follows;  it calculates the zero positions of the encoders
  //these encoder zeros are then used in all subsequent calculations
  //AltEnc0 and AzEnc0 have units of encoder ticks (Radians... DJM)
  if(NxpModelSt.AzIncCwB) {
    NxpModelSt.AzEncZeroD = AzEnc + Phi;
  }
  else {
    NxpModelSt.AzEncZeroD = AzEnc - Phi;
  }

  if(NxpModelSt.AltIncUpB) {
    NxpModelSt.AltEncZeroD = AltEnc - (PI_D / 2 - Theta);
  }
  else {
    NxpModelSt.AltEncZeroD = AltEnc - Theta;
  }

  if (NxpModelSt.AzEncZeroD < 0.0) NxpModelSt.AzEncZeroD += (2 * PI_D);
  if (NxpModelSt.AzEncZeroD > (2 * PI_D)) NxpModelSt.AzEncZeroD -= (2 * PI_D);
  if (NxpModelSt.AltEncZeroD < 0.0) NxpModelSt.AltEncZeroD += (2 * PI_D);
  if (NxpModelSt.AltEncZeroD > (2 * PI_D)) NxpModelSt.AltEncZeroD -= (2 * PI_D);

  //find the mount errors with ErrorWest, ErrorNorth set to zero
  NxpModelSt.ErrorNorthD = 0;
  NxpModelSt.ErrorWestD = 0;
  Nxp_FindErrors();

  //transfer axis pointing errors and recalculate all errors
  NxpModelSt.ErrorNorthD = NxpModelSt.MountErrorsD[2];
  NxpModelSt.ErrorWestD = NxpModelSt.MountErrorsD[3];

  return Nxp_FindErrors();
}



//==============================================================================================================================
// FUNC: Nxp_MdlCalStarErrors()
//==============================================================================================================================
// DESC: this routine calculates the pointing errors to each of the calstars given the current mount model it should be called
//       before displaying the Calstars or the Mount model
//==============================================================================================================================
void Nxp_MdlCalStarErrors(
  void
)
{
  double
    RaRad,
    DecRad;
  double
    RaRad1D,
    DecRad1D;
  Vect3Ty
    R;
  Vect3Ty
    R1;
  double
    LstMarkD,
    JdMarkD;
  double
    AltEnc,
    AzEnc;
  double
    MaxErrorRad = 0.0,
    RMSErrorRad = 0.0;
  double
    MaxErrorAllRad = 0.0,
    RMSErrorAllRad = 0.0;
  double
    N = 0;
  int
    i;

  for (i = 0; i < NxpModelSt.NumCalStarsUT; i++) {
    AltEnc = AlignStarASt[i].AltRadD;
    AzEnc = AlignStarASt[i].AzmRadD;
    LstMarkD = AlignStarASt[i].LstMarkD;
    JdMarkD = AlignStarASt[i].JdMarkD;
    RaRad = AlignStarASt[i].RaRadD;
    DecRad = AlignStarASt[i].DecRadD;

    R = Nxp_FnCtoX(RaRad, DecRad, LstMarkD, JdMarkD);           // find location of star in apparent 3-vector coords
    R = Nxp_FnXtoT(R);                                          // find apparent location of star in T-coords
    R1 = Nxp_FnEnctoT(AltEnc, AzEnc);                           // find estimated apparent location of star in T-coords
                                                                // based on error model and encoder readings

    AlignStarASt[i].DistErrRadD =                               // find error distance in T-coords
      sqrt(
        (R.X - R1.X) * (R.X - R1.X) +
        (R.Y - R1.Y) * (R.Y - R1.Y) +
        (R.Z - R1.Z) * (R.Z - R1.Z)
      );

    Nxp_FnEnctoC(
      LstMarkD,
      JdMarkD,
      AltEnc,
      AzEnc,
      &RaRad1D,
      &DecRad1D
    );

    AlignStarASt[i].RaErrRadD = (RaRad1D - RaRad) * cos(DecRad);  //find RA and Dec error
    AlignStarASt[i].DecErrRadD = DecRad1D - DecRad;

    // calculate RMS error and MaxError
    if (AlignStarASt[i].DistErrRadD > MaxErrorAllRad) {
      MaxErrorAllRad = AlignStarASt[i].DistErrRadD;
    }
    RMSErrorAllRad +=
      AlignStarASt[i].DistErrRadD * AlignStarASt[i].DistErrRadD;

    if (AlignStarASt[i].IsUsedB) {
      if (AlignStarASt[i].DistErrRadD > MaxErrorRad) {
        MaxErrorRad = AlignStarASt[i].DistErrRadD;
          }
      RMSErrorRad +=
        AlignStarASt[i].DistErrRadD * AlignStarASt[i].DistErrRadD;
      N++;
    }
  }

  if (N > 0) {
    RMSErrorRad = sqrt(RMSErrorRad / N);
  }
  if(NxpModelSt.NumCalStarsUT > 0) {
    RMSErrorAllRad =
      sqrt(RMSErrorAllRad / NxpModelSt.NumCalStarsUT);
  }
  NxpModelSt.MaxErrorRadD = MaxErrorRad;
  NxpModelSt.RmsErrorRadD = RMSErrorRad;

  NxpModelSt.MaxErrorAllRadD = MaxErrorAllRad;
  NxpModelSt.RmsErrorAllRadD = RMSErrorAllRad;
}



//==============================================================================================================================
// FUNC: Nxp_MdlUpdate() -- update NXP model
//==============================================================================================================================
// DESC: updates the pointing model
//==============================================================================================================================
void Nxp_MdlUpdate(                                             // RETR: none
  void                                                          // PASS: none
)                                                               //
{
  Nxp_MdlCalibrate();                                           //
  Nxp_MdlCalStarErrors();                                       //
}                                                               //



//==============================================================================================================================
// FUNC: Nxp_MdlReset() -- reset NXP model
//==============================================================================================================================
// DESC: resets all the fields associated with the NXP model, along with the reference point array
//==============================================================================================================================
extern void Nxp_MdlReset(                                       // RETR: none
  void                                                          // PASS: none
)                                                               //
{
  uint32_t                                                           //
    IdxUT;                                                      // index for reference array

  NxpModelSt.PressureD = 1018.0;                                //
  NxpModelSt.TemperatureD = 16.0;                               //
  NxpModelSt.AzEncZeroD = 0.0;                                  //
  NxpModelSt.AltEncZeroD = 0.0;                                 //

  if(NxpModelSt.IsEqAlignB) {                                   //
    NxpModelSt.AltEncZeroD =                                    //
      NxpModelSt.LatRadD < 0 ? 270 * PI_D / 180 : 0;            //
    NxpModelSt.AzIncCwB =                                       //
      NxpModelSt.LatRadD < 0 ? FALSE_D : TRUE_D;                //
    NxpModelSt.AltIncUpB =                                      //
      NxpModelSt.LatRadD < 0 ? FALSE_D : TRUE_D;                //
  }                                                             //
  else {                                                        //
    NxpModelSt.AltEncZeroD = 0.0;                               //
    NxpModelSt.AzIncCwB = TRUE_D;                               //
    NxpModelSt.AltIncUpB = TRUE_D;                              //
  }                                                             //

  NxpModelSt.ErrorNorthD = 0.0;                                 //
  NxpModelSt.ErrorWestD = 0.0;                                  //
  NxpModelSt.InstabilityD = 0.0;                                //
  NxpModelSt.MaxErrorRadD = 0.0;                                //
  NxpModelSt.RmsErrorRadD = 0.0;                                //
  NxpModelSt.MaxErrorAllRadD = 0.0;                             //
  NxpModelSt.RmsErrorAllRadD = 0.0;                             //
  for(IdxUT = 0; IdxUT < NXP_MAXCALSTARS; IdxUT++) {            // for: each of the reference data slots in the aray
    Nxp_MdlCalStarClear(IdxUT);                                 //   clears the reference data
  }                                                             //
  for (IdxUT = 0; IdxUT < NXP_MAXTERMS; IdxUT++) {
    NxpModelSt.MountErrorsD[IdxUT] = 0;
  }
  NxpModelSt.NumCalStarsUT = 0;                                 // clear the reference point count
}                                                               //



//==============================================================================================================================
// FUNC: Nxp_Normalize()
//==============================================================================================================================
// DESC:
//==============================================================================================================================
Vect3Ty Nxp_Normalize(
  Vect3Ty
    P
)
{
  double
    Length = sqrt(P.X * P.X + P.Y * P.Y + P.Z * P.Z);
  if (Length > 0) {
    P.X = P.X / Length;
    P.Y = P.Y / Length;
    P.Z = P.Z / Length;
  }
  return P;
}



//==============================================================================================================================
// FUNC: Nxp_PartialDerivatives()
//==============================================================================================================================
// DESC: precess the Ra and Dec from one epoch Jd0D to another Jd1D
//==============================================================================================================================
void Nxp_PartialDerivatives(
  double
    Theta,
  double
    Phi,
  double
    *P
)
{
  //   given Theta and Phi in the telescope coord system this routine generates the partial
  //   derivatives of the mount errors, P(l,k)
  //   l=1, 2 or 3 (x, y or z)
  //   k=1 phi offset error
  //    =2 theta offset error
  //    =3 azimuth axis pointing error to North
  //    =4 azimuth axis pointing error to West
  //    =5 optical axis non-orthogonality (Cone)
  //    =6 non-orthogonal axes (Hub)
  //   angles are in radians
  //   theta is zero along the az(RA) axis and increases downward
  //   phi is zero pointing North and increases Westward

  double
    CT,
    ST,
    CP,
    SP;

  CT = cos(Theta);
  ST = sin(Theta);
  CP = cos(Phi);
  SP = sin(Phi);

  //phi offset error
  P[0 + 0 * NXP_MAXTERMS] = -SP * ST;
  P[0 + 1 * NXP_MAXTERMS] = CP * ST;
  P[0 + 2 * NXP_MAXTERMS] = 0;
  //theta offset error
  P[1 + 0 * NXP_MAXTERMS] = CP * CT;
  P[1 + 1 * NXP_MAXTERMS] = SP * CT;
  P[1 + 2 * NXP_MAXTERMS] = -ST;
  //azimuth axis pointing error to North
  P[2 + 0 * NXP_MAXTERMS] = CT;
  P[2 + 1 * NXP_MAXTERMS] = 0;
  P[2 + 2 * NXP_MAXTERMS] = -CP * ST;
  //azimuth axis pointing error to West
  P[3 + 0 * NXP_MAXTERMS] = 0;
  P[3 + 1 * NXP_MAXTERMS] = CT;
  P[3 + 2 * NXP_MAXTERMS] = -SP * ST;
  //optical axis non-orthogonality (Cone)
  P[4 + 0 * NXP_MAXTERMS] = -SP;
  P[4 + 1 * NXP_MAXTERMS] = CP;
  P[4 + 2 * NXP_MAXTERMS] = 0;
  //non-orthogonal axes (Hub)
  P[5 + 0 * NXP_MAXTERMS] = -SP * CT;
  P[5 + 1 * NXP_MAXTERMS] = CP * CT;
  P[5 + 2 * NXP_MAXTERMS] = 0;
}



//==============================================================================================================================
// FUNC: Nxp_Precess()
//==============================================================================================================================
// DESC: precess the Ra and Dec from one epoch Jd0D to another Jd1D
//==============================================================================================================================
void Nxp_Precess(
  double
   Jd0D,
 double
   Jd1D,
 double
   RaRad0D,
 double
   DecRad0,
 double
   *RaRad1D,
 double
   *DecRad1D
)
{
  double
    J0,
    J1,
    T,
    dt;
  double
    xi,
    zeta,
    Theta;
  double
    A,
    B,
    C;
  double
    cxi,
    sxi,
    czeta,
    szeta;
  double
    ctheta,
    stheta;
  double
    Xx,
    Yx,
    Zx;
  double
    Xy,
    Yy,
    Zy;
  double
    Xz,
    Yz,
    Zz;
  double
    X0,
    Y0,
    Z0;
  double
    X1,
    Y1,
    Z1;
  double
    Rho;

  J0 = 2000 + (Jd0D - 2451545) / 365.25;                        // calculate Julian Epoch
  J1 = 2000 + (Jd1D - 2451545) / 365.25;

  T = (J0 - 2000) / 100;                                        // Julian Centuries
  dt = (J1 - J0) / 100;

  A = 2306.2181 + 1.39656 * T - 0.000139 * T * T;
  B = 2004.3109 - 0.8533 * T - 0.000217 * T * T;

  xi = A * dt + (0.30188 - 0.000345 * T) * dt * dt;             // xi, zeta and theta are in arcseconds
  zeta = A * dt + (1.09468 + 0.000066 * T) * dt * dt;
  Theta = B * dt + (-0.42665 - 0.000217 * T) * dt * dt;

  C = PI_D / 648000;                                         // turn them into radians
  xi = xi * C;
  zeta = zeta * C;
  Theta = Theta * C;

  cxi = cos(xi);
  sxi = sin(xi);
  czeta = cos(zeta);
  szeta = sin(zeta);
  ctheta = cos(Theta);
  stheta = sin(Theta);

  Xx = cxi * ctheta * czeta - sxi * szeta;
  Yx = -sxi * ctheta * czeta - cxi * szeta;
  Zx = -stheta * czeta;
  Xy = cxi * ctheta * szeta + sxi * czeta;
  Yy = -sxi * ctheta * szeta + cxi * czeta;
  Zy = -stheta * szeta;
  Xz = cxi * stheta;
  Yz = -sxi * stheta;
  Zz = ctheta;

  X0 = cos(RaRad0D) * cos(DecRad0);
  Y0 = sin(RaRad0D) * cos(DecRad0);
  Z0 = sin(DecRad0);

  X1 = Xx * X0 + Yx * Y0 + Zx * Z0;
  Y1 = Xy * X0 + Yy * Y0 + Zy * Z0;
  Z1 = Xz * X0 + Yz * Y0 + Zz * Z0;

  Rho = sqrt(X1 * X1 + Y1 * Y1);
  *DecRad1D = atan2(Z1, Rho);
  *RaRad1D = atan2(Y1, X1);
  if (*RaRad1D < 0) {
    *RaRad1D = *RaRad1D + 2 * PI_D;
  }
  if (*RaRad1D > 2 * PI_D) {
    *RaRad1D = *RaRad1D - 2 * PI_D;
  }
}



//==============================================================================================================================
// FUNC: Nxp_Refract()
//==============================================================================================================================
// DESC:
//==============================================================================================================================
double Nxp_Refract(
  double
    Theta,
  BOOL_D
    ToApparentPositionB
)
{
  // input is zenith distance!
  //  approximate formula for refraction
  //  if ToApparentPositionB is true then routine refracts from no atmosphere to coords with atmosphere
  //  otherwise routine returns coords with no atmosphere
  //  input and output are in radians

  double
    P;                                                          // the pressure in millibars
  double
    Tc;                                                         // the temperature in degrees C
  double
    dEl,
    El;

  P = NxpModelSt.PressureD;
  Tc = NxpModelSt.TemperatureD;

  El = ((PI_DIV_2) - Theta) * (180.0 / PI_D);                     // elevation of star in degrees
  dEl =
    (P * 283.0 / (1010.0 * (273.0 + Tc))) * 1.02 /              // change in elevation in arcminutes
    tan(DEG_TO_RAD_D * (El + 10.3 / (El + 5.11) + 0.0019));
  dEl = (dEl / 60.0) * DEG_TO_RAD_D;                              // change in elevation in radians
  if (ToApparentPositionB) {
    Theta = Theta - dEl;
  }
  else {
    Theta = Theta + dEl;
  }
  return Theta;
}



//==============================================================================================================================
// FUNC: Nxp_RotAboutX()
//==============================================================================================================================
// DESC:
//==============================================================================================================================
Vect3Ty Nxp_RotAboutX(
  Vect3Ty
    P,
  double
    Alfa
)
{
  // rotates 3-vector P about X axis Alfa in radians
  // if alfa is positive Zenith rotates towards positive Y-axis

  double
    tempY = P.Y;
  P.Y = tempY * cos(Alfa) + P.Z * sin(Alfa);
  P.Z = -tempY * sin(Alfa) + P.Z * cos(Alfa);
  return P;
}



//==============================================================================================================================
// FUNC: Nxp_RotAboutY()
//==============================================================================================================================
// DESC:
//==============================================================================================================================
Vect3Ty Nxp_RotAboutY(
  Vect3Ty
    P,
  double
    Alfa
)
{
  // rotates 3-vector P about Y axis Alfa in radians
  // if alfa is positive Z rotates towards positive X-axis

  double
    tempX = P.X;
  P.X = tempX * cos(Alfa) + P.Z * sin(Alfa);
  P.Z = -tempX * sin(Alfa) + P.Z * cos(Alfa);
  return P;
}



//==============================================================================================================================
// FUNC: Nxp_RotAboutZ()
//==============================================================================================================================
// DESC:
//==============================================================================================================================
Vect3Ty Nxp_RotAboutZ(
	Vect3Ty
	  P,
	double
	  Alpha
)
{
	// rotates 3-vector P about Z axis Alpha in radians
	// if Alpha is positive Z rotates towards positive X-axis
	
	double
	x0 = P.X,
	y0 = P.Y;
	P.X = x0 * cos(Alpha) + y0 * sin(Alpha);
	P.Y = -x0 * sin(Alpha) + y0 * cos(Alpha);
	return P;
}



//==============================================================================================================================
// FUNC: Nxp_GetTruePae()
//==============================================================================================================================
// DESC:
//==============================================================================================================================
PaeTy Nxp_GetTruePae(
  void
)
{
	double Lc;
	PaeTy pae;
	
	Vect3Ty polarAxis = {0.0, 0.0, 1.0};
	Vect3Ty phiZeroAxis = {1.0, 0.0, 0.0};
	
	polarAxis = Nxp_FnTtoX(polarAxis);
	phiZeroAxis = Nxp_FnTtoX(phiZeroAxis);
	
	// Rotate in Azimuth to align the polar axis vector in the X-Z plane
	pae.azmRad = atan2(polarAxis.Y, polarAxis.X); // Positive = mount currently rotated too far CCW viewed from zenith
	polarAxis = Nxp_RotAboutZ(polarAxis, -pae.azmRad);
	phiZeroAxis = Nxp_RotAboutZ(phiZeroAxis, -pae.azmRad);
	
	// Undo the latitude adjustment by Un-rotating by co-latitude
	Lc = (PI_D / 2) - NxpModelSt.LatRadD; // get co-latitude of observer in radians
	polarAxis = Nxp_RotAboutY(polarAxis, -Lc);
	phiZeroAxis = Nxp_RotAboutY(phiZeroAxis, -Lc);
	
	// Measure the rotation needed to bring vector into Y-Z plane and Un-rotate
	pae.altRad = atan2(polarAxis.X, polarAxis.Z); // Positive angle = axis currently below NCP and needs to be raised (in Northern hemisphere) or lowered (in Southern hemisphere)
	polarAxis = Nxp_RotAboutY(polarAxis, -pae.altRad);
	phiZeroAxis = Nxp_RotAboutY(phiZeroAxis, -pae.altRad);
	
	// Coordinate system is still probably rolled around Phi
	pae.phiOffsetRad = atan2(phiZeroAxis.Y, phiZeroAxis.X);
	
	return pae;
}



//==============================================================================================================================
// FUNC: Nxp_RotateToZeroPae()
//==============================================================================================================================
// DESC: Will modify the model so that the mount model is polar aligned
//==============================================================================================================================
void Nxp_RotateToZeroPae(
  void
)
{
	PaeTy pae;
	
	pae = Nxp_GetTruePae();
	
	// set the model's PAE terms to zero
	NxpModelSt.ErrorNorthD = NxpModelSt.MountErrorsD[2] = 0;
	NxpModelSt.ErrorWestD = NxpModelSt.MountErrorsD[3] = 0;
	
	// modify the azm zero offset
	if(NxpModelSt.LatRadD > 0)
		NxpModelSt.AzEncZeroD -= pae.phiOffsetRad;
	else
		NxpModelSt.AzEncZeroD += pae.phiOffsetRad;
}


//==============================================================================================================================
// FUNC: Nxp_Init() -- initialize the Point XP model
//==============================================================================================================================
// DESC:
//==============================================================================================================================
void Nxp_Init(                                                        // RETR: none
  double                                                              //
    LonRad,                                                           // PASS: observer longitude in radians E=POS
  double                                                              //
    LatRad,                                                           // PASS: observer latitude in radians N=POS
  bool                                                                //
    isGemB,                                                           // PASS: mount is a GEM
  bool                                                                //
    isEqAlignB,                                                       // PASS: align is eq (GEM / Wedged Alt/Azm)
  double                                                              //
    UtcOffsetDays,                                                    // PASS: UTC offset IN DAYS (departure from orig)
  bool                                                                //
    isDst                                                             // PASS: flag - daylight savings
)                                                                     //
{

  printf("***** NXP_INIT:   Lon: %1.5f  Lat: %1.5f  UtcOffs: %1.5f\n", LonRad, LatRad, UtcOffsetDays );

  PolarAlignGotoB = FALSE_D;                                          //
  PolarAlignStepUT = 0;                                               //

  NxpModelSt.IsEqAlignB = isEqAlignB;                                 // assume Alt/Az align (COSMOS! / CEVO!)
  NxpModelSt.IsGemB = isGemB;                                         //
  NxpModelSt.LatRadD = LatRad;                                        // initialize model starting parameters
  NxpModelSt.LonRadD = LonRad;                                        //
  NxpModelSt.UtcOffsetDaysD = UtcOffsetDays;                          //
  NxpModelSt.IsDstB = isDst;                                          //
  Nxp_MdlReset();                                                     // reset the NXP model
  NxpNumSkyAlignStarsUT = 0;                                          //

  strlcpy ( AlignedStars, "None", ALIGNED_STARS_BUFSIZE );
}



//==============================================================================================================================
// FUNC: Nxp_CalcSeparation() -- calc separation between two spherical coords
//==============================================================================================================================
// DESC: support func for spherical separation calc of telescope and celestial coords
//==============================================================================================================================
double Nxp_CalcSeparation(                                            // RETR: angular separation in radians
  double                                                              //
    a1,                                                               // PASS: coord 1, azm/ra
  double                                                              //
    d1,                                                               // PASS: coord 1, alt/dec
  double                                                              //
    a2,                                                               // PASS: coord 2, azm/ra
  double                                                              //
    d2                                                                // PASS: coord 2, alt/dec
)
{
  double
    aa,
    ad;

  aa = (a1 - a2);
  if (aa > PI_D) aa = 2 * PI_D - aa;
  if (aa < -PI_D) aa = 2 * PI_D + aa;

  ad = (d1 - d2);
  if (ad > PI_D) ad = 2 * PI_D - ad;
  if (ad < -PI_D) ad = 2 * PI_D + ad;

  double
    dist = acos(cos(ad) - cos(d1) * cos(d2) * (1 - cos(aa)));

  if (dist > PI_D)
  {
    dist = (2 * PI_D - dist);
  }
  return dist;
}



//==============================================================================================================================
// FUNC: Nxp_CalcSeparationRaDec() -- calc separation between two celestial coords
//==============================================================================================================================
// DESC: spherical separation calc of celestial coords
//==============================================================================================================================
double Nxp_CalcSeparationRaDec(                                       // RETR: angular separation in radians
  AlignStarTy                                                         //
    StarA,                                                            // PASS: first celestial coord in db entry struct
  AlignStarTy                                                         //
    StarB                                                             // PASS: second celestial coord in db entry struct
)
{
  return Nxp_CalcSeparation(                                          // Call helper function for actual calculation, passing
    StarA.RaRadD,                                                     //   celestial coord params of two references
    StarA.DecRadD,                                                    //
    StarB.RaRadD,                                                     //
    StarB.DecRadD                                                     //
  );                                                                  //
}



//==============================================================================================================================
// FUNC: Nxp_CalcDistAzmAlt() -- calc separation between two telescope coords
//==============================================================================================================================
// DESC: spherical separation calc of telescope coords
//==============================================================================================================================
double Nxp_CalcDistAzmAlt(                                            // RETR: angular separation in radians
  CelestialRefTy                                                      //
    StarA,                                                            // PASS: first telescope coord in align struct
  CelestialRefTy                                                      //
    StarB                                                             // PASS: second telescope coord in align struct
)
{
  return Nxp_CalcSeparation(                                          // Call helper function for actual calculation, passing
    StarA.azm,                                                        //   telescope coord params of two references
    StarA.alt,                                                        //
    StarB.azm,                                                        //
    StarB.alt                                                         //
  );                                                                  //
}



//==============================================================================================================================
// FUNC: Nxp_CalcDistAboveHoriz() -- calc alt above horizon of star/object
//==============================================================================================================================
// DESC: converts celestial coords to horizon coords, and with LST determines angle above horizon
//==============================================================================================================================
double Nxp_CalcDistAboveHoriz(                                        // RETR: distance above horiz in radians
  double                                                              //
    Ra,                                                               // PASS: Ra of star/object
  double                                                              //
    Dec,                                                              // PASS: Dec of star/object
  double                                                              //
    Lst                                                               // PASS: Local Sidereal Time of observation
)
{
  double
    DistFmZentithF;

  DistFmZentithF =                                                    // transforms celestial to horizon (Alt)
    acos(                                                             //
      sin(Dec) * sin(NxpModelSt.LatRadD) +                            //
      cos(Dec) * cos(NxpModelSt.LatRadD) * cos(Ra - Lst)              //
    );                                                                //
  return PI_DIV_2 - DistFmZentithF;                                   //
}



//==============================================================================================================================
// FUNC: Nxp_CalcIsNamedStarUp() -- calc if named star db object is up / used
//==============================================================================================================================
// DESC: calculates altitude above horizon for each named star
//==============================================================================================================================
void Nxp_CalcIsNamedStarUp(
  void
) {
  int s1;

  int TSign;

  for(s1 = 0; s1 < NUM_NAMED_STARS; s1++) {

    if(s1 == Ts1 || s1 == Ts2 || s1 == Ts3){
      TSign = 1;
    }

    if(
      Nxp_CalcDistAboveHoriz(
        NamedStarsDb[s1].RaRadD,
        NamedStarsDb[s1].DecRadD,
        SkyAlignStarASt[0].lst
      ) < 0 && s1 != Ts1
    ) {
      NamedDbStar1[s1].IsUsedB = false;
    }
    else {
      NamedDbStar1[s1].IsUsedB = true;
      NamedDbStar1[s1].RaRadD = NamedStarsDb[s1].RaRadD;
      NamedDbStar1[s1].DecRadD = NamedStarsDb[s1].DecRadD;
      NamedDbStar1[s1].LstMarkD = SkyAlignStarASt[0].lst;
      NamedDbStar1[s1].JdMarkD = SkyAlignStarASt[0].jd;
      TSign = Nxp_GetTSignFromRa(NamedStarsDb[s1].RaRadD, NamedStarsDb[s1].DecRadD);
      Nxp_FnCtoEnc(SkyAlignStarASt[0].lst, NamedDbStar1[s1].JdMarkD, NamedDbStar1[s1].RaRadD, NamedDbStar1[s1].DecRadD, TSign, &NamedDbStar1[s1].AltRadD, &NamedDbStar1[s1].AzmRadD);
    }

    if(
      Nxp_CalcDistAboveHoriz(
        NamedStarsDb[s1].RaRadD,
        NamedStarsDb[s1].DecRadD,
        SkyAlignStarASt[1].lst
    ) < 0 && s1 != Ts2
      ) {
      NamedDbStar2[s1].IsUsedB = false;
    }
    else {
      NamedDbStar2[s1].IsUsedB = true;
      NamedDbStar2[s1].RaRadD = NamedStarsDb[s1].RaRadD;
      NamedDbStar2[s1].DecRadD = NamedStarsDb[s1].DecRadD;
      NamedDbStar2[s1].LstMarkD = SkyAlignStarASt[1].lst;
      NamedDbStar2[s1].JdMarkD = SkyAlignStarASt[1].jd;
      TSign = Nxp_GetTSignFromRa(NamedStarsDb[s1].RaRadD, NamedStarsDb[s1].DecRadD);
      Nxp_FnCtoEnc(SkyAlignStarASt[1].lst, NamedDbStar2[s1].JdMarkD, NamedDbStar2[s1].RaRadD, NamedDbStar2[s1].DecRadD, TSign, &NamedDbStar2[s1].AltRadD, &NamedDbStar2[s1].AzmRadD);
    }

    if(
      Nxp_CalcDistAboveHoriz(
        NamedStarsDb[s1].RaRadD,
        NamedStarsDb[s1].DecRadD,
        SkyAlignStarASt[2].lst
      ) < 0 && s1 != Ts3
    ) {
      NamedDbStar3[s1].IsUsedB = false;
    }
    else {
      NamedDbStar3[s1].IsUsedB = true;
      NamedDbStar3[s1].RaRadD = NamedStarsDb[s1].RaRadD;
      NamedDbStar3[s1].DecRadD = NamedStarsDb[s1].DecRadD;
      NamedDbStar3[s1].LstMarkD = SkyAlignStarASt[2].lst;
      NamedDbStar3[s1].JdMarkD = SkyAlignStarASt[2].jd;
      TSign = Nxp_GetTSignFromRa(NamedStarsDb[s1].RaRadD, NamedStarsDb[s1].DecRadD);
      Nxp_FnCtoEnc(SkyAlignStarASt[2].lst, NamedDbStar3[s1].JdMarkD, NamedDbStar3[s1].RaRadD, NamedDbStar3[s1].DecRadD, TSign, &NamedDbStar3[s1].AltRadD, &NamedDbStar3[s1].AzmRadD);
    }
  }
}



void CorrectRadianRange(double *radval)
{
  if(*radval > PI_D) {
    *radval = (PI_MUL_2 - *radval);
  }
  else if (*radval < -PI_D) {
    *radval = (PI_MUL_2 + *radval);
  }
}

void CalcRiseRunMain(double azm1, double azm2, double alt1, double alt2, double *rise, double *run)
{
  *run = azm1 - azm2;
  CorrectRadianRange(run);

  *rise = alt1 - alt2;
  CorrectRadianRange(rise);
}

void CalcRiseRun(AlignStarTy star1, AlignStarTy star2, double *rise, double *run)
{
  CalcRiseRunMain(star1.AzmRadD, star2.AzmRadD, star1.AltRadD, star2.AltRadD, rise, run);
}

double CalcError(CelestialRefTy star1a, AlignStarTy star2a, CelestialRefTy star1b, AlignStarTy star2b)
{
  double e1, e2, t1, t2;

  t1 = star1a.azm - star2a.AzmRadD;
  CorrectRadianRange(&t1);
  t2 = star1b.azm - star2b.AzmRadD;
  CorrectRadianRange(&t2);
  e1 = t1 - t2;
  CorrectRadianRange(&e1);
  e1 = e1 * e1;

  t1 = star1a.alt - star2a.AltRadD;
  CorrectRadianRange(&t1);
  t2 = star1b.alt - star2b.AltRadD;
  CorrectRadianRange(&t2);
  e2 = t1 - t2;
  CorrectRadianRange(&e1);
  e2 = e2 * e2;
  return e1 + e2;
}

void CalcRiseRun(CelestialRefTy star1, CelestialRefTy star2, double *rise, double *run)
{
  CalcRiseRunMain(star1.azm, star2.azm, star1.alt, star2.alt, rise, run);
}

bool BestOf5(double d1, double d2, double d3, double d4, double d5, double d6, double thresh, double *riserun)
{
  if((fabs(d1) < thresh && fabs(d2) < thresh && fabs(d3) < thresh && fabs(d4) < thresh && fabs(d5) < thresh)) {
    *riserun = fabs(d1) + fabs(d2) + fabs(d3) + fabs(d4) + fabs(d5);
    return true;
  }
  else if((fabs(d1) < thresh && fabs(d2) < thresh && fabs(d3) < thresh && fabs(d4) < thresh && fabs(d6) < thresh)) {
    *riserun = fabs(d1) + fabs(d2) + fabs(d3) + fabs(d4) + fabs(d6);
    return true;
  }
  else if((fabs(d1) < thresh && fabs(d2) < thresh && fabs(d3) < thresh && fabs(d5) < thresh && fabs(d6) < thresh)) {
    *riserun = fabs(d1) + fabs(d2) + fabs(d3) + fabs(d5) + fabs(d6);
    return true;
  }
  else if((fabs(d1) < thresh && fabs(d2) < thresh && fabs(d4) < thresh && fabs(d5) < thresh && fabs(d6) < thresh)) {
    *riserun = fabs(d1) + fabs(d2) + fabs(d4) + fabs(d5) + fabs(d6);
    return true;
  }
  else if((fabs(d1) < thresh && fabs(d3) < thresh && fabs(d4) < thresh && fabs(d5) < thresh && fabs(d6) < thresh)) {
    *riserun = fabs(d1) + fabs(d3) + fabs(d4) + fabs(d5) + fabs(d6);
    return true;
  }
  else if((fabs(d2) < thresh && fabs(d3) < thresh && fabs(d4) < thresh && fabs(d5) < thresh && fabs(d6) < thresh)) {
    *riserun = fabs(d2) + fabs(d3) + fabs(d4) + fabs(d5) + fabs(d6);
    return true;
  }
  else {
    return false;
  }
}

void Nxp_SkyAlignAddStar(CelestialRefTy star)
{
  SkyAlignStarASt[NxpNumSkyAlignStarsUT].ra = star.ra;
  SkyAlignStarASt[NxpNumSkyAlignStarsUT].dec = star.dec;
  SkyAlignStarASt[NxpNumSkyAlignStarsUT].azm = star.azm;
  SkyAlignStarASt[NxpNumSkyAlignStarsUT].alt = star.alt;
  SkyAlignStarASt[NxpNumSkyAlignStarsUT].lst = star.lst;
  SkyAlignStarASt[NxpNumSkyAlignStarsUT].jd = star.jd;
  NxpNumSkyAlignStarsUT++;
}

#ifdef DAN_DBG
void GenSkyAlignFileName(char *fname, int fnamesz)
{
  time_t		timer;
  struct tm	t;

  time ( &timer );
  t = *localtime ( &timer );
  sprintf(
    fname,
    "/%04d%02d%02d_%02d%02d%02d_Stars.dat",
    t.tm_year + 1900,
    t.tm_mon + 1,
    t.tm_mday,
    t.tm_hour,
    t.tm_min,
    t.tm_sec
  );
}
#endif


#define RMS_THRESH      (0.0002908 * 5.5)         // about 5.5'
#define CALCERR_THRESH  (0.0062)                  // threshold for "CalcErr" method (originial skyalign)
#define RISERUN_THRESH1 (0.070)                   // threshold for "RiseRun" method
#define RISERUN_THRESH2 (0.117)
#define DIST_THRESH     (0.015)                   // threshold for "DistErr" method



double Nxp_SkyAlign(void)
{
#ifdef DAN_DBGSKYALIGN
////Align Star 1 -- Alphard (32)
//  SkyAlignStarASt[0].azm = 0.3961846192;
//  SkyAlignStarASt[0].alt = 1.1917573907;
//  SkyAlignStarASt[0].ra = 2.4765654854;
//  SkyAlignStarASt[0].dec = -0.1511213700;
//  SkyAlignStarASt[0].lst = 3.1179523181;
//  SkyAlignStarASt[0].jd = 2456692.9512731500;
////Align Star 2 -- Procyon (163)
//  SkyAlignStarASt[1].azm = 0.9765997320;
//  SkyAlignStarASt[1].alt = 1.0307879031;
//  SkyAlignStarASt[1].ra = 2.0040830402;
//  SkyAlignStarASt[1].dec = 0.0911934530;
//  SkyAlignStarASt[1].lst = 3.1204316387;
//  SkyAlignStarASt[1].jd = 2456692.9516666700;
////Align Star 3 -- Pollux (160)
//  SkyAlignStarASt[2].azm = 1.3522737248;
//  SkyAlignStarASt[2].alt = 1.2611921173;
//  SkyAlignStarASt[2].ra = 2.0303197022;
//  SkyAlignStarASt[2].dec = 0.4891480119;
//  SkyAlignStarASt[2].lst = 3.1226192733;
//  SkyAlignStarASt[2].jd = 2456692.9520138900;


//Align Star 1 -- Segius (187)
SkyAlignStarASt[0].azm = 2.5963823295;
SkyAlignStarASt[0].alt = 5.5907869841;
SkyAlignStarASt[0].ra = 3.8051537452;
SkyAlignStarASt[0].dec = 0.6686016030;
SkyAlignStarASt[0].lst = 1.7831605428;
SkyAlignStarASt[0].jd = 2456713.6815625000;
//Align Star 2 -- Thuban (205)
SkyAlignStarASt[1].azm = 2.3201332118;
SkyAlignStarASt[1].alt = 5.9310513373;
SkyAlignStarASt[1].ra = 3.6843376911;
SkyAlignStarASt[1].dec = 1.1235702500;
SkyAlignStarASt[1].lst = 1.7865149155;
SkyAlignStarASt[1].jd = 2456713.6820949100;
//Align Star 3 -- Tsih (206)
SkyAlignStarASt[2].azm = 1.2714663432;
SkyAlignStarASt[2].alt = 6.0491064389;
SkyAlignStarASt[2].ra = 0.2474353281;
SkyAlignStarASt[2].dec = 1.0597057441;
SkyAlignStarASt[2].lst = 1.7905255781;
SkyAlignStarASt[2].jd = 2456713.6827314800;




  NxpModelSt.NumCalStarsUT = 3;
#endif
  Ts1 = 6;
  Ts2 = 10;
  Ts3 = 17;

  CelestialRefTy
    CalStar1 = SkyAlignStarASt[0],
    CalStar2 = SkyAlignStarASt[1],
    CalStar3 = SkyAlignStarASt[2];
  double
    dist1,
    dist2,
    dist3,
    rise1 = 0.0, run1 = 0.0,
    rise2 = 0.0, run2 = 0.0,
    rise3 = 0.0, run3 = 0.0;

  #ifdef DAN_DBG
  char
    path[1024],
    fname[32];
  //getdocsdir ( path, sizeof ( path ) );
  strlcpy(path, sneakypath, sizeof(path));
  GenSkyAlignFileName(fname, 1024);
  strlcat ( path, fname, sizeof ( path ) );
  FILE *fp = fopen ( path, "wb" );
  if ( fp != NULL ) {
    fwrite (&SkyAlignStarASt, sizeof ( CelestialRefTy ), 3, fp );
  }
  fclose ( fp );
  #endif

  dist1 = Nxp_CalcDistAzmAlt(CalStar1, CalStar2);
  dist2 = Nxp_CalcDistAzmAlt(CalStar2, CalStar3);
  dist3 = Nxp_CalcDistAzmAlt(CalStar3, CalStar1);
  CalcRiseRun(CalStar1, CalStar2, &rise1, &run1);
  CalcRiseRun(CalStar2, CalStar3, &rise2, &run2);
  CalcRiseRun(CalStar3, CalStar1, &rise3, &run3);

  Nxp_MdlReset();
  Nxp_CalcIsNamedStarUp();
  bool matched = false;
#ifdef DAN_DBGPRINT
  printf("bydist, byriserun, byname, ");
  printf("star1, star2, star3, rmserr, , e12, e13, e23, calcerr, ,");
  printf("riseerr1, riseerr2, riseerr3, runerr1, runerr2, runerr3, riserunerr, , ");
  printf("disterr1, disterr2, disterr3, disterr, , ");
  printf("s1, name1, ");
  printf("ra1, dec1, azm1, alt1, lst1, jd1, ");
  printf("s2, name1, ");
  printf("ra2, dec2, azm2, alt2, lst2, jd2," );
  printf("s3, name3, ");
  printf("ra3, dec3, azm3, alt3, lst3, jd3\n");
#endif

  int
    best1 = 0, best2 = 0, best3 = 0;
  double
    besterr = 99999,
    bestcalcerr = 99999,
    bestdisterr = 99999,
    bestriserunerr = 99999;
  bool
    firstmatch = true,
    secondmatch = false;

  int s1, s2, s3;
  for(s1 = 0; !matched && s1 < NUM_NAMED_STARS; s1++) {
    if(NamedDbStar1[s1].IsUsedB) {
      for(s2 = 0; !matched && s2 < NUM_NAMED_STARS; s2++) {
        if(NamedDbStar2[s2].IsUsedB) {
          for(s3 = 0 ; !matched && s3 < NUM_NAMED_STARS; s3++) {
            if(NamedDbStar3[s3].IsUsedB && s1 != s2 && s2 != s3 && s1 != s3) {
              AlignStarTy
                star1 = NamedDbStar1[s1],
                star2 = NamedDbStar2[s2],
                star3 = NamedDbStar3[s3];
              bool
                bydist = false,
                byriserun = false,
                byname = false,
                bycalcerr = false;

              if(s1 == Ts1 && s2 == Ts2 && s3 == Ts3) {
                byname = true;
              }

              double
                dist1b = Nxp_CalcSeparationRaDec(star1, star2),
                dist2b = Nxp_CalcSeparationRaDec(star2, star3),
                dist3b = Nxp_CalcSeparationRaDec(star3, star1),
                rise1b = 0.0, run1b = 0.0,
                rise2b = 0.0, run2b = 0.0,
                rise3b = 0.0, run3b = 0.0;

              CalcRiseRun(star1, star2, &rise1b, &run1b);
              CalcRiseRun(star2, star3, &rise2b, &run2b);
              CalcRiseRun(star3, star1, &rise3b, &run3b);

              double
                calcerr = 9999,
                disterr1 = dist1 - dist1b,
                disterr2 = dist2 - dist2b,
                disterr3 = dist3 - dist3b,
                riseerr1 = rise1 - rise1b,
                riseerr2 = rise2 - rise2b,
                riseerr3 = rise3 - rise3b,
                runerr1 = fabs(run1) - fabs(run1b),
                runerr2 = fabs(run2) - fabs(run2b),
                runerr3 = fabs(run3) - fabs(run3b),
                riserunerr,
                disterr = fabs(disterr1) + fabs(disterr2) + fabs(disterr3);


              if(
                fabs(disterr1) < DIST_THRESH &&
                fabs(disterr2) < DIST_THRESH &&
                fabs(disterr3) < DIST_THRESH
              ) {
                bydist = false;  //DJM: disabled for now
              }

              if(
                BestOf5(riseerr1, riseerr2, riseerr3, runerr1, runerr2, runerr3, RISERUN_THRESH1, &riserunerr) &&
                riserunerr < RISERUN_THRESH2
              ) {
                byriserun = true;
              }

              double e12, e13, e23;
              e12 = CalcError(CalStar1, star1, CalStar2, star2);
              e13 = CalcError(CalStar1, star1, CalStar3, star3);
              e23 = CalcError(CalStar2, star2, CalStar3, star3);
              calcerr = fabs(e12) + fabs(e13) + fabs(e23);

              if(fabs(e12) < CALCERR_THRESH && fabs(e13) < CALCERR_THRESH && fabs(e23) < CALCERR_THRESH) {
                bycalcerr = true;
              }

              if (byriserun || bycalcerr || byname) {
                SkyAlignStarASt[0].ra = star1.RaRadD;
                SkyAlignStarASt[0].dec = star1.DecRadD;
                SkyAlignStarASt[1].ra = star2.RaRadD;
                SkyAlignStarASt[1].dec = star2.DecRadD;
                SkyAlignStarASt[2].ra = star3.RaRadD;
                SkyAlignStarASt[2].dec = star3.DecRadD;
                Nxp_MdlReset();
                Nxp_MdlCalStarAdd(SkyAlignStarASt[0]);
                Nxp_MdlCalStarAdd(SkyAlignStarASt[1]);
                Nxp_MdlCalStarAdd(SkyAlignStarASt[2]);
                Nxp_MdlUpdate();
#ifdef DAN_DBGPRINT
                printf("%s, %s, %s %s, ",bydist == true ? "d" : "-", byriserun == true ? "r" : "-", bycalcerr == true ? "c" : "-", byname == true ? "n" : "-");
                printf("%s, %s, %s, %1.5f, , %1.5f, %1.5f, %1.5f, %1.5f, , ", NamedStarsDb[s1].Name, NamedStarsDb[s2].Name, NamedStarsDb[s3].Name, NxpModelSt.RmsErrorRadD, e12, e13, e23, calcerr);
                printf("%1.5f, %1.5f, %1.5f, %1.5f, %1.5f, %1.5f, %1.5f, , ", riseerr1, riseerr2, riseerr3, runerr1, runerr2, runerr3, riserunerr);
                printf("%1.5f, %1.5f, %1.5f, %1.5f, , ", disterr1, disterr2, disterr3, disterr);
                printf("%03d, %s, ", s1, NamedStarsDb[s1].Name);
                printf("%1.8f, %1.8f, %1.8f, %1.8f, %1.8f, %1.4f, , ", SkyAlignStarASt[0].ra, SkyAlignStarASt[0].dec, SkyAlignStarASt[0].azm, SkyAlignStarASt[0].alt, SkyAlignStarASt[0].lst, SkyAlignStarASt[0].jd );
                printf("%03d, %s, ", s2, NamedStarsDb[s2].Name);
                printf("%1.8f, %1.8f, %1.8f, %1.8f, %1.8f, %1.4f, , ", SkyAlignStarASt[1].ra, SkyAlignStarASt[1].dec, SkyAlignStarASt[1].azm, SkyAlignStarASt[1].alt, SkyAlignStarASt[1].lst, SkyAlignStarASt[1].jd );
                printf("%03d, %s, ", s3, NamedStarsDb[s3].Name);
                printf("%1.8f, %1.8f, %1.8f, %1.8f, %1.8f, %1.4f, ,\n", SkyAlignStarASt[2].ra, SkyAlignStarASt[2].dec, SkyAlignStarASt[2].azm, SkyAlignStarASt[2].alt, SkyAlignStarASt[2].lst, SkyAlignStarASt[2].jd );
#endif
                if(NxpModelSt.RmsErrorRadD < 0.0002908 * 5.5) { // approx 5.5'
                  if(NxpModelSt.RmsErrorRadD < besterr || secondmatch) {
                    bool newbest = false;
                    int Score = 0;
                    if (calcerr < bestcalcerr) Score++;
                    if (disterr < bestdisterr) Score++;
                    if (riserunerr < bestriserunerr) Score++;
                    if (Score > 1) {
                      newbest = true;
                    }
                    if(newbest) {
                      besterr = NxpModelSt.RmsErrorRadD;
                      bestcalcerr = calcerr;
                      bestdisterr = disterr;
                      bestriserunerr = riserunerr;
                      best1 = s1;
                      best2 = s2;
                      best3 = s3;
                      sprintf(AlignedStars, "%s, %s, %s, %0.2f'", NamedStarsDb[s1].Name, NamedStarsDb[s2].Name, NamedStarsDb[s3].Name, NxpModelSt.RmsErrorRadD * 180 * 60 / PI_D);
                    }
                    if(firstmatch) {
                      firstmatch = false;
                      secondmatch = true;
                    }
                    else if(secondmatch) {
                      secondmatch = false;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (besterr < 0.0002908 * 5.5) { // ~10' {
    //Console.WriteLine("^^^^^^^ MATCHED! ^^^^^^^");
    SkyAlignStarASt[0].ra = NamedStarsDb[best1].RaRadD;
    SkyAlignStarASt[0].dec = NamedStarsDb[best1].DecRadD;
    SkyAlignStarASt[1].ra = NamedStarsDb[best2].RaRadD;
    SkyAlignStarASt[1].dec = NamedStarsDb[best2].DecRadD;
    SkyAlignStarASt[2].ra = NamedStarsDb[best3].RaRadD;
    SkyAlignStarASt[2].dec = NamedStarsDb[best3].DecRadD;
    Nxp_MdlReset();
    Nxp_MdlCalStarAdd(SkyAlignStarASt[0]);
    Nxp_MdlCalStarAdd(SkyAlignStarASt[1]);
    Nxp_MdlCalStarAdd(SkyAlignStarASt[2]);
    Nxp_MdlUpdate();
  }
  else {
    //Nxp_Init();
  }

  return besterr;
}



NameStarDbTy NamedStarsDb[NUM_NAMED_STARS] = {
  { "Achernar", 0.4263621200, -0.9989682861, 0.50 },
  { "Acrux", 3.2576483222, -1.1012882140, 1.30 },
  { "Adhara", 1.8265961453, -0.5056605720, 1.50 },
  { "Aldebaran", 1.2039281180, 0.2881393150, 0.90 },
  { "Alderamin", 5.5788576875, 1.0923239120, 2.40 },
  { "Algieba", 2.7051381675, 0.3463024120, 2.60 },
  { "Algol", 0.8210377871, 0.7148091950, 2.10 },
  { "Alhena", 1.7353459690, 0.2862194530, 1.90 },
  { "Alioth", 3.3773342757, 0.9766813040, 1.80 },
  { "Alkaid", 3.6108244230, 0.8606800319, 1.90 },
  { "Almach", 0.5406157360, 0.7387929270, 2.30 },
  { "Alnair", 5.7955097709, -0.8196261060, 1.70 },
  { "Alnilam", 1.4670059600, -0.0209779850, 1.70 },
  { "Alnitak", 1.4868372630, -0.0339079660, 2.10 },
  { "Alpha Centauri", 3.8379877518, -1.0617806499, 1.30 },
  { "Alphard", 2.4765654854, -0.1511213700, 2.00 },
  { "Alphecca", 4.0788724228, 0.4660075310, 2.20 },
  { "Alpheratz", 0.0365995540, 0.5077258791, 2.10 },
  { "Alsuhail", 2.3910879853, -0.7580401270, 2.20 },
  { "Altair", 5.1957710067, 0.1547816160, 0.80 },
  { "Aludra", 1.9377299854, -0.5114347030, 2.50 },
  { "Ankaa", 0.1146812220, -0.7383810300, 2.40 },
  { "Antares", 4.3171024481, -0.4613245550, 1.00 },
  { "Arcturus", 3.7335297960, 0.3347977839, 0.00 },
  { "Arneb", 1.4518085050, -0.3110563611, 2.60 },
  { "Ascella", 4.9855853682, -0.5215093261, 2.60 },
  { "Atria", 4.4011313249, -1.2047620950, 1.90 },
  { "Avior", 2.1926265960, -1.0386404930, 1.90 },
  { "Bellatrix", 1.4186559761, 0.1108234620, 1.60 },
  { "Betelgeuse", 1.5497302031, 0.1292756650, 0.50 },
  { "Canopus", 1.6753066421, -0.9197157940, -0.70 },
  { "Capella", 1.3818208020, 0.8028174220, 0.10 },
  { "Caph", 0.0400465800, 1.0323573080, 2.30 },
  { "Castor", 1.9835666949, 0.5565564099, 2.00 },
  { "Deneb", 5.4167689599, 0.7902900300, 1.30 },
  { "Denebola", 3.0938578985, 0.2543285059, 2.10 },
  { "Diphda", 0.1901972551, -0.3139265550, 2.00 },
  { "Dschubba", 4.1902431914, -0.3948225660, 2.30 },
  { "Dubhe", 2.8960597344, 1.0777553580, 1.80 },
  { "El Nath", 1.4237174311, 0.4992950659, 1.70 },
  { "Eltanin", 4.6975842250, 0.8986505419, 2.20 },
  { "Enif", 5.6905893029, 0.1723512640, 2.40 },
  { "Fomalhaut", 6.0111408367, -0.5170052130, 1.20 },
  { "Gacrux", 3.2775756189, -0.9968157139, 1.60 },
  { "Gemma", 4.0783457697, 0.4662597650, 2.20 },
  { "Gienah", 3.2105637023, -0.3061647851, 2.60 },
  { "Graffias", 4.2125135926, -0.3456720580, 2.60 },
  { "Hadar", 3.6818724135, -1.0537085020, 0.60 },
  { "Hamal", 0.5548968921, 0.4094978759, 2.00 },
  { "Izar", 3.8614842468, 0.4725333510, 2.70 },
  { "Kaus Australis", 4.8178592271, -0.6001265181, 1.90 },
  { "Kochab", 3.8864337285, 1.2942585060, 2.10 },
  { "Kraz", 3.2916342461, -0.4083488670, 2.70 },
  { "Markab", 6.0421640641, 0.2653822580, 2.50 },
  { "Menkalinan", 1.5687368379, 0.7844818661, 1.90 },
  { "Menkar", 0.7953465400, 0.0713790210, 2.50 },
  { "Menkent", 3.6943515177, -0.6347762490, 2.10 },
  { "Merak", 2.8878305070, 0.9840602660, 2.40 },
  { "Miaplacidus", 2.4137903555, -1.2167949760, 1.70 },
  { "Mimosa", 3.3498104334, -1.0417628870, 1.30 },
  { "Mintaka", 1.4486538221, -0.0052359880, 2.20 },
  { "Mirach", 0.3042632490, 0.6216958790, 2.10 },
  { "Mirfak", 0.8915272721, 0.8702406550, 1.80 },
  { "Mirzam", 1.6698437620, -0.3133884121, 2.00 },
  { "Mizar", 3.5077845473, 0.9586270370, 2.30 },
  { "Muphrid", 3.6420003668, 0.3211018940, 2.70 },
  { "Naos", 2.1100376158, -0.6981898780, 2.30 },
  { "Navi", 0.2474353281, 1.0597057441, 2.50 },
  { "Nunki", 4.8039969495, -0.5205976819, 2.70 },
  { "Peacock", 5.3478982676, -0.9902125509, 1.90 },
  { "Phact", 1.4819939750, -0.5947063979, 2.60 },
  { "Phad", 3.1146709499, 0.9371495970, 2.40 },
  { "Polaris", 0.6624048110, 1.5579536119, 2.00 },
  { "Pollux", 2.0303197022, 0.4891480119, 1.10 },
  { "Procyon", 2.0040830402, 0.0911934530, 0.40 },
  { "Rasalhague", 4.6030222861, 0.2192133541, 2.10 },
  { "Regor", 2.1359906618, -0.8261806900, 1.80 },
  { "Regulus", 2.6545236192, 0.2088673330, 1.40 },
  { "Rigel", 1.3724309310, -0.1431460880, 0.10 },
  { "Rukbah", 0.3744473010, 1.0513040200, 2.70 },
  { "Sabik", 4.4958721600, -0.2744480800, 2.40 },
  { "Sadr", 5.3329757025, 0.7026113789, 2.20 },
  { "Saiph", 1.5173761621, -0.1687683940, 2.10 },
  { "Scheat", 6.0378574642, 0.4901370321, 2.40 },
  { "Schedar", 0.1767494930, 0.9867605810, 2.20 },
  { "Scutulum", 2.4307636825, -1.0345488590, 2.30 },
  { "Shaula", 4.5972321563, -0.6475849290, 1.60 },
  { "Sheratan", 0.5002113640, 0.3631689831, 2.60 },
  { "Sirius", 1.7677916395, -0.2917512740, -1.50 },
  { "Spica", 3.5133172410, -0.1948028880, 1.00 },
  { "Suhail", 2.3916576171, -0.7582606351, 2.20 },
  { "Tarazed", 5.1760356962, 0.1852376110, 2.70 },
  { "Unukalhai", 4.1201464053, 0.1121470040, 2.70 },
  { "Vega", 4.8735614102, 0.6769018071, 0.00 },
  { "Wezen", 1.8692127223, -0.4606505670, 1.80 },
  { "Yed Prior", 4.2513820751, -0.0644803170, 2.70 },
  { "Zosma", 2.9413510285, 0.3582046850, 2.60 },
  { "Zubeneshamali", 4.0011978535, -0.1637651160, 2.60 },
  { "Mercury", 0.0, 0.0, 0.0 },
  { "Venus", 0.0, 0.0, 0.0 },
  { "Mars", 0.0, 0.0, 0.0 },
  { "Jupiter", 0.0, 0.0 },
  { "Saturn", 0.0, 0.0 }
};
