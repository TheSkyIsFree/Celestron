/*******************************************************************************************************************************
** FILE: lisasphr.c
** DATE: 2010.02.23 DJM
** DESC: Implements spherical projection routines used in the LISA routines
**
** NOTE: Copywrite (c) 2009 by Celestron.
*******************************************************************************************************************************/
#define __SKYSPHR_C__
#include "lisasphr.h"                                          // interface to LISA spherical projection routines



//==============================================================================================================================
// global variable declarations and subsystem function prototypes
//==============================================================================================================================
ProjectionType
  Projection;



//==============================================================================================================================
// static variable declarations
//==============================================================================================================================



/*******************************************************************************************************************************
********                                       static function implementations                                          ********
*******************************************************************************************************************************/



/*******************************************************************************************************************************
********                                      global function implementations                                           ********
*******************************************************************************************************************************/

//==============================================================================================================================
// FUNC: LisaSphr_SetupProjection()
//==============================================================================================================================
// DESC:
//==============================================================================================================================
void LisaSphr_SetupProjection(
  double
    RaRad,
  double
    DecRad,
  double
    PScale
)
{
  Projection.RaRad = RaRad;
  Projection.DecRad = DecRad;
  Projection.CosRa = cos(RaRad);
  Projection.SinRa = sin(RaRad);
  Projection.CosDec = cos(DecRad);
  Projection.SinDec = sin(DecRad);
  Projection.PScale = PScale;
}



//==============================================================================================================================
// FUNC: LisaSphr_Project()
//==============================================================================================================================
// DESC: given the rotation elements of the projection this routine performs the rotation of a field point at (RaRad, DecRad)
//       and returns the (x,y) vector location of the projected field point rotation carries the center of the field to the
//       +X-axis (X=1,Y=0,Z=0) if OrRad=0, then North is the positive z-direction, otherwise field is rotated by OrRad in a CCW
//       direction the routine returns the Z-value of the projection
//==============================================================================================================================
double LisaSphr_Project(
  double
    RaRad,
  double
    DecRad,
  double
    *X,
  double
    *Y
)
{
  Vect3Type
    R,
    R1,
    R2;
  double
    C;

  C = cos(DecRad);                                              // convert field point to (x,y,z) coords
  R.X = cos(RaRad) * C;
  R.Y = sin(RaRad) * C;
  R.Z = sin(DecRad);

  R1.X = Projection.CosRa * R.X + Projection.SinRa * R.Y;       // rotate field point about Z-axis so that projection is at Y=0
  R1.Y = -Projection.SinRa * R.X + Projection.CosRa * R.Y;
  R1.Z = R.Z;

  R2.X = Projection.CosDec * R1.X + Projection.SinDec * R1.Z;
  R2.Y = R1.Y;
  R2.Z = -Projection.SinDec * R1.X + Projection.CosDec * R1.Z;

  *X = -R2.Y * Projection.PScale;                               // the above are the three-vector coors of the field point that
  *Y = R2.Z * Projection.PScale;                                // project onto 2D surface

  return R2.X;                                                  // the Z-value of projection should be -1 to 1 (is also fit val)
}



//==============================================================================================================================
// FUNC: LisaSphr_DeProject()
//==============================================================================================================================
// DESC: deprojects from catalog plate coordinates (X,Y) to (RA,Dec) projection must be set up before calling
//==============================================================================================================================
void LisaSphr_DeProject(
  double
    X,
  double
    Y,
  double
    *RaRad,
  double
    *DecRad
)
{
  Vect3Type
    R1,
    R2,
    R3;
  double
    Rsq,
    Rho,
    Theta,
    Phi;

  R1.X = 0;
  R1.Y = -X / Projection.PScale;
  R1.Z = Y / Projection.PScale;

  Rsq = 1 - (R1.Y * R1.Y + R1.Z * R1.Z);
  if(Rsq >= 0) {
    R1.X = sqrt(Rsq);
  }

  R2.X = Projection.CosDec * R1.X - Projection.SinDec * R1.Z;
  R2.Y = R1.Y;
  R2.Z = Projection.SinDec * R1.X + Projection.CosDec * R1.Z;

  R3.X = Projection.CosRa * R2.X - Projection.SinRa * R2.Y;
  R3.Y = Projection.SinRa * R2.X + Projection.CosRa * R2.Y;
  R3.Z = R2.Z;

  Rho = sqrt(R3.X * R3.X + R3.Y * R3.Y);
  Theta = atan2(R3.Z, Rho);
  *DecRad = Theta;
  Phi = atan2(R3.Y, R3.X);
  if(Phi < 0) Phi += 2 * PI_D;
  if(Phi > 2 * PI_D) Phi -= 2 * PI_D;
  *RaRad = Phi;
}



