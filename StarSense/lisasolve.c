/*******************************************************************************************************************************
** FILE: lisasolve.c
** DATE: 2010.02.17 DJM
** DESC: Implements the LISA for NexAlign
**
** NOTE: Copywrite (c) 2009 by Celestron.
*******************************************************************************************************************************/
#define LISASOLVE_
#include "lisasolve.h"                                          // LISA functions and data types interface



//==============================================================================================================================
// global variable declarations and subsystem function prototypes
//==============================================================================================================================
LisaPlateStarTy
  Plate[100];                                                   // 40 * 100 = 4000 -- list of extracted stars

int32_t
  LisaCtrXRefSL,
  LisaCtrYRefSL;

uint32_t
  PlateCountUT;

double
  LisaPlateScale,
  LisaPlateRaRad,
  LisaPlateDecRad;

int
  NumCatalog;


//==============================================================================================================================
// static variable declarations
//==============================================================================================================================
#define MAX_CAT_LISA 50000
LisaApmTy
  CatLisa[MAX_CAT_LISA];                                                     //

#define MAX_QUAD_LISA 50000
LisaQuadTy
  QuadLisa[MAX_QUAD_LISA];                                                    //

LisaStarTy
  PlateList[100],                                               // 40 * 100 = 4000  -- list of matched stars from plate
  CatList[100],                                                 // 40 * 100 = 4000  -- list of matched catalog stars in field?
  Catalog[100];                                                 // 40 * 100 = 4000  -- list of potential catalog stars in field?
LisaMatchListTy
  MatchList[100];                                               // 16 * 100 = 1600

LisaAffTransformTy
  Aff;

LisaGeoTransformTy
  Geo;

int
  NumCatalog,
  NumPlate,
  NumMatch,
  NumSavedStars,
  NCatLisa,
  NQuadLisa;

double
  FitError;

BOOL_D
  FoundFit;

const double
  XSizeRad = 6.8 * 3.14159265358 / 180,
  YSizeRad = 5.1 * 3.14159265358 / 180;
  //MaxError = 4.0;



/*******************************************************************************************************************************
********                                       static function implementations                                          ********
*******************************************************************************************************************************/

//==============================================================================================================================
// FUNC: (static) LisaSolve_SwapStars()
//==============================================================================================================================
// DESC: swaps two stars within the star list
//==============================================================================================================================
static void LisaSolve_SwapStars(LisaStarTy *A, LisaStarTy *B)
{
  LisaStarTy
    T;
  memcpy(&T, A, sizeof(LisaStarTy));
  memcpy(A, B, sizeof(LisaStarTy));
  memcpy(B, &T, sizeof(LisaStarTy));
}



//==============================================================================================================================
// FUNC: (static) LisaSolve_SortList()
//==============================================================================================================================
// DESC: sorts the list of stars based on brightness
//==============================================================================================================================
static void LisaSolve_SortList(LisaStarTy *StarList, uint32_t N)
{
  uint32_t
    Split,
    High,
    Low;

  Split = N / 2;
  while(Split) {
    High = Split;
    while(High <= N) {
      Low = High - Split;
      if(StarList[Low].Intensity > StarList[High].Intensity) {
        LisaSolve_SwapStars(&StarList[Low], &StarList[High]);
        if(Split == 1) {
          if(Low > 0) High = Low - 1;
        }
      }
      High++;
    }
    Split = (float)Split / 1.3;
  }

  Low = 1;
  High = N;
  while(Low < High) {
    LisaSolve_SwapStars(&StarList[Low], &StarList[High]);
    Low++; High--;
  }
}



//==============================================================================================================================
// FUNC: (static) LisaSolve_ScanAPMFile()
//==============================================================================================================================
// DESC: Sets the pointer to APM star database, calculates entries
//==============================================================================================================================
void LisaSolve_ScanAPMFile(
	FILE
		*StarFileP,
	const uint32_t
		StarFileOffsetUL,
	const uint32_t
		StarFileSizeUL
)
{
	int i = 0;
	unsigned char starBuffer[10]; // 4 + 4 + 2
	NCatLisa = (int)(StarFileSizeUL / 10);
	fseek ( StarFileP, StarFileOffsetUL, SEEK_SET );
	while (i < MAX_CAT_LISA && i < NCatLisa && fread(starBuffer, 1, 10, StarFileP) > 0)
	{
		memcpy(&CatLisa[i], starBuffer, 10);
		i++;
	}
	fclose(StarFileP);
}



//==============================================================================================================================
// FUNC: (static) LisaSolve_ScanQuadFile()
//==============================================================================================================================
// DESC: Sets the pointer to QUAD star database, calculates entries
//==============================================================================================================================
void LisaSolve_ScanQuadFile(
	FILE
		*QuadFileP,
	const uint32_t
		QuadFileOffsetUL,
	const uint32_t
		QuadFileSizeUL
)
{
	int i = 0;
	unsigned char quadBuffer[40]; // 4 * 4 + 4 * 6
	NQuadLisa = (int)((QuadFileSizeUL) / 40);
	fseek ( QuadFileP, QuadFileOffsetUL, SEEK_SET );
	while (i < MAX_QUAD_LISA && i < NQuadLisa && fread(quadBuffer, 1, 40, QuadFileP) > 0)
	{
		memcpy(&QuadLisa[i], quadBuffer, 40);
		i++;
	}
	fclose(QuadFileP);
}



//==============================================================================================================================
// FUNC: (static) LisaSolve_QuadDistance()
//==============================================================================================================================
// DESC: use LISA to solve plate
//==============================================================================================================================
static double LisaSolve_QuadDistance(
  LisaQuadTy
    Q1,
  LisaQuadTy
    Q2
)
{
  int
    i;

  double
    QuadDistance;

  QuadDistance = 0;
  for(i=1; i<=5; i++) {
    QuadDistance += pow((Q1.D[i] - Q2.D[i]), 2);
  }
  QuadDistance = sqrt(QuadDistance);
  return QuadDistance;
}



//==============================================================================================================================
// FUNC: (static) LisaSolve_FindStartQ()
//==============================================================================================================================
// DESC: finds the starting index in Q[] where D0 is greater than Q[Index - 1] and D0 is <= to Q[Index] and return index. This
//       gives us the starting point in quad database. NOTE: uses binary search algorithm
//==============================================================================================================================
static int LisaSolve_FindStartQ(
  LisaQuadTy
    *Q,
  int
    N,
  double
    D0
)
{
  int
    I0,
    I1,
    I2;

  I0 = 0;
  I1 = 0;
  I2 = N - 1;

  if(D0 <= Q[0].D[0]) {
    return 1;
  }
  if(D0 > Q[N-1].D[0]) {
    return N;
  }
  while(I0 + 1 < I2) {
    I1 = (I0 + I2) / 2;
    if (D0 >= Q[I1].D[0]) I0 = I1;
    if (D0 < Q[I1].D[0]) I2 = I1;
  }
  return I1;
}



//==============================================================================================================================
// FUNC: (static) LisaSolve_PlateDistance()
//==============================================================================================================================
// DESC: use LISA to solve plate
//==============================================================================================================================
static double LisaSolve_PlateDistance(
  LisaPlateStarTy
    *P1,
  LisaPlateStarTy
    *P2
)
{
  double
    PlateDistance;
  PlateDistance = sqrt(
    pow((P1->Xcen - P2->Xcen), 2) +
    pow((P1->Ycen - P2->Ycen), 2)
  );
  return PlateDistance;
}



//==============================================================================================================================
// FUNC: (static) LisaSolve_MakePlateQuad()
//==============================================================================================================================
// DESC: use LISA to solve plate
//==============================================================================================================================
static void LisaSolve_MakePlateQuad(
  LisaPlateStarTy
    *P1,
  LisaPlateStarTy
    *P2,
  LisaPlateStarTy
    *P3,
  LisaPlateStarTy
    *P4,
  LisaQuadTy
    *MakePlateQuad
)
{
  double
    D[6];

  D[0] = LisaSolve_PlateDistance(P1, P2);
  D[1] = LisaSolve_PlateDistance(P1, P3);
  D[2] = LisaSolve_PlateDistance(P1, P4);
  D[3] = LisaSolve_PlateDistance(P2, P3);
  D[4] = LisaSolve_PlateDistance(P2, P4);
  D[5] = LisaSolve_PlateDistance(P3, P4);

  {
    double
      Tmp;
    int
      i,j;
    for(i = 1; i < 6; i++) {
      Tmp = D[i];
      j = i - 1;
      while(Tmp > D[j] && j >= 0) {
        D[j + 1] = D[j];
        j = j - 1;
      }
      D[j + 1] = Tmp;
    }
  }

  //MakePlateQuad->D[0] = D[0];
  MakePlateQuad->D[1] = D[1] / D[0];
  MakePlateQuad->D[2] = D[2] / D[0];
  MakePlateQuad->D[3] = D[3] / D[0];
  MakePlateQuad->D[4] = D[4] / D[0];
  MakePlateQuad->D[5] = D[5] / D[0];
  MakePlateQuad->D[0] =
    MakePlateQuad->D[1] +
    MakePlateQuad->D[2] +
    MakePlateQuad->D[3] +
    MakePlateQuad->D[4] +
    MakePlateQuad->D[5];
}



//==============================================================================================================================
// FUNC: (static) LisaSolve_MakePlate()
//==============================================================================================================================
// DESC: use LISA to solve plate
//'   routine creates a catalog, Cat(), of all catalog stars that are within RadiusRad of
//'   (RACenterRad, DecCenterRad).  The catalog has NumCat members and is of StarListType
//'   projected catalog is centered on RACenterRad,DecCenterRad
//'   the plate scale is 1 arcsecond per pixel
//'   maximum magnitude is MaxMag
//'   .Intensity is set to Vmag
//'   the routine returns the scale used in the projection, PScale = ASPERRAD
//==============================================================================================================================
static void LisaSolve_MakePlate(
  double
    RACenterRad,
  double
    DecCenterRad,
  double
    RadiusRad,
  double
    MaxMag,
  LisaStarTy
    *Catalog,
  int
    *NumCatalogPtr
)
{
  uint32_t
    OldPcntUL,
    PcntUL;
  double
    RadiusAS,
    PScale,
    RARad,
    DecRad,
    Vmag,
    X,
    Y,
    Dist,
    Z=-1;
  int
    i;

  PScale = ASPERRAD;
  RadiusAS = RadiusRad * PScale;
  RadiusAS *= RadiusAS;

  LisaSphr_SetupProjection(RACenterRad, DecCenterRad, PScale);  // setup for projection based on RA and DEC of center based on
  NumCatalog = 0;                                               // the matched QUAD
  X = 0;
  Y = 0;
  OldPcntUL = 9999;

  for(i = 0; i < NCatLisa && NumCatalog < 30; i++) {            // for: each star in image, up to 8 total matched
    PcntUL = (NumCatalog * 100 + 50) / 30;                      //
    if(PcntUL != OldPcntUL) {                                   //
      OldPcntUL = PcntUL;                                       //
      ////Lcd_GotoXY(11, 1);                                      //
      //Lcd_WrStrLine(1, "Solving... ", LCD_CLR_LINE);            //
      //Lcd_WrDec(PcntUL, 3, '%');                                //
    }                                                           //
    //if(                                                         //
    //  (*Keypad_KeyValPUD & KEYPAD_BACK) ||                      //
    //  LisaAbortSolveB == TRUE                                   //
    // ) {                                                        //   check for abort from keypad (this is most heavily exucuted
    //  LisaAbortSolveB = TRUE;                                   //
    //  break;                                                    //
    //}                                                           //
    Vmag = (double)CatLisa[i].Mag / 100;                        //   only compare stars from catalog below brightness threshold
    if(Vmag <= MaxMag) {                                        //   if: within mag range
      RARad = CatLisa[i].RARad;                                 //     attempt to match image plate start to projected plate star
      DecRad = CatLisa[i].DecRad;                               //
      Z = LisaSphr_Project(RARad, DecRad, &X, &Y);              //
      Dist = X * X + Y * Y;                                     //
      if(Dist <= RadiusAS && Z >= 0 ) {                         //     if: distance is within threshold
        Catalog[NumCatalog].Xcen = X;                           //       create catalog of matched stars
        Catalog[NumCatalog].Ycen = Y;                           //
        Catalog[NumCatalog].Intensity = Vmag;                   //
        Catalog[NumCatalog].RARad = RARad;                      //
        Catalog[NumCatalog].DecRad = DecRad;                    //
        NumCatalog++;                                           //
      }                                                         //
    }                                                           //
  }                                                             //
  if(LisaAbortSolveB) return;                           //
  LisaSolve_SortList(Catalog, (uint32_t)NumCatalog);                // Sort Catalog by brightness
  *NumCatalogPtr = NumCatalog;                                  //
}



/*******************************************************************************************************************************
********                                      global function implementations                                           ********
*******************************************************************************************************************************/

//==============================================================================================================================
// FUNC: LisaSolve_SphericalProject()
//==============================================================================================================================
// DESC: use LISA to project Ra, Dec to x,y coords
//==============================================================================================================================
void LisaSolve_SphericalProject(                                 //
  double                                                        //
    RaRadF,                                                     //
  double                                                        //
    DecRadF,                                                    //
  double                                                        //
    *Xcen,                                                      //
  double                                                        //
    *Ycen                                                       //
)                                                               //
{
  double
    U1, V1, U2, V2,
    TempRaF, TempDecF;

  LisaSphr_Project(RaRadF, DecRadF, &U1, &V1);                  //
  LisaAff_Transform(Aff, U1, V1, Xcen, Ycen);                   //
  LisaAff_TransformInv(Aff, *Xcen, *Ycen, &U2, &V2);            //
  LisaSphr_DeProject(U2, V2, &TempRaF, &TempDecF);              //

  if(TempRaF - RaRadF > 0.01 || TempDecF - DecRadF > 0.1) {
    FoundFit = 0;                                           //
  }
}



//==============================================================================================================================
// FUNC: LisaSolve_Solve()
//==============================================================================================================================
// DESC: use LISA to solve plate
//==============================================================================================================================
BOOL_D LisaSolve_Solve(                                           //
  void                                                          //
)                                                               //
{
  //uint8_t                                                         //
  //  SBufAUB[40];                                                //
  double // 142 * 8 = 1136                                      //
    DistError[120],                                             //
    X=0, Y=0,                                                   //
    U=0, V=0,                                                   //
    //Icen, Jcen,                                                 //
    U1 = 0, V1 = 0,                                             //
    Dist,                                                       //
    RARad, DecRad,                                              //
    RARad1, DecRad1,                                            //
    RadiusRad, Norm, MaxMag,                                    //
    PScale,                                                     //
    //Px,Py,                                                      //
    MaxDistDeg,                                                 //
    MaxDist,                                                    //
    MaxPlateDist,                                               //
    QDist, Eps;                                                 //
  int // 15 * 4 =                                               //
    iCat,                                                       //
    iPlate,                                                     //
    i, j,                                                       //
    k, k1, k2, k3, k4,                                          //
    NStars,                                                     //
    NumSearch,                                                  //
    PlateIndex[4];                                              //
  BOOL_D                                                         //
    UsedPlateStar[120];                                         //
  LisaQuadTy                                                    //
    QCat,                                                       //
    QPlate;                                                     //

  RadiusRad =                                                   //
    sqrt(XSizeRad * XSizeRad + YSizeRad * YSizeRad) / 2;        // search radius
  Norm = (RadiusRad * (180 / PI_D) * 60) / 14;                    // significance of 14?
  Norm = 2.5 * log10(Norm * Norm);                              //
  MaxMag = 15 - Norm;                                           // set the magnitude limit for catalog
  NumSearch = 0;                                                // reset LISA parameters
  NumPlate = 0;                                                 //
  //Px = 640;                                                     //
  //Py = 480;                                                     //
  NumPlate = PlateCountUT;                                      //

  //PScale = XSizeRad / Px;
  //PScale = PScale * ASPERRAD;                                 // the approximate plate scale in arcseconds/pix
  //PScale = ASperRAd / PScale;                                 // the approximate plate scale in pixels/radian
  //PScale = 640 / XSizeRad;
        PScale = 1280 / XSizeRad;
  MaxDistDeg = 1.5;                                             // maximum quad distance in degrees (a catalog constant)

  MaxDist = MaxDistDeg * DEG_TO_RAD_D;                            // in radians
  MaxPlateDist = MaxDist * PScale;                              // in pixels

  Eps = 0.01;                                                   // normalized K-space distance
  NumCatalog = 0;
  NumSavedStars = 0;
  FoundFit = 0;

  for(i = 0; i < NumPlate && !FoundFit; i++) {
    NStars = 0;
    j = 0;
    PlateIndex[0] = i;
    //find the 3 brightest stars around this star within MaxDist
    while(NStars < 3 && j < NumPlate) {
      Dist = LisaSolve_PlateDistance(&Plate[i], &Plate[j]);
      if(i != j && Dist < MaxPlateDist) {
        NStars++;
        PlateIndex[NStars] = j;
      }
      j = j + 1;
    }
    if (NStars == 3) {
      LisaSolve_MakePlateQuad(
        &Plate[PlateIndex[0]],
        &Plate[PlateIndex[1]],
        &Plate[PlateIndex[2]],
        &Plate[PlateIndex[3]],
        &QPlate
      );
      k =
/**/    LisaSolve_FindStartQ(
          QuadLisa,
          NQuadLisa,
          QPlate.D[0] - 2 * Eps
        );
      while(
        k < NQuadLisa &&
        QuadLisa[k].D[0] < QPlate.D[0] + 2 * Eps &&
        !FoundFit
      ) {
        QCat = QuadLisa[k];
/**/    QDist = LisaSolve_QuadDistance(QPlate, QCat);
        if(QDist < Eps) {
          NumSearch++;
          //================= Make plate fromthe four catalog stars =======================
          RARad1 = CatLisa[QCat.S[0] - 1].RARad;
          DecRad1 = CatLisa[QCat.S[0] - 1].DecRad;
/**/      LisaSphr_SetupProjection(RARad1, DecRad1, ASPERRAD);
          NumMatch = 4;
#if 0
          // "Plate Stars   : {0:g}, {1:g}, {2:g}, {3:g}\r", PlateIndex[0], PlateIndex[1], PlateIndex[2], PlateIndex[3]));
          // "Matched Quad  : {0:g} Qdist = {1:0.000000}\r", k, Qdist));
          // "Quad Star 1   : RA  {0:g}\r",GetRAStringFromRad(RARad1)));
          // "              : DEC {0:g}\r",GetDecStringFromRad(DecRad1)));
          Serial_WriteString("LISA found Quad match.\r");
          Serial_WriteString("  Plate Stars   : ");
          Serial_WriteDec(PlateIndex[0], 0, ',' << 8 | ' ');
          Serial_WriteDec(PlateIndex[1], 0, ',' << 8 | ' ');
          Serial_WriteDec(PlateIndex[2], 0, ',' << 8 | ' ');
          Serial_WriteDec(PlateIndex[3], 0, '\r');

          Serial_WriteString("  Matched Quad  : ");
          Serial_WriteDec(k, 0, ' ' << 8 | ' ');
          Serial_WriteString("  Qdist = ");
          sprintf((char *)SBufAUB, "%0.6f\r", QDist);
          Serial_WriteString(SBufAUB);
          Serial_WriteString("  Ra            :  ");
          Serial_WriteHms(RARad1, '\r');
          Serial_WriteString("  Dec           : ");
          Serial_WriteDms(DecRad1, TRUE, '\r');
#endif
          if(LisaAbortSolveB) return 0;
          //============== Setup Catalog on Plate Lists for Affine Fit ====================
          for(iCat = 0; iCat < 4; iCat++) {
            X = 0; Y = 0;
            RARad = CatLisa[QuadLisa[k].S[iCat] - 1].RARad;
            DecRad = CatLisa[QuadLisa[k].S[iCat] - 1].DecRad;
/**/        LisaSphr_Project(RARad, DecRad, &X, &Y);
            CatList[iCat].Xcen = X;
            CatList[iCat].Ycen = Y;
            CatList[iCat].RARad = RARad;
            CatList[iCat].DecRad = DecRad;
            MatchList[iCat].StarList2Index = iCat;
            MatchList[iCat].Used = 1;
            PlateList[iCat].Xcen = Plate[PlateIndex[iCat]].Xcen;
            PlateList[iCat].Ycen = Plate[PlateIndex[iCat]].Ycen;
          }
          //============== Match the four cat stars to the four plate stars ===============
          FoundFit = 0;
          for(k1 = 0; k1 < 4 && !FoundFit; k1++) {
            MatchList[0].StarList1Index = k1;
            for(k2 = 0; k2 < 4 && !FoundFit; k2++) {
              if(k2 != k1) {
                MatchList[1].StarList1Index = k2;
                for(k3 = 0; k3 < 4 && !FoundFit; k3++) {
                  if(k3 != k2 && k3 != k1) {
                    MatchList[2].StarList1Index = k3;
                    for(k4 = 0; k4 < 4 && !FoundFit; k4++) {
                      if(k4 != k3 && k4 != k2 && k4 != k1) {
                        MatchList[3].StarList1Index = k4;
                        NumMatch = 4;
/**/                    FitError =
                          LisaAff_Fit(
                            CatList,
                            PlateList,
                            MatchList,
                            NumMatch,
                            &Aff,
                            DistError
                          );
                        if(FitError < 4.0) { //MaxError) {
                          FoundFit = 1;
/**/                      Geo = LisaAff_AffToGeo(Aff);
                          if(fabs(fabs(Geo.Beta) - 1) > 0.5) FoundFit = 0;
                          if(fabs(Geo.Gamma) > 0.05) FoundFit = 0;
                          if(FoundFit) {
#if 0
                            //mUpdateImageStatusDel(String.Format("LISA Found 4-star solution:\r"));
                            //mUpdateImageStatusDel(String.Format("  FitError: {0:0.000}\r", FitError));
                            //mUpdateImageStatusDel(String.Format("  Geometry: {0:0.000} {1:0.000} {2:0.000}\r", Geo.Alfa, Geo.Beta, Geo.Gamma));
                            Serial_WriteString("LISA Found 4-star solution.\r");
                            Serial_WriteString("  FitError      : ");
                            sprintf((char *)SBufAUB, "%0.3f\r", FitError);
                            Serial_WriteString(SBufAUB);
                            Serial_WriteString("  Geometry      : ");
                            sprintf((char *)SBufAUB, "%0.3f %0.3f %0.3f\r", Geo.Alfa, Geo.Beta, Geo.Gamma);
                            Serial_WriteString(SBufAUB);
#endif
/**/                        LisaSolve_MakePlate(
                              RARad1,
                              DecRad1,
                              RadiusRad,
                              MaxMag,
                              &Catalog[0],
                              &NumCatalog
                            );
                            if(LisaAbortSolveB) return 0;

                            for(iPlate = 0; iPlate < NumPlate; iPlate++) {
                              UsedPlateStar[iPlate] = 0;
                            }
                            for(iPlate = 0; iPlate < 4; iPlate++) {
                              UsedPlateStar[PlateIndex[iPlate]] = 1;
                            }
                            for(iPlate = 0; iPlate < NumPlate; iPlate++) {
                              if(!UsedPlateStar[iPlate]) {
                                U = Plate[iPlate].Xcen;
                                V = Plate[iPlate].Ycen;
                                for(iCat = 0; iCat< NumCatalog; iCat++) {
                                  X = Catalog[iCat].Xcen;
                                  Y = Catalog[iCat].Ycen;
/**/                              LisaAff_Transform(Aff, X, Y, &U1, &V1);
                                  Dist = sqrt(pow(U-U1, 2) + pow(V-V1, 2));
                                  if(Dist < 4.0) { //MaxError) {
                                    memcpy(&CatList[NumMatch], &Catalog[iCat], sizeof(LisaStarTy));
                                    //memcpy(&PlateList[NumMatch], &Plate[iPlate], sizeof(LisaStarTy));
                                    PlateList[NumMatch].Xcen = Plate[iPlate].Xcen;
                                    PlateList[NumMatch].Ycen = Plate[iPlate].Ycen;
                                    PlateList[NumMatch].Intensity = 0;
                                    PlateList[NumMatch].RARad = 0;
                                    PlateList[NumMatch].DecRad = 0;
                                    MatchList->StarList1Index = NumMatch;
                                    MatchList->StarList2Index = NumMatch;
                                    MatchList->Used = 1;
                                    NumMatch++;
                                  }
                                }
                              }
                            }
                            if(NumMatch > 4) {
                              //mUpdateImageStatusDel(String.Format("  NumMatch: {0:g}  FitError: {1:0.000}\r", NumMatch, FitError));
/**/                          FitError =
                                LisaAff_Fit (
                                  CatList,
                                  PlateList,
                                  MatchList,
                                  NumMatch,
                                  &Aff,
                                  DistError
                                );
                              //Icen = (Px - 1) / 2 + 64;
                              //Jcen = (Py - 1) / 2 + 4;
                              //if(!AlignIsD3CamB) {
                                LisaAff_TransformInv(
                                  Aff,
                                  (double)(LisaCtrXRefSL + 64),
                                  (double)(LisaCtrYRefSL + 4),
                                  &U,
                                  &V
                                );
                              /*}
                              else {
                                LisaAff_TransformInv(
                                  Aff,
                                  (double)LisaCtrXRefSL,
                                  (double)LisaCtrYRefSL,
                                  &U,
                                  &V
                                );
                              }*/
                              LisaSphr_DeProject(
                                U,
                                V,
                                &LisaPlateRaRad,
                                &LisaPlateDecRad
                              );
                              //LisaSolve_MakePlate(
                              //  LisaPlateRaRad,
                              //  LisaPlateDecRad,
                              //  RadiusRad,
                              //  MaxMag,
                              //  &SavedStars,
                              //  &NumSavedStars
                              //);
                              //mUpdateImageStatusDel(String.Format("Plate Center: RA  {0:g}\r", GetRAStringFromRad(LisaPlateRaRad)));
                              //mUpdateImageStatusDel(String.Format("              Dec {0:g}\r", GetDecStringFromRad(LisaPlateDecRad)));
                              //Serial_WriteString(" Reference Ctr : RA   ");
                              //Serial_WriteHms(LisaPlateRaRad, '\r');
                              //Serial_WriteString("                 Dec ");
                              //Serial_WriteDms(LisaPlateDecRad, TRUE, '\r');
                              LisaPlateScale = 1280 / Geo.Alfa;
                              FoundFit = 1;
                            }
                            else {
                              FoundFit = 0;
                              LisaPlateRaRad = 0.0;
                              LisaPlateDecRad = 0.0;
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
        }
        k = k + 1;
      }
    }
  }
  return FoundFit;
}
