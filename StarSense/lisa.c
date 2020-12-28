/*******************************************************************************************************************************
** FILE: lisa.c -- Lost In Space Algorithm implementation for the NexAlign telescope alignment
** DATE: 2009.03.23 DJM
**
** NOTE: Copywrite (c) 2009 by Celestron.


===========================================
              LISA OVERVIEW
            by: Daniel Goodwin
===========================================


============
   BRIEF:
============
Takes a brightness-organized list of stars pulled from an image (plate-list), selects four to create a quad, compares that quad to a database of quads to find a matching pair, and calculates Ra/Dec coordinate of the specified image position.



============
   TERMS:
============
Quad
- 4 stars (technically the index of the star in the associated star database)
- 5 normalized distances
- 1 normalized distance sum (This term will be referred to as NDS)
- Grouping of four stars which are within a maximum distance of one another
- They are arranged brightest to dimmest
- The six distances between each of the four stars are calculated, then arranged in descending order, then normalized by the longest distance
- The sum of the lower five normalized distances is then calculated and stored in the first distance term (the NDS)
	- Since the first term is always 1.0, therefore it is used for this purpose instead

Quad Database
- A progenerated list of quads which was created from an associated star database and arranged in ascending order of the quad's NDS value

Star Database
- List of stars, organized by star magnitude, from brightest to dimmest



============
   STEPS:
============
1. Pass parameters
	- Pass a plate-list (list of stars from an image) where the stars are arranged brightest to dimmest
	- Image FOV and image size (or if you're Medley, hardcode it in)
2. Create a plate-quad
	- Select the first star in the plate-list and find the next three brightest stars in the plate-list that are within a specified distance and create a plate-quad
	- This distance is based on the characteristics of the quads in the  quad database (StarSense's quad database has quads where all stars are within 3.0 degrees of each other)
3. Find a similar quad to compare
	- Used the plate-quad's NDS to search the quad database for a similar quad
	- Compare the five individual distances between the database-quad and the plate-quad
		- Only progress forward if the error is below a threshold
4. Perform an affine fit
	- Spherically project the database-quad's stars onto a imaginary plate
	- Take the plate-quad and the database-quad, and perform an affine fit between the two
	- Theoretically the stars should be in the same order BUT just in case, all possible arrangements of the four plate-quad stars are affine-fit tested
	- If the error in the affine transfer is low and the amount of skew (gamma term) and stretch (beta term) is below a threshold then continue
5. Test to see if other stars in plate list appear in star database
	- Create a spherical projection and create a sub-catalogue of stars which should appear on the plate
	- Compare the sub-catalogue of stars to the plate list and see if any of them match
	- If enough stars appear then continue
	- If the above does not check out, then select then next star in the plate list and jump back to step 2
6. Calculate Ra/Dec
	- The Ra/Dec value of the image 'centroid' (which for StarSense is calibrated by the user) is calculated



============
  Summary:
============
- In order for a quad pulled from the quad database to be matched to a quad created from the plate-list:
	- The Distances must all matches within a specified threshold
	- An affine fit between the stars in the plate-quad and the database-quad must have minimal skew and stretch
	- A specified number of other stars in the plate-list must be found in the star database



*******************************************************************************************************************************/
#define __LISA_C__
#include "lisa.h"                                              // file system support for image processing  access



//==============================================================================================================================
// global variable declarations and subsystem function prototypes
//==============================================================================================================================
volatile BOOL_D
  LisaHaveSolutionB,
  LisaAbortSolveB;



//==============================================================================================================================
// FUNC: Lisa_SphericalProject()
//==============================================================================================================================
// DESC: use LISA to project Ra, Dec to x,y coords
//==============================================================================================================================
void Lisa_SphericalProject(                                      // RETR: none
  double
    RaRadF,
  double
    DecRadF,
  double
    *Xcen,
  double
    *Ycen
)
{
  LisaSolve_SphericalProject(RaRadF, DecRadF, Xcen, Ycen);
}



//==============================================================================================================================
// FUNC: Lisa_AisSolvePlate()
//==============================================================================================================================
// DESC: sends aux command to AIS and reports the AIS status
//==============================================================================================================================
extern BOOL_D Lisa_AisSolvePlate(                               // RETR: none
	float                                                         //
		*CoordList,                                                 //
	uint32_t                                                          //
		CoordListCountUT                                            //
)
{
	int count = 0;
	for (int i = 0; i < CoordListCountUT; i++)
	{
		Plate[count].Xcen = CoordList[i];
		Plate[count].Ycen = CoordList[i + 1];
		count++;
	}
	PlateCountUT = CoordListCountUT;

	return LisaSolve_Solve();
}



//==============================================================================================================================
// FUNC: Lisa_AisSetCenterCoords()
//==============================================================================================================================
// DESC: sents the default center coordinates for LISA
//==============================================================================================================================
extern void Lisa_AisSetCenterCoords(                            // RETR: none
	int32_t                                                          //
		CtrXRefSL,                                                  //
	int32_t                                                          //
		CtrYRefSL                                                   //
)
{
	LisaCtrXRefSL = CtrXRefSL;
	LisaCtrYRefSL = CtrYRefSL;
}



//==============================================================================================================================
// FUNC: Lisa_GetCenterInRads()
//==============================================================================================================================
// DESC: gets the center reference celestial coords in Rads from last plate solve
//==============================================================================================================================
extern void Lisa_GetCenterInRads(                               // RETR: none
	double                                                        //
		*RaRad,                                                     // PASS: -> recvs Ra coord
	double                                                        //
		*DecRad                                                     // PASS: -> recvs Dec coord
)
{
	*RaRad = LisaPlateRaRad;
	*DecRad = LisaPlateDecRad;
}



//==============================================================================================================================
// FUNC: Lisa_Init
//==============================================================================================================================
// DESC:
//==============================================================================================================================
extern BOOL_D Lisa_Init(                                          // RETR: none
	FILE
		*QuadFileP,
	const uint32_t
		QuadFileOffsetUL,
	const uint32_t
		QuadFileSizeUL,
	FILE
		*StarFileP,
	const uint32_t
		StarFileOffsetUL,
	const uint32_t
		StarFileSizeUL
)
{
	LisaSolve_ScanQuadFile(QuadFileP, QuadFileOffsetUL, QuadFileSizeUL);
	LisaSolve_ScanAPMFile(StarFileP, StarFileOffsetUL, StarFileSizeUL);
	return TRUE_D;
}
