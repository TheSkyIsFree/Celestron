/*******************************************************************************************************************************
** FILE: align.c -- starsense alignment routines
** DATE: 2015.10.09 DJM
** DESC: Implements the align routines for starsense
*******************************************************************************************************************************/
#define __ALIGN_C__
#include "DmTypes.hpp"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "align.hpp"                                                // NexPoint model interface
#include "target.h"


//==============================================================================================================================
// global variable declarations and subsystem function prototypes
//==============================================================================================================================
uint32_t                                                             //
  AlignLastStepUT;                                              // Step/Last Step being processed
uint32_t                                                             //
  AlignStepUT;                                                  // the ss auto current alignment step

BOOL_D                                                          //
  AlignIsEqAlignB;                                              //
BOOL_D                                                          //
  AlignIsGemB;                                                  //
BOOL_D                                                          //
  AlignUseUserAutoAlignB;                                       //




//==============================================================================================================================
// static variable declarations
//==============================================================================================================================
uint32_t                                                             //
  AlignLastFailedStepUT,                                        //
  FirstCapB;                                                    // flags the first capture
uint8_t                                                           //
  AlignStepZoneUB[10];                                          //

const int32_t 
  AlignZoneTblAzmAlt[9 * 2 * 3] = {
//Zone1
//  Azm  Alt       offset
  5,   28,    // 0    0
  5,   33,    // 0    5
  5,   38,    // 0    10
  30,  33,    // 25   5
  30,  38,    // 25   10
  40,  33,    // 35   5
  40,  38,    // 35   10
  55,  33,    // 50   5
  55,  38,    // 50   10
//Zone2
//  Azm  Alt       offset
  105, 44,    // 0    0
  130, 44,    // 25   0
  155, 44,    // 50   0
  180, 44,    // 75   0
  105, 54,    // 0    10
  130, 54,    // 25   10
  155, 54,    // 50   10
  180, 54,    // 75   10
  105, 34,    // 0    -10
//Zone3
//  Azm  Alt       offset
  205, 62,    // 0    0
  235, 62,    // 30   0
  265, 62,    // 60   0
  295, 62,    // 90   0
  325, 62,    // 120  0
  295, 52,    // 90   -10
  265, 52,    // 60   -10
  235, 52,    // 30   -10
  205, 52     // 0    -10
};



const int32_t 
  AlignZoneTblWedged[9 * 2 * 3] = {
//Zone1 Azm     Alt
  95,     40,   //   0     0
  110,     40,   //  15     0
  125,     40,   //  30     0
  110,     45,   //  15     5
  125,     45,   //  30     5
  140,     45,   //  45     5
  110,     50,   //  15    10
  125,     50,   //  30    10
  140,     50,   //  45    10
//Zone 2           //
  165,    -15,   //   0     0
  155,    -15,   // -10     0
  145,    -15,   // -20     0
  170,    -10,   //   5     5
  165,     -5,   //   0    10
  165,      0,   //   0    15
  155,     -5,   // -10    10
  145,     -5,   // -20    10
  145,      0,   // -20    15
//Zone 3           //
  240,     55,   //   0     0
  230,     55,   // -10     0
  220,     55,   // -20     0
  230,     60,   // -10     5
  240,     60,   //   0     5
  230,     50,   // -10    -5
  220,     45,   // -20   -10
  230,     45,   // -10   -10
  240,     45    //   0   -10
};



const int32_t 
  AlignZoneTblEq[9 * 2 * 4] = {
  // Zone 1
  // Azm  Alt         offset
  45,  135,   // 0       0
  45,  140,   // 0       5
  45,  145,   // 0       10
  70,  140,   // 25      5
  70,  145,   // 25      10
  80,  140,   // 35      5
  80,  145,   // 35      10
  95,  140,   // 50      5
  95,  145,   // 50      10

  // Zone 2
  // Azm  Alt         offset
  20,  150,   // 0       0
  45,  150,   // 25      0
  70,  150,   // 50      0
  95,  150,   // 75      0
  20,  160,   // 0       10
  45,  160,   // 25      10
  70,  160,   // 50      10
  95,  160,   // 75      10
  20,  140,   // 0       -10

  // Zone 3
  // Azm  Alt         offset
  130, 60,    // 0       0
  135, 60,    // 30      0
  145, 50,    // 60      0
  115, 50,    // 90      0
  105, 40,    // 120     0
  115, 40,    // 90      -10
  125, 30,    // 60      -10
  135, 30,    // 30      -10
  145, 20,    // 0       -10

  // Zone 4
  // Azm  Alt         offset
  155, 50,    // 30      0
  165, 50,    // 0       0
  145, 50,    // 60      0
  125, 50,    // 90      0
  115, 50,    // 120     0
  105, 40,    // 90      -10
  115, 40,    // 60      -10
  125, 40,    // 30      -10
  135, 40     // 0       -10


};


/*******************************************************************************************************************************
********                                       static function implementations                                          ********
*******************************************************************************************************************************/

/*******************************************************************************************************************************
********                                      global function implementations                                           ********
*******************************************************************************************************************************/

//==============================================================================================================================
// FUNC: Align_AutoPositionInit()
//==============================================================================================================================
// DESC: initializes the auto position move
//==============================================================================================================================
void Align_AutoPositionInit(                                    // RETR: none
  BOOL_D                                                        //
    IsGemB,                                                     // PASS: is telescope GEM or Alt/Azm (incl. wedged)
  BOOL_D                                                        //
    IsEqAlignB                                                  // PASS: is an eq alignment
)                                                               //
{
  AlignIsGemB = IsGemB;                                         //
  AlignIsEqAlignB = IsEqAlignB;                                 //
  memset(&AlignStepZoneUB, 0, sizeof(AlignStepZoneUB));         //
  AlignStepUT = 1;                                              // Initialize the SsAuto step, and set state for
  AlignLastFailedStepUT = 0x0;                                  // clear this before we start
  AlignLastStepUT = 0xff;                                       // clear these
  FirstCapB = TRUE_D;                                           //
}



//==============================================================================================================================
// FUNC: Align_AutoPositionClearCurrentStep()
//==============================================================================================================================
// DESC: clear the current step
//==============================================================================================================================
void Align_AutoPositionClearCurrentStep(                        // RETR: none
  void                                                          // PASS: none
)                                                               //
{
  if(AlignLastStepUT > 0 || AlignLastStepUT < 5) {              //
    AlignStepZoneUB[AlignLastStepUT - 1] |= 0x80;               //   bit7 indicates zone is complete
    if(AlignLastFailedStepUT == AlignLastStepUT) {              //   if: this step was last failed step, then clear last failed
      AlignLastFailedStepUT = 0x00;                             //     step by setting to 0 (not a valid step)
    }                                                           //
  }
}



//==============================================================================================================================
// FUNC: Align_AutoPositionFailCurrentStep()
//==============================================================================================================================
// DESC: clear the current step
//==============================================================================================================================
void Align_AutoPositionFailCurrentStep(                         // RETR: none
  void                                                          // PASS: none
)
{
  AlignLastFailedStepUT = AlignLastStepUT;                      //
}


//==============================================================================================================================
// FUNC: Align_AutoPositionCurrentStep()
//==============================================================================================================================
// DESC: Gets the current align step -- sets the last align step before starting image capture
//==============================================================================================================================
uint32_t  Align_AutoPositionCurrentStep(                             // RETR: the current align step
  BOOL_D                                                        //
    SetLastStepB                                                // PASS: flag wether to set the last step value
)
{
  if(SetLastStepB == TRUE_D) {
    AlignLastStepUT = AlignStepUT;
  }
  return AlignStepUT;
}



//==============================================================================================================================
// FUNC: Align_AutoPositionNextMove()
//==============================================================================================================================
// DESC: performs the automatic move to "next" alignment position, based on the passed StepUT parameter
//==============================================================================================================================
BOOL_D Align_AutoPositionNextMove(                              // RETR: none
  int32_t                                                           //
    *Axis1PosSM,                                                // PASS: -> axis 1 (azm / ra) goto position
  int32_t                                                           //
    *Axis2PosSM                                                 // PASS: -> axis 2 (alt / dec) goto position
)                                                               //
{
  uint32_t                                                           //
    CntUT = 0,                                                      //
    IdxUT,                                                      //
    TblBaseUT,                                                  //
    TblOffsUT;                                                  //

  if(AlignUseUserAutoAlignB) {                                  // if: user auto align, use saved step count
    //CntUT = EepromConfigSt.AisSt.AlignZoneUserZonesUB;        //
  }                                                             //
  else {                                                        // else: use the default 3 for alt/azm, 4 for GEM
    CntUT = AlignIsGemB ? 4 : 3;                                //
  }                                                             //
  for(IdxUT = 0; IdxUT < CntUT; IdxUT++) {                      // for: each alignment reference step (find next zone)
    if(                                                         //   if: (find the next align zone in step to do)
      !(AlignStepZoneUB[IdxUT] & 0x80) &&                       //     this step does not have a solution AND
      (                                                         //
        IdxUT != AlignLastFailedStepUT - 1 ||                   //     this step is not same as current failed step (don't redo
        AlignLastFailedStepUT != AlignLastStepUT                //        previous failed step unless it's only one remaining)
      ) &&                                                      //
      (                                                         //     this step is not current failed step except if it's last
        IdxUT != AlignLastStepUT - 1 ||                         //       (only redo AlignLastStep if it's the last step)
        (AlignStepZoneUB[IdxUT+1] & 0x80) ||                    //
        IdxUT == (CntUT - 1)                                    //
      )                                                         //
    ) {                                                         //
      AlignStepUT = IdxUT + 1;                                  //
      break;                                                    //
    }                                                           //
  }                                                             //

  if(AlignUseUserAutoAlignB) {                                  // if: user starsense align (not implemented)
//    uint32_t                                                         //
//      IdxUT = AlignStepUT - 1;                                  //
//    AzmOffsSM =                                                 //
//      EepromConfigSt.AisSt.AlignZoneUserASW[2 * IdxUT];         //
//    AltOffsSM =                                                 //
//      EepromConfigSt.AisSt.AlignZoneUserASW[2 * IdxUT + 1];     //
  }                                                             //
  else {                                                        // else: normal auto starsense align
    if(AlignStepZoneUB[AlignStepUT - 1] >= 9) {                 //   if: we search all 9 zones...
			//printf("StarSense Step %d, %d\n", AlignStepUT, AlignStepZoneUB[AlignStepUT - 1]);
      return FALSE_D;                                           //
    }                                                           //
    TblBaseUT = 9 * 2 * (AlignStepUT - 1);                      //     base index for step within x2 array
    TblOffsUT = AlignStepZoneUB[AlignStepUT - 1] * 2;           //     offset for zone within step
    if(AlignIsGemB) {                                           //     if: (use GEM table)
      *Axis1PosSM =                                             //
        AlignZoneTblEq[TblBaseUT + TblOffsUT];                  //
      *Axis2PosSM =                                             //
        AlignZoneTblEq[TblBaseUT + TblOffsUT + 1];              //
    }                                                           //
    else if(AlignIsEqAlignB) {                                  //     else if: (use EQ Wedge table)
      *Axis1PosSM =                                             //
        AlignZoneTblWedged[TblBaseUT + TblOffsUT];              //
      *Axis2PosSM =                                             //
        AlignZoneTblWedged[TblBaseUT + TblOffsUT + 1];          //
    }                                                           //
    else {                                                      //     else: (use alt/azm table)
      *Axis1PosSM =                                             //
        AlignZoneTblAzmAlt[TblBaseUT + TblOffsUT];              //
      *Axis2PosSM =                                             //
        AlignZoneTblAzmAlt[TblBaseUT + TblOffsUT + 1];          //
    }                                                           //
  }                                                             //
  AlignStepZoneUB[AlignStepUT - 1]++;                           // move to next zone in anticipation of do-over
  return TRUE_D;                                                //
}

/*******************************************************************************************************************************
********                                       static function implementations                                          ********
*******************************************************************************************************************************/


