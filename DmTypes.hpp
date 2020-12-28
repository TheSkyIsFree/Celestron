#ifndef DMTYPES_H
#define DMTYPES_H
//
//  DmTypes.h
//  Created by Danyal Medley on 10/7/11.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// #defines
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define BUFFSIZE                    32                          // the TCP port we're talkin' thru

#define NULLDATA                    0x0                         // NULL data
#define NULLPTR                     ((void *)0x0)               // NULL pointer

#include <stdint.h>


typedef unsigned int    BOOL_D;           /* true or false                      */


/*----------------------------------------------------------------------------**
 **  "C" Language Conventional Constants                                       **
 **----------------------------------------------------------------------------*/
#define TRUE_D                    1       /* Numeric Value of Logical True      */
#define FALSE_D                   0       /* Numeric Value of Logical False     */


#define PI_D                            3.14159265359                 //
#define PI_MUL_2                        6.28318530716                 //
#define PI_DIV_2                        1.57079632679                 //
#define PI_DIV_3                        1.0471975512                  //
#define PI_DIV_4                        0.78539816339                 //
#define PI_DIV_180                      1.74532925199E-2              //
#define PI_DIV_12                       2.61799387798E-1              //
#define C12_DIV_PI                      3.81971863422                 //
#define C180_DIV_PI                     57.2957795133                 //

#define RAD_TO_DEG_D                    (180 / PI_D)                  //
#define RAD_TO_ARCMIN_D                 (180 / PI_D * 60)             //
#define RAD_TO_ARCSEC_D                 (180 / PI_D * 60 * 60)        //

#define DEG_TO_RAD_D                    (PI_D / 180)                  //
#define RAD_TO_BAM                      2670176.85772                 //
#define BAM_TO_RAD                      3.74507028829E-7              //
#define RAD_TO_BAM16                    10430.3783505                 //
#define BAM16_TO_RAD                    9.58737992426E-5              //

#define SECS_IN_SOLDAY                  86400.0                       //
#define SECS_IN_SIDDAY                  86164.0905                    //
#define SOLDAY_IN_SIDDAY                0.99726956637                 //
#define SIDDAY_IN_SOLDAY                1.00273790935                 //

#define SECS_15                         1.08785234286E-3              // 15 seconds in lst radians
#define SECS_1                          7.25234895243E-5              // 1 second in lst radians
#define SOLAR_SECS_1                    (SECS_1 * 0.99722222222)      // 1 second in solar radians
#define LUNAR_SECS_1                    (SECS_1 * 0.92519867042)      // 1 second in lunar radians

#define RAD_2P5DEG                      4.36332312998E-2              // 2.5 degrees in radians
#define RAD_2P0DEG                      3.49065850398E-2              // 2.0 degrees in radians

#define SIDEREAL_RATE                   15.0412 //(arcsec/sec) based on (apparent sidereal day) =23h56m4s
#define LUNAR_RATE                      14.4966 //(arcsec/sec)  based on (apparent lunar day) = 24h50m.
#define SOLAR_RATE                      15.0000 //based on (apparent solar day)=24h

#define MOTOR_NUMMODELS                  26

#define RESP_NO                         FALSE
#define RESP_YES                        TRUE





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// enumerated types
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef enum                                                    // WiFly Setting for Authentication
{                                                               //
  kOpen,                                                        // "OPEN" no authentication
  kWep,                                                         // "WEP" authentication (uses key)
  kWpa                                                          // "WPA" authentication (used passphrase)
} tWiFlyAuth;



typedef enum                                                    // WiFly Setting for DHCP
{                                                               //
  kDhcpDisabled,                                                //   "OFF" setting
  kDhcpEnabled                                                  //   "ON" setting
} tWiFlyDhcp;



typedef enum                                                    // WiFly Wireless Connection Mode
{                                                               //
  kDirectConnect,                                               //   Direct connection to SkyQ Link
  kAccessPoint                                                  //   Infrastructure connection through WiFi router
} tWiFlyConnectMode;



typedef enum                                                    // WiFly Serial Interface Mode
{                                                               //
  kCommand,                                                     //   Command / Configuration mode (configure WiFly module)
  kAux                                                          //   AUX interface mode -- serial pass-thru to AUX bus
} tWiFlyInterfaceMode;



typedef enum
{
  kPos,
  kNeg
} tDirection;



typedef enum
{
  kMin,
  kMax
} tLimit;



typedef enum                                                    // Get settings state variable enums
{                                                               //
  kGetSettingsInit,                                             //   initialized state for get settings state
  kGetSettingsSendCmd,                                          //   send the '$$$' to enter command mode
  kGetSettingsLoadConfig,                                       //   load the config to access saved settings
  kGetSettingsSendGete,                                         //   "get everything" to parse the settings
  kGetSettingsSendExit,                                         //   exit the command mode
  kGetSettingsDone,                                             //   done getting the settings
  kGetSettingsError                                             //   stuck in an error state
} tWiFlyGetSettingsState;



typedef enum                                                    // Set setting state typedef enums
{                                                               //
  kSetSettingsInit,                                             //   initialized state
  kSetSettingsSendCmd,                                          //   send the '$$$' to enter command mode
  kSetSettingsSsid,                                             //   set ssid
  kSetSettingsAuth,                                             //   set auth mode
  kSetSettingsKeyPhrase,                                        //   set key/phrase
  kSetSettingsDhcp,                                             //   set DHCP mode
  kSetSettingsIp,                                               //   set IP
  kSetSettingsNm,                                               //   set Network Mask
  kSetSettingsGw,                                               //   set Gateway
  kSetSettingsSave,                                             //   save configuration
  kSetSettingsSendExit,                                         //   exit command mode
  kSetSettingsDone,                                             //   done setting / save config
  kSetSettingsError,                                            //   stuck in an error
} tWiFlySetSettingsState;



typedef enum                                                    // Wifly command "type" for get/set settings operations
{                                                               //
  kWiFlyCmd_Enter = 0x100,                                      //   enter command mode
  kWiFlyCmd_LoadConfig = 0x101,                                 //   load WiFly's saved config to be read by get settings
  kWiFlyCmd_GetAllSettings = 0x102,                             //   get all settings
  kWiFlyCmd_SetSetting = 0x103,                                 //   save setting
  kWiFlyCmd_SaveConfig = 0x104,                                 //   save configuration to WiFly
  kWiFlyCmd_Exit = 0x1ff                                        //   exit command mode
} tWiFlyCmd;



typedef enum                                                    // get telescope settings state enum
{                                                               //
  kAzmGetVer,                                                   //   get AZM mc version
  kAzmGetModel,                                                 //   get telescepe model
  kAzmGetBacklash,                                              //   get AZM backlash settings
  kAzmGetCordwrap,                                              //   get telescope cordwrap setting
  kAltGetVer,                                                   //   get ALT mc version
  kAltGetBacklash,                                              //   get ALT backlash setting
  kTelescopeGetSettingsDone,                                    //   done getting telescope settings
  kTelescepeGetSettingsError,                                   //   error getting settings
} tMcStatusState;                                               //



typedef enum                                                    // Goto operation state enums
{                                                               //
  kGotoIdle,                                                    //   goto op in idle state
  kGotoFast,                                                    //   send fast goto command
  kGotoWaitFast,                                                //   wait for fast goto to complete (checking goto over)
  kGotoSlow,
  kGotoWaitSlow,                                                //   wait for slow goto to complete (checking goto over)
  kTracking,
  kReset                                                        //   reset the state machine
} tMcGotoState;                                                 //   Goto oepration state variable



typedef enum                                                    // MC address enums
{                                                               //
  kAddrAzm = 0x10,                                              //   AZM address 16 / 0x11
  kAddrAlt = 0x11,                                              //   ALT address 17 / 0x12
	kAddrFocuser = 0x12,
	kAddrPhone = 0x20,
  kAddrWiFi = 0xb5,                                             //   WiFi "device"
  kAddrPwr = 0xb6,                                              //   Power charge "device"
  kAddrDcp = 0xb7,                                              //   USB charge port
  kAddrAis = 0xb4,                                              //   AIS Camera
  kAddrLmp = 0xbf,                                              //   "lamps"
} tAuxAddr;                                                     //



typedef enum
{                                                               // Motor Controller command enums
  kMtr_GetEncPosition = 0x01,                                   //   get position command
  kMtr_GotoFast = 0x02,                                         //   goto fast command
  kMtr_PositionSet = 0x04,                                      //   set motor position
  kMtr_GetModel = 0x05,                                         //   get telescope model
  kMtr_TrackPositive = 0x06,                                    //   set positive tracking rate
  kMtr_TrackNegative = 0x07,                                    //   set negative tracking rate
  kMtr_MoveSwitch = 0x0b,                                       //   send axis to switch position
  kMtr_BacklashPosSet = 0x10,                                   //   get positive backlash setting
  kMtr_BacklashNegSet = 0x11,                                   //   get negative backlash setting
  kMtr_MoveSwitchOver = 0x12,                                   //   check if move switch is over
  kMtr_GotoOver = 0x13,                                         //   check goto over
  kMtr_GotoSlow = 0x17,                                         //   slow goto command
  kMtr_RaLimitMinSet = 0x1a,                                    //   set ra slew limit min
  kMtr_RaLimitMaxSet = 0x1b,                                    //   set ra slew limit max
  kMtr_RaLimitMinGet = 0x1c,                                    //   get ra slew limit min
  kMtr_RaLimitMaxGet = 0x1d,                                    //   get ra slew limit max
  kMtr_RaLimitEnaGet = 0x1e,                                    //   get the enabled setting for ra limit
  kMtr_RaLimitEnaSet = 0x1f,                                    //   enable ra slew limits
  kMtr_CustomRate9Set = 0x20,                                   //   set the custom rate9 rate
  kMtr_CustomRate9Get = 0x21,                                   //   get the custom rate9 rate
  kMtr_CustomRate9EnaSet = 0x22,                                //   enable the custom rate9 feature
  kMtr_CustomRate9EnaGet = 0x23,                                //   get the enabled setting for custom rate9
  kMtr_PmSlew = 0x24,                                           //   positive slew command
  kMtr_NmSlew = 0x25,                                           //   negative direction slew
  kMtr_CordwrapOn = 0x38,                                       //   turn on cordwrap
  kMtr_CordwrapOff = 0x39,                                      //   turn off cordwrap
  kMtr_CordwrapPosSet = 0x3a,                                   //   set the cordwrap position
  kMtr_CordwrapEnaGet = 0x3b,                                   //   get the enabled setting for cordwrap
  kMtr_CordwrapPosGet = 0x3c,                                   //   get the cordwrap position
  kMtr_BacklashPosGet = 0x40,                                   //   get positive backlash setting
  kMtr_BacklashNegGet = 0x41,                                   //   get negative backlash setting
  kMtr_GuideRateSet = 0x46,                                     //   set the autoguider rate
  kMtr_GuideRateGet = 0x47,                                     //   get the autoguider rate
  kMtr_UnrecognizedCommand = 0xf0,								//	 when a device does not recognize a command sent to it, this will be returned
  kMtr_ApproachDirGet = 0xfc,                                   //   gets the approach direction setting
  kMtr_ApproachDirSet = 0xfd,                                   //   sets the approach direction

  kAis_StartCapture = 0x90,                                     // start capture
  kAis_GetStatus = 0x91,                                        // get capture status
  kAis_GetPlateList = 0x92,                                     // get the plate list (extended packet response)
  kAis_SetGain = 0x93,                                          // set the CCD gain -- old camera cmd
  kAis_ProcessImage = 0x94,                                     // command to start processing image (opt idx param new cam)
  kAis_ProcessImagePrgs = 0x95,                                 // process captured image (opt idx param new cam)
  kAis_ProcessImageResult = 0x96,                               // gets the result of image processing -- diagnostic
  kAis_GetBrightest = 0x97,                                     // gets the location of brightest star
  kAis_BumpUsb = 0x98,                                          // cycle USB interface to allow PC to resync
  kAis_BootSwap = 0x98,                                         // force swap update of starsense program
  kAis_Reset = 0x9f,                                            // reset the state of the AIS camera
	kAis_SetCenterRef = 0x3E,                                     // stores the center reference value in camera (two signed int32 values, x,y, and then 0)
	kAis_GetCenterRef = 0x3F,                                     // reads out the center reference value in camera (two signed int32 values, x,y)


                                                                // "Lamp" Controller command enums
  kLmp_LampBrightness = 0x10,                                   //

                                                                // WiFi device command enums
  kWiFi_Status = 0x10,                                          //
  kWiFi_Reset = 0x12,                                           //

                                                                // DCP device command enums
  kDcp_Status = 0x10,                                           //

                                                                // Power Controller command enums
  kPwr_Status = 0x10,                                           //
  kPwr_CurrentLimit = 0x18,                                     //

                                                                // All device command enums
  kCmd_GetVer = 0xfe,                                           //   get version command
} tAuxCommand;                                                  //



typedef enum                                                    // Motor slew speed enums
{                                                               //
  kRate0 = 0,                                                   //
  kRate1 = 1,                                                   //
  kRate2 = 2,                                                   //
  kRate3 = 3,                                                   //
  kRate4 = 4,                                                   //
  kRate5 = 5,                                                   //
  kRate6 = 6,                                                   //
  kRate7 = 7,                                                   //
  kRate8 = 8,                                                   //
  kRate9 = 9                                                    //
} tSlewRate;                                                    //



typedef enum                                                    // Motor slew diretion enums
{                                                               //
  kStop,                                                        //
  kMoveNeg,                                                     //
  kMovePos                                                      //
} tMcSlewDir;                                                   //



typedef enum                                                    // Telescope models
{                                                               //
  kUndef,                                                       //   Undefined
  kNexStarGps,                                                  //   NexStar GPS
  kNexStarGpsSa,                                                //   NexStar GPS
  kNexStari,                                                    //   NexStar i
  kNexStarSe,                                                   //   NexStar SE
  kCge,                                                         //   CGE Series
  kAdvancedGt,                                                  //   Advanced GT
  kSlt,                                                         //   SLT Series
  kLegend,                                                      //   Legend
  kCpc,                                                         //   CPC Series
  kNexStarGt,                                                   //   NexStar GT
  kNexStarSe45,                                                 //   NexStar SE 4/5
  kNexStarSe68,                                                 //   NexStar SE 6/8
  kCgePro,                                                      //   CGE Pro
  kCgem,                                                        //   CGEM
  kLcm,                                                         //   LCM Series
  kSkyProdigy,                                                  //   SkyProdigy
  kCpcDelux,                                                    //   CPC Deluxe
  kNexStarGt16,                                                 //   NexStar GT.
  kStarSeekerGt,                                                //   NexStar GT for Orion
  kAdvancedVx,                                                  //   Advanced VX
  kCosmosGt,                                                    //   Cosmos GT 90
  kEvolution,                                                   //   Celestron Evolution (CEVO)
} tTelescopeModel;



typedef enum                                                    // Alignment mode enum
{                                                               //
  kAlignModeAzmAlt,                                             //   Azm Alt mode
  kAlignModeEqN,                                                //   Equatorial North
  kAlignModeEqS,                                                //   Equatorial South
} tAlignMode;                                                   //   alignment mode type



typedef enum                                                    // Tracking rate enums -- actual MC settings
{                                                               //
  kTrackingSidereal = 0,                                        //
  kTrackingSolar = 1,                                           //
  kTrackingLunar = 2                                            //
} tTrackingRate;                                                //



typedef enum
{
  kAlignSkyAlign = 0,                                           // explicit, matches settings.h values
  kAlignManual = 1,
  kAlignEquatorial = 2,
  kAlignStarSense = 3
} tAlignType;



typedef enum
{
  kChargePortAuto,
  kChargePortOn
} tDcpChargeMode;



typedef enum
{
  kLvlLow,
  kLvlMed,
  kLvlHigh
} tBattLvl;



typedef enum
{
  kLampTray,
  kLampWiFi,
  kLampPower,
  kLampMount,
  kNumEnd
} tLamp;



typedef enum
{
  kDischarging,
  kCharging,
  kCharged,
  kFault
} tBattStatus;


typedef enum   // Make these explicit since they are matched in the Java code
{
  kAisNoErr = 0,
  kAisCommErr = 1,
  kAisTooFewStars = 2,
  kAisNoSolve = 3,
  kNoMoreGotos = 4
} tAlignError;

typedef enum  // Make these explicit since they are matched in the Java code
{
  kAlignIdle = 0,                                // waiting to start alignment
  kStartAlignment = 1,                           // start the alignment -- initiate parameters for align process
  kReStartAlignment = 2,                         // restarts an alignment -- after a calibrate, does not init position
  kAlignInProgress = 3,                          // alignment is in progress -- the ais / goto state handlers are runing
  kAlignDone = 4,                                // alignment is done -- the ais / goto state handlers exited, align complete
  kAlignFailed = 5,                              // alignment failed -- arises when pointing exceeds predefined threshold
  kAlignAddManualReference = 6,                  // add an alignment reference
  kAlignAddManualReferenceDone = 7,              // done with add reference
  kAlignAddManualReferenceFailed = 8,            // add manual reference failed RMS pointing threshold
  kAlignAddReferenceInit = 9,                    // init for add an alignment reference
  kAlignAddReference = 10,                       // add an alignment reference
  kAlignAddReferenceDone = 11,                   // done with add reference
  kAlignAddReferenceFailed = 12,                 // add alignment failed RMS pointing threshold
  kAlignCalibrateInit = 13,                      // init for calibrate the starsense camera
  kAlignCalibrate = 14,                          // calibrate the starsense camera
  kAlignCalibrateDone = 15,                      // done with calibrate
  kAlignError = 16                               // alignment failed -- there was a problem during alignment
} tStarSenseAlignState;

typedef enum  // Make these explicit since they are matched in the Java code
{
  kAisIdle = 0,                                  // waiting for goto operation to complete, or AIS capture not yet initiated
  kWaitStartCapture = 1,                         // start capture initiated -- wait for any goto op to complete
  kStartCapture = 2,                             // start the capture process
  kCheckCaptureDone = 3,                         // check if the capture is completed, get status
  kProcessImg = 4,                               // start image processing
  kCheckProcessingDone = 5,                      // check if image processing is done
  kGetPlate = 6,                                 // get the plate when image processing is done
  kSolvePlate= 7,                                // solve the plate when get plate is done
  kSolved = 8,                                   // solved message with delay
  kAisError= 9                                   // indicates an error in process -- error is captured by AlignImgSystemError
} tAisCameraState;


typedef enum   // Make these explicit since they are matched in the Java code
{
  kAlignGotoIdle = 0,                            // waiting for ais image capture to complete, or alignment not initiated
  kGotoAlignPosition = 1,                        // determine next goto position, and slew
  kCheckGotoAlignPositionDone = 2,               // check if goto operation is over
  kAlignGotoError = 3                            // indicates an error in goto state
} tAlignGotoState;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// typedef structs
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef struct CelestialRefTg {
  double
    ra;                                                               //  Right Ascention of point
  double
    dec;                                                              //  Declination of point
  double
    alt;                                                              //  Altitude of point
  double
    azm;                                                              //  Azimuth of point
  double
    lst;                                                              // "days" since power on at time of star astroment
  double
    jd;                                                               // Julian
} CelestialRefTy;




typedef struct sMcVersionTg                                     // Version information
{                                                               //
  int VerMaj;                                                   //
  int VerMin;                                                   //
  int VerBld;                                                   //
} sMcVersion;



typedef struct
{
  tTrackingRate trackingRate;
  bool trackingEnabled;
  bool trackingUpdateEnabled;
  bool trackRaOnly;
  double nextUpdateTick;
} sTrackingSettings;



typedef struct
{
  tStarSenseAlignState
    State;
  tAlignGotoState
    AlignGotoState;
  tAisCameraState
    AisCameraState;
  tAisCameraState
    AisContinueState;
  double
    NextUpdateTick;
  int
    Xcen;
  int
    Ycen;
  double
    skyObjRa;
  double
    skyObjDec;
  bool
    isManualAlign;
  tAlignError
    alignError;
} sStarSenseAlignState;



typedef struct
{
  float battVolt;
  tBattLvl battLvl;
  tBattStatus battStatus;
  double nextUpdateTick;
} sCevoStatus;



typedef struct sMountSettingsTg                                 // Motor controller settings
{                                                               //
  //On connect :: group 0 -- read on connect to mount
  sMcVersion VersionInfo;
  tTelescopeModel telescopeModel;                               //    telescope model
  bool haveGroup0;

  //Common settings :: group 1                                  //
  tDirection approachDirectionAzm;                              //
  tDirection approachDirectionAlt;                              //
  int backlashAzmPos;                                           //   backlash settings for each axis / dir
  //int backlashAzmNeg;                                           //
  int backlashAltPos;                                           //
  //int backlashAltNeg;                                           //
  bool useGroup1;                                               //
  bool haveGroup1;                                              //

  //Has Custom Rate 9 :: group 2
  bool customRate9AzmEnabled;                                   //
  int customRate9Azm;                                           //   expressed in milidegrees
  int customRate9MaxValue;                                      //   expressed in milidegrees
  //bool customRate9AltEnabled;                                 //
  //int customRate9Alt;                                         //   must enforce alt setting = azm setting
  bool useGroup2;                                               //
  bool haveGroup2;

  //Can Do Wedge | Is Gem :: group 3                            //
  int guideRateRa;                                              //   guider port rates
  int guideRateDec;                                             //
  bool useGroup3;                                               //
  bool haveGroup3;                                              //

  //Has Ra Limits :: group 4
  bool raLimitsEnabled;                                         //
  int raLimitMax;                                               //
  int raLimitMin;                                               //
  bool useGroup4;                                               //
  bool haveGroup4;                                              //

  // group 5 and 6
  int mountLighting;                                            //

  //Is Cevo :: group 5
  int trayLighting;                                             //
  tDcpChargeMode chargeMode;                                    //
  float maxCurrent;                                             //
  bool useGroup5;                                               //
  bool haveGroup5;                                              //
                                                                //
  //Has Mount Lighting :: group 6
  bool useGroup6;                                               //
  bool haveGroup6;                                              //

} sMountSettings;                                               //



typedef struct sUserSettingsTg {
//  bool                                                                //
//    trackingUpdateEnabled;                                                  // default tracking on/off after alignment (on, off)
//  tTrackingRate                                                       //
//    trackingRate;                                                     // default tracking rate after alignment (sidereal, solar, lunar)
//  tAlignType                                                          //
//    alignType;                                                        // alignment mode
//  int                                                                 //
//    slewLimitMin;                                                     // lower alt limit for alt/azm scope 0-90 deg
//  int                                                                 //
//    slewLimitMax;                                      \
//    cordwrapEnabledB;
  tTelescopeModel                                                     // READ ONLY -- telescope model
    telescopeModel;                                                    //

  // user settings for common settings
  bool                                                                //
    reverseLeftRight;                                                 // (motor setting*) sets approach, button dir for lf/rt
  bool                                                                //
    reverseUpDown;                                                    // (motor setting*) sets approach, button dir for up/dn
  int                                                                 //
    backlashAzm;                                                      // (motor setting) azm(ra) backlash 0-99
  int                                                                 //
    backlashAlt;                                                      // (motor setting) alt(dec) backlash 0-99

  // user settings for group 2
  float                                                               //
    fastSlewRate;                                                     // (motor setting) the default max rate (rate9) 0.25 - max*
  bool                                                                //
    fastSlewRateEnabled;                                              // (motor setting) is custom max rate enabled
  float                                                               //
    fastSlewRateMax;                                                  //

  // user settings for group 3
  int                                                                 //
    guideRateRa;                                                      //
  int                                                                 //
    guideRateDec;                                                     //

  // user settings for group 4                                        //
  int                                                                 //
    raLimitMin;                                                       // west limit (North EQ)
  int                                                                 //
    raLimitMax;                                                       // east limit (North EQ)
  bool                                                                //
    raLimitsEnabled;                                                   // use the software ra limits

  // user settings for group 5
  int                                                                 //
    cevoTrayLighting;                                                 // (motor setting) CEVO lighting setting 0-10
  //int                                                                 //
  //  cevoMountlighting;                                                // (motor setting) CEVO mount lighting setting 0-10
  tDcpChargeMode                                                      //
    cevoChargePortMode;                                               // (motor setting) USB charge port auto or ON (verride)
  float                                                               //
    cevoMaxCurrent;                                                   //

  // user settings for group 5 & 6
  int                                                                 //
    mountlighting;                                                    // (motor setting) mount lighting setting 0-10


} sUserSettings;



typedef struct sTelescopeConfigTg {
  bool
    telIsGemB,                                                  //
    telIsCevoB,                                                 //
    telHasLightControlsB,                                       //
    telHasSwitchesB,                                            //
    telHasPecB,                                                 //
    telHasRaLimitsB,                                            //
    telCanDoStarSenseB,                                         //
    telCanDoWedgeB,                                             //
    telHasCustomRate9B,                                         //
    telIsEqAlignB;                                              //
} sTelescopeConfig;



typedef struct DmsTg {
    uint32_t
      DegreesUM;
    uint32_t
      MinutesUT;
    uint32_t
      SecondsUT;
    BOOL_D
      MinusB;
} DmsTy;



typedef struct HmsTg {
    uint32_t
      HoursUT;
    uint32_t
      MinutesUT;
    uint32_t
      SecondsUT;
    uint32_t
      TenthsUT;
} HmsTy;



#define AIS_STA_RDY             0x01                            // status flags from AIS camera
#define AIS_STA_IMG             0x02                            //
#define AIS_STA_PROC            0x04                            //
#define AIS_STA_PLATE           0x08                            //
#define AIS_STA_ERROR           0x10                            //
#define AIS_STA_TIMEDOUT        0xff                            // (not from camera, but indicates a time-out on AUX bus)

#define AIS_STA_RDYIMG          (AIS_STA_RDY | AIS_STA_IMG)     //
#define AIS_STA_RDYPLATE        (AIS_STA_RDY | AIS_STA_PLATE)   //



#endif
