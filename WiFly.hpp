//
//  WiFly.h
//  Created by Danyal Medley on 10/7/11.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef WIFLY_HPP
#define WIFLY_HPP

#include "AuxRecv.hpp"
#include "AuxSend.hpp"

class WiFly
{
public:
  AuxRecv                                                       //
    *mRecvResp;                                                 // the last successfully parsed AUX packet
  CTelescope *mTelescope;

  WiFly(CTelescope *pTelescope);
  
  
  virtual ~WiFly(void);


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  #pragma mark - AUX Operations
  uint8_t *ParseAux(
    unsigned char
      *recvBuff,
    int
      len,
    tAuxAddr
      DstAxis,
    tAuxCommand
      Cmd
  );

  void ParseAuxEx(
    uint8_t                                                             //
      *recvBuff,                                                      //
    int                                                               //
      len                                                             //
  );

  //============================================================================================================================
  // AuxPacketMaster -- does the packet exchange with AUX device
  //
  //   formwats an AUX packet, sends to AUX bus, the response is captured and parsed into (new AuxRecv) *mRecvResp
  //============================================================================================================================
  void AuxPacketMaster(                                               // RETR: NONE
    tAuxAddr                                                          //
      Axis,                                                           // PASS: enumerated address
    tAuxCommand                                                       //
      Cmd,                                                            // PASS: AUX command
    uint8_t                                                             //
      *data,                                                          // PASS: -> data buffer (or NULL)
    int                                                               //
      dlen                                                            // PASS: -> data length (or 0)
  );
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  #pragma mark - AUX Device common commands
  bool AuxDevGetVer(
    tAuxAddr
      Addr,
    sMcVersion
      *VerInfo
  );


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  #pragma mark - AIS Camera Specific Commands
  bool AuxAisGetPlate(                                                // RETR: result of operation -- have response
    float                                                             //
      *PlateArray,                                                    // PASS: -> 2xn array for plate coords
    int                                                               //
      *PlateCount                                                     // PASS: -> Return plate count
  );


  bool AuxAisProcessImg(                                               // RETR: result of operation -- have response
    int                                                                //
      alignIdx                                                         // PASS: alignment reference index 0-6
  );

  bool AuxAisStartCapture(                                             // RETR: result of operation -- have response
    int                                                                //
      alignIdx                                                         // PASS: alignment reference index 0-6
  );

  bool AuxAisReset(                                                    // RETR: result of operation -- have response
    void
  );
	
	bool AuxAisSaveCalibratedCentroid(int x, int y);
	bool AuxAisReadCalibratedCentroid(int* x, int* y);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  #pragma mark - AUX Motor Controller Specific Commands
  bool AuxMcCmdGotoFast(
    tAuxAddr
      Addr,
    int
      EncPositionBAM
  );

  bool AuxMcCmdGotoSlow(
    tAuxAddr
      Addr,
    int
      EncPositionBAM
  );

  bool AuxMcCmdMoveSwitch(
    tAuxAddr
      Addr
  );

  bool AuxMcCmdSlew(
    tAuxAddr
      Addr,
    tSlewRate
      SlewRate,
    tAuxCommand
      SlewCommand
  );

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  #pragma mark - AUX Motor Controller Querry Commands
  bool AuxMcQryGotoIsOver(                                            // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address
    bool                                                              //
      *isOver                                                         // PASS: -> to accept goto is over setting
  );                                                                  //



  bool AuxMcQryModel(                                                 // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address
    tTelescopeModel                                                   //
      *model                                                          // PASS: -> to accept model number
  ) ;



  bool AuxMcQryMoveSwitchIsOver(                                      // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address
    bool                                                              //
      *isOver                                                         // PASS: -> to accept move switch is over setting
  );                                                                  //


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  #pragma mark - AUX Motor Controller Get / Set Commands
  bool AuxMcApproachDirGet(                                           // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address
    tDirection                                                        //
      *dir                                                            // PASS: -> to accept approach direction value
  );                                                                  //



  bool AuxMcApproachDirSet(                                           // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address
    tDirection                                                        //
      dir                                                             // PASS: approach direction
  );                                                                  //



  bool AuxMcBacklashGet(                                              // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address
    tDirection                                                        //
      dir,                                                            // PASS: direction
    int                                                               //
      *value                                                          // PASS: -> to accept blacklash setting result
  );                                                                  //



  bool AuxMcBacklashSet(                                              // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address
    tDirection                                                        //
      dir,                                                            // PASS: direction
    int                                                               //
      value                                                           // PASS: backlash value (0-99)
  );                                                                  //



  bool AuxMcCordwrapPosGet(                                           // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address
    int                                                               //
      *position                                                       // PASS: -> cordwrap result
  );                                                                  //



  bool AuxMcCordwrapPosSet(                                           // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address (axis address)
    int                                                               //
      position                                                        // PASS: cordwrap position
  );



  bool AuxMcCordwrapEnaGet(                                           // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address
    bool                                                              //
      *enabled                                                        // PASS: -> to accept cordwrap enabled setting
  );                                                                  //



  bool AuxMcCordwrapEnaSet(                                           // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address
    bool                                                              //
      enabled                                                         // PASS: cordwrap enabled setting
  );                                                                  //



  bool AuxMcCustomRate9EnaGet(                                        // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address
    bool                                                              //
      *enabled                                                        // PASS: -> to accept custom rate enable value
  );                                                                  //



  bool AuxMcCustomRate9EnaSet(                                        // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address
    bool                                                              //
      enabled                                                         // PASS: custom rate enable
  );                                                                  //



  bool AuxMcCustomRate9Get(                                           // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address
    int                                                               //
      *value,                                                         // PASS: -> to accept custome rate9 value
    int                                                               //
      *maxValue                                                       // PASS: the max rate9 value
  );                                                                  //



  bool AuxMcCustomRate9Set(                                           // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address
    int                                                               //
      value                                                           // PASS: custom rate9 value
  );                                                                  //



  bool AuxMcEncPositionGet(                                           // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address
    int                                                               //
      *encPositionBAM                                                 // PASS: -> int BAM position to accept result
  );                                                                  //



  bool AuxMcEncPositionSet(                                           // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address (axis address)
    int                                                               //
      encPositionBAM                                                  // PASS: encoder position
  );



  bool AuxMcGuideRateGet(                                             // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address
    int                                                               //
      *value                                                          // PASS: -> to accept blacklash setting result
  );                                                                  //



  bool AuxMcGuideRateSet(                                             // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address
    int                                                               //
      value                                                           // PASS: backlash value
  );                                                                  //



  bool AuxMcRaLimitGet(                                               // RETR: result of operation -- have response
    tLimit                                                            //
      Limit,                                                          // PASS: which limit (min or max) to get
    int                                                               //
      *value                                                          // PASS: -> to accept blacklash setting result
  );                                                                  //



  bool AuxMcRaLimitSet(                                               // RETR: result of operation -- have response
    tLimit                                                            //
      Limit,                                                          // PASS: which limit (min or max) to get
    int                                                               //
      value                                                           // PASS: backlash value
  );                                                                  //



  bool AuxMcRaLimitEnaGet(                                            // RETR: result of operation -- have response
    bool                                                              //
      *enabled                                                        // PASS: -> to accept blacklash setting result
  );                                                                  //



  bool AuxMcRaLimitEnaSet(                                            // RETR: result of operation -- have response
    bool                                                              //
      enabled                                                         // PASS: backlash value
  );                                                                  //



  bool AuxMcTrackingSet(                                              // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address
    double                                                            //
      rateAsec                                                        // PASS: tracking rate in arcsec/s
  );



  bool AuxMcTrackingSet(                                              // RETR: result of operation -- have response
    tAuxAddr                                                          //
      addr,                                                           // PASS: device address
    tDirection                                                        //
      dir,                                                            // PASS: direction
    tTrackingRate                                                     //
      mode                                                            // PASS: tracking mode (enumerated value)
  );

  bool AuxCevoBattStatusGet(                                          // RETR: result of operation -- have response
    float                                                             //
      *battVolts,                                                     // PASS: -> rcv battery voltage
    tBattLvl                                                          //
      *battLvl,                                                       // PASS: -> rcv battery level LOW, MED, HIGH
    tBattStatus                                                       //
      *battStatus                                                     // PASS: -> rcv battery status Charging, Charged,
  );                                                                  //       Discharging, Fault

  bool AuxLightingGet(                                            // RETR: result of operation -- have response
    tLamp                                                             //
    lamp,                                                           // PASS: which 'lamp'
    int                                                               //
    *value                                                          // PASS: -> rcv lamp brightness setting
  );                                                                  //

  bool AuxLightingSet(                                            // RETR: result of operation -- have response
    tLamp                                                             //
    lamp,                                                           // PASS: which 'lamp'
    int                                                               //
    value                                                           // PASS: lamp brightness 0-10
  );                                                                  //

  bool AuxCevoDcpChargeModeGet(                                       // RETR: result of operation -- have response
    tDcpChargeMode                                                    //
      *mode                                                           // PASS: -> rec charge mode setting
  );                                                                  //

  bool AuxCevoDcpChargeModeSet(                                       // RETR: result of operation -- have response
    tDcpChargeMode                                                    //
      mode                                                            // PASS: charge mode: AUTO or ON
  );                                                                  //

  bool AuxCevoPwrMaxCurrentGet(                                       // RETR: result of operation -- have response
    float                                                             //
      *maxCurrent                                                     // PASS: -> rcv max current setting
  );                                                                  //

  bool AuxCevoPwrMaxCurrentSet(                                       // RETR: result of operation -- have response
    float                                                             //
      maxCurrent                                                      // PASS: max current setting 2.0 - 5.0
  );                                                                  //

};


#endif //ifdef WIFLY_HPP

