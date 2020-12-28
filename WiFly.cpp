 //
//  WiFly.cpp
//  Created by Danyal Medley on 10/7/13.
//
#define WIFLY_CPP
#include "DmTypes.hpp"
#include "CTelescope.h"
#include "AuxSend.hpp"
#include "AuxRecv.hpp"
#include "WiFly.hpp"

//#define PI 3.1415926535897932

const uint8_t lightingLevelLUT[11] = {
  0, 8, 12, 17, 25, 37, 55, 80, 119, 174, 255
};

uint8_t
  exBuffer[1024];

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma mark - AUX Operations
//==============================================================================================================================
// WiFlyParseAux --
//==============================================================================================================================
uint8_t  *WiFly::ParseAux (                                                //
  uint8_t                                                               //
    *recvBuff,                                                        //
  int                                                                 //
    len,                                                              //
  tAuxAddr                                                            //
    addr,                                                             //
  tAuxCommand                                                         //
    cmd                                                               //
)
{
  uint8_t                                                               //
    *dirPtr = recvBuff;                                               // set up a scanning pointer to receive buffer
  uint8_t                                                               //
    *pktdata;                                                         // points to start of packet

  mRecvResp = NULLPTR_AUXRECV;                                        // always reset this -- or it uses prev. packet
  while(1) {                                                          // while: scanning the recv data
    pktdata = ((uint8_t *)0x0);                                          //   init packet pointer
    long remaining = (recvBuff + len - dirPtr);                       //   set remaining count for scanning
    while ((recvBuff + len - dirPtr) > 0) {                           //   while more to scan
      char somChar = *dirPtr++;                                       //     check for start of message byte 0x3b
      if(somChar == '\x3b') {                                         //     if: we found it, next byte is the length
                                                                      //       of the packet. If the full packet is
        char lenChar = *dirPtr;                                       //       cont'd in buffer, then set pkt ptr and
        if(lenChar + 3 <= remaining) {                                //       break out of loop
          pktdata = dirPtr;                                           //
          break;                                                      //
        }                                                             //
      }                                                               //
    }                                                                 //
    if(pktdata == NULLPTR) {                                          //   if: we didn't get full packet, or not found
      break;                                                          //     exit from loop, and wait for more data.
    }                                                                 //     this is the exit mech. from while loop
    else {                                                            //   else: (found a full packet)
      mRecvResp = new AuxRecv(pktdata);                               //     allocate an AuxRecv object to accept the
      if (mRecvResp->wasValidated) {                                  //     if: valid packet
        if (                                                          //
					mRecvResp->mDstAddr != kAddrPhone ||
          mRecvResp->mSrcAddr != addr ||                              //
          mRecvResp->mCmd != cmd                                      //
        ) {                                                           //       if: correct address
          dirPtr += mRecvResp->mLen + 2;                              //
          mRecvResp = NULLPTR_AUXRECV;                                //
        }                                                             //
        else {                                                        //
          dirPtr += mRecvResp->mLen + 2;                              //
          break;                                                      //       NOTE: we're not building buffers of packets
        }                                                             //         anymore, so exit after first successful
      }                                                               //         decode. It should be the correct resp.
      else {                                                          //     else: return a null packet
        mRecvResp = NULLPTR_AUXRECV;                                  //
      }                                                               //
    }                                                                 //
  }                                                                   // loop continues while there is data to scan
  return dirPtr;
}



//==============================================================================================================================
// WiFlyParseAux --
//==============================================================================================================================
void WiFly::ParseAuxEx (                                              //
  uint8_t                                                               //
    *recvBuff,                                                        //
  int                                                                 //
    len                                                               //
)
{
  uint8_t                                                               //
    *dirPtr = recvBuff;                                               // set up a scanning pointer to= receive buffer
  int                                                                 //
    pktLen;                                                           //
  long                                                                //
    remaining,
    recvBuffLast;

  recvBuffLast = (long)recvBuff + len;                                //
  pktLen = 0;                                                         //
  while(1) {                                                          // while: scanning the recv data
    remaining = (recvBuff + len - dirPtr);                            //   set remaining count for scanning
    while (remaining > 0) {                                           //   while more to scan
      char somChar = *dirPtr++;                                       //     check for start of EX message byte 0x3c
      if(somChar == '\x3c') {                                         //     if: we found it, next byte is the length
        if(remaining >= 4) {                                          //
          pktLen = *dirPtr++ << 24;                                   //
          pktLen += *dirPtr++ << 16;                                  //
          pktLen += *dirPtr++ << 8;                                   //
          pktLen += *dirPtr++;                                        //
          break;                                                      //
        }                                                             //
      }                                                               //
    }                                                                 //
    if(pktLen == 0) {                                                 //   if: we didn't get full packet, or not found
      break;                                                          //     exit from loop, and wait for more data.
    }                                                                 //     this is the exit mech. from while loop
    else {                                                            //   else: (found a full packet)
      uint8_t                                                           //
        *exBuffPtr = exBuffer;                                        //
      memset(exBuffer, 0x00, sizeof(exBuffer));                       //
      pktLen++;                                                       //   makes sure to copy last bite
      while(pktLen-- && (long)dirPtr < recvBuffLast) {                //
        *exBuffPtr++ = *dirPtr;                                       //
        if(*dirPtr++ == '\x3b') {                                     //
          if(*dirPtr < 3) {                                           //
            dirPtr+= 2;                                               //
          }                                                           //
        }                                                             //
      }                                                               //
      break;                                                          //
    }                                                                 //
  }                                                                   // loop continues while there is data to scan
}



//==============================================================================================================================
// AuxPacketMaster --
//   formwats an AUX packet, sends to AUX bus, the response is captured, and parsed into (new AuxRecv) *mRecvResp
//==============================================================================================================================
void WiFly::AuxPacketMaster(                                          // RETR: NONE
  tAuxAddr                                                            //
    addr,                                                             // PASS: enumerated address
  tAuxCommand                                                         //
    cmd,                                                              // PASS: AUX command
  uint8_t                                                               //
    *data,                                                            // PASS: -> data buffer (or NULL)
  int                                                                 //
    dlen                                                              // PASS: -> data length (or 0)
)
{
  uint8_t
    recvBuff[1024] = { 0 };//,
    //*dirPtr;
  memset(recvBuff, 0xff, 1024);
  AuxSend msg = AuxSend((uint8_t)cmd, dlen, data, (uint8_t)kAddrPhone, (uint8_t)addr);

  int
	timeout = (mTelescope->moduleType == MODULE_TYPE_TCELL) ? 3000 : 150; // ms  the TCell modules experienced periodic long pauses in communication
  if(addr == kAddrAis)
  {
	/* 
			The SS camera response time is fairly long, especially if it is
			any of the following commands which can take over 6 seconds to complete
	 */
		//if(cmd == kAis_GetPlateList || cmd == kAis_StartCapture || cmd == kAis_ProcessImage)
		if(addr == kAddrAis) // starsense camera protocal sucks, wait a while for all commands to respond
			timeout = 6500; // ms
		else
			timeout = 3000;
  }
	
	int err = mTelescope->CelAuxCmdResp(&msg, &mRecvResp, exBuffer, timeout);
	if(err != TEL_NO_ERROR)
		mRecvResp = NULLPTR_AUXRECV;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma mark - AUX Device common commands

//==============================================================================================================================
// AuxDevGetVer --
//
//   get version information for device (not just motor control, but any device)
//==============================================================================================================================
bool WiFly::AuxDevGetVer(                                             // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  sMcVersion                                                          //
    *verInfo                                                          // PASS: -> to version info object to receive result
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response

  AuxPacketMaster(addr, kCmd_GetVer, 0x0, 0);                         // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       // if: we have a valid response
    uint8_t *data = mRecvResp->mDataPtr;                                //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 2) {                                    //   if: an old 2-byte version
      verInfo->VerMaj = (int)data[0];                                 //     record the two-byte version info
      verInfo->VerMin = (int)data[1];                                 //
      haveResponse = true;                                            //     indicate that we have result
    }                                                                 //
    else if(mRecvResp->mDataLen == 4) {                               //   else if: newer 4-byte version (HcPro, cam)
      verInfo->VerMaj = (int)data[0];                                 //     recrod the four-byte version info
      verInfo->VerMin = (int)data[1];                                 //
      verInfo->VerBld = (int)data[2] * 256 + (int)data[3];            //
      haveResponse = true;                                            //     indicate that we have result
    }                                                                 //
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma mark - AIS Camera Specific Commands
bool WiFly::AuxAisGetPlate(                                           // RETR: result of operation -- have response
  float                                                               //
    *PlateArray,                                                      // PASS: -> 2xn array for plate coords
  int                                                                 //
    *PlateCount                                                       // PASS: -> Return plate count
)
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response
  int
    thePlateCount = 0;

  *PlateCount = 0;                                                    //
  PlateArray[0] = 0.0;                                                //
  PlateArray[1] = 0.0;                                                //
  AuxPacketMaster(kAddrAis, kAis_GetPlateList, 0x0, 0);                //


#if false
  printf("mLen:     %02x\r\n", mRecvResp->mLen);
  printf("mSrcAddr: %02x\r\n", mRecvResp->mSrcAddr);
  printf("mDstAddr: %02x\r\n", mRecvResp->mDstAddr);
  printf("mCmd:     %02x\r\n", mRecvResp->mCmd);
  printf("mDataLen: %02x\r\n", mRecvResp->mDataLen);
  for(int i=0; i < mRecvResp->mDataLen; i++)
  {
    printf("  Data[%02d] %02x\r\n", i, mRecvResp->mDataPtr[i]);
  }
  printf("mCsum:    %02x\r\n", mRecvResp->mCsum);
#endif

  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       // if: we have a valid response
    uint8_t *data = mRecvResp->mDataPtr;                                //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 2) {                                    //   if: an old 2-byte version
      if((data[0] & AIS_STA_RDYPLATE) == AIS_STA_RDYPLATE) {          //
        thePlateCount = (int)data[1];                                   //     record the two-byte version info
        haveResponse = true;                                          //     indicate that we have result
      }                                                               //
    }                                                                 //
  }                                                                   //
  if(haveResponse) {
    memcpy(PlateArray, exBuffer, thePlateCount * 8);
    *PlateCount = thePlateCount;                                   //     record the two-byte version info
  }
  return haveResponse;
}



bool WiFly::AuxAisProcessImg(                                         // RETR: result of operation -- have response
  int                                                                 //
    alignIdx                                                          // PASS: alignment reference index 0-6
)
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response
  uint8_t                                                               //
    buff[1];                                                          // transmit data for BAM position (24-bit pos)

  if(alignIdx > 6) {                                                  // Only enough image slots for 7 (0-6), so image idx 7-n
    alignIdx = 6;                                                     //   uses the last
  }                                                                   //

  buff[0] = (uint8_t)(alignIdx);                                        // set the align reference image index

  AuxPacketMaster(kAddrAis, kAis_ProcessImage, &buff[0], 1);          // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       // if: we have a valid response
    haveResponse = true;                                              //   flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}



bool WiFly::AuxAisStartCapture(                                       // RETR: result of operation -- have response
  int                                                                 //
    alignIdx                                                          // PASS: alignment reference index 0-6
)
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response
  uint8_t                                                               //
    buff[3];                                                          // transmit data for BAM position (24-bit pos)

  if(alignIdx > 6) {                                                  // Only enough image slots for 7 (0-6)
    alignIdx = 6;                                                     //
  }                                                                   //

  buff[0] = 0x03;                                                     // always 1000ms
  buff[1] = 0xe8;                                                     //
  buff[2] = (uint8_t)(alignIdx);                                        // set the align reference image index

  AuxPacketMaster(kAddrAis, kAis_StartCapture, &buff[0], 3);          // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       // if: we have a valid response
    haveResponse = true;                                              //   flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}



bool WiFly::AuxAisReset(
    void
)
{
  bool                                                                //
      haveResponse = false;                                           // flags if we have valid response

  AuxPacketMaster(kAddrAis, kAis_Reset, 0x0, 0);                      // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       // if: we have a valid response
    haveResponse = true;                                              //   flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}



bool WiFly::AuxAisSaveCalibratedCentroid(int x, int y) {
	uint8_t data[9];
	data[3] = (x >> 24) & 0xFF;
	data[2] = (x >> 16) & 0xFF;
	data[1] = (x >> 8) & 0xFF;
	data[0] = (x >> 0) & 0xFF;
	data[7] = (y >> 24) & 0xFF;
	data[6] = (y >> 16) & 0xFF;
	data[5] = (y >> 8) & 0xFF;
	data[4] = (y >> 0) & 0xFF;
	data[8] = 0; // profile number... SkyPortal does not support multiple profiles so just keep this value as 0
	
	AuxPacketMaster(kAddrAis, kAis_SetCenterRef, data, 9);                      // get a packet from the packet master send proc
	return mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated;
}



bool WiFly::AuxAisReadCalibratedCentroid(int* x, int* y) {
	uint8_t profileNumber = 0;
	AuxPacketMaster(kAddrAis, kAis_GetCenterRef, &profileNumber, 1);                      // get a packet from the packet master send proc
	if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)
	{
		uint8_t *data = mRecvResp->mDataPtr;                                //   setup ptr to data in packet
		if(mRecvResp->mDataLen == 8) {
			*x = (data[3] << 24) | (data[2] << 16) | (data[1] << 8) | data[0];
			*y = (data[7] << 24) | (data[6] << 16) | (data[5] << 8) | data[4];
			return true;
		}
		else return false;
	}
	else return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma mark - AUX Motor Controller Specific Commands
//==============================================================================================================================
// AuxMcCmdGotoFast --
//
//   Sends the AUX command for fast goto operation
//==============================================================================================================================
bool WiFly::AuxMcCmdGotoFast(                                         // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address (axis address)
  int                                                                //
    encPositionBAM                                                    // PASS: encoder position
)
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response
  uint8_t                                                               //
    buff[3];                                                          // transmit data for BAM position (24-bit pos)

  buff[2] = (uint8_t)(encPositionBAM % 256);                            //   in big-endian order by converting the int
  encPositionBAM /= 256;                                              //   position into a three-byte value
  buff[1] = (uint8_t)(encPositionBAM % 256);                            //
  encPositionBAM /= 256;                                              //
  buff[0] = (uint8_t)(encPositionBAM % 256);                            //

  AuxPacketMaster(addr, kMtr_GotoFast, &buff[0], 3);                  // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       // if: we have a valid response
    haveResponse = true;                                              //   flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}



//==============================================================================================================================
// AuxMcCmdGotoSlow --
//
//   Sends the AUX command for slow goto operation
//==============================================================================================================================
bool WiFly::AuxMcCmdGotoSlow(                                         // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address (axis address)
  int                                                                 //
    encPositionBAM                                                    // PASS: encoder position
)
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response
  uint8_t                                                               //
    buff[3];                                                          // transmit data for BAM position (24-bit pos)

  buff[2] = (uint8_t)(encPositionBAM % 256);                            //   in big-endian order by converting the int
  encPositionBAM /= 256;                                              //   position into a three-byte value
  buff[1] = (uint8_t)(encPositionBAM % 256);                            //
  encPositionBAM /= 256;                                              //
  buff[0] = (uint8_t)(encPositionBAM % 256);                            //

  AuxPacketMaster(addr, kMtr_GotoSlow, &buff[0], 3);                  // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       // if: we have a valid response
    haveResponse = true;                                              //   flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}



//==============================================================================================================================
// AuxMcCmdMoveSwitch --
//
//   Command axis to move to switch position
//==============================================================================================================================
bool WiFly::AuxMcCmdMoveSwitch(                                       // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr                                                              // PASS: device address
)
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response

  AuxPacketMaster(addr, kMtr_MoveSwitch, 0x0, 0);                     // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       // if: we have a valid response
    haveResponse = true;                                              //   flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}



//==============================================================================================================================
// AuxMcCmdSlew --
//   Interface to MC for slew command
//==============================================================================================================================
bool WiFly::AuxMcCmdSlew(                                             // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  tSlewRate                                                           //
    slewRate,                                                         // PASS: enumerated slew speed
  tAuxCommand                                                         //
    slewCommand                                                       // PASS: enumerated slew command
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response

  AuxPacketMaster(addr, slewCommand, (uint8_t *)&slewRate, 1);          // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       // if: we have a valid response
    haveResponse = true;                                              //   flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma mark - AUX Motor Controller Querry Commands
//==============================================================================================================================
// AuxMcGotoIsOver --
//   get axis position
//==============================================================================================================================
bool WiFly::AuxMcQryGotoIsOver(                                       // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  bool                                                                //
    *isOver                                                           // PASS: -> to accept goto is over setting
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response

  AuxPacketMaster(addr, kMtr_GotoOver, 0x0, 0);                       // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    unsigned char *data = mRecvResp->mDataPtr;                        //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 1)                                      //   if: correct response length
    {                                                                 //
      *isOver = data[0] == 0x0 ? false : true;                        //     set the return value to response
      haveResponse = true;                                            //     indicate that we have result
    }                                                                 //
  }                                                                   //
  return haveResponse;                                                //
}



//==============================================================================================================================
// AuxMcGetModel --
//   get telescope model information for device (not just motor control, but any device)
//==============================================================================================================================
bool WiFly::AuxMcQryModel(                                            // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  tTelescopeModel                                                     //
    *model                                                            // PASS: -> to accept model number
)
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response

  AuxPacketMaster(addr, kMtr_GetModel, 0x0, 0);                       // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    unsigned char *data = mRecvResp->mDataPtr;                        //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 2)                                      //   if: correct response length
    {                                                                 //
      *model = (tTelescopeModel)data[0];                              //     set the return value to response and
      haveResponse = true;                                            //     indicate that we have result
    }                                                                 //
  }                                                                   //
  return haveResponse;                                                //
}



//==============================================================================================================================
// AuxMcGetMoveSwitchOverState --
//   check if move to switch operation is over
//==============================================================================================================================
bool WiFly::AuxMcQryMoveSwitchIsOver(                                 // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  bool                                                                //
    *isOver                                                           // PASS: -> to accept move switch is over setting
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response

  AuxPacketMaster(addr, kMtr_MoveSwitchOver, 0x0, 0);                 // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    unsigned char *data = mRecvResp->mDataPtr;                        //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 1)                                      //   if: correct response length
    {                                                                 //
      *isOver = data[0] == 0x0 ? false : true;                        //     set the return value to response
      haveResponse = true;                                            //     indicate that we have result
    }                                                                 //
  }                                                                   //
  return haveResponse;                                                //
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma mark - AUX Motor Controller Get / Set Commands
//==============================================================================================================================
// AuxMcApproachDirGet --
//   get goto approach direction
//==============================================================================================================================
bool WiFly::AuxMcApproachDirGet(                                      // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  tDirection                                                          //
    *dir                                                              // PASS: -> to accept approach direction value
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response

  AuxPacketMaster(addr, kMtr_ApproachDirGet, 0x0, 0);                 // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    unsigned char *data = mRecvResp->mDataPtr;                        //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 1)                                      //   if: correct response length
    {                                                                 //
      *dir = data[0] == 0 ? kPos : kNeg;                              //     set the return value to response
      haveResponse = true;                                            //     indicate that we have result
    }                                                                 //
  }                                                                   //
  return haveResponse;                                                //
}



//==============================================================================================================================
// AuxMcApproachDirSet --
//   set goto approach direction
//==============================================================================================================================
bool WiFly::AuxMcApproachDirSet(                                      // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  tDirection                                                          //
    dir                                                               // PASS: approach direction
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid result
  uint8_t                                                               //
    data;                                                             // hold the single byte of data

  data = dir == kPos ? 0x0 : 0x1;                                     // determine data for specified direction
  AuxPacketMaster(addr, kMtr_ApproachDirSet, &data, 1);               // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    haveResponse = true;                                              //   flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}



//==============================================================================================================================
// AuxMcBacklashGet --
//   get telescope backlash setting  (0-99)
//==============================================================================================================================
bool WiFly::AuxMcBacklashGet(                                         // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  tDirection                                                          //
    dir,                                                              // PASS: direction
  int                                                                 //
    *value                                                            // PASS: -> to accept blacklash setting result
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response
  tAuxCommand                                                         //
    cmd;                                                              // hold the command -- pos or neg to get backlash

  cmd = dir == kPos ? kMtr_BacklashPosGet : kMtr_BacklashNegGet;      // prep the command based on flag
  AuxPacketMaster(addr, cmd, 0x0, 0);                                 // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       // if: we have a valid response
    unsigned char *data = mRecvResp->mDataPtr;                        //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 1) {                                    //   if: correct response length
      *value = (int)data[0]; //(int)((double)data[0] * (99.0 / 255.0) + 0.5);         //     set the return value to response
      haveResponse = true;                                            //     indicate that we have result
    }                                                                 //
  }                                                                   //
  return haveResponse;                                                //
}




//==============================================================================================================================
// AuxMcBacklashSet --
//   set telescope backlash setting  (0-99)
//==============================================================================================================================
bool WiFly::AuxMcBacklashSet(                                         // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  tDirection                                                          //
    dir,                                                              // PASS: direction
  int                                                                 //
    value                                                             // PASS: backlash value (0-99)
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid result
  tAuxCommand                                                         //
    cmd;                                                              // hold the command -- pos or neg to set backlash
  uint8_t                                                               //
    data;                                                             // hold the single byte of data

  if(value <= 99) {                                                   // if: valid range for value
    cmd = dir == kPos ? kMtr_BacklashPosSet : kMtr_BacklashNegSet;    //   prep the command based on flag
    data = (uint8_t)value; //((double)value *(255.0 / 99.0) + 0.5);              //   scale the 0-99 value to 0-255
    AuxPacketMaster(addr, cmd, &data, 1);                             //   get a packet from the packet master send proc
    if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {     //   if: we have a valid response
      haveResponse = true;                                            //     flag that we have the result
    }                                                                 //
  }
  return haveResponse;                                                // return the result of operation
}




//==============================================================================================================================
// AuxMcCordwrapPosGet --
//   get coordwrap position
//==============================================================================================================================
bool WiFly::AuxMcCordwrapPosGet(                                      // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  int                                                                 //
    *position                                                         // PASS: -> to accept cordwrap position
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response

  AuxPacketMaster(addr, kMtr_CordwrapPosGet, 0x0, 0);                 // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       // if: we have a valid response
    unsigned char *data = mRecvResp->mDataPtr;                        //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 3) {                                    //   if: correct response format length
      int                                                             //
        BAM;                                                          //
      BAM = (int)data[0];                                             //
      BAM = BAM * 256 + (int)data[1];                                 //
      BAM = BAM * 256 + (int)data[2];                                 //
      *position = BAM;                                                //
      haveResponse = true;                                            //   indicate that we have result
    }                                                                 //
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}



//==============================================================================================================================
// AuxMcCordwrapPosSet
//   set the cordwrap position
//==============================================================================================================================
bool WiFly::AuxMcCordwrapPosSet(                                      // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address (axis address)
  int                                                                 //
    position                                                          // PASS: cordwrap position
)
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response
  uint8_t                                                               //
    buff[3];                                                          // transmit data for BAM position (24-bit pos)

  buff[2] = (uint8_t)(position % 256);                                  // in big-endian order by converting the int
  position /= 256;                                                    //   position into a three-byte value
  buff[1] = (uint8_t)(position % 256);                                  //
  position /= 256;                                                    //
  buff[0] = (uint8_t)(position % 256);                                  //

  AuxPacketMaster(addr, kMtr_CordwrapPosSet, &buff[0], 3);            // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       // if: we have a valid response
    haveResponse = true;                                              //   flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}



//==============================================================================================================================
// AuxMcCordwrapEnaGet --
//   get cordwrap enabled setting
//==============================================================================================================================
bool WiFly::AuxMcCordwrapEnaGet(                                      // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  bool                                                                //
    *enabled                                                          // PASS: -> to accept cordwrap enabled setting
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response

  AuxPacketMaster(addr, kMtr_CordwrapEnaGet, 0x0, 0);                 // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    unsigned char *data = mRecvResp->mDataPtr;                        //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 1)                                      //   if: correct response length
    {                                                                 //
      *enabled = data[0] == 0x0 ? false : true;                       //     set the return value to response
      haveResponse = true;                                            //     indicate that we have result
    }                                                                 //
  }                                                                   //
  return haveResponse;                                                //
}



//==============================================================================================================================
// AuxMcCordwrapEnaSet --
//   enable/disable cordwrap
//==============================================================================================================================
bool WiFly::AuxMcCordwrapEnaSet(                                      // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  bool                                                                //
    enabled                                                           // PASS: cordwrap enabled setting
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid result
  tAuxCommand                                                         //
    cmd;                                                              // the cordwrap command based on setting

  cmd = enabled == false ? kMtr_CordwrapOff : kMtr_CordwrapOn;        // determine which cordwrap command to use
  AuxPacketMaster(addr, cmd, 0x0, 0);                                 // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    haveResponse = true;                                              //   flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}



 //==============================================================================================================================
// AuxMcCustomRate9Get --
//   get custom rate9 and max rate values
//==============================================================================================================================
bool WiFly::AuxMcCustomRate9Get(                                      // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  int                                                                 //
    *value,                                                           // PASS: -> to accept custome rate9 value
  int                                                                 //
    *maxValue                                                         // PASS: the max rate9 value
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response

  AuxPacketMaster(addr, kMtr_CustomRate9Get, 0x0, 0);                 // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       // if: we have a valid response
    uint8_t *data = mRecvResp->mDataPtr;                                //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 4) {                                    //   if: an original 1-byte response
      *value = (int)data[0]*256 + data[1];                            //     get the rate9 value
      if(maxValue != NULLPTR) {                                       //     if: valid ptr for maxvalue...
        *maxValue = (int)data[2] * 256 + (int)data[3];                //       get that value
      }                                                               //
      haveResponse = true;                                            //     indicate that we have result
    }                                                                 //
  }                                                                   //
  return haveResponse;                                                //
}



//==============================================================================================================================
// AuxMcCustomRate9Set
//   set the custom rate9 rate
//==============================================================================================================================
bool WiFly::AuxMcCustomRate9Set(                                      // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  int                                                                 //
    value                                                             // PASS: custom rate9 value
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid result
  uint8_t                                                               //
    data[2];                                                          // hold the single byte of data

  data[1] = (uint8_t)(value % 256);                                     //
  data[0] = (uint8_t)(value / 256);                                     //
  AuxPacketMaster(addr, kMtr_CustomRate9Set, data, 2);                // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    haveResponse = true;                                              //   flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation

}



//==============================================================================================================================
// AuxMcCustomRate9EnaGet --
//   get custom rate9 enabled setting
//==============================================================================================================================
bool WiFly::AuxMcCustomRate9EnaGet(                                   // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  bool                                                                //
    *enabled                                                          // PASS: -> to accept custom rate enable value
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response

  AuxPacketMaster(addr, kMtr_CustomRate9EnaGet , 0x0, 0);             // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    unsigned char *data = mRecvResp->mDataPtr;                        //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 1)                                      //   if: correct response length
    {                                                                 //
      *enabled = data[0] == 0x0 ? false : true;                       //     set the return value to response
      haveResponse = true;                                            //     indicate that we have result
    }                                                                 //
  }                                                                   //
  return haveResponse;                                                //
}



//==============================================================================================================================
// AuxMcCustomRate9EnaSet --
//   enable/disable custom rate9
//==============================================================================================================================
bool WiFly::AuxMcCustomRate9EnaSet(                                   // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  bool                                                                //
    enabled                                                           // PASS: custom rate enable
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid result
  uint8_t                                                               //
    data;                                                             // hold the single byte of data

  data = enabled == false ? 0x0 : 0x1;                                // scale the 0-99 value to 0-255
  AuxPacketMaster(addr, kMtr_CustomRate9EnaSet, &data, 1);            // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    haveResponse = true;                                              //   flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}



 //==============================================================================================================================
// AuxMcEncPositionGet --
//   get axis position
//==============================================================================================================================
bool WiFly::AuxMcEncPositionGet(                                      // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  int                                                                //
    *encPositionBAM                                                   // PASS: -> long BAM position to accept result
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response

  AuxPacketMaster(addr, kMtr_GetEncPosition, 0x0, 0);                 // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       // if: we have a valid response
    unsigned char *data = mRecvResp->mDataPtr;                        //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 3) {                                    //   if: correct response format length
      int                                                             //
        BAM;                                                          //
      BAM = (int)data[0];                                             //
      BAM = BAM * 256 + (int)data[1];                                 //
      BAM = BAM * 256 + (int)data[2];                                 //
      *encPositionBAM = BAM;                                          //
      haveResponse = true;                                            //   indicate that we have result
    }                                                                 //
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}



//==============================================================================================================================
// AuxMcEncPositionSet --
//   get axis position
//==============================================================================================================================
bool WiFly::AuxMcEncPositionSet(                                      // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address (axis address)
  int                                                                 //
    encPositionBAM                                                    // PASS: encoder position
)
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response
  uint8_t                                                               //
    buff[3];                                                          // transmit data for BAM position (24-bit pos)

  buff[2] = (uint8_t)(encPositionBAM % 256);                            // in big-endian order by converting the int
  encPositionBAM /= 256;                                              //   position into a three-byte value
  buff[1] = (uint8_t)(encPositionBAM % 256);                            //
  encPositionBAM /= 256;                                              //
  buff[0] = (uint8_t)(encPositionBAM % 256);                            //

  AuxPacketMaster(addr, kMtr_PositionSet, &buff[0], 3);               // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       // if: we have a valid response
    haveResponse = true;                                              //   flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}



//==============================================================================================================================
// AuxMcGuideRateGet --
//   get telescope guide rate (0-99)
//==============================================================================================================================
bool WiFly::AuxMcGuideRateGet(                                        // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  int                                                                 //
    *value                                                            // PASS: -> to accept guide rate value (0-99)
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response

  AuxPacketMaster(addr, kMtr_GuideRateGet, 0x0, 0);                   // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    unsigned char *data = mRecvResp->mDataPtr;                        //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 1)                                      //   if: correct response length
    {                                                                 //
      *value = (int)((double)data[0] * (99.0 / 255.0) + 0.5);         //     set the return value to response
      haveResponse = true;                                            //     indicate that we have result
    }                                                                 //
  }                                                                   //
  return haveResponse;                                                //
}



//==============================================================================================================================
// AuxMcGuideRateSet --
//   set telescope guide rate (0-99)
//==============================================================================================================================
bool WiFly::AuxMcGuideRateSet(                                        // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  int                                                                 //
    value                                                             // PASS: backlash value
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid result
  uint8_t                                                               //
    data;                                                             // hold the single byte of data

  data = (uint8_t)((double)value *(255.0 / 99.0) + 0.5);                // scale the 0-99 value to 0-255
  AuxPacketMaster(addr, kMtr_GuideRateSet, &data, 1);                 // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    haveResponse = true;                                              //   flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}




//==============================================================================================================================
// AuxMcRaLimitGet --
//   get software RA limit position
//==============================================================================================================================
bool WiFly::AuxMcRaLimitGet(                                          // RETR: result of operation -- have response
  tLimit                                                              //
    Limit,                                                            // PASS: which limit (min or max) to get
  int                                                                 //
    *value                                                            // PASS: -> to accept limit value
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response
  tAuxCommand
    cmd;

  cmd = Limit == kMin ? kMtr_RaLimitMinGet : kMtr_RaLimitMaxGet;      // determine which command to send
  AuxPacketMaster(kAddrAzm, cmd, 0x0, 0);                             // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       // if: we have a valid response
    uint8_t *data = mRecvResp->mDataPtr;                                //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 1) {                                    //   if: an original 1-byte response
      *value = (int)data[0] * 2;                                      //
      haveResponse = true;                                            //     indicate that we have result
    }                                                                 //
    else if(mRecvResp->mDataLen == 4) {                               //   else if: newer 4-byte version (HcPro, cam)
      float tmp = (float)data[1] + (float)data[0] * 256;              //
      float mcscale = (float)data[3] + (float)data[2] * 256;          //
      *value = (int)(tmp * 90 / mcscale + 0.5);                       //
      haveResponse = true;                                            //     indicate that we have result
    }                                                                 //
    if(Limit == kMax)  {                                              //
      *value = 180 - *value;                                          //
    }                                                                 //
    else if(*value > 270) *value -= 360;                              //
  }                                                                   //
  return haveResponse;                                                //
}



//==============================================================================================================================
// AuxMcRaLimitSet --
//   set telescope RA limit
//==============================================================================================================================
bool WiFly::AuxMcRaLimitSet(                                          // RETR: result of operation -- have response
  tLimit                                                              //
    Limit,                                                            // PASS: which limit (min or max) to get
  int                                                                 //
    value                                                             // PASS: the limit value
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response
  tAuxCommand                                                         //
    cmd;                                                              //
  int
    mcscale = 0;

  cmd = Limit == kMin ? kMtr_RaLimitMinGet : kMtr_RaLimitMaxGet;      // determine which command to send
  AuxPacketMaster(kAddrAzm, cmd, 0x0, 0);                             // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       // if: we have a valid response
    uint8_t *data = mRecvResp->mDataPtr;                                //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 1) {                                    //   if: an original 1-byte response
      mcscale = 45 * 256;
    }                                                                 //
    else if(mRecvResp->mDataLen == 4) {                               //   else if: newer 4-byte version (HcPro, cam)
      mcscale = (int)data[3] + (int)data[2] * 256;                    //
    }                                                                 //
  }                                                                   //
  cmd = Limit == kMin ? kMtr_RaLimitMinSet : kMtr_RaLimitMaxSet;      // determine which command to send
  int tmp = value;                                                    //
  if(Limit == kMax) {                                                 //
    tmp = 180 - tmp;                                                  //
  }                                                                   //
  else if(tmp < 0) {                                                  //
    tmp += 360;                                                       //
  }                                                                   //
  mcscale /= 90;
  tmp = tmp * mcscale;                                                    //
  uint8_t data[2];                                                      //
  data[0] = (uint8_t)(tmp / 256);                                       //
  data[1] = (uint8_t)(tmp % 256);                                       //
  AuxPacketMaster(kAddrAzm, cmd, data, 2);                            // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    haveResponse = true;                                              //   flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}




//==============================================================================================================================
// AuxMcRaLimitEnaGet --
//   get the current ra limit enabled setting
//==============================================================================================================================
bool WiFly::AuxMcRaLimitEnaGet(                                       // RETR: result of operation -- have response
  bool                                                                //
    *enabled                                                          // PASS: backlash value
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response

  AuxPacketMaster(kAddrAzm, kMtr_RaLimitEnaGet, 0x0, 0);              // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    unsigned char *data = mRecvResp->mDataPtr;                        //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 1)                                      //   if: correct response length
    {                                                                 //
      *enabled = data[0] == 0x0 ? false : true;                       //     set the return value to response
      haveResponse = true;                                            //     indicate that we have result
    }                                                                 //
  }                                                                   //
  return haveResponse;                                                //
}




//==============================================================================================================================
// AuxMcRaLimitEnaSet --
//   enable/disable telescope ra slew limits
//==============================================================================================================================
bool WiFly::AuxMcRaLimitEnaSet(                                       // RETR: result of operation -- have response
  bool                                                                //
    enabled                                                           // PASS: -> to accept guide rate value (0-99)
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid result
  uint8_t                                                               //
    data;                                                             // hold the single byte of data

  data = enabled == false ? 0x0 : 0x1;                                // 0x0=disabled, 0x1=enabled based on setting
  AuxPacketMaster(kAddrAzm, kMtr_RaLimitEnaSet, &data, 1);            // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    haveResponse = true;                                              //   flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}



//==============================================================================================================================
// AuxMcTrackingSet --
//   Set the tracking rate to the motor controller based on ArcSec/S rate value
//==============================================================================================================================
bool WiFly::AuxMcTrackingSet(                                         // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  double                                                              //
    rateAsec                                                          // PASS: tracking rate in arcsec/s
)
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response
  uint8_t                                                               //
    buff[3];                                                          // transmit data for BAM position (24-bit pos)
  tAuxCommand                                                         //
    cmd;                                                              // determine command based on sign of rate
  int                                                                 //
    rate;                                                             // int value of rate * 1024 rounded

  cmd = rateAsec < 0 ? kMtr_TrackNegative : kMtr_TrackPositive;       // set command based on pos or neg tracking
  rate = (int)(fabs(rateAsec) * 1024 + 0.5);                          //
  buff[2] = (uint8_t)(rate % 256);                                      // in big-endian order by converting the int
  rate /= 256;                                                        //   position into a three-byte value
  buff[1] = (uint8_t)(rate % 256);                                      //
  rate /= 256;                                                        //
  buff[0] = (uint8_t)(rate % 256);                                      //

  AuxPacketMaster(addr, cmd, &buff[0], 3);                            // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       // if: we have a valid response
    haveResponse = true;                                              //   flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}



//==============================================================================================================================
// AuxMcTrackingSet --
//   Set the tracking rate to the motor controller based on tracking mode (sidereal, solar, lunar)
//==============================================================================================================================
bool WiFly::AuxMcTrackingSet(                                         // RETR: result of operation -- have response
  tAuxAddr                                                            //
    addr,                                                             // PASS: device address
  tDirection                                                          //
    dir,                                                              // PASS: direction
  tTrackingRate                                                       //
    mode                                                              // PASS: tracking mode (enumerated value)
)
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response
  tAuxCommand                                                         //
    cmd;                                                              // command based on direction
  uint8_t                                                               //
    data[2];                                                          // tracking mode data value (enumerated value)


  cmd = dir == kPos ? kMtr_TrackPositive : kMtr_TrackNegative;        // set the tracking command based on direction setting
  data[0] = 255;                                                      //
  data[1] = 255-(uint8_t)mode;                                          // set tracking mode to enum value
  AuxPacketMaster(addr, cmd, data, 2);                                // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    haveResponse = true;                                              //   flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}



//==============================================================================================================================
// AuxCevoBattStatusGet --
//   get CEVO battery status
//==============================================================================================================================
bool WiFly::AuxCevoBattStatusGet(                                     // RETR: result of operation -- have response
  float                                                               //
    *battVolts,                                                       // PASS: -> rcv battery voltage
  tBattLvl                                                            //
    *battLvl,                                                         // PASS: -> rcv battery level LOW, MED, HIGH
  tBattStatus                                                         //
    *battStatus                                                       // PASS: -> rcv battery status Charging, Charged,
)                                                                     //       Discharging, Fault
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid response
  int
    centiVolts;

  AuxPacketMaster(kAddrPwr, kPwr_Status, 0x0, 0);                     // get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    unsigned char *data = mRecvResp->mDataPtr;                        //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 6)                                      //   if: correct response length
    {                                                                 //
      *battStatus = (tBattStatus)data[0];                             //     set the return value to response
      *battLvl = (tBattLvl)data[1];                                   //
      centiVolts =                                                    //
        data[2] << 24 | data[3] << 16 | data[4] << 8 | data[5];       //
      *battVolts = (float)centiVolts / 1000000.0;                     //
      haveResponse = true;                                            //     indicate that we have result
    }                                                                 //
  }                                                                   //
  return haveResponse;                                                //
}



//==============================================================================================================================
// AuxLightingGet --
//   get telescope lighting 0-10,
//==============================================================================================================================
bool WiFly::AuxLightingGet(                                           // RETR: result of operation -- have response
  tLamp                                                               //
    lamp,                                                             // PASS: which 'lamp'
  int                                                                 //
    *value                                                            // PASS: -> rcv tray lamp brightness setting
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid result

  if(lamp >= kNumEnd) {                                               //
    return haveResponse;                                              //
  }                                                                   //
  AuxPacketMaster(kAddrLmp, kLmp_LampBrightness, (uint8_t *)&lamp, 1);  //   get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    unsigned char *data = mRecvResp->mDataPtr;                        //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 1)                                      //   if: correct response length
    {                                                                 //
      int c;                                                          //
      c = 0;                                                          //
      while(++c < 11) {                                               //
        if(lightingLevelLUT[c] > (int)data[0]) {                      //
          break;                                                      //
        };                                                            //
      }                                                               //
      *value = c -1;                                                  //
      haveResponse = true;                                            //     indicate that we have result
    }                                                                 //
  }                                                                   //
  return haveResponse;                                                //
}




//==============================================================================================================================
// AuxLightingSet --
//   set telescope lighting 0-10 for one of the three 'lamps' WiFi, Logo, Tray
//==============================================================================================================================
bool WiFly::AuxLightingSet(                                           // RETR: result of operation -- have response
  tLamp                                                               //
    lamp,                                                             // PASS: which 'lamp'
  int                                                                 //
    value                                                             // PASS: lamp led brightness 0-10
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid result
  uint8_t                                                               //
    data[2];                                                          // hold the target "lamp" and brightness level

  if(lamp >= kNumEnd) {                                               //
    return haveResponse;                                              //
  }                                                                   //
  if(value >= 0 && value <= 10) {                                     // if: valid range for value
    data[0] = (uint8_t)lamp;                                            //   select the 'lamp' lighting
    data[1] = (uint8_t)(lightingLevelLUT[value]);                       //   set lighting level
    AuxPacketMaster(kAddrLmp, kLmp_LampBrightness, data, 2);          //   get a packet from the packet master send proc
    if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {     //   if: we have a valid response
      haveResponse = true;                                            //     flag that we have the result
    }                                                                 //
  }
  return haveResponse;                                                // return the result of operation
}




//==============================================================================================================================
// AuxLightingGet --
//   get telescope lighting 0-10,
//==============================================================================================================================
bool WiFly::AuxCevoDcpChargeModeGet(                                  // RETR: result of operation -- have response
  tDcpChargeMode                                                      //
    *mode                                                             // PASS: -> rec charge mode setting
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid result

  AuxPacketMaster(kAddrDcp, kDcp_Status, 0x0, 0);                     //   get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    unsigned char *data = mRecvResp->mDataPtr;                        //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 1)                                      //   if: correct response length
    {                                                                 //
      *mode = data[0] == 0x0 ? kChargePortAuto : kChargePortOn;       //     set the return value to response
      haveResponse = true;                                            //     indicate that we have result
    }                                                                 //
  }                                                                   //
  return haveResponse;                                                //
}




 //==============================================================================================================================
// AuxCevoDcpChargeModeSet --
//   set telescope USB charge mode AUTO or ON
//==============================================================================================================================
bool WiFly::AuxCevoDcpChargeModeSet(                                  // RETR: result of operation -- have response
  tDcpChargeMode                                                      //
    mode                                                              // PASS: charge mode: AUTO or ON
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid result
  uint8_t                                                               //
    data[1];                                                          // hold the target "lamp" and brightness level

  data[0] = mode == kChargePortAuto ? 0x0 : 0x1;                      //   select the logo lighting
  AuxPacketMaster(kAddrDcp, kDcp_Status, data, 1);                    //   get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       //   if: we have a valid response
    haveResponse = true;                                              //     flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}



//==============================================================================================================================
// AuxLightingGet --
//   get telescope max input current setting
//==============================================================================================================================
bool WiFly::AuxCevoPwrMaxCurrentGet(                                  // RETR: result of operation -- have response
  float                                                               //
    *maxCurrent                                                       // PASS: -> rcv max current setting
)                                                                     //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid result

  AuxPacketMaster(kAddrPwr, kPwr_CurrentLimit, 0x0, 0);               //   get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated)         // if: we have a valid response
  {                                                                   //
    unsigned char *data = mRecvResp->mDataPtr;                        //   setup ptr to data in packet
    if(mRecvResp->mDataLen == 2)                                      //   if: correct response length
    {                                                                 //
      *maxCurrent = (float)((data[0] << 8) + data[1]) / 1000;           //     set the return value to response
      haveResponse = true;                                            //     indicate that we have result
    }                                                                 //
  }                                                                   //
  return haveResponse;                                                //
}




 //==============================================================================================================================
// AuxCevoPwrMaxCurrentSet --
//   set max input current for CEVO
//==============================================================================================================================
bool WiFly::AuxCevoPwrMaxCurrentSet(                                  // RETR: result of operation -- have response
   float                                                              //
     maxCurrent                                                       // PASS: max current setting 2.0 - 5.0
 )                                                                    //
{
  bool                                                                //
    haveResponse = false;                                             // flags if we have valid result
  uint8_t                                                               //
    data[2];                                                          // hold the target "lamp" and brightness level
  uint16_t 
    maxCurrentVal;

  if(maxCurrent < 2.0) maxCurrent = 2.0;                              //
  if(maxCurrent > 5.0) maxCurrent = 5.0;                              //
  maxCurrentVal = maxCurrent * 1000;                                  //   select the logo lighting
  data[0] = (uint8_t)(maxCurrentVal >> 8);                              //
  data[1] = (uint8_t)maxCurrentVal;                                     //
  AuxPacketMaster(kAddrPwr, kPwr_CurrentLimit, data, 2);              //   get a packet from the packet master send proc
  if(mRecvResp != NULLPTR_AUXRECV && mRecvResp->wasValidated) {       //   if: we have a valid response
    haveResponse = true;                                              //     flag that we have the result
  }                                                                   //
  return haveResponse;                                                // return the result of operation
}





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//#pragma mark - Telescope Operations
//
////==============================================================================================================================
//// TelescopeInitAzmAlt --
////
////   Enables or disables the updating of telescope position / status
////==============================================================================================================================
//- (void) TelescopeInitAzmAlt
//{
//    [self AuxMcCordwrapPosSet :kAddrAzm :180 * DEGTOBAM];
//    [self AuxMcCordwrapEnaSet :kAddrAzm :true];
//    [self AuxMcCordwrapEnaSet :kAddrAlt :false];
//    [self AuxMcPositionSet :kAddrAzm :0 * DEGTOBAM];
//    [self AuxMcPositionSet :kAddrAlt :0 * DEGTOBAM];
//    [self AuxMcTrackingSet :kAddrAzm :0];
//    [self AuxMcTrackingSet :kAddrAlt :0];
//}
//
//
//
////==============================================================================================================================
//// TelescopeInitEquatorial --
////
////   Enables or disables the updating of telescope position / status
////==============================================================================================================================
//- (void) TelescopeInitEquatorial
//{
//    [self AuxMcCordwrapPosSet :kAddrAzm :(int)(270.0 * DEGTOBAM)];
//    [self AuxMcCordwrapPosSet :kAddrAlt :(int)(270.0 * DEGTOBAM)];
//    [self AuxMcCordwrapEnaSet :kAddrAzm :true];
//    [self AuxMcCordwrapEnaSet :kAddrAlt :true];
//    [self AuxMcPositionSet :kAddrAzm :(int)(90.0 * DEGTOBAM)];
//    [self AuxMcPositionSet :kAddrAlt :(int)(90.0 * DEGTOBAM)];
//}
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//#pragma mark - AUX Operations
//
////==============================================================================================================================
//// WiFlyPacketMasterSend
////
////   Master WiFly Module send procedure
////==============================================================================================================================
//-(bool) WiFlyPacketMasterSend :(unsigned char *)cmdDataPtr :(tWiFlyCmd)cmd :(NSString *)valstr
//{
//    DLog(@"WiFlyMasterSend");
//
//    int len = 0;
//
//    AuxSend *msg = [[AuxSend alloc] init];
//    msg.mCmd = cmd;
//    [msg commandString :cmdDataPtr];
//    //[self TcpConnectRequest];
//    //[self TcpTransaction :msg :mRecvBuff :&len :valstr];
//
//    bool haveResponse = [self WiFlyParseGetSettings :mRecvBuff :msg.mCmd];
//
//    msg = nil;
//    return haveResponse;
//}
//
//
////==============================================================================================================================
//// TelescopeRestoreSettingsFromMc
////
////   Retrieve / Restore the TelescopeSettings from motor controllor
////==============================================================================================================================
//- (void) TelescopeRestoreSettingsFromMc
//{
//    if(telescopeSettings.hasCustomRate9 == true)
//    {
//        telescopeSettings.customRate9AzmEnabled = mMotorControlSettings.customRate9AzmEnabled;
//        telescopeSettings.customRate9Azm = mMotorControlSettings.customRate9Azm;
//        telescopeSettings.customRate9AltEnabled = mMotorControlSettings.customRate9AltEnabled;
//        telescopeSettings.customRate9Alt = mMotorControlSettings.customRate9Alt;
//        telescopeSettings.customRate9MaxValue = mMotorControlSettings.customRate9MaxValue;
//    }
//    if(telescopeSettings.hasRaLimits == true)
//    {
//        telescopeSettings.raLimitsEnabled = mMotorControlSettings.raLimitsEnabled;
//        telescopeSettings.raLimitEast = mMotorControlSettings.raLimitEast;
//        telescopeSettings.raLimitWest = mMotorControlSettings.raLimitWest;
//    }
//    if(telescopeSettings.canDoEquatorial == true)
//    {
//        telescopeSettings.guideRateRa = mMotorControlSettings.guideRateRa;
//        telescopeSettings.guideRateDec = mMotorControlSettings.guideRateDec;
//    }
//    telescopeSettings.backlashAzmPos = mMotorControlSettings.backlashAzmPos;
//    telescopeSettings.backlashAzmNeg = mMotorControlSettings.backlashAzmNeg;
//    telescopeSettings.backlashAltPos = mMotorControlSettings.backlashAltPos;
//    telescopeSettings.backlashAltNeg = mMotorControlSettings.backlashAltNeg;
//    telescopeSettings.approachAzmPos = mMotorControlSettings.approachAzmPos;
//    telescopeSettings.approachAltPos = mMotorControlSettings.approachAltPos;
//    telescopeSettings.cordwrapEnabled = mMotorControlSettings.cordwrapEnabled;
//}
//
////==============================================================================================================================
//// TelescopeCommonSettingsGet
////
////   Update the MC settings from the TelescopeSettings before writing to MC
////==============================================================================================================================
//- (void) TelescopeUpdateMcSettings
//{
//    if(telescopeSettings.hasCustomRate9 == true)
//    {
//        mMotorControlSettings.customRate9AzmEnabled = telescopeSettings.customRate9AzmEnabled;
//        mMotorControlSettings.customRate9Azm = telescopeSettings.customRate9Azm;
//        mMotorControlSettings.customRate9AltEnabled = telescopeSettings.customRate9AltEnabled;
//        mMotorControlSettings.customRate9Alt = telescopeSettings.customRate9Alt;
//        mMotorControlSettings.customRate9MaxValue = telescopeSettings.customRate9MaxValue;
//    }
//    if(telescopeSettings.hasRaLimits == true)
//    {
//        mMotorControlSettings.raLimitsEnabled = telescopeSettings.raLimitsEnabled;
//        mMotorControlSettings.raLimitEast = telescopeSettings.raLimitEast;
//        mMotorControlSettings.raLimitWest = telescopeSettings.raLimitWest;
//    }
//    if(telescopeSettings.canDoEquatorial == true)
//    {
//        mMotorControlSettings.guideRateRa = telescopeSettings.guideRateRa;
//        mMotorControlSettings.guideRateDec = telescopeSettings.guideRateDec;
//    }
//    mMotorControlSettings.backlashAzmPos = telescopeSettings.backlashAzmPos;
//    mMotorControlSettings.backlashAzmNeg = telescopeSettings.backlashAzmNeg;
//    mMotorControlSettings.backlashAltPos = telescopeSettings.backlashAltPos;
//    mMotorControlSettings.backlashAltNeg = telescopeSettings.backlashAltNeg;
//    mMotorControlSettings.approachAzmPos = telescopeSettings.approachAzmPos;
//    mMotorControlSettings.approachAltPos = telescopeSettings.approachAltPos;
//
//    // Fixes crashing bug
//    telescopeSettings.slewLimitAltMaxDeg = 90;
//    telescopeSettings.slewLimitAltMinDeg = 0;
//
//}
//
//
//
////==============================================================================================================================
//// TelescopeCommonSettingsGet
////
////   Used for getting the common telescope settings
////==============================================================================================================================
//- (bool) TelescopeCommonSettingsGet
//{
//    bool gotTheSettings = false;
//    bool t1 = false, t2 = false, t3 = false, t4 = false, t5 = false, t6 = false, t7 = false;
//    int attempts = 5;
//    while (attempts--) {
//        if(!t1) t1 = t1 | [self AuxMcBacklashGet :kAddrAzm :true :&mMotorControlSettings.backlashAzmPos];
//        if(!t2) t2 = t2 | [self AuxMcBacklashGet :kAddrAzm :false :&mMotorControlSettings.backlashAzmNeg];
//        if(!t3) t3 = t3 | [self AuxMcBacklashGet :kAddrAlt :true :&mMotorControlSettings.backlashAltPos];
//        if(!t4) t4 = t4 | [self AuxMcBacklashGet :kAddrAlt :false :&mMotorControlSettings.backlashAltNeg];
//        if(!t5) t5 = t5 | [self AuxMcApproachDirGet :kAddrAzm :&mMotorControlSettings.approachAzmPos];
//        if(!t6) t6 = t6 | [self AuxMcApproachDirGet :kAddrAlt :&mMotorControlSettings.approachAltPos];
//        if(!t7) t7 = t7 | [self AuxMcCordwrapEnaGet :kAddrAzm :&mMotorControlSettings.cordwrapEnabled];
//        if(
//           t1 && t2 && t3 && t4 && t5 && t6 && t7
//           )
//        {
//            gotTheSettings = true;
//            break;
//        }
//    }
//    return gotTheSettings;
//}
//
//
//
////==============================================================================================================================
//// TelescopeCommonSettingsSet
////
////   Used for saving the common telescope settings
////==============================================================================================================================
//- (bool) TelescopeCommonSettingsSet
//{
//    bool savedTheSettings = false;
//    bool t1 = false, t2 = false, t3 = false, t4 = false, t5 = false, t6 = false;
//    int attempts = 5;
//    while (attempts--) {
//        if(!t1) t1 = t1 | [self AuxMcBacklashSet :kAddrAzm :true :mMotorControlSettings.backlashAzmPos];
//        if(!t2) t2 = t2 | [self AuxMcBacklashSet :kAddrAzm :false :mMotorControlSettings.backlashAzmNeg];
//        if(!t3) t3 = t3 | [self AuxMcBacklashSet :kAddrAlt :true :mMotorControlSettings.backlashAltPos];
//        if(!t4) t4 = t4 | [self AuxMcBacklashSet :kAddrAlt :false :mMotorControlSettings.backlashAltNeg];
//        if(!t5) t5 = t5 | [self AuxMcApproachDirSet :kAddrAzm :mMotorControlSettings.approachAzmPos];
//        if(!t6) t6 = t6 | [self AuxMcApproachDirSet :kAddrAlt :mMotorControlSettings.approachAltPos];
//        if(
//           t1 && t2 && t3 && t4 && t5 && t6
//           )
//        {
//            savedTheSettings = true;
//            break;
//        }
//    }
//    return savedTheSettings;
//}
//
//
//
////==============================================================================================================================
//// TelescopeCustomRate9SettingsGet
////
////   Gets the custom rate 9 settings from MC
////==============================================================================================================================
//- (bool) TelescopeCustomRate9SettingsGet
//{
//    bool gotTheSettings = false;
//    bool t1 = false, t2 = false, t3 = false, t4 = false;
//    int attempts = 5;
//    while (attempts--) {
//        if(!t1) t1 = [self AuxMcCustomRate9Get :kAddrAzm :&mMotorControlSettings.customRate9Azm :&mMotorControlSettings.customRate9MaxValue];
//        if(!t2) t2 = [self AuxMcCustomRate9EnaGet :kAddrAzm :&mMotorControlSettings.customRate9AzmEnabled];
//        if(!t3) t3 = [self AuxMcCustomRate9Get :kAddrAlt :&mMotorControlSettings.customRate9Alt :nil];
//        if(!t4) t4 = [self AuxMcCustomRate9EnaGet :kAddrAlt :&mMotorControlSettings.customRate9AltEnabled];
//        if(
//           t1 && t2 && t3 && t4
//           )
//        {
//            gotTheSettings = true;
//            break;
//        }
//    }
//    return gotTheSettings;
//}
//
//
//
////==============================================================================================================================
//// TelescopeCustomRate9SettingsSet
////
////   Sends the custom rate 9 settings to MC
////==============================================================================================================================
//- (bool) TelescopeCustomRate9SettingsSet
//{
//    bool savedTheSettings = false;
//    bool t1 = false, t2 = false, t3 = false, t4 = false;
//    int attempts = 5;
//    while (attempts--) {
//        if(!t1) t1 = t1 | [self AuxMcCustomRate9Set :kAddrAzm :mMotorControlSettings.customRate9Azm];
//        if(!t2) t2 = t2 | [self AuxMcCustomRate9EnaSet :kAddrAzm :mMotorControlSettings.customRate9AzmEnabled];
//        if(!t3) t3 = t3 | [self AuxMcCustomRate9Set  :kAddrAlt :mMotorControlSettings.customRate9Alt];
//        if(!t4) t4 = t4 | [self AuxMcCustomRate9EnaSet :kAddrAlt :mMotorControlSettings.customRate9AltEnabled];
//        if(
//           t1 && t2 && t3 && t4
//           )
//        {
//            savedTheSettings = true;
//            break;
//        }
//    }
//    return savedTheSettings;
//}
//
//
//
////==============================================================================================================================
//// TelescopeGuideRatesGet
////
////   Gets the guide rate settings from MC
////==============================================================================================================================
//- (bool) TelescopeGuideRatesGet
//{
//    bool gotTheSettings = false;
//    bool t1 = false, t2 = false;
//    int attempts = 5;
//    while (attempts--) {
//        if(!t1) t1 = [self AuxMcGuideRateGet :kAddrAzm :&mMotorControlSettings.guideRateRa];
//        if(!t2) t2 = [self AuxMcGuideRateGet :kAddrAlt :&mMotorControlSettings.guideRateDec];
//        if(
//           t1 && t2
//           )
//        {
//            gotTheSettings = true;
//            break;
//        }
//    }
//    return gotTheSettings;
//}
//
//
//
////==============================================================================================================================
//// TelescopeGuideRatesSet
////
////   Sends the guide rate settings to MC
////==============================================================================================================================
//- (bool) TelescopeGuideRatesSet
//{
//    bool savedTheSettings = false;
//    bool t1 = false, t2 = false;
//    int attempts = 5;
//    while (attempts--) {
//        if(!t1) t1 = t1 | [self AuxMcGuideRateSet :kAddrAzm :mMotorControlSettings.guideRateRa];
//        if(!t2) t2 = t1 | [self AuxMcGuideRateSet :kAddrAlt :mMotorControlSettings.guideRateDec];
//        if(
//           t1 && t2
//           )
//        {
//            savedTheSettings = true;
//            break;
//        }
//    }
//    return savedTheSettings;
//}
//
//
//
////==============================================================================================================================
//// TelescopeRaLimitSettingsGet
////
////   Gets the ra limit settings from MC
////==============================================================================================================================
//- (bool) TelescopeRaLimitSettingsGet
//{
//    bool gotTheSettings = false;
//    bool t1 = false, t2 = false, t3 = false;
//    int attempts = 5;
//    while (attempts--) {
//        if(!t1) t1 = t1 | [self AuxMcRaLimitEnaGet :kAddrAzm :&mMotorControlSettings.raLimitsEnabled];
//        if(!t2) t2 = t2 | [self AuxMcRaLimitGet :kAddrAzm :true :&mMotorControlSettings.raLimitWest];
//        if(!t3) t3 = t3 | [self AuxMcRaLimitGet :kAddrAzm :false :&mMotorControlSettings.raLimitEast];
//        if(
//           t1 && t2 && t3
//           )
//        {
//            gotTheSettings = true;
//            break;
//        }
//        else {
//            gotTheSettings = false;
//        }
//    }
//    return gotTheSettings;
//}
//
//
//
////==============================================================================================================================
//// TelescopeRaLimitSettingsSet
////
////   Sends the ra limit settings to MC
////==============================================================================================================================
//- (bool) TelescopeRaLimitSettingsSet
//{
//    bool gotTheSettings = false;
//    bool t1 = false, t2 = false, t3 = false;
//    int attempts = 5;
//    while (attempts--) {
//        if(!t1) t1 = t1 | [self AuxMcRaLimitEnaSet :kAddrAzm :mMotorControlSettings.raLimitsEnabled];
//        if(!t2) t2 = t2 | [self AuxMcRaLimitSet :kAddrAzm :true :mMotorControlSettings.raLimitWest];
//        if(!t3) t3 = t3 | [self AuxMcRaLimitSet :kAddrAzm :false :mMotorControlSettings.raLimitEast];
//        if(
//           t1 && t2 && t3
//           )
//        {
//            gotTheSettings = true;
//            break;
//        }
//        else {
//            gotTheSettings = false;
//        }
//    }
//    return gotTheSettings;
//}
//
//
//
////==============================================================================================================================
//// handlerGetSettings
////
////   Used for getting the WiFly settings in command mode
////==============================================================================================================================
//- (void)handlerGetSettings
//{
//    DLog(@"**** handlerGetSettings ****");
//    mGetSettingsHandlerRunning = true;
//    mHaveSettingsWiFly = false;
//    mHaveSettingsTelescope = false;
//    @autoreleasepool {
//        mWiFlyGetSettingsState = kGetSettingsSendCmd;
//        mWiFlyParsedResponse = false;
//        while (
//               mWiFlyGetSettingsState != kGetSettingsDone &&
//               mWiFlyGetSettingsState != kGetSettingsError
//               )
//        {
//            DLog(@"  WiFlyGetSettings State:");
//            switch (mWiFlyGetSettingsState)
//            {
//                case kGetSettingsSendCmd:
//                {
//                    DLog(@"    Enter CMD mode");
//                    if([self WiFlyPacketMasterSend :(unsigned char *)"$$$" :kWiFlyCmd_Enter :@"CMD"])
//                    {
//                        mWiFlyGetSettingsState = kGetSettingsLoadConfig;
//                    }
//                }
//                    break;
//
//                case kGetSettingsLoadConfig:
//                {
//                    DLog(@"     Load config to WiFly module");
//                    if([self WiFlyPacketMasterSend :(unsigned char *)"load config\r\n" :kWiFlyCmd_LoadConfig :@"AOK"])
//                    {
//                        mWiFlyGetSettingsState = kGetSettingsSendGete;
//                    }
//                }
//
//                case kGetSettingsSendGete:
//                {
//
//                    DLog(@"    Get settings");
//                    if([self WiFlyPacketMasterSend :(unsigned char *)"get e\r\n" :kWiFlyCmd_GetAllSettings : @">"]) //@"<2.28.2>"])
//                    {
//                        mWiFlyGetSettingsState = kGetSettingsSendExit;
//                    }
//                }
//                    break;
//
//                case kGetSettingsSendExit:
//                {
//                    DLog(@"    Exit CMD mode");
//                    if([self WiFlyPacketMasterSend :(unsigned char *)"exit\r\n" :kWiFlyCmd_Exit :@"EXIT"])
//                    {
//                        mWiFlyGetSettingsState = kGetSettingsDone;
//                    }
//                }
//                    break;
//                default:
//                    break;
//            }
//        }
//        int attempts = 5;
//        while (attempts--)
//        {
//            if(
//               [self AuxDevGetVer :kAddrAzm :&mMotorControlSettings.AxisAzm] &&
//               [self AuxMcGetModel :kAddrAzm :&mMotorControlSettings.Model] &&
//               [self AuxDevGetVer :kAddrAlt :&mMotorControlSettings.AxisAlt] &&
//               [self AuxMcApproachDirGet :kAddrAzm :&telescopeSettings.approachAzmPos] &&
//               [self AuxMcApproachDirGet :kAddrAlt :&telescopeSettings.approachAltPos]
//               )
//            {
//                // NOTE: until we have a way of setting defaults, these much match for the UP / RIGHT centering to work
//                telescopeSettings.buttonDirAzmReversed = telescopeSettings.approachAzmPos;
//                telescopeSettings.buttonDirAltReversed = telescopeSettings.approachAltPos;
//                mHaveSettingsTelescope = true;
//                break;
//            }
//        }
//        mGetSettingsHandlerRunning = false;
//            //  NSLog(@"**** handlerGetSettings -- exit ****");
//    } //[pool release];
//}
//
//
//
//
////==============================================================================================================================
//// handlerWiFlySetSettings
////
////   Updating the WiFly settings in command mode
////==============================================================================================================================
//- (void)handlerWiFlySetSettings
//{
//    DLog(@"**** handlerWiFlySetSettings ****");
//    @autoreleasepool {
//
//        NSString *command;
//        mWiFlySetSettingsState = kSetSettingsSendCmd;
//        mWiFlyParsedResponse = false;
//        while (
//               mWiFlySetSettingsState != kSetSettingsDone &&
//               mWiFlySetSettingsState != kSetSettingsError
//               )
//        {
//            DLog(@"  SetSettings State:");
//            switch (mWiFlySetSettingsState)
//            {
//                case kSetSettingsSendCmd:
//                    DLog(@"    Enter CMD mode");
//                    [self WiFlyPacketMasterSend :(unsigned char *)"$$$" :kWiFlyCmd_Enter :@"CMD"];
//                    break;
//
//                case kSetSettingsSsid:
//                    command = [NSString stringWithFormat:@"set w s %s\r\n", [mWiFlySsidSetting UTF8String]];
//                    [self WiFlyPacketMasterSend :(unsigned char *)[command UTF8String] :kWiFlyCmd_SetSetting :@"AOK"];
//                    break;
//
//                case kSetSettingsAuth:
//                    command = [NSString stringWithFormat:@"set w a %d\r\n", (int)mWiFlyAuthSetting];
//                    [self WiFlyPacketMasterSend :(unsigned char *)[command UTF8String] :kWiFlyCmd_SetSetting :@"AOK"];
//                    break;
//
//                case kSetSettingsKeyPhrase:
//                    if (mWiFlyAuthSetting == kWep) {
//                        command = [NSString stringWithFormat:@"set w k %s\r\n", [mWiFlyKeyPhrase UTF8String]];
//                    }
//                    else
//                    {
//                        command = [NSString stringWithFormat :@"set w p %s\r\n", [mWiFlyKeyPhrase UTF8String]];
//                    }
//                    DLog(@"    command: %s", [command UTF8String]);
//                    [self WiFlyPacketMasterSend :(unsigned char *)[command UTF8String] :kWiFlyCmd_SetSetting :@"AOK"];
//                    break;
//
//                case kSetSettingsDhcp:
//                    command = [NSString stringWithFormat:@"set ip d %d\r\n", (int)mWiFlyDhcpSetting];
//                    DLog(@"    command: %s", [command UTF8String]);
//                    [self WiFlyPacketMasterSend :(unsigned char *)[command UTF8String] :kWiFlyCmd_SetSetting :@"AOK"];
//                    break;
//
//                case kSetSettingsIp:
//                    command = [NSString stringWithFormat:@"set ip a %s\r\n", [mWiFlyIpSetting UTF8String]];
//                    DLog(@"    command: %s", [command UTF8String]);
//                    [self WiFlyPacketMasterSend :(unsigned char *)[command UTF8String] :kWiFlyCmd_SetSetting :@"AOK"];
//                    break;
//
//                case kSetSettingsNm:
//                    command = [NSString stringWithFormat:@"set ip n %s\r\n", [mWiFlyNmSetting UTF8String]];
//                    DLog(@"    command: %s", [command UTF8String]);
//                    [self WiFlyPacketMasterSend :(unsigned char *)[command UTF8String] :kWiFlyCmd_SetSetting :@"AOK"];
//                    break;
//
//                case kSetSettingsGw:
//                    command = [NSString stringWithFormat:@"set ip g %s\r\n", [mWiFlyGwSetting UTF8String]];
//                    DLog(@"    command: %s", [command UTF8String]);
//                    [self WiFlyPacketMasterSend :(unsigned char *)[command UTF8String] :kWiFlyCmd_SetSetting :@"AOK"];
//                    break;
//
//                case kSetSettingsSave:
//                    command = @"save\r\n";
//                    DLog(@"    command: %s", [command UTF8String]);
//                    [self WiFlyPacketMasterSend :(unsigned char *)[command UTF8String] :kWiFlyCmd_SetSetting :@">"];
//                    break;
//
//                case kSetSettingsSendExit:
//                    DLog(@"    Exit CMD mode");
//                    [self WiFlyPacketMasterSend :(unsigned char *)"exit\r\n" :kWiFlyCmd_Exit :@"EXIT"];
//                    break;
//                default:
//                    break;
//            }
//            int sleepCount = 0;
//            while(mWiFlyParsedResponse == false && mWiFlySetSettingsState != kSetSettingsError)
//            {
//                [NSThread sleepForTimeInterval:0.5];
//                if(++sleepCount == 10) {
//                    //[self TcpDisconnect];
//                    mWiFlySetSettingsState = kSetSettingsError;
//                }
//            }
//            if(mWiFlyParsedResponse == true)
//            {
//                mWiFlyParsedResponse = false;
//                mWiFlySetSettingsState++;
//            }
//        }
//        DLog(@"**** handlerWiFlySetSettings -- exit ****");
//    } //[pool release];
//}
//
//
//
//
////==============================================================================================================================
//// WiFlyParseGetSettings --
////
////   Handles the WiFly "get settings" process
////==============================================================================================================================
//- (bool) WiFlyParseGetSettings :(unsigned char *)recvBuff :(tWiFlyCmd)Cmd
//{
//    NSString *Response;
//    mWiFlyParsedResponse = false;
//    DLog(@"  Parse Command Mode Response: ");
//    switch (Cmd) {
//        case kWiFlyCmd_Enter:
//        {
//            DLog(@"    -Enter");
//            Response = [NSString stringWithUTF8String :(const char *)recvBuff];
//            if (Response != nil)
//            {
//                mWiFlyParsedResponse = true;
//            }
//        }
//            break;
//
//        case kWiFlyCmd_LoadConfig:
//        {
//            DLog(@"    -Load Config");
//            Response = [NSString stringWithUTF8String :(const char *)recvBuff];
//            if (Response != nil)
//            {
//                mWiFlyParsedResponse = true;
//            }
//        }
//            break;
//
//        case kWiFlyCmd_GetAllSettings:
//        {
//            DLog(@"    -Get Everything");
//            Response = [NSString stringWithUTF8String: (const char *)(recvBuff+8)];
//
//            NSArray *Entries = [Response componentsSeparatedByString :@"\r\n"]; // Query - need to alloc array before use? Confirmed - no.
//
//
//            if (Entries.count < 5)
//            {
//                DLog(@"      !Parse Error");
//                break;
//            }
//            for (NSString *entry in Entries)
//            {
//                DLog(@"      Entry: %s",[entry UTF8String]);
//
//                NSArray *params = [entry componentsSeparatedByString :@"="];
//
//
//
//
//                if(params.count >= 2)
//                {
//                    NSString *entryId = [params objectAtIndex:0];
//                    NSString *entryValue = [params objectAtIndex:1];
//                    if([entryId isEqualToString :@"SSID"])
//                    {
//                        DLog(@"      FOUND %s=%s", [entryId UTF8String], [entryValue UTF8String]);
//                        mWiFlySsidSetting = entryValue;
//                    }
//                    else if([entryId isEqualToString :@"Auth"])
//                    {
//                        DLog(@"      FOUND %s=%s", [entryId UTF8String], [entryValue UTF8String]);
//                        if([entryValue isEqualToString :@"OPEN"])
//                        {
//                            mWiFlyAuthSetting = kOpen;
//                        }
//                        else if([entryValue isEqualToString :@"WEP"])
//                        {
//                            mWiFlyAuthSetting = kWep;
//                        }
//                        else
//                        {
//                            mWiFlyAuthSetting = kWpa;
//                        }
//                    }
//                    else if(
//                            [entryId isEqualToString :@"Passphrase"] ||
//                            [entryId isEqualToString :@"Key"]
//                            )
//                    {
//                        DLog(@"      FOUND %s=%s", [entryId UTF8String], [entryValue UTF8String]);
//                        mWiFlyKeyPhrase = [entryValue stringByReplacingOccurrencesOfString:@" " withString:@""];
//                    }
//                    else if([entryId isEqualToString :@"DHCP"])
//                    {
//                        DLog(@"      FOUND %s=%s", [entryId UTF8String], [entryValue UTF8String]);
//                        if([entryValue isEqualToString :@"ON"])
//                        {
//                            mWiFlyDhcpSetting = kDhcpEnabled;
//                        }
//                        else
//                        {
//                            mWiFlyDhcpSetting = kDhcpDisabled;
//                        }
//                    }
//                    else if([entryId isEqualToString :@"IP"])
//                    {
//                        NSArray *combinedIpandPort = [entryValue componentsSeparatedByString:@":"];
//
//
//                        NSString *entryIp = [combinedIpandPort objectAtIndex:0];
//                        DLog(@"      FOUND %s=%s", [entryId UTF8String], [entryIp UTF8String]);
//                        mWiFlyIpSetting = entryIp;
//                    }
//                    else if([entryId isEqualToString :@"NM"])
//                    {
//                        DLog(@"      FOUND %s=%s", [entryId UTF8String], [entryValue UTF8String]);
//                        mWiFlyNmSetting = entryValue;
//                    }
//                    else if([entryId isEqualToString :@"GW"])
//                    {
//                        DLog(@"      FOUND %s=%s", [entryId UTF8String], [entryValue UTF8String]);
//                        mWiFlyGwSetting = entryValue;
//                    }
//                }
//            }
//            mWiFlyParsedResponse = true;
//        }
//            break;
//
//        case kWiFlyCmd_SetSetting:                                             // case: save the settings config to WiFly
//        {
//            DLog(@"    Set");
//            Response = [NSString stringWithUTF8String: (const char *)recvBuff];
//            DLog(@"      %s", [Response UTF8String]);
//            mWiFlyParsedResponse = true;
//        }
//            break;
//
//        case kWiFlyCmd_SaveConfig:
//        {
//            DLog(@"     Save Config");
//            Response = [NSString stringWithUTF8String: (const char *)recvBuff];
//            DLog(@"       %s", [Response UTF8String]);
//            mWiFlyParsedResponse = true;
//        }
//
//        case kWiFlyCmd_Exit:                                                    // case: exit WiFly command mode
//        {
//            DLog(@"    Exit");
//            mWiFlyParsedResponse = true;
//            mHaveSettingsWiFly = true;
//        }
//            break;
//    }
//    return mWiFlyParsedResponse;
//}
//
//
//
////==============================================================================================================================
//// WiFlyUpdateSettings --
////
////   updates the WiFly connectivity settings
////==============================================================================================================================
//- (BOOL) WiFlyUpdateSettings :(NSString*)ssid :(tWiFlyAuth)security :(NSString*)keyphrase :(tWiFlyDhcp)dhcp :(NSString*)ip :(NSString*)nm :(NSString*)gw
//{
//    int sleepCount;
//
//    mWiFlySsidSetting = [NSString stringWithString:ssid];
//    mWiFlyAuthSetting = security;
//    mWiFlyKeyPhrase = [NSString stringWithString:keyphrase];
//    mWiFlyDhcpSetting = dhcp;
//    mWiFlyIpSetting = [NSString stringWithString:ip];
//    mWiFlyNmSetting = [NSString stringWithString:nm];
//    mWiFlyGwSetting = [NSString stringWithString:gw];
//    mWiFlySetSettingsState = kSetSettingsInit;
//    [NSThread detachNewThreadSelector:@selector(handlerWiFlySetSettings) toTarget:self withObject:nil];
//    sleepCount = 0;
//    while (
//           mWiFlySetSettingsState != kSetSettingsDone &&
//           mWiFlySetSettingsState != kSetSettingsError
//           ) {
//        [NSThread sleepForTimeInterval:0.50];
//        if (sleepCount++ > 20) {
//            mWiFlySetSettingsState = kSetSettingsError;
//            break;
//        }
//
//    }
//    return mWiFlySetSettingsState == kSetSettingsDone ? true : false;
//}
//
//
//
//
////==============================================================================================================================
//// WiFlyGetSettings --
////
////   Get settings that were parsed from WiFly module
////==============================================================================================================================
//- (void) WiFlyGetSettings
//{
//    ipAddress = [NSString stringWithString:mWiFlyIpSetting]; // Bug fix 1.5.10d. Precviously this was set to mWiFlySsidSetting.
//    networkName = [NSString stringWithString:mWiFlySsidSetting];
//    networkMask = [NSString stringWithString:mWiFlyNmSetting];
//    networkGateway = [NSString stringWithString:mWiFlyGwSetting];
//    password = [NSString stringWithString:mWiFlyKeyPhrase];
//    networkSecurity = mWiFlyAuthSetting;
//
//    //  NSLog(@"WiFly - network is %@ / %@",networkName,mWiFlySsidSetting);
//}
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//#pragma mark - Reachability
//
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//#pragma mark - Initialization / Life Cycle
////==============================================================================================================================
//// initWithSettings --
////
////   Initializes the WiFly interface with specified settings
//// this returned 'nil' on fail, but this meant some memory and observers were left dangling
//// so i changed it to return self, but to set a flag - error - which we could act on.
////==============================================================================================================================
//- (id) initWithSettings :(NSString *) Ip : (tWiFlyConnectMode) ConnectMode : (bool *)error
//{
//    DLog(@"**** WiFly initwithsettings ****");
//
//    self = [super init];
//    if (self)
//    {
//
//        *error = false;
//
//        // Initialize WiFly
//        if(ConnectMode == kDirectConnect)
//        {
//            mWiFlyIP = @"1.2.3.4";
//            DLog(@"  WiFly Direct Access - 1.2.3.4");
//        }
//        else
//        {
//            mWiFlyIP = Ip;
//            DLog(@"  WiFly IP Address : %s", [mWiFlyIP UTF8String]);
//        }
//
//        // Initialize WiFly Module settings
//        mWiFlyGetSettingsState = kGetSettingsInit;               // Initialize state for getting settings from WiFly Module
//
//        // Initialize Aux Interface objects
//        mAuxRecvPkts = [[NSMutableArray alloc] init]; // crash init with size?
//
//        // Initialize WiFi "Reachability"
//        //mWiFiReachability = [[Reachability reachabilityForLocalWiFi] retain];
//        //mWiFiReachability = [Reachability reachabilityForLocalWiFi] ; // when does this go to nil? // crash suspect //
//        //[mWiFiReachability startNotifier];
//        //[self updateInterfaceWithReachability: mWiFiReachability];
//        //[[NSNotificationCenter defaultCenter] addObserver: self selector: @selector(reachabilityChanged:) name: kReachabilityChangedNotification object: nil]; // remember to cancel this
//
//        // Initialize WiFly connectivity
//        mTcpMaintainConnection = false;
//        mTcpHandlerRunning = false;
//        mIsConnected = false;
//
//        mGetSettingsHandlerRunning = false;
//        mWiFlySsidSetting = @"";
//        mWiFlyAuthSetting = kOpen;
//        mWiFlyKeyPhrase = @"";
//        mWiFlyDhcpSetting = kDhcpDisabled;
//        mWiFlyIpSetting = @"";
//        mWiFlyNmSetting = @"";
//        mWiFlyGwSetting = @"";
//
//        mHaveSettingsWiFly = false;
//        mHaveSettingsTelescope = false;
//
//        //dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
//        //dispatch_async(queue, ^(void) { [self handlerGetSettings]; });
//
//        [NSThread detachNewThreadSelector:@selector(handlerGetSettings) toTarget:self withObject:nil];
//
//        //while (mTcpHandlerRunning == false) {                           // wait for the TCP handler to start -- we will us this
//        //
//        //    [NSThread sleepForTimeInterval: 0.01];                      //   signal as a sign that things are not working...
//        //}
//        //while (mGetSettingsHandlerRunning == true) {                    // get the WiFly settings to initiate the connection
//
//        //    [NSThread sleepForTimeInterval: 0.01];                      //
//        //    if(mTcpHandlerRunning == false)
//        //    {
//
//                //NSLog(@"WiFly init failed!");
//
//        //        *error = true;
//
//        //        [[NSNotificationCenter defaultCenter] removeObserver:self name:kReachabilityChangedNotification object:nil];// We need to remove the observer, but still takes tcpDisconnect 30 seconds to appear..
//
//        //       [mWiFiReachability stopNotifier];
//
//        //        mAuxRecvPkts = nil;
//
//
//        //        return self;                 //   if: the tcp handler exits, things failed, return
//        //    }
//        //}
//    }
//    return self;
//}
//
//
//
////==============================================================================================================================
//// init --
////
////   Initializes the WiFfly interface with default settings
////==============================================================================================================================
//- (id) init  : (bool *)error; {
//    return [self initWithSettings:@"" :kDirectConnect : error];
//}
//
//
//
////==============================================================================================================================
//// dealloc --
////==============================================================================================================================
//- (void)dealloc {
//    [self CeaseAndDesist];
//}
//
//
//
//-(void) CeaseAndDesist
//{
//        //NSLog(@"CeaseAndDesist 1");
//
//    //**djm not sure this is necessary, but seemed reasonable
//    //[[NSNotificationCenter defaultCenter] removeObserver :self];
//    //[mWiFiReachability stopNotifier];
//    //mWiFiReachability = nil;
//    //**djm
//
//    //[self TcpDisconnect];
//    mWiFlyIpSetting = nil;
//    mWiFlySsidSetting = nil;
//    mWiFlyKeyPhrase = nil;
//    mWiFlyNmSetting = nil;
//    mWiFlyGwSetting = nil;
//
//    mAuxRecvPkts = nil;
//
//    //NSLog(@"CeaseAndDesist 2");
//}
//
//@end
//#endif



WiFly::WiFly(
  CTelescope
    *pTelescope
)
{
  mTelescope = pTelescope;
}



WiFly::~WiFly()
{
}


