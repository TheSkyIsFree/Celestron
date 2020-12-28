//
//  AuxRecv.cpp
//  Created by Danyal Medley on 9/11/13.
//
#define AUXRECV_CPP
//#include "DmTypes.hpp"
#include "AuxRecv.hpp"


AuxRecv::AuxRecv(uint8_t *Packet)
{
  int i;

  mLen = Packet[0];                                             // map the packet data into member vars
  mSrcAddr = Packet[1];                                         //
  mDstAddr = Packet[2];                                         //
  mCmd = Packet[3];                                             //
  mDataLen = 0;                                                 //
  mCsum = 0;
	
  if (mLen > RECVBUFFSZ) {                                      // if: mLen or invalid packet can be bogus, prevent overrun
    wasValidated = false;                                       //
  }                                                             //
  else {                                                        // else: valid length
    for(i = 0; i <= mLen; i++) {                                //   for: calculate csum
      mCsum += buff[i] = Packet[i];                             //
    }                                                           //
    if(((Packet[i] + mCsum) & 0xff) == 0x00) {                  //   if: valid checksum
      wasValidated = true;                                      //     flag that packet was validated
      if(mLen > 3) {                                            //     if: packet returned data, calc length and set ptr
        mDataLen = mLen - 3;                                    //
        mDataPtr = (unsigned char*)&Packet[4];                  //
      }                                                         //
    }
    else {                                                      //   else: checksum invalid
      wasValidated = false;                                     //
    }                                                           //
  }                                                             //
}



AuxRecv::~AuxRecv()
{
}
