//
//  AuxSend.cpp
//  Created by Danyal Medley on 9/11/13.
//
#define AUXSEND_CPP
//#include "DmTypes.hpp"
#include "AuxSend.hpp"


AuxSend::AuxSend(
  uint8_t
    cmd,
  uint8_t
    len,
  uint8_t
    *data,
  uint8_t
    srcaddr,
  uint8_t
    dstaddr
)
{
  int i;
  
  mLen = len + 3;
  mSrcAddr = srcaddr;
  mDstAddr = dstaddr;
  mCmd = cmd;
  mDataLen = len;
  
  buff[0] = 0x3b;
  buff[1] = mLen;
  buff[2] = mSrcAddr;
  buff[3] = mDstAddr;
  buff[4] = mCmd;
  if(data == 0x0) len = 0;
  for(i = 5; i < 5 + mDataLen; i++) {
    buff[i] = *data++;
  }
  mCsum = 0;
  for(i = 1; i < 5 + mDataLen; i++) {
    mCsum += buff[i];
  }
  buff[i] = 0 - mCsum;
  mValid = true;
}



AuxSend::AuxSend(
  uint8_t
    *string
)
{
  uint8_t *buffPtr = string;
  mLen = 0;
  while(1)
  {
    *buffPtr = *string;
    mLen++;
    if(*buffPtr == '\x0')
    {
      mValid = true;
      break;
    }
    else if(mLen >= 40) {
      mValid = false;
      break;
    }
    buffPtr++;
    string++;
  }
}



AuxSend::~AuxSend(void)
{
}
