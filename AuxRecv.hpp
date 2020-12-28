//
//  AuxRecv.hpp
//  Created by Danyal Medley on 9/11/2013
//
#ifndef AUXRECV_HPP
#define AUXRECV_HPP
#include "DmTypes.hpp"

#define RECVBUFFSZ                40                            //

class AuxRecv
{
public:
  uint8_t                                                         //
    mLen,                                                       //
    mSrcAddr,                                                   //
    mDstAddr,                                                   //
    mCmd,                                                       //
    mDataLen,                                                   //
    mCsum,                                                      //
    buff[RECVBUFFSZ],                                           //
    *mDataPtr;                                                  //
  bool                                                          //
    wasValidated;                                               //

  AuxRecv(                                                      //
    uint8_t                                                       //
      *Packet                                                   //
  );                                                            //

  virtual ~AuxRecv(void);                                       //
};

#define NULLPTR_AUXRECV         ((AuxRecv *)0x0)

#endif
