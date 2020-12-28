//
//  AuxSend.hpp
//  Created by Danyal Medley on 9/11/13.
//
#ifndef AUXSEND_HPP
#define AUXSEND_HPP
#include "DmTypes.hpp"


class AuxSend
{
public:
  uint8_t
    mLen,
    mSrcAddr,
    mDstAddr,
    mCmd,
    mDataLen,
    mCsum,
    buff[40];
  bool
    mValid;

  AuxSend(
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
  );

  AuxSend(
    uint8_t
      *data
  );

  virtual ~AuxSend(void);
};

#endif

