#pragma once
/*
 *      Copyright (C) 2010-2012 Team XBMC
 *      http://xbmc.org
 *
 *  This Program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2, or (at your option)
 *  any later version.
 *
 *  This Program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with XBMC; see the file COPYING.  If not, write to
 *  the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 *  http://www.gnu.org/copyleft/gpl.html
 *
 *  Calculations derived from the BS2B project by Boris Mikhaylov
 *  http://bs2b.sourceforge.net/
 *
 *  Recoded and ported to XBMC AudioEngine by DDDamian May 2012
 */

#include "Interfaces/AEDSP.h"
//#include <samplerate.h>

/**
 * CAEDSPHeadphones - An implmentation of HRTF (Head-Related Transfer Function)
 *                    or crossfeed circuit for binaural headphone listening
 */
class CAEDSPHeadphonesHRTF
{
public:
  CAEDSPHeadphonesHRTF();
  virtual ~CAEDSPHeadphonesHRTF();

  /** Must call Initialize first to set DSP parameters */
  virtual bool Initialize(const CAEChannelInfo& channels, const unsigned int sampleRate);

  virtual void DeInitialize();

  /** Channels always two - downmix first, sample rate can be defined in call
   *    if not defined default of 44100hz returned */
  virtual void GetOutputFormat(CAEChannelInfo& channels, unsigned int& sampleRate);

  virtual unsigned int Process(float *data, unsigned int samples);

  virtual float *GetOutput(unsigned int& samples);

  virtual double GetDelay();

private:

  typedef struct
  {
    double level;                  /* Crossfeed level */
    uint32_t srate;                /* Sample rate (Hz) */
    double a0_lo, b1_lo;           /* Lowpass IIR filter coefficients */
    double a0_hi, a1_hi, b1_hi;    /* Highboost IIR filter coefficients */
    double gain;                   /* Global gain against overloading */
    /* Buffer of last filtered sample: [0] 1st channel, [1] 2nd channel */
    struct {double asis[ 2 ], lo[ 2 ], hi[ 2 ]; } lfs;
  } t_hrtfd;

  typedef t_hrtfd* t_hrtfdp;
  t_hrtfdp         hrtfdp;

  float*           pReturnBuffer;  /* Store float buffer data for GetOutput() */
  unsigned int     iReturnSamples; /* Store number of samples for GetOutput() */

  //void         DSPCrossFeed(float *sample, int n, int m_SampleRate);
  //void         DSPCrossFeedInit(int sampleRate);
};

