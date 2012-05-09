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
 *  Added to XBMC AudioEngine by DDDamian May 2012
 */

#include "Interfaces/AEDSP.h"

#define FILTER_SECTIONS 2

/* dont use double precision on embedded devices */
#if defined(__arm__) || defined(TARGET_DARWIN_IOS)
  #define LFP_TYPE float
#else
  #define LFP_DOUBLE
  #define LFP_TYPE double
#endif

/**
 * CAEDSPLowPassFilter - A fourth-order (24db/octave) Butterworth low-pass resonant filter
 */
class CAEDSPLowPassFilter
{
public:
  CAEDSPLowPassFilter();
  virtual ~CAEDSPLowPassFilter();

  /** Must call Initialize first to set DSP parameters */
  virtual bool Initialize(const CAEChannelInfo& channels, const unsigned int sampleRate);

  virtual void DeInitialize();

  virtual void GetOutputFormat(CAEChannelInfo& channels, unsigned int& sampleRate);

  virtual unsigned int Process(float *data, unsigned int samples);

  virtual float *GetOutput(unsigned int& samples);

  virtual double GetDelay();

  virtual void OnSettingChange(std::string setting);

private:

  typedef struct
  {
    unsigned int length;     /* size of filter */
    LFP_TYPE history[2 * FILTER_SECTIONS *     sizeof(LFP_TYPE)];   /* pointer to history in filter */
    LFP_TYPE coef   [4 * FILTER_SECTIONS + 1 * sizeof(LFP_TYPE)];   /* pointer to coefficients of filter */
  } FILTER;

  typedef struct
  {
    LFP_TYPE a0, a1, a2;     /* numerator coefficients */
    LFP_TYPE b0, b1, b2;     /* denominator coefficients */
  } BIQUAD;

  FILTER  m_iirF;            /* the history and its coefficients */
  FILTER* m_iir;             /* pointer to above */

  BIQUAD  coefs[FILTER_SECTIONS];  /* filter coefficients */

  void  CAEDSPLowPassFilter::szxform   (LFP_TYPE *a0, LFP_TYPE *a1, LFP_TYPE *a2,   /* numerator coefficients */
                                        LFP_TYPE *b0, LFP_TYPE *b1, LFP_TYPE *b2,   /* denominator coefficients */
                                        LFP_TYPE fc,           /* Filter cutoff frequency */
                                        LFP_TYPE fs,           /* sampling rate */
                                        LFP_TYPE *k,           /* overall gain factor */
                                        LFP_TYPE *coef);

  void CAEDSPLowPassFilter::bilinear   (LFP_TYPE a0, LFP_TYPE a1, LFP_TYPE a2,    /* numerator coefficients */
                                        LFP_TYPE b0, LFP_TYPE b1, LFP_TYPE b2,    /* denominator coefficients */
                                        LFP_TYPE *k,           /* overall gain factor */
                                        LFP_TYPE fs,           /* sampling rate */
                                        LFP_TYPE *coef);       /* pointer to 4 iir coefficients */

  void  CAEDSPLowPassFilter::prewarp   (LFP_TYPE *a0, LFP_TYPE *a1, LFP_TYPE *a2, LFP_TYPE fc, LFP_TYPE fs);

  float*           m_returnBuffer;  /* store float buffer data for GetOutput() */
  unsigned int     m_returnSamples; /* store number of samples for GetOutput() */
  unsigned int     m_channelCount;  /* number of channels passed by engine */
  unsigned int     m_lfeChannel;    /* channel number of the lfe channel */

};

