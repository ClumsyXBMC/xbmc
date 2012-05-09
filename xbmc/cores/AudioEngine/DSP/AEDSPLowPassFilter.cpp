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

#include "AEDSPLowPassFilter.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "utils/log.h"

CAEDSPLowPassFilter::CAEDSPLowPassFilter() :
  m_returnBuffer  (NULL ),
  m_returnSamples (0    ),
  m_iir             (&m_iirF)
{
  memset(&m_iirF, 0, sizeof(m_iirF));
}

CAEDSPLowPassFilter::~CAEDSPLowPassFilter()
{
  DeInitialize();
}

void CAEDSPLowPassFilter::DeInitialize()
{
  memset(&m_iirF, 0, sizeof(m_iirF));
  if (m_iir->history) free(m_iir->history);
  if (m_iir->coef)    free(m_iir->coef);
}

bool CAEDSPLowPassFilter::Initialize(const CAEChannelInfo& channels, const unsigned int sampleRate)
{
  LFP_TYPE   *coef;
  LFP_TYPE   fs, fc;    /* Sampling frequency, cutoff frequency */
  LFP_TYPE   Q;         /* Resonance > 1.0 < 1000 */
  unsigned nInd;
  LFP_TYPE   a0, a1, a2, b0, b1, b2;
  LFP_TYPE   k;         /* overall gain factor */

  m_iir->length = FILTER_SECTIONS;   /* Number of filter sections */

  /* Store number of channels for Process() */
  m_channelCount = channels.Count();

  /* Check that we have an LFE channel and get its number */
  m_lfeChannel = 0;

  for (unsigned int c = 0; c < m_channelCount; c++)
  {
    if (channels[c] == AE_CH_LFE)
    {
      m_lfeChannel = c;
      break;
    }
  }

  if (!m_lfeChannel) return false;  /* no LFE channel present */

  /* Set up history and coefficient storage memory */
  memset(m_iir->history, 0, 2 * m_iir->length * sizeof(LFP_TYPE));
  memset(m_iir->coef, 0, (4 * m_iir->length + 1) * sizeof(LFP_TYPE));

  /* Setup filter s-domain coefficients */
  /* Section 1 */
  coefs[0].a0 = 1.0;
  coefs[0].a1 = 0;
  coefs[0].a2 = 0;
  coefs[0].b0 = 1.0;
  coefs[0].b1 = 0.765367;
  coefs[0].b2 = 1.0;

  /* Section 2 */
  coefs[1].a0 = 1.0;
  coefs[1].a1 = 0;
  coefs[1].a2 = 0;
  coefs[1].b0 = 1.0;
  coefs[1].b1 = 1.847759;
  coefs[1].b2 = 1.0;

  k = 1.0;              /* Set overall filter gain */
  coef = m_iir->coef + 1; /* Skip k, or gain */

  Q = 1;                /* Resonance */
  fc = 140;             /* Filter cutoff (Hz) */
  fs = sampleRate;      /* Sampling frequency (Hz) */

  /* Compute z-domain coefficients for each biquad section
     for new Cutoff Frequency and Resonance */
  for (nInd = 0; nInd < m_iir->length; nInd++)
  {
    a0 = coefs[nInd].a0;
    a1 = coefs[nInd].a1;
    a2 = coefs[nInd].a2;

    b0 = coefs[nInd].b0;
    b1 = coefs[nInd].b1 / Q;  /* Divide by resonance or Q */
    b2 = coefs[nInd].b2;

    szxform(&a0, &a1, &a2, &b0, &b1, &b2, fc, fs, &k, coef);

    coef += 4;  /* Point to next filter section */
  }

  /* Update overall filter gain in coef array */
  m_iir->coef[0] = k;

  return true;
}

void CAEDSPLowPassFilter::GetOutputFormat(CAEChannelInfo& channels, unsigned int& sampleRate)
{
  return;
}

unsigned int CAEDSPLowPassFilter::Process(float *data, unsigned int samples)
{
  float* pSampleBuffer = data;  /* incoming buffer of samples */
  float* pOffset;               /* storage for current sample address */

  LFP_TYPE *hist1_ptr, *hist2_ptr, *coef_ptr;
  LFP_TYPE new_hist, history1, history2;
  LFP_TYPE output;

  unsigned int frames = samples / m_channelCount;

  for (unsigned int c = 0; c < frames; c++)
  {
    coef_ptr = m_iir->coef;  /* pointer to cycle coefficients */
    
    /* calculate where our next LFE sample is */
    pOffset = pSampleBuffer + (c * m_channelCount + m_lfeChannel);

    /* 1st number of coefficients array is overall input scale factor or filter gain */
    output = (LFP_TYPE)pOffset[0] * (*coef_ptr++);
    
    hist1_ptr = m_iir->history;  /* first history */
    hist2_ptr = hist1_ptr + 1; /* next history */

    for (unsigned int i = 0 ; i < m_iir->length; i++)
    {
      history1 = *hist1_ptr;   /* history values */
      history2 = *hist2_ptr;

      output =   output - history1 * (*coef_ptr++);
      new_hist = output - history2 * (*coef_ptr++);  /* poles */

      output = new_hist + history1 * (*coef_ptr++);
      output = output   + history2 * (*coef_ptr++);  /* zeros */

      *hist2_ptr++ = *hist1_ptr;
      *hist1_ptr++ = new_hist;
      hist1_ptr++;
      hist2_ptr++;
    }

    pOffset[0] = (float)output;
  }

  return samples;
}

float *CAEDSPLowPassFilter::GetOutput(unsigned int& samples)
{
  samples = m_returnSamples;

  return m_returnBuffer;
}

double CAEDSPLowPassFilter::GetDelay()
{
  return 0.0;
}

void CAEDSPLowPassFilter::OnSettingChange(std::string setting)
{
  /* TODO: re-init based on GUI setting */
}

void CAEDSPLowPassFilter::szxform(LFP_TYPE *a0, LFP_TYPE *a1, LFP_TYPE *a2, /* numerator coefficients */
                                  LFP_TYPE *b0, LFP_TYPE *b1, LFP_TYPE *b2, /* denominator coefficients */
                                  LFP_TYPE fc,         /* Filter cutoff frequency */
                                  LFP_TYPE fs,         /* sampling rate */
                                  LFP_TYPE *k,         /* overall gain factor */
                                  LFP_TYPE *coef)      /* pointer to 4 iir coefficients */
{
  /* Calculate a1 and a2 and overwrite the original values */
  prewarp(a0, a1, a2, fc, fs);
  prewarp(b0, b1, b2, fc, fs);
  bilinear(*a0, *a1, *a2, *b0, *b1, *b2, k, fs, coef);
}

void CAEDSPLowPassFilter::bilinear(LFP_TYPE a0, LFP_TYPE a1, LFP_TYPE a2,    /* numerator coefficients */
                                   LFP_TYPE b0, LFP_TYPE b1, LFP_TYPE b2,    /* denominator coefficients */
                                   LFP_TYPE *k,        /* overall gain factor */
                                   LFP_TYPE fs,        /* sampling rate */
                                   LFP_TYPE *coef)     /* pointer to 4 iir coefficients */
{
  LFP_TYPE ad, bd;

  ad = 4. * a2 * fs * fs + 2. * a1 * fs + a0;
  bd = 4. * b2 * fs * fs + 2. * b1* fs + b0;

  *k *= ad/bd;

  *coef++ = (2. * b0 - 8. * b2 * fs * fs) / bd;
  *coef++ = (4. * b2 * fs * fs - 2. * b1 * fs + b0) / bd;
  *coef++ = (2. * a0 - 8. * a2 * fs * fs) / ad;
  *coef   = (4. * a2 * fs * fs - 2. * a1 * fs + a0) / ad;
}

void CAEDSPLowPassFilter::prewarp(LFP_TYPE *a0, LFP_TYPE *a1, LFP_TYPE *a2, LFP_TYPE fc, LFP_TYPE fs)
{
    LFP_TYPE wp, pi;

    pi = 4.0 * atan(1.0);
    wp = 2.0 * fs * tan(pi * fc / fs);

    *a2 = (*a2) / (wp * wp);
    *a1 = (*a1) / wp;
}
