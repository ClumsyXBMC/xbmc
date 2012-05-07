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

#include "AEDSPHeadphonesHRTF.h"

#include "utils/MathUtils.h"
#include "utils/log.h"
#include "settings/AdvancedSettings.h"

/* Minimum/maximum sample rate (Hz) */
#define hrtf_MINSRATE 2000
#define hrtf_MAXSRATE 384000

/* Minimum/maximum cut frequency (Hz) */
#define hrtf_MINFCUT 300
#define hrtf_MAXFCUT 2000

/* Minimum/maximum feed level (dB * 10 @ low frequencies) */
#define hrtf_MINFEED 10   /* 1 dB  */
#define hrtf_MAXFEED 150  /* 15 dB */

/* Minimum/maximum delays (uSec) */
#define hrtf_MINDELAY  90  /*  90 uS */
#define hrtf_MAXDELAY 620  /* 620 uS */

/* Minimum/maximum gains (scale) */
#define hrtf_MINGAIN 0.7
#define hrtf_MAXGAIN 1.2

/* Default sample rate (Hz) */
#define hrtf_DEFAULT_SRATE   44100

/* Lowpass filter */
#define lo_filter(in, out_1) \
  (hrtfdp->a0_lo * in + hrtfdp->b1_lo * out_1)

/* Highboost filter */
#define hi_filter(in, in_1, out_1) \
  (hrtfdp->a0_hi * in + hrtfdp->a1_hi * in_1 + hrtfdp->b1_hi * out_1)

/* Define value for PI */
#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

struct hrtfModel
{
  std::string szModel;
  uint32_t    uiCutFreq;
  double      dFeedLvl;
};

static const hrtfModel hrtfModels[] =     { {"DEFAULT", 700,  3.5},
                                            {"CMOY",    700,  6.0},
                                            {"JMEIER",  650,  9.5},
                                            {"WIDE",   1200,  1.1},
                                            {"NARROW",  600,  1.1} };


CAEDSPHeadphonesHRTF::CAEDSPHeadphonesHRTF()
{
  hrtfdp = NULL;
  if((hrtfdp = (t_hrtfdp)malloc(sizeof(t_hrtfd))) != NULL)
  {
    memset(hrtfdp, 0, sizeof(t_hrtfd));
  }
}

CAEDSPHeadphonesHRTF::~CAEDSPHeadphonesHRTF()
{
  DeInitialize();
}

void CAEDSPHeadphonesHRTF::DeInitialize()
{
  if (hrtfdp)
  {
    free(hrtfdp);
    hrtfdp = NULL;
  }
}

bool CAEDSPHeadphonesHRTF::Initialize(const CAEChannelInfo& channels, const unsigned int sampleRate)
{
  if (channels.Count() != 2 || sampleRate < hrtf_MINSRATE || sampleRate > hrtf_MAXSRATE)
    return false;

  double Fc_lo; /* Lowpass filter cut frequency (Hz) */
  double Fc_hi; /* Highboost filter cut frequency (Hz) */
  double G_lo;  /* Lowpass filter gain (multiplier) */
  double G_hi;  /* Highboost filter gain (multiplier) */
  double GB_lo; /* Lowpass filter gain (dB) */
  double GB_hi; /* Highboost filter gain (dB) (0 dB is high) */
  double level; /* Feeding level (dB) (level = GB_lo - GB_hi) */
  double x;

  /* Get advancedsettings.xml parameters if set */
  std::string dspHRTFModel     = g_advancedSettings.dspHRTFModel.ToUpper().c_str();
  int         dspHRTFCutFreq   = g_advancedSettings.dspHRTFCutFreq;
  double      dspHRTFFeedLvl   = g_advancedSettings.dspHRTFFeedLvl;
  double      dspHRTFGain      = g_advancedSettings.dspHRTFGain;

  bool extSettings = false;
  bool extModel    = false;

  /* Determine if we use a standard model */
  if (dspHRTFModel != "")
  {
    for (int j = 0; j < sizeof(hrtfModels)/sizeof(hrtfModel); j++)
    {
      if (dspHRTFModel == hrtfModels[j].szModel)
      {
        Fc_lo = hrtfModels[j].uiCutFreq;
        level = hrtfModels[j].dFeedLvl;
        extModel = true;
        break;
      }
    }
    if (!extModel)
      CLog::Log(LOGERROR, __FUNCTION__": Invalid Model selected for Headphones DSP");
  }
  /* No standard model - check for settings */
  else
  {
    if (dspHRTFCutFreq >= hrtf_MINFCUT  && dspHRTFCutFreq <= hrtf_MINFCUT  &&
        dspHRTFFeedLvl >= hrtf_MINFEED  && dspHRTFFeedLvl <= hrtf_MAXFEED)
    {
      Fc_lo = dspHRTFCutFreq;
      level = dspHRTFFeedLvl;
      extSettings = true;
    }
    else
      CLog::Log(LOGERROR, __FUNCTION__": Invalid Settings selected for Headphones DSP");
  }

  if (!extSettings && !extModel)
  {
    Fc_lo = 700.0;
    level = 4.5;
  }

  hrtfdp->srate = sampleRate;
  hrtfdp->level = level;

  GB_lo = level * -5.0 / 6.0 - 3.0;
  GB_hi = level / 6.0 - 3.0;

  G_lo  = pow(10, GB_lo / 20.0);
  G_hi  = 1.0 - pow(10, GB_hi / 20.0);
  Fc_hi = Fc_lo * pow(2.0, (GB_lo - 20.0 * log10(G_hi )) / 12.0);

  x = exp(-2.0 * M_PI * Fc_lo / (double)hrtfdp->srate);
  hrtfdp->b1_lo = x;
  hrtfdp->a0_lo = G_lo * (1.0 - x);

  x = exp(-2.0 * M_PI * Fc_hi / (double)hrtfdp->srate);
  hrtfdp->b1_hi = x;
  hrtfdp->a0_hi = 1.0 - G_hi * (1.0 - x);
  hrtfdp->a1_hi = -x;

  if (dspHRTFGain < hrtf_MINGAIN || dspHRTFGain > hrtf_MAXGAIN)
  {
    hrtfdp->gain = 1.0 / ((1.0 - G_hi + G_lo) * 0.9);
  }
  else
  {
    hrtfdp->gain = dspHRTFGain;
  }

  return true;
}

void CAEDSPHeadphonesHRTF::GetOutputFormat(CAEChannelInfo& channels, unsigned int& sampleRate)
{
  if((sampleRate > hrtf_MAXSRATE) || (sampleRate < hrtf_MINSRATE || sampleRate == NULL))
  {
    hrtfdp->srate = hrtf_DEFAULT_SRATE;
    sampleRate    = hrtf_DEFAULT_SRATE;
  }
  else
  {
    hrtfdp->srate = sampleRate;
  }

  channels = AE_CH_LAYOUT_2_0;

  return;
}

unsigned int CAEDSPHeadphonesHRTF::Process(float *data, unsigned int samples)
{
  double       sample_d      [2];
  float*       pSampleBuf  = data;
  unsigned int count       = samples/2;
  
  pReturnBuffer            = pSampleBuf;
  iReturnSamples           = samples;

  if (count > 0)
  {
    while(count--)
    {
      sample_d[0] = (double)pSampleBuf[0];
      sample_d[1] = (double)pSampleBuf[1];

      /* Lowpass filter */
      hrtfdp->lfs.lo[0] = lo_filter(sample_d[0], hrtfdp->lfs.lo[0]);
      hrtfdp->lfs.lo[1] = lo_filter(sample_d[1], hrtfdp->lfs.lo[1]);

      /* Highboost filter */
      hrtfdp->lfs.hi[0] =
        hi_filter(sample_d[0], hrtfdp->lfs.asis[0], hrtfdp->lfs.hi[0]);
      hrtfdp->lfs.hi[1] =
        hi_filter(sample_d[1], hrtfdp->lfs.asis[1], hrtfdp->lfs.hi[1]);
      hrtfdp->lfs.asis[0] = sample_d[0];
      hrtfdp->lfs.asis[1] = sample_d[1];

      /* Crossfeed */
      sample_d[0] = hrtfdp->lfs.hi[0] + hrtfdp->lfs.lo[1];
      sample_d[1] = hrtfdp->lfs.hi[1] + hrtfdp->lfs.lo[0];

      /* Bass boost requires allpass attenuation */
      sample_d[0] *= hrtfdp->gain;
      sample_d[1] *= hrtfdp->gain;

      pSampleBuf[0] = (float)sample_d[0];
      pSampleBuf[1] = (float)sample_d[1];

      pSampleBuf += 2;
    }
  }

  return samples;
}

float *CAEDSPHeadphonesHRTF::GetOutput(unsigned int& samples)
{
  samples = iReturnSamples;

  return pReturnBuffer;
}

double CAEDSPHeadphonesHRTF::GetDelay()
{
  return 0.0;
}
