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

#include <string.h>
#include "utils/MathUtils.h"
#include "utils/log.h"
#include "settings/AdvancedSettings.h"

/* Minimum/maximum sample rate (Hz) */
#define HRTF_MINSRATE 2000
#define HRTF_MAXSRATE 384000

/* Minimum/maximum cut frequency (Hz) */
#define HRTF_MINFCUT 300
#define HRTF_MAXFCUT 2000

/* Minimum/maximum feed level (dB * 10 @ low frequencies) */
#define HRTF_MINFEED 10   /* 1 dB  */
#define HRTF_MAXFEED 150  /* 15 dB */

/* Minimum/maximum delays (uSec) */
#define HRTF_MINDELAY  90  /*  90 uS */
#define HRTF_MAXDELAY 620  /* 620 uS */

/* Minimum/maximum gains (scale) */
#define HRTF_MINGAIN 0.7
#define HRTF_MAXGAIN 1.2

/* Default sample rate (Hz) */
#define HRTF_DEFAULT_SRATE   44100

/* Lowpass filter */
#define lo_filter(in, out_1) \
  (m_hrtfd.a0_lo * in + m_hrtfd.b1_lo * out_1)

/* Highboost filter */
#define hi_filter(in, in_1, out_1) \
  (m_hrtfd.a0_hi * in + m_hrtfd.a1_hi * in_1 + m_hrtfd.b1_hi * out_1)

/* Define value for PI */
#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

struct hrtfModel
{
  const char*        name;
  const unsigned int cutFreq;
  const HRFT_TYPE    feedLvl;
};

static const hrtfModel hrtfModels[] =     { {"DEFAULT", 700,  3.5},
                                            {"CMOY",    700,  6.0},
                                            {"JMEIER",  650,  9.5},
                                            {"WIDE",   1200,  1.1},
                                            {"NARROW",  600,  1.1},
                                            {NULL, 0, 0.0} };


CAEDSPHeadphonesHRTF::CAEDSPHeadphonesHRTF() :
  m_returnBuffer (NULL),
  m_returnSamples(0   )
{
  memset(&m_hrtfd, 0, sizeof(m_hrtfd));
}

CAEDSPHeadphonesHRTF::~CAEDSPHeadphonesHRTF()
{
  DeInitialize();
}

void CAEDSPHeadphonesHRTF::DeInitialize()
{
  memset(&m_hrtfd, 0, sizeof(m_hrtfd));
}

bool CAEDSPHeadphonesHRTF::Initialize(const CAEChannelInfo& channels, const unsigned int sampleRate)
{
  if (channels.Count() != 2 || !channels.ContainsChannels(AE_CH_LAYOUT_2_0))
    return false;

  if (sampleRate < HRTF_MINSRATE || sampleRate > HRTF_MAXSRATE)
       m_sampleRate = HRTF_DEFAULT_SRATE;
  else m_sampleRate = sampleRate;
  m_channels = channels;

  HRFT_TYPE Fc_lo; /* Lowpass filter cut frequency (Hz) */
  HRFT_TYPE Fc_hi; /* Highboost filter cut frequency (Hz) */
  HRFT_TYPE G_lo;  /* Lowpass filter gain (multiplier) */
  HRFT_TYPE G_hi;  /* Highboost filter gain (multiplier) */
  HRFT_TYPE GB_lo; /* Lowpass filter gain (dB) */
  HRFT_TYPE GB_hi; /* Highboost filter gain (dB) (0 dB is high) */
  HRFT_TYPE level; /* Feeding level (dB) (level = GB_lo - GB_hi) */
  HRFT_TYPE x;

  /* Get advancedsettings.xml parameters if set */
  std::string dspHRTFModel   = g_advancedSettings.dspHRTFModel.ToUpper();
  int         dspHRTFCutFreq = g_advancedSettings.dspHRTFCutFreq;
  HRFT_TYPE   dspHRTFFeedLvl = g_advancedSettings.dspHRTFFeedLvl;
  HRFT_TYPE   dspHRTFGain    = g_advancedSettings.dspHRTFGain;


  /* set the defaults */
  Fc_lo = hrtfModels[0].cutFreq;
  level = hrtfModels[0].feedLvl;

  if (dspHRTFModel.empty())
  {
    /* No standard model - check for settings */
    if (dspHRTFCutFreq >= HRTF_MINFCUT  && dspHRTFCutFreq <= HRTF_MINFCUT  &&
        dspHRTFFeedLvl >= HRTF_MINFEED  && dspHRTFFeedLvl <= HRTF_MAXFEED)
    {
      Fc_lo = dspHRTFCutFreq;
      level = dspHRTFFeedLvl;
    }
    else
      CLog::Log(LOGERROR, "CAEDSPHeadphonesHRTF::Initialize - Invalid Settings selected for Headphones DSP");
  }
  else
  {
    /* Determine if we use a standard model */
    bool found = false;
    for(const hrtfModel *model = &hrtfModels[0]; model->name != NULL; ++model)
      if (dspHRTFModel == model->name)
      {
        Fc_lo = model->cutFreq;
        level = model->feedLvl;
        found = true;
        break;
      }

    if (!found)
      CLog::Log(LOGERROR, "CAEDSPHeadphonesHRTF::Initialize - Invalid Model selected for Headphones DSP");
  }

  m_hrtfd.srate = m_sampleRate;
  m_hrtfd.level = level;

  GB_lo = level * -5.0 / 6.0 - 3.0;
  GB_hi = level / 6.0 - 3.0;

  G_lo  = pow(10, GB_lo / 20.0);
  G_hi  = 1.0 - pow(10, GB_hi / 20.0);
  Fc_hi = Fc_lo * pow(2.0, (GB_lo - 20.0 * log10(G_hi )) / 12.0);

  x = exp(-2.0 * M_PI * Fc_lo / (HRFT_TYPE)m_hrtfd.srate);
  m_hrtfd.b1_lo = x;
  m_hrtfd.a0_lo = G_lo * (1.0 - x);

  x = exp(-2.0 * M_PI * Fc_hi / (HRFT_TYPE)m_hrtfd.srate);
  m_hrtfd.b1_hi = x;
  m_hrtfd.a0_hi = 1.0 - G_hi * (1.0 - x);
  m_hrtfd.a1_hi = -x;

  if (dspHRTFGain < HRTF_MINGAIN || dspHRTFGain > HRTF_MAXGAIN)
  {
    m_hrtfd.gain = 1.0 / ((1.0 - G_hi + G_lo) * 0.9);
  }
  else
  {
    m_hrtfd.gain = dspHRTFGain;
  }

  return true;
}

void CAEDSPHeadphonesHRTF::GetOutputFormat(CAEChannelInfo& channels, unsigned int& sampleRate)
{
  sampleRate = m_sampleRate;
  channels   = m_channels;
  return;
}

unsigned int CAEDSPHeadphonesHRTF::Process(float *data, unsigned int samples)
{
  m_returnBuffer  = data;
  m_returnSamples = samples;

  for(unsigned int i = 0; i < samples; i += 2, data += 2)
  {
#ifdef HRTF_DOUBLE
    HRFT_TYPE sample_d[2];

    sample_d[0] = (HRFT_TYPE)data[0];
    sample_d[1] = (HRFT_TYPE)data[1];

    /* Lowpass filter */
    m_hrtfd.lfs.lo[0] = lo_filter(sample_d[0], m_hrtfd.lfs.lo[0]);
    m_hrtfd.lfs.lo[1] = lo_filter(sample_d[1], m_hrtfd.lfs.lo[1]);

    /* Highboost filter */
    m_hrtfd.lfs.hi[0] =
      hi_filter(sample_d[0], m_hrtfd.lfs.asis[0], m_hrtfd.lfs.hi[0]);
    m_hrtfd.lfs.hi[1] =
      hi_filter(sample_d[1], m_hrtfd.lfs.asis[1], m_hrtfd.lfs.hi[1]);
    m_hrtfd.lfs.asis[0] = sample_d[0];
    m_hrtfd.lfs.asis[1] = sample_d[1];

    /* Crossfeed */
    sample_d[0] = m_hrtfd.lfs.hi[0] + m_hrtfd.lfs.lo[1];
    sample_d[1] = m_hrtfd.lfs.hi[1] + m_hrtfd.lfs.lo[0];

    /* Bass boost requires allpass attenuation */
    sample_d[0] *= m_hrtfd.gain;
    sample_d[1] *= m_hrtfd.gain;

    data[0] = (float)sample_d[0];
    data[1] = (float)sample_d[1];
#else
    /* Lowpass filter */
    m_hrtfd.lfs.lo[0] = lo_filter(data[0], m_hrtfd.lfs.lo[0]);
    m_hrtfd.lfs.lo[1] = lo_filter(data[1], m_hrtfd.lfs.lo[1]);

    /* Highboost filter */
    m_hrtfd.lfs.hi[0] =
      hi_filter(data[0], m_hrtfd.lfs.asis[0], m_hrtfd.lfs.hi[0]);
    m_hrtfd.lfs.hi[1] =
      hi_filter(data[1], m_hrtfd.lfs.asis[1], m_hrtfd.lfs.hi[1]);
    m_hrtfd.lfs.asis[0] = data[0];
    m_hrtfd.lfs.asis[1] = data[1];

    /* Crossfeed */
    data[0] = m_hrtfd.lfs.hi[0] + m_hrtfd.lfs.lo[1];
    data[1] = m_hrtfd.lfs.hi[1] + m_hrtfd.lfs.lo[0];

    /* Bass boost requires allpass attenuation */
    data[0] *= m_hrtfd.gain;
    data[1] *= m_hrtfd.gain;
#endif
  }

  return samples;
}

float *CAEDSPHeadphonesHRTF::GetOutput(unsigned int& samples)
{
  samples = m_returnSamples;
  return m_returnBuffer;
}

double CAEDSPHeadphonesHRTF::GetDelay()
{
  return 0.0;
}

void CAEDSPHeadphonesHRTF::OnSettingChange(std::string setting)
{
  /* TODO: re-init based on GUI setting */
}

