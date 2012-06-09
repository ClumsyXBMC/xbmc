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
 *
 *  Added to XBMC AudioEngine by DDDamian June 2012
 */


#include "AEDSPDRCompressor.h"
#include "settings/AdvancedSettings.h"
#include "settings/GUIsettings.h"
#include "settings/Settings.h"


///////////////////////////////////////////////////////////////////////////


CAEDSPDRCompressor::CAEDSPDRCompressor() :
  m_returnBuffer (NULL),
  m_returnSamples(0   ),
  m_channelFL    (0   ),
  m_channelFR    (0   ),
  m_channelCount (0   )
{
}


CAEDSPDRCompressor::~CAEDSPDRCompressor()
{
  DeInitialize();
}

void CAEDSPDRCompressor::DeInitialize()
{
}

///////////////////////////////////////////////////////////////////////////


bool CAEDSPDRCompressor::Initialize(const CAEChannelInfo& channels, const unsigned int sampleRate)
{
  m_sampleRate   = sampleRate;
  m_channels     = channels;
  m_channelCount = channels.Count();
  m_rampCoeff    = SinglePoleCoeff(m_sampleRate, 0.05);

  if (!m_channelCount)
    return false;

  SetParameterDefaults();

  /* Get our FL & FR channel offsets */
  /* FR = FL if mono                 */
  m_channelFL = m_channelFR = 0;

  for (unsigned int c = 0; c < m_channelCount; c++)
  {
    if (m_channels[c] == AE_CH_FL)
      m_channelFL = c;
    if (m_channels[c] == AE_CH_FR)
      m_channelFR = c;
  }

  if (m_channelCount == 1)
    m_channelFR = m_channelFL;  /* must be mono */

  return true;
}

void CAEDSPDRCompressor::OnSettingChange(std::string setting)
{
  Initialize(m_channels, m_sampleRate);
}

void CAEDSPDRCompressor::SetParameterDefaults()
{
  m_rampCoeff = SinglePoleCoeff(m_sampleRate, 0.05);
  m_cvRelease = 0.0;
  m_cvAttack = 0.0;
  m_cvSmooth = 0.0;

  m_logThreshold[0] = 0.0;
  m_logThreshold[1] = 0.0;
  m_logKneeWidth = 0.0;
  m_logGain[0] = 0.0;
  m_logGain[1] = 0.0;
  m_slope = 0.0;

  m_crestPeak = MINVAL;
  m_crestRms = MINVAL;

  m_attackTime = 0.0;
  m_attackCoeff = 0.0;

  m_releaseTime = 0.0;
  m_releaseCoeff = 0.0;

  m_metaKneeMult = 0.0;
  m_metaMaxAttackTime = 0.0;
  m_metaMaxReleaseTime = 0.0;
  m_metaCrestTime = 0.0;
  m_metaCrestCoeff = 0.0;
  m_metaAdaptTime = 0.0;
  m_metaAdaptCoeff = 0.0;

  m_logThreshold[kTarget] = 0.0f;
  m_logGain[kTarget]      = 0.0f;
  m_attackTime            = LOG2PARAM(0.005, MINATTACKTIME, MAXATTACKTIME);
  m_releaseTime           = LOG2PARAM(0.05, MINRELEASETIME, MAXRELEASETIME);

  m_bAutoKnee      = true;
  m_bAutoGain      = true;
  m_bAutoAttack    = true;
  m_bAutoRelease   = true;
  m_bNoClipping    = true;

  /* Get advancedsettings.xml parameters if set */
  m_bAutoGain    = g_advancedSettings.dspDRCAutoGain;
  m_bAutoAttack  = g_advancedSettings.dspDRCAutoAttack;
  m_bAutoRelease = g_advancedSettings.dspDRCAutoRelease;
  m_bAutoKnee    = g_advancedSettings.dspDRCAutoKnee;
  m_attackTime   = 0.005; //PARAM2LOG(g_advancedSettings.dspDRCAttackTime, MINATTACKTIME, MAXATTACKTIME);
  m_releaseTime  = 0.200; //PARAM2LOG(g_advancedSettings.dspDRCReleaseTime, MINRELEASETIME, MAXRELEASETIME);
  m_slope        = -0.7; //-(1.0 - (1.0 / g_advancedSettings.dspDRCRatio));
  m_logGain[kTarget]      = 0.0; //20.0 * log(abs(g_advancedSettings.dspDRCGain) / MAXGAINDB);
  m_logThreshold[kTarget] = -20.00; //log(pow(10.0, g_advancedSettings.dspDRCThreshold * -MAXGAINDB / 20.0));
  m_logKneeWidth          = -log(pow(10.0, g_advancedSettings.dspDRCKneeWidth * MAXGAINDB / 20.0));

  m_metaKneeMult       = 2.0; //LIN2PARAM(2.0, MINKNEEMULT, MAXKNEEMULT);
  m_metaMaxAttackTime  = 0.8; //LOG2PARAM(0.08, MINATTACKTIME, MAXATTACKTIME);
  m_metaMaxReleaseTime = 1.0; //LOG2PARAM(1.0, MINRELEASETIME, MAXRELEASETIME);
  m_metaCrestTime      = 0.2; //LOG2PARAM(0.2, MINCRESTTIME, MAXCRESTTIME);
  m_metaAdaptTime      = 2.5; //LOG2PARAM(2.5, MINADAPTTIME, MAXADAPTTIME);

  m_attackCoeff    = SinglePoleCoeff(m_sampleRate, m_attackTime);
  m_releaseCoeff   = SinglePoleCoeff(m_sampleRate, m_releaseTime);
  m_metaCrestCoeff = SinglePoleCoeff(m_sampleRate, m_metaCrestTime);
  m_metaAdaptCoeff = SinglePoleCoeff(m_sampleRate, m_metaAdaptTime);
}

unsigned int CAEDSPDRCompressor::Process(float *data, unsigned int samples)
{
  if (!g_guiSettings.GetBool("audiooutput.dtshdpassthrough"))
    return samples;
  
  float *pCurrentFrame = data;

  m_returnBuffer  = data;
  m_returnSamples = samples;

  m_slope = (double)g_settings.m_currentVideoSettings.m_VolumeAmplification / 60.0;

  unsigned int frames = samples / m_channelCount; 

  for (unsigned int c = 0; c < frames; c++)
  {
    // Process
    pCurrentFrame = data + (c * m_channelCount);

    float inL = pCurrentFrame[m_channelFL];
    float inR = pCurrentFrame[m_channelFR];

    float absL = ABS(inL);
    float absR = ABS(inR);
    float maxLR = MAX(absL, absR);

    float mult = (float) ProcessSidechain(maxLR);

    for (unsigned int d = 0; d < m_channelCount; d++)
    {
      float out = pCurrentFrame[d] * mult;
      pCurrentFrame[d] = out;
    }

    // Ramp parameters
    m_logThreshold[kCurrent] = MIX(m_logThreshold[kTarget], m_logThreshold[kCurrent], m_rampCoeff);
    m_logGain[kCurrent]      = MIX(m_logGain[kTarget], m_logGain[kCurrent], m_rampCoeff);
  }

  // Kill denormals
  m_cvAttack               = KillDenormal(m_cvAttack);
  m_cvRelease              = KillDenormal(m_cvRelease);
  m_cvSmooth               = KillDenormal(m_cvSmooth);
  m_crestRms               = KillDenormal(m_crestRms);
  m_crestPeak              = KillDenormal(m_crestPeak);
  m_logThreshold[kCurrent] = KillDenormal(m_logThreshold[kCurrent]);
  m_logGain[kCurrent]      = KillDenormal(m_logGain[kCurrent]);

  return m_returnSamples;
}

double CAEDSPDRCompressor::GetDelay()
{
  return 0.0;
}

float *CAEDSPDRCompressor::GetOutput(unsigned int& samples)
{
  samples = m_returnSamples;
  return m_returnBuffer;
}

void CAEDSPDRCompressor::GetOutputFormat(CAEChannelInfo& channels, unsigned int& sampleRate)
{
  sampleRate = m_sampleRate;
  channels   = m_channels;
  return;
}

double CAEDSPDRCompressor::ProcessSidechain(double inAbs)
{
  // Square of crest factor
  double inSquare = MAX(SQUARE(inAbs), MINVAL);
  m_crestRms = MIX(inSquare, m_crestRms, m_metaCrestCoeff);
  m_crestPeak = MAX(MIX(inSquare, m_crestPeak, m_metaCrestCoeff), inSquare);
  double crestSquare = m_crestPeak / m_crestRms;

  // Attack and release coefficients
  double myAttackCoeff = m_attackCoeff;
  double myAttackTime = m_attackTime;
  if (m_bAutoAttack)
  {
    myAttackTime = 2.0 * m_metaMaxAttackTime / crestSquare;
    myAttackCoeff = SinglePoleCoeff(m_sampleRate, myAttackTime);
  }

  double myReleaseCoeff = m_releaseCoeff;
  if (m_bAutoRelease)
  {
    double myReleaseTime = 2.0 * m_metaMaxReleaseTime / crestSquare;
    myReleaseCoeff = SinglePoleCoeff(m_sampleRate, myReleaseTime - myAttackTime);
  }

  // Log conversion and overshoot
  double logIn = log(MAX(inAbs, MINVAL));
  double logOvershoot = logIn - m_logThreshold[kCurrent];

  // Set ratio/slope
  double mySlope = m_slope;
  if (m_bAutoKnee)
    mySlope = -1.0;

  // Set estimate for average CV
  double cvEstimate = m_logThreshold[kCurrent] * -mySlope / 2.0;

  // Set knee width
  double myLogWidth = m_logKneeWidth;
  if (m_bAutoKnee)
    myLogWidth = MAX(-(m_cvSmooth + cvEstimate) * m_metaKneeMult, 0.0);

  // Soft knee rectification
  double cv = 0.0;
  if (logOvershoot >= myLogWidth / 2.0)
  {
    cv = logOvershoot;
  }
  else if (logOvershoot > -myLogWidth / 2.0 && logOvershoot < myLogWidth / 2.0)
  {
    cv = 1.0 / (2.0 * myLogWidth) * SQUARE(logOvershoot + myLogWidth / 2.0);
  }

  // Multiply by negative slope for positive CV
  cv *= -mySlope;

  // Release and Attack
  m_cvRelease = MAX(cv, MIX(cv, m_cvRelease, myReleaseCoeff));
  m_cvAttack = MIX(m_cvRelease, m_cvAttack, myAttackCoeff);

  // Invert CV again
  cv = -m_cvAttack;

  // Smooth CV
  m_cvSmooth = MIX(cv - cvEstimate, m_cvSmooth, m_metaAdaptCoeff);

  // Make-up gain
  if (m_bAutoGain)
  {
    // Check for clipping
    if (m_bNoClipping && logIn + cv - (m_cvSmooth + cvEstimate) > MAXCLIPLOG)
    {
      m_cvSmooth = logIn + cv - cvEstimate - MAXCLIPLOG;
    }

    // Apply automatic gain
    cv -= m_cvSmooth + cvEstimate;
  }
  else
  {
    // Apply static gain
    cv += m_logGain[kCurrent];
  }

  // VCA
  return exp(cv);
}

// Calculate one pole filter coefficient
double CAEDSPDRCompressor::SinglePoleCoeff(unsigned int m_sampleRate, double tau)
{
  if (tau > 0.0)
  {
    return 1.0 - exp(-1.0 / (tau * (double)m_sampleRate));
  }
  return 1.0;
}

// Kill denormal number
double CAEDSPDRCompressor::KillDenormal(double value)
{
  static const double denormal = 1E-18;
  value += denormal;
  value -= denormal;
  return value;
}
