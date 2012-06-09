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

#include <math.h>
#include "Interfaces/AEDSP.h"

// Macros used in calcs
#define MIN(a, b) (((a)<(b))?(a):(b))
#define MAX(a, b) (((a)>(b))?(a):(b))
#define CLAMP(x, min, max) ((x)>(max)?(max):(((x)<(min))?(min):(x)))
#define ABS(x) ((x)<0?(-(x)):(x))
#define SQUARE(x) ((x)*(x))
#define MIX(x0, y1, coeff) (((x0)-(y1))*(coeff)+(y1))

#define PARAM2LOG(x, mn, mx) (exp(log(mn) + (x) * (log(mx) - log(mn))))
#define LOG2PARAM(x, mn, mx) ((log(x) - log(mn)) / (log(mx) - log(mn)))

#define PARAM2LIN(x, mn, mx) ((mn) + (x) * ((mx) - (mn)))
#define LIN2PARAM(x, mn, mx) (((x) - (mn)) / ((mx) - (mn)))

// Constants
static const double MINVAL         = 1.0e-6;
static const double MINKNEEMULT    = 0.0;
static const double MAXKNEEMULT    = 4.0;
static const double MINATTACKTIME  = 0.0001;
static const double MAXATTACKTIME  = 0.2;
static const double MINRELEASETIME = 0.005;
static const double MAXRELEASETIME = 2.0;
static const double MAXGAINDB      = 60.0;
static const double MINADAPTTIME   = MAXRELEASETIME / 2.0;
static const double MAXADAPTTIME   = MAXRELEASETIME * 2.0;
static const double MINCRESTTIME   = MINRELEASETIME;
static const double MAXCRESTTIME   = MAXRELEASETIME;
static const double MINMETERLOG    = log(pow(10.0, -MAXGAINDB / 20.0));
static const double MAXCLIPLOG     = log(pow(10.0, -0.01 / 20.0));

static const int kCurrent = 0;
static const int kTarget  = 1;

///////////////////////////////////////////////////////////////////////////

class CAEDSPDRCompressor : public IAEDSP
{
public:
  CAEDSPDRCompressor();
  virtual ~CAEDSPDRCompressor();

  /** Must call Initialize first to set DSP parameters */
  virtual bool Initialize(const CAEChannelInfo& channels, const unsigned int sampleRate);

  virtual void DeInitialize();

  virtual void GetOutputFormat(CAEChannelInfo& channels, unsigned int& sampleRate);

  virtual unsigned int Process(float *data, unsigned int samples);

  virtual float *GetOutput(unsigned int& samples);

  virtual double GetDelay();

  virtual void   OnSettingChange(std::string setting);

//////////////////////////////////////////////////////////////////////////////////////////

private:
  double ProcessSidechain(double absIn);
  double SinglePoleCoeff(unsigned int m_sampleRate, double tau);
  double KillDenormal(double value);

  void   SetParameterDefaults();

  double m_rampCoeff;

  /* Use auto adjust */
  bool   m_bAutoKnee;
  bool   m_bAutoGain;
  bool   m_bAutoAttack;
  bool   m_bAutoRelease;

  double m_cvRelease;
  double m_cvAttack;
  double m_cvSmooth;

  double m_logThreshold[2];
  double m_logKneeWidth;
  double m_logGain[2];
  double m_slope;

  double m_crestPeak;
  double m_crestRms;

  double m_attackTime;
  double m_attackCoeff;

  double m_releaseTime;
  double m_releaseCoeff;

  double m_metaKneeMult;
  double m_metaMaxAttackTime;
  double m_metaMaxReleaseTime;
  double m_metaCrestTime;
  double m_metaCrestCoeff;
  double m_metaAdaptTime;
  double m_metaAdaptCoeff;
  bool   m_bNoClipping;

  unsigned int   m_sampleRate;
  unsigned int   m_channelCount;
  unsigned int   m_channelFL;
  unsigned int   m_channelFR;
  CAEChannelInfo m_channels;

  unsigned int m_returnSamples;
  float*       m_returnBuffer;
};
