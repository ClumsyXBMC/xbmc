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
 */

#include "Interfaces/AEDSP.h"
#include <samplerate.h>

/**
 * CAEDSPResample - A basic resampler that uses SSRC
 */
class CAEDSPResample : public IAEDSP
{
public:
  CAEDSPResample();
  virtual ~CAEDSPResample();

  virtual bool Initialize(const CAEChannelInfo& channels, const unsigned int sampleRate);

  void   SetSampleRate(const unsigned int sampleRate);
  void   SetRatio(const double ratio);
  double GetRatio();

  virtual void DeInitialize();

  virtual void GetOutputFormat(CAEChannelInfo& channels, unsigned int& sampleRate);

  virtual unsigned int Process(float *data, unsigned int samples);

  virtual float *GetOutput(unsigned int& samples);

  virtual double GetDelay();

private:
  CAEChannelInfo m_channels;
  unsigned int   m_sampleRate, m_outputRate;
  float*         m_buffer;

  SRC_STATE*     m_ssrc;
  unsigned int   m_ssrcDataSamples;
  SRC_DATA       m_ssrcData;
};

