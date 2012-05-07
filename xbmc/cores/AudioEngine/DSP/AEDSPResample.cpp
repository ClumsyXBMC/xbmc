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

#include "system.h"
#include "AEDSPResample.h"
#include <string.h>

#include "utils/MathUtils.h"

CAEDSPResample::CAEDSPResample()
{
  m_ssrc = NULL;
  memset(&m_ssrcData, 0, sizeof(m_ssrcData));
}

CAEDSPResample::~CAEDSPResample()
{
  DeInitialize();
}

void CAEDSPResample::DeInitialize()
{
  if (m_ssrc)
  {
    src_delete(m_ssrc);
    m_ssrc = NULL;
  }

  if (m_ssrcData.data_out)
    _aligned_free(m_ssrcData.data_out);

  memset(&m_ssrcData, 0, sizeof(m_ssrcData));
  m_ssrcDataSamples = 0;

  m_sampleRate = 0;
  m_outputRate = 0;
  m_channels.Reset();
}

bool CAEDSPResample::Initialize(const CAEChannelInfo& channels, const unsigned int sampleRate)
{
  DeInitialize();

  int err;
  m_ssrc = src_new(SRC_SINC_BEST_QUALITY, channels.Count(), &err);

  m_channels   = channels;
  m_sampleRate = sampleRate;
  return true;
}

void CAEDSPResample::SetSampleRate(const unsigned int sampleRate)
{
  m_outputRate         = sampleRate;
  m_ssrcData.src_ratio = (double)m_outputRate / (double)m_sampleRate;
}

void CAEDSPResample::SetRatio(const double ratio)
{
  m_ssrcData.src_ratio = ratio;
  m_outputRate = (unsigned int)std::ceil(ratio * m_sampleRate);
}

double CAEDSPResample::GetRatio()
{
  return m_ssrcData.src_ratio;
}

void CAEDSPResample::GetOutputFormat(CAEChannelInfo& channels, unsigned int& sampleRate)
{
  channels   = m_channels;
  sampleRate = m_outputRate;
}

unsigned int CAEDSPResample::Process(float *data, unsigned int samples)
{
  unsigned int needSamples = samples * (unsigned int)std::ceil(m_ssrcData.src_ratio);
  if (m_ssrcDataSamples < needSamples)
  {
    if (m_ssrcData.data_out)
      _aligned_free(m_ssrcData.data_out);

    m_ssrcDataSamples        = needSamples;
    m_ssrcData.data_out      = (float*)_aligned_malloc(m_ssrcDataSamples * sizeof(float), 16);
    m_ssrcData.output_frames = m_ssrcDataSamples / m_channels.Count();
  }

  m_ssrcData.data_in      = data;
  m_ssrcData.input_frames = samples / m_channels.Count();
  if (src_process(m_ssrc, &m_ssrcData) != 0)
    return 0;

  return m_ssrcData.input_frames_used * m_channels.Count();
}

float *CAEDSPResample::GetOutput(unsigned int& samples)
{
  samples = m_ssrcData.output_frames_gen * m_channels.Count();
  return m_ssrcData.data_out;
}

double CAEDSPResample::GetDelay()
{
  return 0.0;
}
