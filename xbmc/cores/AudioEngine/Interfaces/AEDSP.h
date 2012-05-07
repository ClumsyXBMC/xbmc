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

#include "Utils/AEChannelInfo.h"
#include <vector>

class IAEDSP;
typedef std::vector<IAEDSP*> AEDSPList;

/**
 * IAEDSP Interface
 */
class IAEDSP
{
public:
  IAEDSP() {}
  virtual ~IAEDSP() {}

  /**
   * Initialize the DSP
   * @param channels The channel layout the DSP is going to be provided with
   * @param sampleRate The sample rate the data provided to the DSP is at
   * @return True if Initialize was successfull
   */
  virtual bool Initialize(const CAEChannelInfo& channels, const unsigned int sampleRate) = 0;

  /**
   * DeInitlize the DSP
   */
  virtual void DeInitialize() = 0;

  /**
   * Get the output format of the DSP
   * @param channels The channel layout the DSP is going to return
   * @param sampleRate The sample rate the DSP is going to return
   */
  virtual void GetOutputFormat(CAEChannelInfo& channels, unsigned int& sampleRate) = 0;

  /**
   * Process PCM float data
   * @param data The interleaved PCM float samples
   * @param samples The number of samples
   * @return The number of samples consumed
   */
  virtual unsigned int Process(float *data, unsigned int samples) = 0;

  /**
   * Get the output buffer
   * @param samples The number of samples in the buffer returned
   * @return the buffer containing LPCM float samples
   */
  virtual float *GetOutput(unsigned int& samples) = 0;

  /**
   * Returns the delay caused by any buffering inside the DSP
   * @return The delay in seconds that the DSP intoduces to the stream
   */
  virtual double GetDelay() = 0;

  /**
   * Called by IAE when a setting change occurs
   * @note This method is called in a thread safe manner, it is safe to re-initialize if required
   * @param The name of the setting that changed
   */
  virtual void OnSettingChange(std::string setting) {};

};

