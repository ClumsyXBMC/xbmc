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

#include "AEAudioFormat.h"
#include "Interfaces/AEStream.h"

class ISoftAEStream : public IAEStream
{
protected:
  friend class CSoftAE;
  friend class CSoftAEPCMStream;
  friend class CSoftAERAWStream;

  ISoftAEStream(enum AEDataFormat format, unsigned int sampleRate, unsigned int encodedSamplerate, CAEChannelInfo channelLayout, unsigned int options) :
    m_parent(NULL),
    m_slave (NULL) {}

  virtual ~ISoftAEStream() {}

  virtual void Initialize() = 0;
  virtual void InitializeRemap() = 0;
  virtual void Destroy() = 0;
  virtual uint8_t* GetFrame() = 0;

  virtual bool IsPaused   () = 0;
  virtual void SetPaused  (const bool paused) = 0;
  virtual bool IsDestroyed() = 0;
  virtual bool IsValid    () = 0;
  virtual const bool IsRaw() const = 0;
  virtual const CAEChannelInfo& GetChannelLayout() const = 0;

  /* parent/slave stream */
  ISoftAEStream *m_parent;
  ISoftAEStream *m_slave;

private:
};

