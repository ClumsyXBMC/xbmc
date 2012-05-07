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
#include "threads/SingleLock.h"
#include "utils/log.h"
#include "utils/MathUtils.h"

#include "AEFactory.h"
#include "Utils/AEUtil.h"

#include "SoftAE.h"
#include "SoftAERAWStream.h"

/* typecast AE to CSoftAE */
#define AE (*((CSoftAE*)CAEFactory::AE))

using namespace std;

CSoftAERAWStream::CSoftAERAWStream(enum AEDataFormat dataFormat, unsigned int sampleRate, unsigned int encodedSampleRate, CAEChannelInfo channelLayout, unsigned int options) :
  ISoftAEStream(dataFormat, sampleRate, encodedSampleRate, channelLayout, options),
  m_valid           (false),
  m_delete          (false),
  m_refillBuffer    (0    ),
  m_inputBuffer     (NULL ),
  m_framesBuffered  (0    ),
  m_packet          (NULL ),
  m_draining        (false)
{
  m_initDataFormat        = dataFormat;
  m_initSampleRate        = sampleRate;
  m_initEncodedSampleRate = encodedSampleRate;
  m_initChannelLayout     = channelLayout;
  m_chLayoutCount         = channelLayout.Count();
  m_paused                = (options & AESTREAM_PAUSED) != 0;
  m_autoStart             = (options & AESTREAM_AUTOSTART) != 0;

  if (m_autoStart)
    m_paused = true;

  ASSERT(m_initChannelLayout.Count());
}

void CSoftAERAWStream::InitializeRemap()
{
}

void CSoftAERAWStream::Initialize()
{
  CExclusiveLock lock(m_lock);
  if (m_valid)
    InternalFlush();

  enum AEDataFormat useDataFormat;

  /* we are raw, which means we need to work in the output format */
  useDataFormat       = AE.GetSinkDataFormat();
  m_initChannelLayout = AE.GetSinkChLayout  ();
  m_samplesPerFrame   = m_initChannelLayout.Count();

  m_bytesPerSample  = (CAEUtil::DataFormatToBits(useDataFormat) >> 3);
  m_bytesPerFrame   = m_bytesPerSample * m_initChannelLayout.Count();

  m_aeChannelLayout = AE.GetChannelLayout();
  m_waterLevel      = m_initSampleRate;
  m_refillBuffer    = m_waterLevel;

  m_format.m_dataFormat    = useDataFormat;
  m_format.m_sampleRate    = m_initSampleRate;
  m_format.m_encodedRate   = m_initEncodedSampleRate;
  m_format.m_channelLayout = m_initChannelLayout;
  m_format.m_frames        = m_initSampleRate / 2;
  m_format.m_frameSamples  = m_format.m_frames * m_initChannelLayout.Count();
  m_format.m_frameSize     = m_bytesPerFrame;
  m_chLayoutCount          = m_format.m_channelLayout.Count();

  if (!m_inputBuffer)
    m_inputBuffer = new CAEBuffer();

  m_inputBuffer->Alloc(m_format.m_frames * m_format.m_channelLayout.Count() * m_bytesPerFrame);

  m_packet      = NULL;
  m_valid       = true;
}

void CSoftAERAWStream::Destroy()
{
  CExclusiveLock lock(m_lock);
  m_valid       = false;
  m_delete      = true;
}

CSoftAERAWStream::~CSoftAERAWStream()
{
  CExclusiveLock lock(m_lock);

  InternalFlush();
  CLog::Log(LOGDEBUG, "CSoftAERAWStream::~CSoftAERAWStream - Destructed");
}

inline unsigned int CSoftAERAWStream::GetSpace()
{
  if (!m_valid || m_draining)
    return 0;

  if (m_framesBuffered >= m_waterLevel)
    return 0;

  return ((m_waterLevel - m_framesBuffered) * m_bytesPerFrame) + m_inputBuffer->Free();
}

unsigned int CSoftAERAWStream::AddData(void *data, unsigned int size)
{
  if (!m_valid || size == 0 || data == NULL)
    return 0;

  /* if the stream is draining */
  if (m_draining)
  {
    /* if the stream has finished draining, cork it */
    if (m_packet && !m_packet->Used() && m_outBuffer.empty())
      m_draining = false;
    else
      return 0;
  }

  /* honour the size reported by GetSpace() */
  size = std::min(size, GetSpace());
  if (size == 0)
    return 0;

  unsigned int used = 0;
  while(size > 0)
  {
    size_t copy = std::min(m_inputBuffer->Free(), (size_t)size);
    if (copy > 0)
    {
      m_inputBuffer->Push(data, copy);
      used += copy;
      size -= copy;
      data  = (uint8_t*)data + copy;
    }

    if (m_inputBuffer->Free() == 0)
    {
      int frames = m_inputBuffer->Used() / m_bytesPerFrame;

      CExclusiveLock lock(m_lock);
      m_outBuffer.push_back(m_inputBuffer);
      m_inputBuffer = new CAEBuffer();
      m_inputBuffer->Alloc(m_format.m_frames * m_format.m_channelLayout.Count() * m_bytesPerFrame);
      m_framesBuffered += frames;
      m_refillBuffer    = std::max(0, (int)m_refillBuffer - (int)frames);
      lock.Leave();

      /* if the stream is flagged to autoStart when the buffer is full, then do it */
      if (m_autoStart && m_framesBuffered >= m_waterLevel)
        Resume();
    }
  }

  return used;
}

uint8_t* CSoftAERAWStream::GetFrame()
{
  CExclusiveLock lock(m_lock);

  /* if we have been deleted or are refilling but not draining */
  if (!m_valid || m_delete || (m_refillBuffer && !m_draining))
    return NULL;

  /* if the packet is empty, advance to the next one */
  if (!m_packet || m_packet->CursorEnd())
  {
    delete m_packet;
    m_packet = NULL;

    /* no more packets, return null */
    if (m_outBuffer.empty())
    {
      if (m_draining)
        return NULL;
      else
      {
        /* underrun, we need to refill our buffers */
        CLog::Log(LOGDEBUG, "CSoftAERAWStream::GetFrame - Underrun");
        ASSERT(m_waterLevel > m_framesBuffered);
        m_refillBuffer = m_waterLevel - m_framesBuffered;
        return NULL;
      }
    }

    /* get the next packet */
    m_packet = m_outBuffer.front();
    m_outBuffer.pop_front();
  }

  /* fetch one frame of data */
  uint8_t *ret = (uint8_t*)m_packet->CursorRead(m_bytesPerFrame);

  --m_framesBuffered;
  return ret;
}

double CSoftAERAWStream::GetDelay()
{
  if (m_delete)
    return 0.0;

  int frames = m_inputBuffer->Used() / m_bytesPerFrame;
  frames += m_framesBuffered;

  double delay = AE.GetDelay();
  delay += (double)frames / (double)m_initSampleRate;

  return delay;
}

double CSoftAERAWStream::GetCacheTime()
{
  if (m_delete)
    return 0.0;

  int frames = m_inputBuffer->Free() / m_format.m_frameSize;
  frames += m_waterLevel - m_refillBuffer;

  double time = AE.GetCacheTime();
  time += (double)frames / (double)m_initSampleRate;

  return time;
}

double CSoftAERAWStream::GetCacheTotal()
{
  if (m_delete)
    return 0.0;

  int frames = m_inputBuffer->Size() / m_format.m_frameSize;
  frames += m_waterLevel;

  double total = AE.GetCacheTotal();
  total += (double)frames / (double)m_initSampleRate;

  return total;
}

void CSoftAERAWStream::Pause()
{
  if (m_paused)
    return;
  m_paused = true;
  AE.PauseStream(this);
}

void CSoftAERAWStream::Resume()
{
  if (!m_paused)
    return;
  m_paused    = false;
  m_autoStart = false;
  AE.ResumeStream(this);
}

void CSoftAERAWStream::Drain()
{
  CSharedLock lock(m_lock);
  m_draining = true;
}

bool CSoftAERAWStream::IsDrained()
{
  CSharedLock lock(m_lock);
  return (m_draining && !m_packet && m_outBuffer.empty());
}

void CSoftAERAWStream::Flush()
{
  CLog::Log(LOGDEBUG, "CSoftAERAWStream::Flush");
  CExclusiveLock lock(m_lock);
  InternalFlush();

  /* internal flush does not do this as these samples are still valid if we are re-initializing */
  m_inputBuffer->Empty();
}

void CSoftAERAWStream::InternalFlush()
{
  /* invalidate any incoming samples */
  m_inputBuffer->Empty();

  /*
    clear the current buffered packet, we cant delete the data as it may be
    in use by the AE thread, so we just seek to the end of the buffer
  */
  if (m_packet)
    m_packet->CursorSeek(m_packet->Used());

  /* clear any other buffered packets */
  while (!m_outBuffer.empty())
  {
    CAEBuffer *p = m_outBuffer.front();
    m_outBuffer.pop_front();
    delete p;
  }

  /* reset our counts */
  m_framesBuffered = 0;
  m_refillBuffer   = m_waterLevel;
}

void CSoftAERAWStream::RegisterSlave(IAEStream *slave)
{
  ISoftAEStream *s = (ISoftAEStream*)slave;

  CSharedLock lock(m_lock);
  m_slave = s;
  s->m_parent = this;
}

