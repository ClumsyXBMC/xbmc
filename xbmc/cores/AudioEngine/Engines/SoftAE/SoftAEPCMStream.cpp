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
#include "SoftAEPCMStream.h"

/* typecast AE to CSoftAE */
#define AE (*((CSoftAE*)CAEFactory::AE))

using namespace std;

CSoftAEPCMStream::CSoftAEPCMStream(enum AEDataFormat dataFormat, unsigned int sampleRate, unsigned int encodedSampleRate, CAEChannelInfo channelLayout, unsigned int options) :
  ISoftAEStream(dataFormat, sampleRate, encodedSampleRate, channelLayout, options),
  m_convertBuffer   (NULL ),
  m_valid           (false),
  m_delete          (false),
  m_volume          (1.0f ),
  m_rgain           (1.0f ),
  m_refillBuffer    (0    ),
  m_convertFn       (NULL ),
  m_framesBuffered  (0    ),
  m_newPacket       (NULL ),
  m_packet          (NULL ),
  m_vizPacketPos    (NULL ),
  m_draining        (false),
  m_vizBufferSamples(0    ),
  m_audioCallback   (NULL ),
  m_fadeRunning     (false)
{
  m_resample = NULL;

  m_initDataFormat        = dataFormat;
  m_initSampleRate        = sampleRate;
  m_initEncodedSampleRate = encodedSampleRate;
  m_initChannelLayout     = channelLayout;
  m_chLayoutCount         = channelLayout.Count();
  m_forceResample         = (options & AESTREAM_FORCE_RESAMPLE) != 0;
  m_paused                = (options & AESTREAM_PAUSED) != 0;
  m_autoStart             = (options & AESTREAM_AUTOSTART) != 0;

  if (m_autoStart)
    m_paused = true;

  ASSERT(m_initChannelLayout.Count());
}

void CSoftAEPCMStream::InitializeRemap()
{
  CExclusiveLock lock(m_lock);

  /* re-init the remappers */
  m_remap   .Initialize(m_initChannelLayout, AE.GetChannelLayout()           , false, false, AE.GetStdChLayout());
  m_vizRemap.Initialize(m_initChannelLayout, CAEChannelInfo(AE_CH_LAYOUT_2_0), false, true);

  /*
  if the layout has changed we need to drop data that was already remapped
  */
  if (AE.GetChannelLayout() != m_aeChannelLayout)
  {
    InternalFlush();
    m_aeChannelLayout = AE.GetChannelLayout();
    m_aeBytesPerFrame = AE.GetChannelLayout().Count() * sizeof(float);
  }
}

void CSoftAEPCMStream::Initialize()
{
  CExclusiveLock lock(m_lock);
  if (m_valid)
  {
    InternalFlush();
    delete m_newPacket;

    if (m_convert)
      _aligned_free(m_convertBuffer);

    if (m_resample)
    {
      delete m_resample;
      m_resample = NULL;
    }
  }

  if (!m_initChannelLayout.Count())
  {
    m_valid = false;
    return;
  }

  m_bytesPerSample  = (CAEUtil::DataFormatToBits(m_initDataFormat) >> 3);
  m_bytesPerFrame   = m_bytesPerSample * m_initChannelLayout.Count();

  m_aeChannelLayout = AE.GetChannelLayout();
  m_aeBytesPerFrame = AE.GetChannelLayout().Count() * sizeof(float);
  m_waterLevel      = m_initSampleRate / 2;
  m_refillBuffer    = m_waterLevel;

  m_format.m_dataFormat    = m_initDataFormat;
  m_format.m_sampleRate    = m_initSampleRate;
  m_format.m_encodedRate   = m_initEncodedSampleRate;
  m_format.m_channelLayout = m_initChannelLayout;
  m_format.m_frames        = m_initSampleRate / 8;
  m_format.m_frameSamples  = m_format.m_frames * m_initChannelLayout.Count();
  m_format.m_frameSize     = m_bytesPerFrame;

  m_newPacket = new PPacket();
  if (
    !m_remap   .Initialize(m_initChannelLayout, m_aeChannelLayout               , false, false, AE.GetStdChLayout()) ||
    !m_vizRemap.Initialize(m_initChannelLayout, CAEChannelInfo(AE_CH_LAYOUT_2_0), false, true))
  {
    m_valid = false;
    return;
  }

  m_newPacket->data.Alloc(m_format.m_frameSamples * sizeof(float));

  m_packet = NULL;

  m_inputBuffer.Alloc(m_format.m_frames * m_format.m_frameSize);

  if (m_forceResample || m_initSampleRate != AE.GetSampleRate())
  {
    CLog::Log(LOGDEBUG, "CSoftAEPCMStream::Initialize - Resampling %u to %u", m_initSampleRate, AE.GetSampleRate());
    m_resample = new CAEDSPResample();
    m_resample->Initialize(m_initChannelLayout, m_initSampleRate);
    m_resample->SetSampleRate(AE.GetSampleRate());
  }

  m_convert = m_initDataFormat != AE_FMT_FLOAT;

  /* if we need to convert, set it up */
  if (m_convert)
  {
    /* get the conversion function and allocate a buffer for the data */
    CLog::Log(LOGDEBUG, "CSoftAEPCMStream::CSoftAEPCMStream - Converting from %s to AE_FMT_FLOAT", CAEUtil::DataFormatToStr(m_initDataFormat));
    m_convertFn = CAEConvert::ToFloat(m_initDataFormat);
    if (m_convertFn)
      m_convertBuffer = (float*)_aligned_malloc(m_format.m_frameSamples * sizeof(float), 16);
    else
      m_valid         = false;
  }
  else
    m_convertBuffer = (float*)m_inputBuffer.Raw(m_format.m_frames * m_format.m_frameSize);

  m_chLayoutCount = m_format.m_channelLayout.Count();
  m_valid = true;
}

void CSoftAEPCMStream::Destroy()
{
  CExclusiveLock lock(m_lock);
  m_valid       = false;
  m_delete      = true;
}

CSoftAEPCMStream::~CSoftAEPCMStream()
{
  CExclusiveLock lock(m_lock);

  InternalFlush();
  if (m_convert)
    _aligned_free(m_convertBuffer);

  delete m_resample;
  CLog::Log(LOGDEBUG, "CSoftAEPCMStream::~CSoftAEPCMStream - Destructed");
}

/* inlined as AddData uses it */
inline unsigned int CSoftAEPCMStream::GetSpace()
{
  if (!m_valid || m_draining)
    return 0;

  if (m_framesBuffered >= m_waterLevel)
    return 0;

  return m_inputBuffer.Free() + (std::max(0U, (m_waterLevel - m_framesBuffered)) * m_format.m_frameSize);
}

unsigned int CSoftAEPCMStream::AddData(void *data, unsigned int size)
{
  CExclusiveLock lock(m_lock);
  if (!m_valid || size == 0 || data == NULL)
    return 0;

  /* if the stream is draining */
  if (m_draining)
  {
    /* if the stream has finished draining, cork it */
    if (m_packet && !m_packet->data.Used() && m_outBuffer.empty())
      m_draining = false;
    else
      return 0;
  }

  /* dont ever take more then GetSpace advertises */
  size = std::min(size, GetSpace());
  if (size == 0)
    return 0;

  unsigned int taken = 0;
  while(size)
  {
    unsigned int copy = std::min((unsigned int)m_inputBuffer.Free(), size);
    if (copy > 0)
    {
      m_inputBuffer.Push(data, copy);
      size  -= copy;
      taken += copy;
      data   = (uint8_t*)data + copy;
    }

    if (m_inputBuffer.Free() == 0)
    {
      unsigned int consumed = ProcessFrameBuffer();
      m_inputBuffer.Shift(NULL, consumed);
    }
  }

  return taken;
}

unsigned int CSoftAEPCMStream::ProcessFrameBuffer()
{
  uint8_t     *data;
  unsigned int frames, consumed, sampleSize;

  /* convert the data if we need to */
  unsigned int samples;
  if (m_convert)
  {
    data       = (uint8_t*)m_convertBuffer;
    samples    = m_convertFn(
      (uint8_t*)m_inputBuffer.Raw(m_inputBuffer.Used()),
      m_inputBuffer.Used() / m_bytesPerSample,
      m_convertBuffer
    );
    sampleSize = sizeof(float);
  }
  else
  {
    data       = (uint8_t*)m_inputBuffer.Raw(m_inputBuffer.Used());
    samples    = m_inputBuffer.Used() / m_bytesPerSample;
    sampleSize = m_bytesPerSample;
  }

  if (samples == 0)
    return 0;

  /* resample it if we need to */
  if (m_resample)
  {
    consumed = m_resample->Process((float*)data, samples) * m_bytesPerSample;
    data = (uint8_t*)m_resample->GetOutput(samples);
    m_resample->GetOutput(samples);
    if (!samples)
      return consumed;

    frames = samples / m_chLayoutCount;
  }
  else
  {
    data     = (uint8_t*)m_convertBuffer;
    frames   = samples / m_chLayoutCount;
    consumed = frames * m_bytesPerFrame;
  }

  if (m_refillBuffer)
  {
    if (frames > m_refillBuffer)
      m_refillBuffer = 0;
    else
      m_refillBuffer -= frames;
  }

  /* buffer the data */
  m_framesBuffered += frames;

  size_t remaining = samples * sampleSize;
  while (remaining)
  {
    size_t copy = std::min(m_newPacket->data.Free(), remaining);
    m_newPacket->data.Push(data, copy);
    data      += copy;
    remaining -= copy;

    /* wait till we have a full packet, or no more data before processing the packet */
    if (!m_draining && m_newPacket->data.Free() > 0)
      continue;

    /* make a new packet for downmix/remap */
    PPacket *pkt = new PPacket();

    /* downmix/remap the data */
    size_t frames = m_newPacket->data.Used() / m_chLayoutCount / sizeof(float);
    size_t used   = frames * m_aeChannelLayout.Count() * sizeof(float);
    pkt->data.Alloc(used);
    m_remap.Remap(
      (float*)m_newPacket->data.Raw (m_newPacket->data.Used()),
      (float*)pkt        ->data.Take(used),
      frames
    );

    /* downmix for the viz if we have one */
    if (m_audioCallback)
    {
      size_t vizUsed = frames * 2 * sizeof(float);
      pkt->vizData.Alloc(vizUsed);
      m_vizRemap.Remap(
        (float*)m_newPacket->data   .Raw (m_newPacket->data.Used()),
        (float*)pkt        ->vizData.Take(vizUsed),
        frames
      );
    }

    /* add the packet to the output */
    m_outBuffer.push_back(pkt);
    m_newPacket->data.Empty();
  }

  /* if the stream is flagged to autoStart when the buffer is full, then do it */
  if (m_autoStart && m_framesBuffered >= m_waterLevel)
    Resume();

  return consumed;
}

uint8_t* CSoftAEPCMStream::GetFrame()
{
  CExclusiveLock lock(m_lock);

  /* if we are fading, this runs even if we have underrun as it is time based */
  if (m_fadeRunning)
  {
    m_volume += m_fadeStep;
    m_volume = std::min(1.0f, std::max(0.0f, m_volume));
    if (m_fadeDirUp)
    {
      if (m_volume >= m_fadeTarget)
        m_fadeRunning = false;
    }
    else
    {
      if (m_volume <= m_fadeTarget)
        m_fadeRunning = false;
    }
  }

  /* if we have been deleted or are refilling but not draining */
  if (!m_valid || m_delete || (m_refillBuffer && !m_draining))
    return NULL;

  /* if the packet is empty, advance to the next one */
  if (!m_packet || m_packet->data.CursorEnd())
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
        CLog::Log(LOGDEBUG, "CSoftAEPCMStream::GetFrame - Underrun");
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
  uint8_t *ret = (uint8_t*)m_packet->data.CursorRead(m_aeBytesPerFrame);

  /* we have a frame, if we have a viz we need to hand the data to it */
  if (m_audioCallback && !m_packet->vizData.CursorEnd())
  {
    float *vizData = (float*)m_packet->vizData.CursorRead(2 * sizeof(float));
    memcpy(m_vizBuffer + m_vizBufferSamples, vizData, 2 * sizeof(float));
    m_vizBufferSamples += 2;
    if (m_vizBufferSamples == 512)
    {
      m_audioCallback->OnAudioData(m_vizBuffer, 512);
      m_vizBufferSamples = 0;
    }
  }

  --m_framesBuffered;
  return ret;
}

double CSoftAEPCMStream::GetDelay()
{
  if (m_delete)
    return 0.0;

  double delay = AE.GetDelay();
  delay += (double)(m_inputBuffer.Used() / m_format.m_frameSize) / (double)m_format.m_sampleRate;
  delay += (double)m_framesBuffered                              / (double)AE.GetSampleRate();

  if (m_resample)
    delay += m_resample->GetDelay();

  return delay;
}

double CSoftAEPCMStream::GetCacheTime()
{
  if (m_delete)
    return 0.0;

  double time;
  time  = (double)(m_inputBuffer.Free() / m_format.m_frameSize) / (double)m_format.m_sampleRate;
  time += (double)(m_waterLevel - m_refillBuffer)               / (double)AE.GetSampleRate();
  time += AE.GetCacheTime();
  return time;
}

double CSoftAEPCMStream::GetCacheTotal()
{
  if (m_delete)
    return 0.0;

  double total;
  total  = (double)(m_inputBuffer.Size() / m_format.m_frameSize) / (double)m_format.m_sampleRate;
  total += (double)m_waterLevel                                  / (double)AE.GetSampleRate();
  total += AE.GetCacheTotal();
  return total;
}

void CSoftAEPCMStream::Pause()
{
  if (m_paused)
    return;
  m_paused = true;
  AE.PauseStream(this);
}

void CSoftAEPCMStream::Resume()
{
  if (!m_paused)
    return;
  m_paused    = false;
  m_autoStart = false;
  AE.ResumeStream(this);
}

void CSoftAEPCMStream::Drain()
{
  CSharedLock lock(m_lock);
  m_draining = true;
}

bool CSoftAEPCMStream::IsDrained()
{
  CSharedLock lock(m_lock);
  return (m_draining && !m_packet && m_outBuffer.empty());
}

void CSoftAEPCMStream::Flush()
{
  CLog::Log(LOGDEBUG, "CSoftAEPCMStream::Flush");
  CExclusiveLock lock(m_lock);
  InternalFlush();

  /* internal flush does not do this as these samples are still valid if we are re-initializing */
  m_inputBuffer.Empty();
}

void CSoftAEPCMStream::InternalFlush()
{
  /* invalidate any incoming samples */
  m_newPacket->data.Empty();

  /*
    clear the current buffered packet, we cant delete the data as it may be
    in use by the AE thread, so we just seek to the end of the buffer
  */
  if (m_packet)
    m_packet->data.CursorSeek(m_packet->data.Used());

  /* clear any other buffered packets */
  while (!m_outBuffer.empty())
  {
    PPacket *p = m_outBuffer.front();
    m_outBuffer.pop_front();
    delete p;
  }

  /* reset our counts */
  m_framesBuffered = 0;
  m_refillBuffer   = m_waterLevel;
}

double CSoftAEPCMStream::GetResampleRatio()
{
  if (!m_resample)
    return 1.0f;

  CSharedLock lock(m_lock);
  return m_resample->GetRatio();
}

bool CSoftAEPCMStream::SetResampleRatio(double ratio)
{
  if (!m_resample)
    return false;

  CSharedLock lock(m_lock);
  m_resample->SetRatio(ratio);
  return true;
}

void CSoftAEPCMStream::RegisterAudioCallback(IAudioCallback* pCallback)
{
  CExclusiveLock lock(m_lock);
  m_vizBufferSamples = 0;
  m_audioCallback = pCallback;
  if (m_audioCallback)
    m_audioCallback->OnInitialize(2, m_initSampleRate, 32);
}

void CSoftAEPCMStream::UnRegisterAudioCallback()
{
  CExclusiveLock lock(m_lock);
  m_audioCallback = NULL;
  m_vizBufferSamples = 0;
}

void CSoftAEPCMStream::FadeVolume(float from, float target, unsigned int time)
{
  CExclusiveLock lock(m_lock);
  float delta   = target - from;
  m_fadeDirUp   = target > from;
  m_fadeTarget  = target;
  m_fadeStep    = delta / (((float)AE.GetSampleRate() / 1000.0f) * (float)time);
  m_fadeRunning = true;
}

bool CSoftAEPCMStream::IsFading()
{
  CSharedLock lock(m_lock);
  return m_fadeRunning;
}

void CSoftAEPCMStream::RegisterSlave(IAEStream *slave)
{
  ISoftAEStream *s = (ISoftAEStream*)slave;

  CSharedLock lock(m_lock);
  m_slave = s;
  s->m_parent = this;
}
