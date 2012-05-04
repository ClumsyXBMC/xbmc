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

#include <list>

#include "threads/SharedSection.h"

#include "AEAudioFormat.h"
#include "ISoftAEStream.h"
#include "Utils/AEConvert.h"
#include "Utils/AERemap.h"
#include "Utils/AEBuffer.h"

#include "DSP/AEDSPResample.h"

class CSoftAEPCMStream : public ISoftAEStream
{
protected:
  friend class CSoftAE;
  CSoftAEPCMStream(enum AEDataFormat format, unsigned int sampleRate, unsigned int encodedSamplerate, CAEChannelInfo channelLayout, unsigned int options);
  virtual ~CSoftAEPCMStream();

  virtual void Initialize();
  virtual void InitializeRemap();
  virtual void Destroy();
  virtual uint8_t* GetFrame();

  virtual bool IsPaused   () { return m_paused; }
  virtual void SetPaused  (const bool paused) { m_paused = paused; }
  virtual bool IsDestroyed() { return m_delete; }
  virtual bool IsValid    () { return m_valid;  }
  virtual const bool IsRaw() const { return false; } 
  virtual const CAEChannelInfo& GetChannelLayout() const { return m_initChannelLayout; }

public:
  virtual unsigned int      GetSpace        ();
  virtual unsigned int      AddData         (void *data, unsigned int size);
  virtual double            GetDelay        ();
  virtual bool              IsBuffering     () { return m_refillBuffer > 0; }
  virtual double            GetCacheTime    ();
  virtual double            GetCacheTotal   ();

  virtual void              Pause           ();
  virtual void              Resume          ();
  virtual void              Drain           ();
  virtual bool              IsDraining      () { return m_draining;    }
  virtual bool              IsDrained       ();
  virtual void              Flush           ();

  virtual float             GetVolume       ()             { return m_volume; }
  virtual float             GetReplayGain   ()             { return m_rgain ; }
  virtual void              SetVolume       (float volume) { m_volume = std::max( 0.0f, std::min(1.0f, volume)); }
  virtual void              SetReplayGain   (float factor) { m_rgain  = std::max(-1.0f, std::max(1.0f, factor)); }

  virtual const unsigned int      GetFrameSize   () const  { return m_format.m_frameSize; }
  virtual const unsigned int      GetChannelCount() const  { return m_initChannelLayout.Count(); }
  
  virtual const unsigned int      GetSampleRate  () const  { return m_initSampleRate; }
  virtual const unsigned int GetEncodedSampleRate() const  { return m_initEncodedSampleRate; }
  virtual const enum AEDataFormat GetDataFormat  () const  { return m_initDataFormat; }
  
  virtual double            GetResampleRatio();
  virtual bool              SetResampleRatio(double ratio);
  virtual void              RegisterAudioCallback(IAudioCallback* pCallback);
  virtual void              UnRegisterAudioCallback();
  virtual void              FadeVolume(float from, float to, unsigned int time);
  virtual bool              IsFading();
  virtual void              RegisterSlave(IAEStream *stream);
private:
  void InternalFlush();
  void CheckResampleBuffers();

  CSharedSection    m_lock;
  enum AEDataFormat m_initDataFormat;
  unsigned int      m_initSampleRate;
  unsigned int      m_initEncodedSampleRate;
  CAEChannelInfo    m_initChannelLayout;
  unsigned int      m_chLayoutCount;
  
  typedef struct
  {
    CAEBuffer data;
    CAEBuffer vizData;
  } PPacket;

  AEAudioFormat m_format;

  bool                    m_forceResample; /* true if we are to force resample even when the rates match */
  CAEDSPResample         *m_resample;      /* not null if the audio is to be resampled */
  bool                    m_convert;       /* true if the bitspersample needs converting */
  float                  *m_convertBuffer; /* buffer for converted data */
  bool                    m_valid;         /* true if the stream is valid */
  bool                    m_delete;        /* true if CSoftAE is to free this object */
  CAERemap                m_remap;         /* the remapper */
  float                   m_volume;        /* the volume level */
  float                   m_rgain;         /* replay gain level */
  unsigned int            m_waterLevel;    /* the fill level to fall below before calling the data callback */
  unsigned int            m_refillBuffer;  /* how many frames that need to be buffered before we return any frames */

  CAEConvert::AEConvertToFn m_convertFn;

  CAEBuffer           m_inputBuffer;
  unsigned int        m_bytesPerSample;
  unsigned int        m_bytesPerFrame;
  CAEChannelInfo      m_aeChannelLayout;
  unsigned int        m_aeBytesPerFrame;
  unsigned int        m_framesBuffered;
  std::list<PPacket*> m_outBuffer;
  unsigned int        ProcessFrameBuffer();
  PPacket            *m_newPacket;
  PPacket            *m_packet;
  uint8_t            *m_packetPos;
  float              *m_vizPacketPos;
  bool                m_paused;
  bool                m_autoStart;
  bool                m_draining;

  /* vizualization internals */
  CAERemap           m_vizRemap;
  float              m_vizBuffer[512];
  unsigned int       m_vizBufferSamples;
  IAudioCallback    *m_audioCallback;

  /* fade values */
  bool               m_fadeRunning;
  bool               m_fadeDirUp;
  float              m_fadeStep;
  float              m_fadeTarget;
  unsigned int       m_fadeTime;
};

