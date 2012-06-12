/**
 * @file Streams.cpp
 * @author Pekka Mikkola <pjmikkol@cs.helsinki.fi>
 * @author Dominik Kempa <dominik.kempa@cs.helsinki.fi>
 *
 * @section LICENSE
 *
 * This file is part of bwtc.
 *
 * bwtc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * bwtc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with bwtc.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @section DESCRIPTION
 *
 * Implementation for InStream and OutStream.
 */

#include <cassert>

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <algorithm>

#include "globaldefs.hpp"
#include "Streams.hpp"

#include "Profiling.hpp"

namespace bwtc {

OutStream::OutStream(const std::string& file_name)
    :m_name(file_name) {
  m_buffer = new byte[kBufferSize];
  if (!m_buffer) {
    fprintf(stderr,"Memory allocation error!\n");
    exit(1);
  }
  if (m_name != "") {
    m_fileptr = fopen(m_name.c_str(), "w");
    if (!m_fileptr) {
      perror(m_name.c_str());
      exit(1);
    }
  } else {
    m_fileptr = stdout;
  }
  m_filled = 0;
  assert(m_fileptr);
}

OutStream::~OutStream() {
  PROFILE("OutStream::~OutStream()");
  if (m_filled > 0) {
    flush();
  }
  if (m_fileptr != stdout) {
    fclose(m_fileptr);
  }
  delete[] m_buffer;
}

void OutStream::writeByte(byte b) {
  m_buffer[m_filled++] = b;
  if (m_filled == kBufferSize) {
    fwrite(m_buffer, 1, m_filled, m_fileptr);
    m_filled = 0;
  }
}

void OutStream::flush() {
  if (m_filled > 0) {
    fwrite(m_buffer, 1, m_filled, m_fileptr);
    m_filled = 0;
  }
  fflush(m_fileptr);
}

long int OutStream::getPos() {
  flush();
  return ftell(m_fileptr);
}

void OutStream::writeBlock(byte *begin, byte *end) {
  uint64 size = end - begin;
  if (m_filled + size < kBufferSize) {
    std::copy(begin, end, m_buffer + m_filled);
    m_filled += size;
  } else {
    flush();
    fwrite(begin, 1, size, m_fileptr);
  }
}

void OutStream::write48bits(uint64 to_written, long int position) {
  assert((to_written & (((uint64)0xFFFF) << 48)) == 0);
  flush();
  long int current = ftell(m_fileptr);
  // fseeking with negative offset from SEEK_CUR might be faster.
  fseek(m_fileptr, position, SEEK_SET);
  for(int i = 5; i >= 0; --i) {
    byte b = 0xFF & (to_written >> i*8);
    fputc(b, m_fileptr);
  }
  fseek(m_fileptr, current, SEEK_SET);
}

uint64 InStream::read48bits() {
  uint64 result = 0;
  for(int i = 0; i < 6; ++i) {
    result <<= 8;
    result |= static_cast<byte>(fetchByte());
  }
  return result;
}

InStream::InStream(const std::string &file_name) :
    m_name(file_name), m_buffer(0), m_bitsInBuffer(0)
{
  m_bigbuf = new byte[kBigbufSize];
  if (!m_bigbuf) {
    fprintf(stderr,"Memory allocation error!\n");
    exit(1);
  }
  if (m_name != "") {
    m_fileptr = fopen(m_name.c_str(), "r");
    if (!m_fileptr) {
      perror(m_name.c_str());
      exit(1);
    }
  } else {
    m_fileptr = stdin;
  }
  m_bigbuf_left = 0;
  m_bigbuf_pos = 0;
  assert(m_fileptr);
}

InStream::~InStream() {
  if (m_fileptr != stdin) {
    fclose(m_fileptr);
  }
  delete[] m_bigbuf;
}

uint64 InStream::readBlock(byte* to, uint64 max_block_size) {
  assert(m_bitsInBuffer == 0);
  uint64 have_read = 0;
  int32 c;
  while (have_read < max_block_size && (c = fetchByte()) != EOF) {
    to[have_read++] = static_cast<byte>(c);
  }
  return have_read;
}

int32 InStream::fetchByte() {
  if (m_bigbuf_left <= 0) {
    m_bigbuf_pos = 0;
    m_bigbuf_left = fread(m_bigbuf, 1, kBigbufSize, m_fileptr);
    if (m_bigbuf_left <= 0) return EOF;
  }
  --m_bigbuf_left;
  return m_bigbuf[m_bigbuf_pos++];
}

int32 InStream::peekByte() {
  if (m_bigbuf_left <= 0) {
    m_bigbuf_pos = 0;
    m_bigbuf_left = fread(m_bigbuf, 1, kBigbufSize, m_fileptr);
    if (m_bigbuf_left <= 0) return EOF;
  }
  return m_bigbuf[m_bigbuf_pos];
}

} //namespace bwtc

