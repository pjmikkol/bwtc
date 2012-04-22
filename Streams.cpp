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
 * Implementation for InStream, OutStream, RawInStream and RawOutStream.
 */

#include <cassert>

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <algorithm>

#include "globaldefs.hpp"
#include "Streams.hpp"

namespace bwtc {

OutStream::OutStream(std::string file_name) :
    m_name(file_name), m_to(NULL), m_outfile(NULL)
{
  if (m_name != "") {
    m_outfile = new std::ofstream(m_name.c_str());
    m_to = m_outfile;
  } else
    m_to = &std::cout;
  assert(m_to);
}

OutStream::~OutStream() {
  if (m_outfile) {
    m_outfile->close();
    delete m_outfile;
  }
}

void OutStream::flush() {
  m_to->flush();
}

void OutStream::writeByte(byte b) {
  m_to->put(b);
}

std::streampos OutStream::getPos() const {
  return m_to->tellp();
}

void OutStream::writeBlock(std::vector<byte>::const_iterator begin,
                           std::vector<byte>::const_iterator end) {
  std::copy(begin, end, std::ostream_iterator<byte>(*m_to));
}

void OutStream::write48bits(uint64 to_written, std::streampos position) {
  assert((to_written & (((uint64)0xFFFF) << 48)) == 0 );
  std::streampos current = m_to->tellp();
  m_to->seekp(position);
  for(int i = 5; i >= 0; --i) {
    byte b = 0xFF & (to_written >> i*8);
    m_to->put(b);
  }
  m_to->flush();
  m_to->seekp(current);
}

uint64 InStream::read48bits() {
  uint64 result = 0;
  for(int i = 0; i < 6; ++i) {
    result <<= 8;
    result |= readByte();
  }
  return result;
}

InStream::InStream(const std::string& file_name) :
    m_name(file_name), m_from(0), m_infile(0), m_buffer(0),
    m_bitsInBuffer(0)
{
  if (m_name != "") {
    m_infile = new std::ifstream(m_name.c_str());
    m_from = dynamic_cast<std::istream*>(m_infile);
  } else {
    m_from = &std::cin;
  }
  assert(m_from);
}

InStream::~InStream() {
  if (m_infile) {
    m_infile->close();
    delete m_infile;
  }
}

std::streamsize InStream::readBlock(byte* to, std::streamsize max_block_size) {
  if(!*m_from) return 0; /* stream has reached EOF */
  assert(m_bitsInBuffer == 0);
  m_from->read(reinterpret_cast<char*>(to), max_block_size);
  return m_from->gcount();
}

RawOutStream::RawOutStream(std::string file_name)
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

RawOutStream::~RawOutStream() {
  if (m_filled > 0) {
    flush();
  }
  if (m_fileptr != stdout) {
    fclose(m_fileptr);
  }
  delete[] m_buffer;
}

void RawOutStream::writeByte(byte b) {
  m_buffer[m_filled++] = b;
  if (m_filled == kBufferSize) {
    fwrite(m_buffer, 1, m_filled, m_fileptr);
    m_filled = 0;
  }
}

void RawOutStream::flush() {
  if (m_filled > 0) {
    fwrite(m_buffer, 1, m_filled, m_fileptr);
    m_filled = 0;
  }
}

long int RawOutStream::getPos() {
  flush();
  return ftell(m_fileptr);
}

void RawOutStream::writeBlock(byte *begin, byte *end) {
  uint64 size = end - begin;
  if (m_filled + size < kBufferSize) {
    std::copy(begin, end, m_buffer + m_filled);
    m_filled += size;
  } else {
    flush();
    fwrite(begin, 1, size, m_fileptr);
  }
}

void RawOutStream::write48bits(uint64 to_written, long int position) {
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

uint64 RawInStream::read48bits() {
  uint64 result = 0;
  for(int i = 0; i < 6; ++i) {
    result <<= 8;
    result |= static_cast<byte>(fetchByte());
  }
  return result;
}

RawInStream::RawInStream(const std::string &file_name) :
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

RawInStream::~RawInStream() {
  if (m_fileptr != stdin) {
    fclose(m_fileptr);
  }
  delete[] m_bigbuf;
}

uint64 RawInStream::readBlock(byte* to, uint64 max_block_size) {
  assert(m_bitsInBuffer == 0);
  uint64 have_read = 0;
  int32 c;
  while (have_read < max_block_size && (c = fetchByte()) != EOF) {
    to[have_read++] = static_cast<byte>(c);
  }
  return have_read;
}

int32 RawInStream::fetchByte() {
  if (m_bigbuf_left <= 0) {
    m_bigbuf_pos = 0;
    m_bigbuf_left = fread(m_bigbuf, 1, kBigbufSize, m_fileptr);
    if (m_bigbuf_left <= 0) return EOF;
  }
  --m_bigbuf_left;
  return m_bigbuf[m_bigbuf_pos++];
}

int32 RawInStream::peekByte() {
  if (m_bigbuf_left <= 0) {
    m_bigbuf_pos = 0;
    m_bigbuf_left = fread(m_bigbuf, 1, kBigbufSize, m_fileptr);
    if (m_bigbuf_left <= 0) return EOF;
  }
  --m_bigbuf_left;
  return m_bigbuf[m_bigbuf_pos];
}

} //namespace bwtc

