/**
 * @file Streams.cpp
 * @author Pekka Mikkola <pjmikkol@cs.helsinki.fi>
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

InStream::InStream(std::string file_name) :
    m_name(file_name), m_from(NULL), m_infile(NULL)
{
  if (m_name != "") {
    m_infile = new std::ifstream(m_name.c_str());
    m_from = dynamic_cast<std::istream*>(m_infile);
  } else
    m_from = &std::cin;
  assert(m_from);
}

InStream::~InStream() {
  if (m_infile) {
    m_infile->close();
    delete m_infile;
  }
}

std::streamsize InStream::readBlock(byte* to, std::streamsize max_block_size) {
  if(! *m_from) return 0; /* stream has reached EOF */
  m_from->read(reinterpret_cast<char*>(to), max_block_size);
  return m_from->gcount();
}

} //namespace bwtc
