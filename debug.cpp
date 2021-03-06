/* debug
 * print debug messages only if cmake build type "Debug" is set in CMakeLists.txt
 *
 * Copyright (C) 2017 Jan Felix Schmidt <janschmidt@mailbox.org>
 *   
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "debug.hpp"

void polematrix::debug(const std::string& funcname, const std::string& msg)
{
  #ifndef NDEBUG
  std::cout << "DEBUG[" << funcname << "]: " << msg << std::endl;
  #else
  auto a=funcname; auto b=msg; // avoid "unused parameter" warnings
  #endif
}
