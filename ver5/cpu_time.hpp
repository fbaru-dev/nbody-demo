/*
    This file is part of the example codes which have been used
    for the "Code Optmization Workshop".
    
    Copyright (C) 2016  Fabio Baruffa <fbaru-dev@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _CPUTIME_HPP
#define _CPUTIME_HPP

#include <sys/time.h>
#include <sys/types.h>
#include <sys/resource.h>

// Return number of microseconds since 1.1.1970, in a 64 bit integer.

class CPUTime {
private:
    double wctime;
    
    inline double readTime() 
    {
      struct timeval tp;

      gettimeofday(&tp,NULL);
      wctime = (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6;
      return wctime;
    }
public:
    CPUTime() : wctime(0.0) { }
        
    inline double start() { return readTime(); }
    inline double stop()  { return readTime(); }
    
};

#endif
