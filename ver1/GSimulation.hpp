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

#ifndef _GSIMULATION_HPP
#define _GSIMULATION_HPP

#include <random>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>

#include <omp.h>

#include "Particle.hpp"

class GSimulation 
{
public:
  GSimulation();
  ~GSimulation();
  
  void init();
  void set_number_of_particles(int N);
  void set_number_of_steps(int N);
  void start();
  
private:
  Particle *particles;
  
  int       _npart;		//number of particles
  int	    _nsteps;	//number of integration steps
  real_type _tstep;		//time step of the simulation

  int	    _sfreq;		//sample frequency
  
  real_type _kenergy;	//kinetic energy
  
  double _totTime;		//total time of the simulation
  double _totFlops;		//total number of flops 
   
  void init_pos();	
  void init_vel();
  void init_acc();
  void init_mass();
    
  inline void set_npart(const int &N){ _npart = N; }
  inline int get_npart() const {return _npart; }
  
  inline void set_tstep(const real_type &dt){ _tstep = dt; }
  inline real_type get_tstep() const {return _tstep; }
  
  inline void set_nsteps(const int &n){ _nsteps = n; }
  inline int get_nsteps() const {return _nsteps; }
  
  inline void set_sfreq(const int &sf){ _sfreq = sf; }
  inline int get_sfreq() const {return _sfreq; }
  
  void print_header();
  
};

#endif
