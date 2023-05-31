// -*- c++ -*-
#ifndef __NM_RING_H
#define __NM_RING_H

#include <vector>
#include <random>
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// Copyright 2019 Johannes Hausmann, Nicole Voges (drjoe@free.fe, nicole.voges@gmx.com)
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// stepper class with delay (all members public => struct)
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
struct delay_stepper_t
{ // class with all members public = struct ;)
  // ctor/dtor/init.
  delay_stepper_t(); 
  ~delay_stepper_t() {}
  bool init(size_t Size, size_t Delay, double DeltaT, double J0, double J1, double J2, double S2, double Noise, double Iext, size_t offset);
  void set_connectivity(double J0, double J1, double J2, double S2);
  // main stepping routine, call in main for each time step.
  std::vector<double> step(size_t ti);
  // getter for rate vector and past rate vectors
  const std::vector<double>& get_rate(void) const { return vR; }
  const std::vector<double>& get_past1_rate(void) const { return vpast1; }
  const std::vector<double>& get_past3_rate(void) const { return vpast3; }
  // setter for rate vector => needed for coupling two rings.
  void set_other_rate(const std::vector<double>& crvOther) {  vOther = crvOther; }
  // internals: Runge-Kutta for DDE and derivative methods.
  void rk4step(size_t ti, const std::vector<double>& crvRt, const std::vector<double>& crvRp1,
	       const std::vector<double>& crvRp2, const std::vector<double>& crvRp3);
  void derivative(double t, const std::vector<double>& crvR, const std::vector<double>& crvRpast, std::vector<double>& rvdR);
  //
  size_t delay;
  double dt,j0,j1,j2,noise,iext,iext0;    
  std::vector<double> vR,vdR,vj,vjOther,vOther;          // rateVector, derivative of rateVector, coupling consts (i->j)
  std::vector<double> vk1,vk2,vk3,vk4,vtmp1,vtmp2,vtmp3; // standard RungeKutta variables
  std::vector<double> vpast1,vpast2,vpast3,vhalf,viext;  // new RK variables required for integration with delays & extInput
  std::vector<std::vector<double> > vvR,vvdR,vvR2;       // ring buffers for R, dR/dt, and the middle points R2.
  std::uniform_real_distribution<double> unidist;        // will iwann iwo Glvert verwenden, name unidist
};
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
struct stim_interval_t
{ // class with all members public = struct ;)
  size_t beg,end,pos; 
  // constructor, initializes member vars => defines beg and end of stimulus
  // inline => declaration & definition (inside struct) 
  stim_interval_t(size_t Beg=0, size_t End=0, size_t Pos=0) : beg(Beg),end(End),pos(Pos) {} 
};
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
struct multi_stim_interval_t
{ // class with all members public = struct ;)
  size_t beg,end;
  std::vector<size_t> vpos;
  // constructor, initializes member vars => defines beg and end of stimulus
  // inline => declaration & definition (inside struct) 
  multi_stim_interval_t(size_t Beg=0, size_t End=0, const std::vector<size_t>& vPos=std::vector<size_t>()) : beg(Beg),end(End),vpos(vPos) {} 
};
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<double> generate_stimulus(size_t Size, size_t Index, double J0, double J1, int mode=0);

#endif
