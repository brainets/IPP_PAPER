#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <vector>
#include <random>

#include "nmring.h"
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// Copyright 2019 Johannes Hausmann, Nicole Voges (drjoe@free.fe, nicole.voges@gmx.com)
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// Mersenne Twister RNG.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static std::random_device rd;  // seeder.
static std::mt19937 gen(rd()); // gives rnd integer number using seed.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
delay_stepper_t::delay_stepper_t() //constructor definiton 
  : delay(0),dt(0.0),j0(0.0),j1(0.0),j2(0.0),noise(0.0),iext(0.0),vR(),vdR(),vj(),vjOther(),vOther(),vk1(),vk2(),
    vk3(),vk4(),vtmp1(),vtmp2(),vtmp3(),vpast1(),vpast2(),vpast3(),vhalf(),viext(),vvR(),vvdR(),vvR2(),unidist()
{}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool delay_stepper_t::init(size_t Size, size_t Delay, double DeltaT, double J0, double J1, double J2, double S2, double Noise, double Iext, size_t offset)
{
  if( (2>Size) || (0==Delay) || (DeltaT<1.0e-12) )
    return false;  // 1 ring only, with delay, integration time step not too small
  delay=Delay;
  dt=DeltaT;
  j0=J0;
  j1=J1;
  j2=J2;      // for coupling betw rings
  iext=Iext;  // ongoing drive to all nodes
  noise=Noise;
  //fprintf(stderr,"   noise %f \n",Size,noise);
  unidist=std::uniform_real_distribution<double>(-Noise,Noise);// fct object, can be called => provides 1 rd value [-Noise,Noise] per call
  // allocate vector buffers.
  std::vector<double>(Size,0.0).swap(vR);
  std::vector<double>(Size,0.0).swap(vdR);
  std::vector<double>(Size,0.0).swap(vj);
  std::vector<double>(Size,0.0).swap(vjOther);
  std::vector<double>(Size,0.0).swap(vk1);
  std::vector<double>(Size,0.0).swap(vk2);
  std::vector<double>(Size,0.0).swap(vk3);
  std::vector<double>(Size,0.0).swap(vk4);
  std::vector<double>(Size,0.0).swap(vtmp1);
  std::vector<double>(Size,0.0).swap(vtmp2);
  std::vector<double>(Size,0.0).swap(vtmp3);
  std::vector<double>(Size,0.0).swap(vpast1);
  std::vector<double>(Size,0.0).swap(vpast2);
  std::vector<double>(Size,0.0).swap(vpast3);
  std::vector<double>(Size,0.0).swap(vhalf);
  std::vector<double>(Size,0.0).swap(viext); //per unit external input
  // allocate ring buffers.
  std::vector<std::vector<double> >(1+Delay).swap(vvR);
  std::vector<std::vector<double> >(1+Delay).swap(vvdR);
  std::vector<std::vector<double> >(Delay).swap(vvR2);
  for(size_t i=0;i<vvR.size();++i)
    vvR[i].resize(Size,0.0);
  for(size_t i=0;i<vvdR.size();++i)
    vvdR[i].resize(Size,0.0);
  for(size_t i=0;i<vvR2.size();++i)
    vvR2[i].resize(Size,0.0);

  set_connectivity(J0,J1,J2,S2);
  
  std::uniform_real_distribution<double> unidistinit(0.0,0.1); 
  std::vector<double> vinit(Size,0.0);

  for(size_t k=0;k<Size;++k)
    vinit[k]=unidistinit(gen); //init rd value in [0,0.1], different for each k

  // initialize state and derivative ring buffers, and helper variables
  for(size_t k=0;k<Size;++k)
    {
      for(size_t j=0;j<=Delay;++j)
	{
	  vvR[j][k]=vinit[k];
	  if(j!=Delay)
	    vvR2[j][k]=vvR[j][k];  // initial value for time D before 0.
	}
      vR[k]=vpast3[k]=vvR[Delay][k];
      vhalf[k]=vvdR[Delay][k];
    }
  return true;
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void delay_stepper_t::set_connectivity(double J0, double J1, double J2, double S2)
{
  size_t Size = vj.size();
  if(Size==0)
    return;
  j0=J0;
  j1=J1;
  j2=J2;      // coupling betw rings
  double dTheta=2.0*M_PI/Size,sum=0.0; // internal connections => coupling within each ring
  for(size_t j=1;j<Size;++j)           // skip j=0: no self-connection for internal couplings!
    {
      vj[j]=j0+j1*cos(-dTheta*j);
      sum+=vj[j];
    }
  // Gauss  external
  if(0.0<S2)
    {
    size_t SizeHalf=Size/2;
    double sig2 = S2*S2;
    double gaussfac=0.5/sig2;
    double fac = J2;
    for(size_t i=0;i<=SizeHalf;++i)
      {
      double x=exp(-gaussfac*i*i);
      //fprintf(stderr,"i=%zu Size-i=%zu => %zu %zu %f\n",i,Size-i,(i+Index)%Size,((Size-i)+Index)%Size,x);
      vjOther[i]=fac*x;
      vjOther[(Size-i)%Size]=fac*x;
      }
    }
  else // reset external input connections to all zero.
    std::vector<double>(Size,0.0).swap(vjOther);

  iext=0.1*(1.0-sum); // external input consistent with homogenous stationary state R(t)=const=0.1!
  if (iext<0.0)       // no negative currents allowed!
    {
    iext = 0.0;
    fprintf(stderr,"   sum of %zu couplings: %f => iext %f\n",Size,sum,iext);
    }
} 

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void delay_stepper_t::rk4step(size_t ti, const std::vector<double>& crvRt, const std::vector<double>& crvRp1,
			      const std::vector<double>& crvRp2, const std::vector<double>& crvRp3)
{
  double t=dt*ti,half_dt=0.5*dt,sixth_dt=dt/6.0;
  // k1: calculate derivative at t,Rt with past Rt_past1
  derivative(t,crvRt,crvRp1,vk1);
  for(size_t i=0;i<vtmp1.size();++i)
    vtmp1[i]=crvRt[i]+half_dt*vk1[i];
  // k2: calculate derivative at t+dt/2,Rt+k1/2 with past Rt_past2
  derivative(t+half_dt,vtmp1,crvRp2,vk2);
  for(size_t i=0;i<vtmp2.size();++i)
    vtmp2[i]=crvRt[i]+half_dt*vk2[i];
  // k3: calculate derivative at t+dt/2,Rt+k2/2 with past Rt_past2
  derivative(t+half_dt,vtmp2,crvRp2,vk3);
  for(size_t i=0;i<vtmp3.size();++i)
    vtmp3[i]=crvRt[i]+dt*vk3[i];
  // k4: calculate derivative at t+dt,Rt+k3 with past Rt_past3
  derivative(t+dt,vtmp3,crvRp3,vk4);
  // update Rt
  for(size_t i=0;i<vR.size();++i)
    vR[i]=crvRt[i]+sixth_dt*(vk1[i]+double(2)*vk2[i]+double(2)*vk3[i]+vk4[i]);
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// step function called in main
std::vector<double> delay_stepper_t::step(size_t ti)
{
  size_t trb=ti%delay;               // ti => t_ringbuffer 
  size_t t0=1+((delay+trb-1)%delay); // some magic
  if(0!=trb)
    vvR[trb].swap(vR);
  // snapshot.
  vR=vvR[t0];    // rkVariable
  vpast1=vpast3; // new rkVariable due to delay
  vpast2=vhalf;  // new rkVariable due to delay
  vpast3=vvR[trb+1]; // new rkVariable due to delay taken from ring buffer R
  vhalf=vvR2[trb];   // new rkVariable due to delay taken from ring buffer R2 (Rhalf)
  if(1==trb) //true all 100 steps because delay=100. Only update RB0=vvR[0] when RB1=vvR[1] has been calculated
    vvR[0]=vvR.back();
  // update.
  rk4step(ti,vR,vpast1,vpast2,vpast3); //rk4 call
  vvR[trb+1]=vR;                       
  derivative(dt*ti,vR,vpast3,vdR);     //derivative call
  for(size_t i=0;i<vhalf.size();++i)
    vvR2[trb][i]=0.5*(vvR[trb][i]+vR[i])+0.125*dt*(vvdR[trb][i]-vdR[i]);
  vvdR[trb+1].swap(vdR);
  return vvR[trb+1];      //return rate vector for ti+1
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void delay_stepper_t::derivative(double t, const std::vector<double>& crvR,
				 const std::vector<double>& crvRpast,
				 std::vector<double>& rvdR)
{
  size_t n=vdR.size();  // n = total number of nodes in 1 ring
  for(size_t i=0;i<n;++i)
    {               // calculate F(iext,R_0(t-D)...R_n-1(t-D))...      
      size_t k=n-i; // periodic boundary conditions in the coupling matrix
      double fsum=0.0;
      for(size_t j=0;j<n;++j,++k)    //iext, viext = member vars, belong to stepper...
      fsum+=vj[k%n]*crvRpast[j];   // according to periodic boundary conditions in coupling matrix vj
      fsum+=iext*(1.0+unidist(gen))+viext[i]; // small noise in external input for each unit WITH external stimulus viext[i]

      // now coupling between the 2 rings with j2 
      if(!vOther.empty())
	{
	  k=n-i;
	  for(size_t j=0;j<n;++j,++k)
	    fsum += j2*vjOther[k%n]*vOther[j];  //MH-like coupling betw rings //fsum+=j2*vOther[(i+n/2)%n]; from other ring's opposite unit.
	}
    // pseudo-sigmoid: never negative!
      if(fsum<0.0)
      fsum=0.0;
      // derivative.
      rvdR[i]=fsum-crvR[i];
    }
  return;
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// generate external input representing a specific/targeted stimulus to the SE ring
//
// modified joe 20200709: used to be called with J0=-JStim and J1=+JStim. we now
//                        expect J0=JStim and J1=STD for Gauss and Gumbel. That is
//                        to say, J0=stimulus maximum, J1=stimulus sigma.
//
#warning "joe 20200709: modified generate_stimulus!"
#warning "joe 20200709: case zero looks broken, limited to a range of width 5 always!?"
std::vector<double> generate_stimulus(size_t Size, size_t Index, double J0, double J1, int mode)
{
  std::vector<double> vs(Size,0.0);
  if(Index>=Size)
    {
      fprintf(stderr,"error: invalid Index %zu >= %zu\n",Index,Size);
      return vs;
    }
  size_t SizeHalf=Size/2;
  size_t k=Size-Index; // 
  size_t range = 5;
  double sigma2= J1*J1; // joe 20200709: J1 is now stimulus sigma.
  double gaussfac=0.5/sigma2;
  double fac = J0;  // joe 20200709: J0 is now stimulus maximums => gauss normalized to max=J0 at 0.
  //fprintf(stderr,"Index inject =%zu mode=%d fac=%f sigm=%f ampl=%f\n",Index,mode,fac,sigma2,J0);
  const double mu=2.53;
  double dTheta=(2.0*M_PI)/Size,thetaCtr=dTheta*Index;
  double dGumbelBeta=(0.0<J1) ? 1.0/J1 : 1.0;
  double x=0;
  switch(mode)
    {
    case 0:
      if(Index>=range/2)
        k=Index-range/2;
      else
        k=Size-((range/2)-Index);
      for(size_t i=0;i<=range;++i,++k)
        {
          vs[k%Size]=J0*1.;
          //fprintf(stderr,"i=%zu k=%zu => mod=%zu\n",i,k,k%Size);
        }
      break;
    case 1:
      // Gauss   
      for(size_t i=0;i<=SizeHalf;++i)
	{
          x=exp(-gaussfac*i*i);
	  //fprintf(stderr,"i=%zu Size-i=%zu => %zu %zu %f\n",i,Size-i,(i+Index)%Size,((Size-i)+Index)%Size,x);
          vs[(i+Index)%Size]=fac*x;
          vs[(Size+Index-i)%Size]=fac*x;
	}
      break;
    case 2:
      // Arird 2007
      for(size_t i=0;i<Size;++i,++k)
        vs[i]=J0+J1*exp(mu*(cos(thetaCtr-dTheta*i)-1.0));
      break;
    case 3:
      // Gumbel distribution
      fac=fabs(J0)*exp(1.0)/dGumbelBeta;      
      for(size_t i=0;i<SizeHalf;++i,++k)
	{
	//x=dGumbelBeta*(J0*i); // z = (x-µ)/ß and we use µ=0 always.
	x=(0.0<J0) ? dGumbelBeta*(0.5*i) : -dGumbelBeta*(0.5*i); // z = (x-µ)/ß and we use µ=0 always.
	vs[(i+Index)%Size]=fac*dGumbelBeta*exp(-(x+exp(-x))); // pdf=exp(-(z+exp(-z)))/ß
	vs[(Size+Index-i)%Size]=fac*dGumbelBeta*exp(-(-x+exp(x)));
	}
      break;
    case 4:
      for(size_t i=Index;i<Index+SizeHalf/2;++i)
	vs[i%Size]=J0;
    default:
      fprintf(stderr,"choose a mode \n");
    }
  return vs;
}
