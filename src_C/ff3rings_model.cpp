// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// Copyright 2019 Johannes Hausmann, Nicole Voges (drjoe@free.fe, nicole.voges@gmx.com)
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <vector>
#include <random>
#include <getopt.h>

#include "nmring.h"
#include "util.h"
#include "threadpool.h"
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static std::random_device rd;
static std::mt19937 gen(rd());
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// MT tasks and runners for the simulation.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
struct no_stats_task_t
{
  size_t index,tbeg,tend,ringpos,delayff;
  double stim_noiseOFF; 
  std::vector<delay_stepper_t>* pvs1;
  std::vector<delay_stepper_t>* pvs2;
  std::vector<delay_stepper_t>* pvs3;
  std::vector<std::vector<std::vector<double> > >* pvvvRateWinR1; 
  std::vector<std::vector<std::vector<double> > >* pvvvRateWinR2;
  std::mt19937 taskgen;
};
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void run_no_stats_task(no_stats_task_t* pt)
  {
  if( (NULL==pt) || (NULL==pt->pvs1) || (NULL==pt->pvs2)|| (NULL==pt->pvs3) || (NULL==pt->pvvvRateWinR1)  || (NULL==pt->pvvvRateWinR2)|| (pt->index>=pt->pvs1->size()) )
    return;
  size_t tb=pt->tbeg,te=pt->tend,index=pt->index,ringpos=pt->ringpos,delayff=pt->delayff;
  double stim_noiseOFF=pt->stim_noiseOFF;
  std::uniform_real_distribution<float> stim_unidist(-stim_noiseOFF,stim_noiseOFF);
  delay_stepper_t& rStepperR1=(*(pt->pvs1))[index];
  delay_stepper_t& rStepperR2=(*(pt->pvs2))[index];
  delay_stepper_t& rStepperR3=(*(pt->pvs3))[index];
  std::vector<std::vector<std::vector<double> > >& rvvvRateR1=*(pt->pvvvRateWinR1);
  std::vector<std::vector<std::vector<double> > >& rvvvRateR2=*(pt->pvvvRateWinR2);
  //
  size_t windowSize=rvvvRateR1.size(),pastpos=0;
  size_t dim=rStepperR1.viext.size();
  for(size_t t=tb;t<te;++t)
    {
    if(0.0<stim_noiseOFF)
      for(size_t k=0;k<dim;++k)
        rStepperR1.viext[k]=stim_unidist(pt->taskgen); // SET current stimulus to noise.
    pastpos=(ringpos+windowSize-delayff)%windowSize; 
    rStepperR2.set_other_rate(rvvvRateR1[pastpos][index]);
    rStepperR3.set_other_rate(rvvvRateR2[pastpos][index]);
    rvvvRateR1[ringpos][index]=rStepperR1.step(t);
    rvvvRateR2[ringpos][index]=rStepperR2.step(t);
    rStepperR3.step(t);
    ringpos = (ringpos+1)%windowSize;    //update ringpos...
    }
  pt->ringpos=ringpos;
  return;
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
struct stats_task_t
{
  size_t index,tbeg,tend,ringpos,delayff;
  double stim_noiseOFF; 
  std::vector<delay_stepper_t>* pvs1;
  std::vector<delay_stepper_t>* pvs2;
  std::vector<delay_stepper_t>* pvs3;
  std::vector<std::vector<std::vector<double> > >* pvvvRateWinR1;
  std::vector<std::vector<std::vector<double> > >* pvvvStimWinR1;
  std::vector<std::vector<std::vector<double> > >* pvvvRateWinR2;
  std::vector<std::vector<std::vector<double> > >* pvvvRateWinR3;
  std::mt19937 taskgen;
};
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void run_stats_task(stats_task_t* pt)
  {
  if( (NULL==pt) || (NULL==pt->pvs1) || (NULL==pt->pvs2) || (NULL==pt->pvs3) || (NULL==pt->pvvvRateWinR1) || (NULL==pt->pvvvStimWinR1) ||
      (NULL==pt->pvvvRateWinR2) || (NULL==pt->pvvvRateWinR3) || (pt->index>=pt->pvs1->size()) )
    return;
  size_t tb=pt->tbeg,te=pt->tend,index=pt->index,ringpos=pt->ringpos,delayff=pt->delayff;
  double stim_noiseOFF=pt->stim_noiseOFF;
  std::uniform_real_distribution<float> stim_unidist(-stim_noiseOFF,stim_noiseOFF);
  delay_stepper_t& rStepperR1=(*(pt->pvs1))[index];
  delay_stepper_t& rStepperR2=(*(pt->pvs2))[index];
  delay_stepper_t& rStepperR3=(*(pt->pvs3))[index];
  std::vector<std::vector<std::vector<double> > >& rvvvRateR1=*(pt->pvvvRateWinR1);
  std::vector<std::vector<std::vector<double> > >& rvvvStimR1=*(pt->pvvvStimWinR1);
  std::vector<std::vector<std::vector<double> > >& rvvvRateR2=*(pt->pvvvRateWinR2);
  std::vector<std::vector<std::vector<double> > >& rvvvRateR3=*(pt->pvvvRateWinR3);
  //
  size_t windowSize=rvvvRateR1.size(),pastpos=0;
  size_t dim=rStepperR1.viext.size();
  for(size_t t=tb;t<te;++t)
    {
    if(0.0<stim_noiseOFF)
      for(size_t k=0;k<dim;++k)
        rStepperR1.viext[k]=stim_unidist(pt->taskgen); // SET current stimulus to noise.
    pastpos=(ringpos+windowSize-delayff)%windowSize;
    rStepperR2.set_other_rate(rvvvRateR1[pastpos][index]);
    rStepperR3.set_other_rate(rvvvRateR2[pastpos][index]);
    rvvvRateR1[ringpos][index]=rStepperR1.step(t);
    rvvvRateR2[ringpos][index]=rStepperR2.step(t);
    rvvvRateR3[ringpos][index]=rStepperR3.step(t);
    rvvvStimR1[ringpos][index]=rStepperR1.viext; // pass back the noisy stimulus actually used in this time step!
    ringpos = (ringpos+1)%windowSize;    //update ringpos...
    }
  pt->ringpos=ringpos;
  return;
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
struct stats1_task_t
{
  size_t beg,end,step,t,ringpos,delayff;
  double stim_noiseON,stim_noiseOFF;
  std::vector<delay_stepper_t>* pvs1;
  std::vector<delay_stepper_t>* pvs2;
  std::vector<delay_stepper_t>* pvs3;
  std::vector<std::vector<std::vector<double> > >* pvvvRateWinR1;
  std::vector<std::vector<std::vector<double> > >* pvvvStimWinR1;
  std::vector<std::vector<std::vector<double> > >* pvvvRateWinR2;
  std::vector<std::vector<std::vector<double> > >* pvvvRateWinR3;
  std::mt19937 taskgen;
};
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void run_stats1_task(stats1_task_t* pt)
  {
  if( (NULL==pt) || (NULL==pt->pvs1) || (NULL==pt->pvs2) || (NULL==pt->pvs3) || (NULL==pt->pvvvRateWinR1) || (NULL==pt->pvvvStimWinR1) ||
      (NULL==pt->pvvvRateWinR2) || (NULL==pt->pvvvRateWinR3) )
    return;
  size_t beg=pt->beg,end=pt->end,step=pt->step,t=pt->t,windowSize=0,pastpos=0,ringpos=pt->ringpos,delayff=pt->delayff,dim=0;
  bool bStimIsOn=(0.0<pt->stim_noiseON);
  double stim_noise=pt->stim_noiseON;
  if(!bStimIsOn)
    stim_noise=pt->stim_noiseOFF;
  std::uniform_real_distribution<float> stim_unidist(-stim_noise,stim_noise);
  std::vector<delay_stepper_t>& rvStepperR1=*(pt->pvs1);
  std::vector<delay_stepper_t>& rvStepperR2=*(pt->pvs2);
  std::vector<delay_stepper_t>& rvStepperR3=*(pt->pvs3);
  std::vector<std::vector<std::vector<double> > >& rvvvRateR1=*(pt->pvvvRateWinR1);
  std::vector<std::vector<std::vector<double> > >& rvvvStimR1=*(pt->pvvvStimWinR1);
  std::vector<std::vector<std::vector<double> > >& rvvvRateR2=*(pt->pvvvRateWinR2);
  std::vector<std::vector<std::vector<double> > >& rvvvRateR3=*(pt->pvvvRateWinR3);
  //
  windowSize=pt->pvvvRateWinR1->size();
  pastpos=(ringpos+windowSize-delayff)%windowSize;
  for(size_t j=beg;j<end;j+=step)
    {
    dim=rvStepperR1[j].viext.size();
    if(0.0<stim_noise)
      {
      if(bStimIsOn)
        for(size_t k=0;k<dim;++k)
          rvStepperR1[j].viext[k]+=stim_unidist(pt->taskgen); // ADD noise to current stimulus.
      else
        for(size_t k=0;k<dim;++k)
          rvStepperR1[j].viext[k]=stim_unidist(pt->taskgen);  // SET noise as current stimulus.
      }
    rvStepperR2[j].set_other_rate(rvvvRateR1[pastpos][j]);
    rvStepperR3[j].set_other_rate(rvvvRateR2[pastpos][j]);
    rvvvRateR1[ringpos][j]=rvStepperR1[j].step(t);
    rvvvRateR2[ringpos][j]=rvStepperR2[j].step(t);
    rvvvRateR3[ringpos][j]=rvStepperR3[j].step(t);
    rvvvStimR1[ringpos][j]=rvStepperR1[j].viext; // pass back the noisy stimulus actually used in this time step!
    }
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<std::vector<double> > time_compress(std::vector<std::vector<double> >& crvv, size_t nmax)
  {
  std::vector<std::vector<double> > vv; // compressed_times x nodes
  if((crvv.size()%nmax)!=0)
    return vv;
  vv.resize(crvv.size()/nmax);
  double nrm=1.0/nmax;
  for(size_t j=0;j<crvv.size();)
    {
    std::vector<double> v(crvv[j].size(),0.0);
    for(size_t n=0;n<nmax;++n,++j)
      {
      for(size_t k=0;k<v.size();++k)
	v[k]+=crvv[j][k];
      }
    for(size_t k=0;k<v.size();++k)
      v[k]*=nrm;
    vv[(j-nmax)/nmax].swap(v);
    }
  return vv;
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// main driver program.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int main(int argc, char* argv[])
  {
  size_t nsteps=8500, tis=4000, tstimON=5000, tstimOFF=6500; // the real thing: final timings, determined by SB
  //size_t nsteps=850, tis=401, tstimON=500, tstimOFF=650;  // minitest of functionality.
  double JStim=2.0, SStim=8.0; //JStim=stimAmpl default 2 //SStim=stimWidth default 8               
  int MStim=1;                                            // stimulus mode, default=1=>gauss
  
  double stim_noiseON=0.0*JStim, stim_noiseOFF=0.0;

  size_t dim=100, delay=10, delayFF=20; //delayFF new for coupling betw 2 rings
  double dt=0.01;
  double Noise=0.50, Iext=0.0; //defaults: add 50% noise, gaussian stim with sigm2=64
  size_t npos = 4, posStep = dim/npos; 
  
  std::string sBaseDir="3FFres";
  make_dir(sBaseDir);
  
#if 1
  double J0 = -25.0,J1 = 11.0;   // SB
  bool bIsSU = false;
#else
  double J0 = -30.0,J1 = -8.0;  // SU
  bool bIsSU = true;
#endif
 
  double J2=0.35, S2=3.0; //feed forward (ff) coupling betw rings
  
  size_t trials = 100; //40000
  size_t windowSize= 400; 
  //
  // MT stuff: initialize a thread pool with nThreads workers.
  //
  size_t nThreads=std::thread::hardware_concurrency(); //NV schleppi: 4
  thread_pool_t mtpool(nThreads);
  //
  // init parallel steppers with different initial conditions.
  std::vector<delay_stepper_t> vstepperR1(trials);
  std::vector<delay_stepper_t> vstepperR2(trials);
  std::vector<delay_stepper_t> vstepperR3(trials);
  for(size_t i=0;i<trials;++i)
    {
    bool bOk=vstepperR1[i].init(dim,delay,dt,J0/dim,J1/dim,0.0,0.0,Noise,Iext,0);//bottom ring, receives only stimulus
    if(bOk)
      ;
    else
      fprintf(stderr,"stepperR1 %zu: FAILED!\n",i);
    bOk=vstepperR2[i].init(dim,delay,dt,J0/dim,J1/dim,J2,S2,Noise,Iext,0);//2nd ring: receives FF input from R1
    if(bOk)
      ;
    else
      fprintf(stderr,"stepperR2 %zu: FAILED!\n",i);
    bOk=vstepperR3[i].init(dim,delay,dt,J0/dim,J1/dim,J2,S2,Noise,Iext,0);//3nd ring: receives FF input from R2
    if(bOk)
      ;
    else
      fprintf(stderr,"stepperR3 %zu: FAILED!\n",i);
    }
  // init stimuli.
  std::vector<double> vnull(dim,0.0);
  std::vector<std::vector<double> > vvStimBase(npos); // npos x nodes
  for(size_t i=0;i<npos;++i)    
    vvStimBase[i]=generate_stimulus(dim,i*posStep,JStim,SStim,MStim); 
  fprintf(stderr,"generate stim end\n");
  // will save average firing rates to file for t=tis...nsteps, only every 10th step.
  size_t saveStep = 10;
  std::vector<std::vector<std::vector<double> > > vvvRateAvgR1, vvvRateAvgR2, vvvRateAvgR3;
  vvvRateAvgR1.reserve(nsteps-tis);
  vvvRateAvgR2.reserve(nsteps-tis);
  vvvRateAvgR3.reserve(nsteps-tis);
  // stimulus injection matrix init to zero => just once write out as .npy, only every 10th value according to avFR
  std::vector<std::vector<std::vector<double> > > vvStimSave(npos);
  for(size_t np=0;np<npos;++np)
    {
    vvStimSave[np].resize((nsteps-tis)/saveStep);// time x nodes
    reset(vvStimSave[np],vnull);
    }
  //ring buffer.
  std::vector<std::vector<std::vector<double> > > vvvRateWinR1,vvvRateWinR2,vvvRateWinR3; //ring buffers for past rates: windowSize x trials x nodes
  std::vector<std::vector<std::vector<double> > > vvvStimWin; //ring buffer for stim: windowSize x trials x nodes!
  //
  // initialize struct to handle average and variance calculations.
  //
  double avgSpaceR1=0.0,varSpaceR1=0.0,avgSpaceR2=0.0,varSpaceR2=0.0,avgSpaceR3=0.0,varSpaceR3=0.0;
  std::vector<double> vZeroN(dim,0.0); //to reset
  node_avg_t avgFR1,avgFR2,avgFR3;
  avgFR1.init(trials);
  avgFR1.reset(vZeroN);
  avgFR2.init(trials);
  avgFR2.reset(vZeroN);
  avgFR3.init(trials);
  avgFR3.reset(vZeroN);
  std::vector<double> vAvRtR1(dim,0.0);
  std::vector<double> vAvRt2R1(dim,0.0);
  std::vector<double> vAvRtR2(dim,0.0);
  std::vector<double> vAvRt2R2(dim,0.0);
  std::vector<double> vAvRtR3(dim,0.0);
  std::vector<double> vAvRt2R3(dim,0.0);
  
  size_t ringpos=0,ringcur=0; // ringpos: oldest vvRate in ring buffer, initially at zero
  //
  // stimulus intervals. simple case with just one time interval here.
  //
  std::vector<stim_interval_t> vsi(1); // 1 periods with stim
  stim_interval_t firstStim(tstimON,tstimOFF); // first stim defined
  vsi[0]=firstStim;
  size_t nxtstim=0; // stimulus interval index
  //
  // MT tasks.
  //
  vvvRateWinR1.resize(windowSize);
  vvvRateWinR2.resize(windowSize);
  vvvRateWinR3.resize(windowSize);
  vvvStimWin.resize(windowSize);
  for(size_t i=0;i<vvvRateWinR1.size();++i)
    {
    vvvRateWinR1[i].resize(trials);
    vvvRateWinR2[i].resize(trials);
    vvvRateWinR3[i].resize(trials);
    vvvStimWin[i].resize(trials);
    for(size_t j=0;j<vvvRateWinR1[i].size();++j)
      {
      vvvRateWinR1[i][j].resize(dim,0.0);
      vvvRateWinR2[i][j].resize(dim,0.0);
      vvvRateWinR3[i][j].resize(dim,0.0);
      vvvStimWin[i][j].resize(dim,0.0);
      }
    }
  no_stats_task_t transientTask;
  transientTask.tbeg=0;
  transientTask.tend=tis-windowSize;
  transientTask.ringpos=0;
  transientTask.delayff=delayFF;
  transientTask.stim_noiseOFF=stim_noiseOFF;
  transientTask.pvs1=&vstepperR1;
  transientTask.pvs2=&vstepperR2;
  transientTask.pvs3=&vstepperR3;
  transientTask.pvvvRateWinR1=&vvvRateWinR1;
  transientTask.pvvvRateWinR2=&vvvRateWinR2; 
  std::vector<no_stats_task_t>  vTransientTask(trials,transientTask);
  stats_task_t runTask;
  runTask.tbeg=transientTask.tend;
  runTask.tend=tis;
  runTask.ringpos=0;
  runTask.delayff=delayFF;
  runTask.stim_noiseOFF=stim_noiseOFF;
  runTask.pvs1=&vstepperR1;
  runTask.pvs2=&vstepperR2;
  runTask.pvs3=&vstepperR3;
  runTask.pvvvRateWinR1=&vvvRateWinR1;
  runTask.pvvvStimWinR1=&vvvStimWin;
  runTask.pvvvRateWinR2=&vvvRateWinR2;
  runTask.pvvvRateWinR3=&vvvRateWinR3;
  std::vector<stats_task_t>  vRunTask(trials,runTask);
  stats1_task_t runTask1;       // real task: just integrate one time step (ringbuffers are NOT updated, this is done ihn the main loop below!).
  runTask1.end=trials;
  runTask1.step=nThreads;
  runTask1.t=tis;
  runTask1.delayff=delayFF;
  runTask1.ringpos=0;
  runTask1.stim_noiseOFF=stim_noiseOFF;
  runTask1.stim_noiseON=0.0;
  runTask1.pvs1=&vstepperR1;
  runTask1.pvs2=&vstepperR2;
  runTask1.pvs3=&vstepperR3;
  runTask1.pvvvRateWinR1=&vvvRateWinR1;
  runTask1.pvvvStimWinR1=&vvvStimWin;
  runTask1.pvvvRateWinR2=&vvvRateWinR2;
  runTask1.pvvvRateWinR3=&vvvRateWinR3;
  std::vector<stats1_task_t> vRunTask1(nThreads,runTask1);
  for(size_t i=0;i<vRunTask1.size();++i)
    {
    vRunTask1[i].beg=i;
    vRunTask[i].taskgen.seed(rd());
    }
  START(integrate); // timer macros to print out run time: START(name);do something; STOP(name);
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  // run transient task.
  //
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  START(transient);
  for(size_t i=0;i<vTransientTask.size();++i)
    { // transient task: run from 0...tis-windowSize fill only R1,R2 steppers' ringbuffer.
    vTransientTask[i].index=i;
    vTransientTask[i].taskgen.seed(rd());
    mtpool.run_task(std::bind(run_no_stats_task,&vTransientTask[i]));
    }
  mtpool.wait();
  for(size_t i=0;i<vRunTask.size();++i)
    { // transient task: run from tis-windowSize...tis fill all the ringbuffers (R1, R2, R3 and Stim).
    vRunTask[i].ringpos=vTransientTask[i].ringpos;    
    vRunTask[i].index=i;
    vRunTask[i].taskgen.seed(rd());
    mtpool.run_task(std::bind(run_stats_task,&vRunTask[i]));
    }
  mtpool.wait();
  STOPN(transient,tis);
  ringpos=vRunTask[0].ringpos%windowSize;
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  // main integration and evaluation loop.
  //
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  bool bStimIsOn=false,bStoppingStim=false;
  size_t turnOffUpdateStimRingBuffer=0;
  for(size_t ti=tis;ti<=nsteps;++ti)
    {
    //if(ti==nsteps)
    //break; // trick to properly handle average rates above at the end of the main loop.
    if(ti<nsteps)
      {
      //
      // switch stimulus according to time intervals defined in vsi.
      //
      if( (nxtstim<vsi.size()) && (ti==vsi[nxtstim].beg) )      // beg => ON
	{
	avgSpaceR1=0.0,varSpaceR1=0.0, avgSpaceR2=0.0,varSpaceR2=0.0, avgSpaceR3=0.0,varSpaceR3=0.0; // switching => evaluate space avg and var.
	avgFR1.get_avg_var_space(avgSpaceR1,varSpaceR1);
	avgFR2.get_avg_var_space(avgSpaceR2,varSpaceR2);
	avgFR3.get_avg_var_space(avgSpaceR3,varSpaceR3);
	fprintf(stderr,"%zu: avgFRspace=%f +/- %f (cnt %zu), %f +/- %f, %f +/- %f \n",ti,
		avgSpaceR1,sqrt(varSpaceR1),avgFR1.m_count,avgSpaceR2,sqrt(varSpaceR2),avgSpaceR3,sqrt(varSpaceR3));
	avgFR1.reset(vZeroN);
	avgFR2.reset(vZeroN);
	avgFR3.reset(vZeroN);
	// do stimulus ring buffer updates now. also remember when to stop the updates.
	for(size_t i=0;i<vRunTask1.size();++i)
	  {
	  vRunTask1[i].stim_noiseOFF=0.0;
	  vRunTask1[i].stim_noiseON=stim_noiseON; 
	  }
	bStimIsOn=true;
	bStoppingStim=false;
	turnOffUpdateStimRingBuffer=vsi[nxtstim].end+windowSize;
	fprintf(stderr,"bStimIsOn: %s\nbStoppingStimRingBuffer: %s turn off at %zu.\n",
		bStimIsOn?"true":"false",bStoppingStim?"true":"false",
		turnOffUpdateStimRingBuffer);
	}
      else if( (nxtstim<vsi.size()) && (ti==vsi[nxtstim].end) ) // end => OFF
	{
	//
	avgSpaceR1=0.0,varSpaceR1=0.0, avgSpaceR2=0.0,varSpaceR2=0.0, avgSpaceR3=0.0,varSpaceR3=0.0; // switching => evaluate space avg and var.
	avgFR1.get_avg_var_space(avgSpaceR1,varSpaceR1);
	avgFR2.get_avg_var_space(avgSpaceR2,varSpaceR2);
	avgFR3.get_avg_var_space(avgSpaceR3,varSpaceR3);
	fprintf(stderr,"%zu: avgFRspace=%f +/- %f (cnt %zu), %f +/- %f, %f +/- %f \n",ti,
		avgSpaceR1,sqrt(varSpaceR1),avgFR1.m_count,avgSpaceR2,sqrt(varSpaceR2),avgSpaceR3,sqrt(varSpaceR3));
	avgFR1.reset(vZeroN);
	avgFR2.reset(vZeroN);
	avgFR3.reset(vZeroN);
	//
	for(size_t i=0;i<vRunTask1.size();++i)
	  {
	  vRunTask1[i].stim_noiseON=0.0;
	  vRunTask1[i].stim_noiseOFF=stim_noiseOFF; 
	  }
	for(size_t i=0;i<vstepperR1.size();++i)
	  vstepperR1[i].viext = vnull;
	++nxtstim;
	bStimIsOn=false;
	bStoppingStim=true;
	fprintf(stderr,"bStimIsOn: %s\nbStoppingStim: %s turn off at %zu.\n",
		bStimIsOn?"true":"false",bStoppingStim?"true":"false",
		turnOffUpdateStimRingBuffer);
	}
      if(bStoppingStim && (turnOffUpdateStimRingBuffer<=ti) )
	{ // we need to update the stimulus ring buffer windowSize steps after turning off the stimulus...
	bStoppingStim=false;
	fprintf(stderr,"bStimIsOn: %s\nbStoppingStim: %s turn off at %zu.\n",
		bStimIsOn?"true":"false",bStoppingStim?"true":"false",
		turnOffUpdateStimRingBuffer);
	}
      //
      // add noise to stimulus?
      //
      if(bStimIsOn) // true during stimulation
	for(size_t i=0;i<vstepperR1.size();++i)
	  vstepperR1[i].viext = vvStimBase[i%npos];
      //
      // integrate one step, so now vvRate contains the firing rates for time ti.
      //
      for(size_t i=0;i<vRunTask1.size();++i)
	{
	vRunTask1[i].t=ti;
	vRunTask1[i].ringpos=ringpos;
	mtpool.run_task(std::bind(run_stats1_task,&vRunTask1[i]));
	}
      mtpool.wait();
      //
      // update the ring buffers...
      //
      ringcur=ringpos;                       // index of current rates in ring buffer.
      ringpos = (ringpos+1)%windowSize;      // update index ringpos so that it again points to the oldest rates.
      //
      // ...then update averages and variances with the new rates...
      //
      avgFR1.update(vvvRateWinR1[ringcur]);
      avgFR2.update(vvvRateWinR2[ringcur]);
      avgFR3.update(vvvRateWinR3[ringcur]);
      }
    //
    // preparing to write out FR and Stim
    if( (ti>tis) && (0==(ti%saveStep)) )
      { // save FRs to write out data, happens every windowSize time steps ===> should happen every 10th time step
      vvvRateAvgR1.push_back(get_ringbuffer_avg(vvvRateWinR1,ringcur,saveStep));
      vvvRateAvgR2.push_back(get_ringbuffer_avg(vvvRateWinR2,ringcur,saveStep));
      vvvRateAvgR3.push_back(get_ringbuffer_avg(vvvRateWinR3,ringcur,saveStep));	
      fprintf(stderr,"FR out ti=%zu\n",ti);
      }
    // to save stimulus injection matrix to write out just once!
    if( (ti>=tstimON) && (ti<tstimOFF) )
      {
      if( (ti-tis)%saveStep == 0)
	{
	for(size_t np=0;np<npos;++np)
	  vvStimSave[np][(ti-tis)/saveStep]=vvvStimWin[ringcur][np];
	}
      }
    }
  //
  // end of main integration loop, dump run time on stderr.
  //
  STOPN(integrate,(nsteps-tis)); // timer macros to print out run time: START(name);do something; STOP(name);  
  avgSpaceR1=0.0,varSpaceR1=0.0, avgSpaceR2=0.0,varSpaceR2=0.0, avgSpaceR3=0.0,varSpaceR3=0.0; // end of integration => evaluate space avg and var.
  avgFR1.get_avg_var_space(avgSpaceR1,varSpaceR1);
  avgFR2.get_avg_var_space(avgSpaceR2,varSpaceR2);
  avgFR3.get_avg_var_space(avgSpaceR3,varSpaceR3);
  fprintf(stderr,"%zu: avgFRspace=%f +/- %f (cnt %zu), %f +/- %f, %f +/- %f \n",nsteps,
	  avgSpaceR1,sqrt(varSpaceR1),avgFR1.m_count,avgSpaceR2,sqrt(varSpaceR2),avgSpaceR3,sqrt(varSpaceR3));
  avgFR1.reset(vZeroN); // well...are we ever going to reuse it?
  avgFR2.reset(vZeroN);
  avgFR3.reset(vZeroN);
  // write out numpy files: FR and Stim: save FR all trials, all nodes, time averaged with deltaT=10:
  char bufFR1[1024];
  char bufFR2[1024];
  char bufFR3[1024];
  snprintf(bufFR1,1024,"%s/%s_frR1avg10_ff3rings_trials%zu_ws%zu_npos%zu.npy",sBaseDir.c_str(),bIsSU?"su":"sb",trials,windowSize,npos);
  snprintf(bufFR2,1024,"%s/%s_frR2avg10_ff3rings_trials%zu_ws%zu_npos%zu.npy",sBaseDir.c_str(),bIsSU?"su":"sb",trials,windowSize,npos);
  snprintf(bufFR3,1024,"%s/%s_frR3avg10_ff3rings_trials%zu_ws%zu_npos%zu.npy",sBaseDir.c_str(),bIsSU?"su":"sb",trials,windowSize,npos);
  fprintf(stderr,"size1=%zu  %zu  %zu\n",vvvRateAvgR1.size(),vvvRateAvgR1[0].size(),vvvRateAvgR1[0][0].size()); // time, trials, nodes
  npy_save(bufFR1,vvvRateAvgR1); 
  npy_save(bufFR2,vvvRateAvgR2);
  npy_save(bufFR3,vvvRateAvgR3);
  fprintf(stderr,"saved FRs\n");
  // just once if no stimNoise...
  char bufFRstim[1024];
  snprintf(bufFRstim,1024,"%s/%s_stim-fravg10_ff3rings_trials%zu_ws%zu_npos%zu.npy",sBaseDir.c_str(),bIsSU?"su":"sb",trials,windowSize,npos);
  npy_save(bufFRstim,vvStimSave);
  
  fprintf(stderr,"saved FR1 size= %zu %zu %zu\n",vvvRateAvgR1.size(),vvvRateAvgR1[0].size(),vvvRateAvgR1[0][0].size()); // time, trials, nodes
  fprintf(stderr,"saved FR2 size= %zu %zu %zu\n",vvvRateAvgR2.size(),vvvRateAvgR2[0].size(),vvvRateAvgR2[0][0].size()); // time, trials, nodes
  fprintf(stderr,"saved FR3 size= %zu %zu %zu\n",vvvRateAvgR3.size(),vvvRateAvgR3[0].size(),vvvRateAvgR3[0][0].size()); // time, trials, nodes
  fprintf(stderr,"saved Stim size= %zu %zu %zu\n",vvStimSave.size(),vvStimSave[0].size(),vvStimSave[0][0].size()); // npos, time
  
  return 0;
  }
