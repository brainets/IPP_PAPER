// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// Copyright 2019 Johannes Hausmann, Nicole Voges (drjoe@free.fr, nicole.voges@gmx.com)
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
  std::vector<std::vector<std::vector<double> > >* pvvvRateWinR1; 
  std::vector<std::vector<std::vector<double> > >* pvvvRateWinR2;
  std::mt19937 taskgen;
};
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void run_no_stats_task(no_stats_task_t* pt)
  {
  if( (NULL==pt) || (NULL==pt->pvs1) || (NULL==pt->pvs2) || (NULL==pt->pvvvRateWinR1)  || (NULL==pt->pvvvRateWinR2)|| (pt->index>=pt->pvs1->size()) )
    return;
  size_t tb=pt->tbeg,te=pt->tend,index=pt->index,ringpos=pt->ringpos,delayff=pt->delayff;// windowSize=0;//new ws=0!!
  double stim_noiseOFF=pt->stim_noiseOFF;
  std::uniform_real_distribution<float> stim_unidist(-stim_noiseOFF,stim_noiseOFF);
  delay_stepper_t& rStepperR1=(*(pt->pvs1))[index];
  delay_stepper_t& rStepperR2=(*(pt->pvs2))[index];
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
    rStepperR1.set_other_rate(rvvvRateR2[pastpos][index]);
    rvvvRateR1[ringpos][index]=rStepperR1.step(t);
    rvvvRateR2[ringpos][index]=rStepperR2.step(t);
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
  std::vector<std::vector<std::vector<double> > >* pvvvRateWinR1;
  std::vector<std::vector<std::vector<double> > >* pvvvStimWinR1;
  std::vector<std::vector<std::vector<double> > >* pvvvRateWinR2;
  std::mt19937 taskgen;
};
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void run_stats1_task(stats1_task_t* pt)
  {
  if( (NULL==pt) || (NULL==pt->pvs1) || (NULL==pt->pvs2) || (NULL==pt->pvvvRateWinR1) || (NULL==pt->pvvvStimWinR1) ||
      (NULL==pt->pvvvRateWinR2) )
    return;
  size_t beg=pt->beg,end=pt->end,step=pt->step,t=pt->t,windowSize=0,pastpos=0,ringpos=pt->ringpos,delayff=pt->delayff,dim=0;
  bool bStimIsOn=(0.0<pt->stim_noiseON);
  double stim_noise=pt->stim_noiseON;
  if(!bStimIsOn)
    stim_noise=pt->stim_noiseOFF;
  std::uniform_real_distribution<float> stim_unidist(-stim_noise,stim_noise);
  std::vector<delay_stepper_t>& rvStepperR1=*(pt->pvs1);
  std::vector<delay_stepper_t>& rvStepperR2=*(pt->pvs2);
  std::vector<std::vector<std::vector<double> > >& rvvvRateR1=*(pt->pvvvRateWinR1);
  std::vector<std::vector<std::vector<double> > >& rvvvStimR1=*(pt->pvvvStimWinR1);
  std::vector<std::vector<std::vector<double> > >& rvvvRateR2=*(pt->pvvvRateWinR2);
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
    rvStepperR1[j].set_other_rate(rvvvRateR2[pastpos][j]);
    rvvvRateR1[ringpos][j]=rvStepperR1[j].step(t);
    rvvvRateR2[ringpos][j]=rvStepperR2[j].step(t);
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
  bool bIsAt = true;//false; // true or false representing attend-IN and attend-OUT
  size_t trials, nsteps, tis, tswitch, tstimON, tstimOFF, tstimCue;
  if(bIsAt == true)
    trials=2000, nsteps=9500, tis=3400, tswitch=4400, tstimON=4500, tstimOFF=6500, tstimCue=8000; // the real thing: attend IN switch WM to SB shortly before ON
  //trials=200, nsteps=950, tis=340, tswitch=440, tstimON=450, tstimOFF=650, tstimCue=800; // minitest of functionality
  else 
    trials=2000, nsteps=9500, tis=3400, tswitch=9900, tstimON=4500, tstimOFF=6500, tstimCue=8000; // the real thing: attend OUT no switch
  //trials=200, nsteps=950, tis=340, tswitch=990, tstimON=450, tstimOFF=650, tstimCue=800; // minitest of functionality
  //
  double JStim=2.0,SStim=8.0;                                 // stimAmplitude and STD
  int MStim=1;                                                // stimulus mode, default=1=>gauss
  double stim_noiseON=0.0*JStim, stim_noiseOFF=0.0;

  size_t dim=100, delay=10, delayFF=20; //delayFF new for coupling betw 2 rings, MUST have delayff<windowSize
  double dt=0.01;
  double Noise=0.50, Iext=0.0; //defaults: add 50% noise, gaussian stim with sigm2=64 centered at node pos
  double J0sb = -25.0,J1sb = 11.0; 
  double J0su = -30.0,J1su = -8.0; 
  double J2SEfromWM = 0.23, S2SEfromWM = 6.0; 
  double J2WMfromSE = 0.15, S2WMfromSE = 3.0;
  size_t npos = 1;
  size_t cpos = 10;
  trials=100;
  size_t windowSize=400; //MUST have windowSize>delayff
  std::vector<size_t> vPosOne(npos,0); //first injection positions
  vPosOne[0] = 50; 
  std::vector<size_t> vPosTwo(cpos,0); //second injection positions
  for(size_t i=0;i<cpos;++i)
    vPosTwo[i] = 10*i;
  std::string sBaseDir="2FBres";
  std::string s=make_dir(sBaseDir);
  fprintf(stderr,"storing result in base directory <%s>...\n",s.c_str());
  fprintf(stderr,"Parameters: ws=%zu npos=%zu cpos=%zu trials=%zu Attstate=%s\n",windowSize,npos,cpos,trials,bIsAt?"YES":"NO");
  //
  // MT stuff: initialize a thread pool with nThreads workers.
  //
  size_t nThreads=std::thread::hardware_concurrency(); //NV schleppi: 4
  thread_pool_t mtpool(nThreads);
  //
  // init parallel steppers with different initial conditions.
  std::vector<delay_stepper_t> vstepperR1(trials);
  std::vector<delay_stepper_t> vstepperR2(trials);
  for(size_t i=0;i<trials;++i)
    {
    bool bOk=vstepperR1[i].init(dim,delay,dt,J0su/dim,J1su/dim,J2SEfromWM,S2SEfromWM,Noise,Iext,0); //bottom SEring: receives stim and FBinput from R2=WM
    if(bOk)
      ;
    else
      fprintf(stderr,"stepperR1 %zu: FAILED!\n",i);
    bOk=vstepperR2[i].init(dim,delay,dt,J0su/dim,J1su/dim,J2WMfromSE,S2WMfromSE,Noise,Iext,0);     //top WMring: receives FFinput from R1=SE
    if(bOk)
      ;
    else
      fprintf(stderr,"stepperR2 %zu: FAILED!\n",i);
    }
  // init stimuli.
  std::vector<double> vnull(dim,0.0);
  std::vector<std::vector<double> > vvStimBase(dim); // nodes x nodes => ALL possible stimuli (positions)
  for(size_t i=0;i<dim;++i)   
    vvStimBase[i]=generate_stimulus(dim,i,JStim,SStim,MStim);
  fprintf(stderr,"generate stim end\n");
  // will save average firing rates to file for t=tis...nsteps, only every 10th step.
  size_t saveStep = 10;
  std::vector<std::vector<std::vector<double> > > vvvRateAvgR1, vvvRateAvgR2;
  vvvRateAvgR1.reserve(nsteps-tis);
  vvvRateAvgR2.reserve(nsteps-tis);
  std::vector<std::vector<std::vector<double> > > vvStimSave(npos);// stimInjMatrix to zero=> save evry 10th step in .npy
  for(size_t np=0;np<npos;++np)
    {
    vvStimSave[np].resize((nsteps-tis)/saveStep);// time x nodes
    reset(vvStimSave[np],vnull);
    }
  std::vector<std::vector<std::vector<double> > > vvStimSaveCue(cpos);// stimInjMatrix to zero=> save evry 10th step in .npy
  for(size_t np=0;np<cpos;++np)
    {
    vvStimSaveCue[np].resize((nsteps-tis)/saveStep);// time x nodes
    reset(vvStimSaveCue[np],vnull);
    }
  //ring buffer. vvvRsteWinR1,R2 also to be used for global binning initialization. 
  std::vector<std::vector<std::vector<double> > > vvvRateWinR1,vvvRateWinR2; //ring buffers for past rates: windowSize x trials x nodes
  std::vector<std::vector<std::vector<double> > > vvvStimWin; //ring buffer for stim: windowSize x trials x nodes!
  //
  // initialize struct to handle average and variance calculations.
  //
  double avgSpaceR1=0.0,varSpaceR1=0.0,avgSpaceR2=0.0,varSpaceR2=0.0;
  std::vector<double> vZeroN(dim,0.0); //to reset
  node_avg_t avgFR1,avgFR2;
  avgFR1.init(trials);
  avgFR1.reset(vZeroN);
  avgFR2.init(trials);
  avgFR2.reset(vZeroN);
  std::vector<double> vAvRtR1(dim,0.0);
  std::vector<double> vAvRt2R1(dim,0.0);
  std::vector<double> vAvRtR2(dim,0.0);
  std::vector<double> vAvRt2R2(dim,0.0);
  //
  bool bStimIsOn=false,bStoppingStim=false;
  size_t turnOffUpdateStimRingBuffer=0;
  size_t ringpos=0,ringcur=0;
  //
  // stimulus intervals. complex case with two time intervals and multiple positions.
  //
  std::vector<multi_stim_interval_t> vsi(2);
  multi_stim_interval_t firstStim(tstimON,tstimOFF,vPosOne);
  multi_stim_interval_t secondStim(tstimCue,nsteps,vPosTwo);
  vsi[0]=firstStim;
  vsi[1]=secondStim;
  size_t nxtstim=0; // stimulus interval index
  //
  // MT tasks.
  //
  vvvRateWinR1.resize(windowSize);
  vvvRateWinR2.resize(windowSize);
  vvvStimWin.resize(windowSize);
  for(size_t i=0;i<vvvRateWinR1.size();++i)
    {
    vvvRateWinR1[i].resize(trials);
    vvvRateWinR2[i].resize(trials);
    vvvStimWin[i].resize(trials);
    for(size_t j=0;j<vvvRateWinR1[i].size();++j)
      {
      vvvRateWinR1[i][j].resize(dim,0.0);
      vvvRateWinR2[i][j].resize(dim,0.0);
      vvvStimWin[i][j].resize(dim,0.0);
      }
    }
  no_stats_task_t transientTask;
  transientTask.tbeg=0;
  transientTask.tend=tis;
  transientTask.ringpos=0;
  transientTask.delayff=delayFF;
  transientTask.stim_noiseOFF=stim_noiseOFF;
  transientTask.pvs1=&vstepperR1;
  transientTask.pvs2=&vstepperR2;
  transientTask.pvvvRateWinR1=&vvvRateWinR1;
  transientTask.pvvvRateWinR2=&vvvRateWinR2; 
  std::vector<no_stats_task_t>  vTransientTask(trials,transientTask);
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
  runTask1.pvvvRateWinR1=&vvvRateWinR1;
  runTask1.pvvvStimWinR1=&vvvStimWin;
  runTask1.pvvvRateWinR2=&vvvRateWinR2;
  std::vector<stats1_task_t> vRunTask1(nThreads,runTask1);
  for(size_t i=0;i<vRunTask1.size();++i)
    {
    vRunTask1[i].beg=i;
    vRunTask1[i].taskgen.seed(rd());
    }
  START(integrate); // timer macros to print out run time: START(name);do something; STOP(name);
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  // run transient task.
  //
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  START(transient);
  for(size_t i=0;i<vTransientTask.size();++i)
    { // transient task: run from 0...tis-windowSize using/filling only R1 steppers' ringbuffer.
    vTransientTask[i].index=i;
    vTransientTask[i].taskgen.seed(rd());
    mtpool.run_task(std::bind(run_no_stats_task,&vTransientTask[i]));
    }
  mtpool.wait();
  STOPN(transient,tis);
  ringpos=vTransientTask[0].ringpos%windowSize;//ok
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  // main integration and evaluation loop.
  //
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  bStimIsOn=false;
  bStoppingStim=false;
  turnOffUpdateStimRingBuffer=0;
  nxtstim = 0;
  for(size_t ti=tis;ti<=nsteps;++ti)
    {
    if(ti<nsteps)
      {
      if(ti==tswitch)
	{
	avgSpaceR1=0.0,varSpaceR1=0.0, avgSpaceR2=0.0,varSpaceR2=0.0; // switching => evaluate space avg and var.
	avgFR1.get_avg_var_space(avgSpaceR1,varSpaceR1);
	avgFR2.get_avg_var_space(avgSpaceR2,varSpaceR2);
	fprintf(stderr,"%zu: avgFRspace=%f +/- %f (cnt %zu), %f +/- %f \n",ti,avgSpaceR1,sqrt(varSpaceR1),avgFR1.m_count,avgSpaceR2,sqrt(varSpaceR2));
	avgFR1.reset(vZeroN);
	avgFR2.reset(vZeroN);
	if(bIsAt)
	  {
	  fprintf(stderr,"SWITCH at ti=%zu\n",ti);
	  for(size_t wm=0;wm<vstepperR2.size();++wm)
	    vstepperR2[wm].set_connectivity(J0sb/dim,J1sb/dim,J2WMfromSE,S2WMfromSE);
	  }
	else
	  fprintf(stderr,"NO switch \n");
	}
      //
      // switch stimulus according to time intervals defined in vsi.
      //
      if( (nxtstim<vsi.size()) && (ti==vsi[nxtstim].beg) )      // beg => ON
	{
	avgSpaceR1=0.0,varSpaceR1=0.0, avgSpaceR2=0.0,varSpaceR2=0.0; // switching => evaluate space avg and var.
	avgFR1.get_avg_var_space(avgSpaceR1,varSpaceR1);
	avgFR2.get_avg_var_space(avgSpaceR2,varSpaceR2);
	fprintf(stderr,"%zu: avgFRspace=%f +/- %f (cnt %zu), %f +/- %f \n",ti,avgSpaceR1,sqrt(varSpaceR1),avgFR1.m_count,avgSpaceR2,sqrt(varSpaceR2));
	avgFR1.reset(vZeroN);
	avgFR2.reset(vZeroN);
	// do stimulus ring buffer updates now. also remember when to stop the updates.
	for(size_t i=0;i<vRunTask1.size();++i)
	  {
	  vRunTask1[i].stim_noiseOFF=0.0;
	  vRunTask1[i].stim_noiseON=stim_noiseON; 
	  }
	size_t curr_npos = vsi[nxtstim].vpos.size(); 
	bStimIsOn=true;
	bStoppingStim=false;
	
	for(size_t i=0;i<vstepperR1.size();++i)//trials
	  vstepperR1[i].viext = vvStimBase[vsi[nxtstim].vpos[i%curr_npos]];
          
	turnOffUpdateStimRingBuffer=vsi[nxtstim].end+windowSize;
	fprintf(stderr,"bStimIsOn: %s\nbStoppingStimRingBuffer: %s turn off at %zu.\n",
		bStimIsOn?"true":"false",bStoppingStim?"true":"false",
		turnOffUpdateStimRingBuffer);
	}
      else if( (nxtstim<vsi.size()) && (ti==vsi[nxtstim].end) ) // end => OFF
	{
	avgSpaceR1=0.0,varSpaceR1=0.0, avgSpaceR2=0.0,varSpaceR2=0.0; // switching => evaluate space avg and var.
	avgFR1.get_avg_var_space(avgSpaceR1,varSpaceR1);
	avgFR2.get_avg_var_space(avgSpaceR2,varSpaceR2);
	fprintf(stderr,"%zu: avgFRspace=%f +/- %f (cnt %zu), %f +/- %f \n",ti,avgSpaceR1,sqrt(varSpaceR1),avgFR1.m_count,avgSpaceR2,sqrt(varSpaceR2));
	avgFR1.reset(vZeroN);
	avgFR2.reset(vZeroN);
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
      // preparing to write out FR and Stim
      }
    if( (ti>tis) && (0==(ti%saveStep)) )
      { // save FRs to write out data, happens every windowSize time steps ===> should happen every 10th time step
      vvvRateAvgR1.push_back(get_ringbuffer_avg(vvvRateWinR1,ringcur,saveStep));
      vvvRateAvgR2.push_back(get_ringbuffer_avg(vvvRateWinR2,ringcur,saveStep));	
      fprintf(stderr,"ti=%zu\n",ti);
      }
    // to save stimulus injection matrix to write out just once!
    if( (ti>=tstimON) && (ti<tstimOFF) )
      {
      if( (ti-tis)%saveStep == 0)
	{
	size_t td=(ti-tis)/saveStep;
	for(size_t np=0;np<npos;++np)
	  if(td<vvStimSave[np].size())
	    vvStimSave[np][td]=vvvStimWin[ringcur][np];
	}
      }
    else if(ti>=tstimCue)
      {
      if( (ti-tis)%saveStep == 0)
	{
	size_t td=(ti-tis)/saveStep;
	for(size_t np=0;np<cpos;++np)
	  if(td<vvStimSaveCue[np].size())
	    vvStimSaveCue[np][td]=vvvStimWin[ringcur][np];
	}
      }
    }
  //
  // end of main integration loop, dump run time on stderr.
  //
  STOPN(integrate,(nsteps-tis)); // timer macros to print out run time: START(name);do something; STOP(name);  
  
  // write out numpy files: FR and Stim: save FR all trials, all nodes, time averaged with deltaT=10:
  char bufFR1[1024];
  char bufFR2[1024];
  snprintf(bufFR1,1024,"%s/%s_frR1avg10_fb_trials%zu_ws%zu_npos%zu.npy",sBaseDir.c_str(),bIsAt?"in":"out",trials,windowSize,npos);
  snprintf(bufFR2,1024,"%s/%s_frR2avg10_fb_trials%zu_ws%zu_npos%zu.npy",sBaseDir.c_str(),bIsAt?"in":"out",trials,windowSize,npos);
  fprintf(stderr,"size1=%zu  %zu  %zu\n",vvvRateAvgR1.size(),vvvRateAvgR1[0].size(),vvvRateAvgR1[0][0].size()); // time, trials, nodes
  npy_save(bufFR1,vvvRateAvgR1); 
  npy_save(bufFR2,vvvRateAvgR2);
  fprintf(stderr,"saved FRs\n");
  
  //just once if no stimNoise
  char bufFRstim[1024];
  snprintf(bufFRstim,1024,"%s/%s_stim-fravg10_fb_trials%zu_ws%zu_npos%zu.npy",sBaseDir.c_str(),bIsAt?"in":"out",trials,windowSize,npos);
  npy_save(bufFRstim,vvStimSave);
  fprintf(stderr,"saved stim\n");
  
  snprintf(bufFRstim,1024,"%s/%s_stimCue-fravg10_fb_trials%zu_ws%zu_npos%zu.npy",sBaseDir.c_str(),bIsAt?"in":"out",trials,windowSize,npos);
  npy_save(bufFRstim,vvStimSaveCue);
  fprintf(stderr,"saved cue\n");
  //
  fprintf(stderr,"saved FR1 size= %zu %zu %zu\n",vvvRateAvgR1.size(),vvvRateAvgR1[0].size(),vvvRateAvgR1[0][0].size()); // time, trials, nodes
  fprintf(stderr,"saved FR2 size= %zu %zu %zu\n",vvvRateAvgR2.size(),vvvRateAvgR2[0].size(),vvvRateAvgR2[0][0].size()); // time, trials, nodes
  
  fprintf(stderr,"saved Stim size= %zu %zu %zu\n",vvStimSave.size(),vvStimSave[0].size(),vvStimSave[0][0].size()); // npos, time
  fprintf(stderr,"saved Cue  size= %zu %zu %zu\n",vvStimSaveCue.size(),vvStimSaveCue[0].size(),vvStimSaveCue[0][0].size()); // npos, time
  //
  return 0;
  }
