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
  size_t index,tbeg,tend;
  std::vector<delay_stepper_t>* pvs;
  std::mt19937 taskgen;
};
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void run_no_stats_task(no_stats_task_t* pt)
  {
  if( (NULL==pt) || (NULL==pt->pvs) || (pt->index>=pt->pvs->size()) )
    return;
  size_t tb=pt->tbeg,te=pt->tend;
  delay_stepper_t& rStepper=(*(pt->pvs))[pt->index];
  //
  for(size_t t=tb;t<te;++t)
    rStepper.step(t);
  //fprintf(stderr,"stepper %zu: iterated from t=%zu...%zu.\n",pt->index,tb,te);
  return;
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
struct stats_task_t
{
  size_t index,tbeg,tend;
  std::vector<delay_stepper_t>* pvs;
  std::vector<std::vector<std::vector<double> > >* pvvvRate;
  std::vector<std::vector<std::vector<double> > >* pvvvStim;
  std::mt19937 taskgen;
};
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void run_stats_task(stats_task_t* pt)
  {
  if( (NULL==pt) || (NULL==pt->pvs) || (NULL==pt->pvvvRate) || (NULL==pt->pvvvStim) || (pt->index>=pt->pvs->size()) )
    return;
  size_t tb=pt->tbeg,te=pt->tend,index=pt->index;
  delay_stepper_t& rStepper=(*(pt->pvs))[pt->index];
  std::vector<std::vector<std::vector<double> > >& rvvvr=*(pt->pvvvRate);
  std::vector<std::vector<std::vector<double> > >& rvvvs=*(pt->pvvvStim);
  //
  for(size_t t=tb,j=0;t<te;++t,++j)
    {
    rvvvr[j][index]=rStepper.step(t);
    rvvvs[j][index]=rStepper.viext; // fill stim ringbuffer with the noisy stimulus actually used in this time step!
    }
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
struct stats1_task_t
{
  size_t beg,end,step,t;
  std::vector<delay_stepper_t>* pvs;
  std::vector<std::vector<double> >* pvvRate;
  std::vector<std::vector<double> >* pvvStim;
  std::mt19937 taskgen;
};
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void run_stats1_task(stats1_task_t* pt)
  {
  if( (NULL==pt) || (NULL==pt->pvs) || (NULL==pt->pvvRate) || (NULL==pt->pvvStim) )
    return;
  size_t beg=pt->beg,end=pt->end,step=pt->step,t=pt->t;
  std::vector<delay_stepper_t>& rvStepper=*(pt->pvs);
  std::vector<std::vector<double> >& rvvr=*(pt->pvvRate);
  std::vector<std::vector<double> >& rvvs=*(pt->pvvStim);
  for(size_t j=beg;j<end;j+=step)
    {
    rvvr[j]=rvStepper[j].step(t);
    rvvs[j]=rvStepper[j].viext; // pass back the stimulus actually used in this time step!
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
  size_t trials=4000, nsteps=8500, tis=4000, tstimON=5000, tstimOFF=6500;  //the real thing: final timings, determined by SB
  //size_t trials=100, nsteps=850, tis=400, tstimON=500, tstimOFF=650;  // minitest of functionality.
  double JStim=2.0, SStim=8.0;  // stimAmplitude and width
  int MStim=1;                  // stimulus mode, default=1=>gauss
  
  size_t dim=100, delay=10;
  double dt=0.01;
  double Noise=0.50, Iext=0.0; //defaults: add 50% noise, gaussian stim with sigm2=64 centered at node 50  
  std::string sBaseDir="1ringRes";
  
  std::string s=make_dir(sBaseDir);
  fprintf(stderr,"storing result in base directory <%s>...\n",s.c_str());
#if 0
  double J0 = -25.0,J1 = 11.0; // SB
  bool bIsSU = false;
#else
  double J0 = -30.0,J1 = -8.0; // SU
  bool bIsSU = true;
#endif
  size_t npos = 4;     
  size_t posStep = dim/npos;
  size_t windowSize= 400;  //memory for rates = ring buffer size
  fprintf(stderr,"Parameters: ws=%zu npos=%zu trials=%zu state=%s\n",windowSize,npos,trials,bIsSU?"su":"sb");
  // MT stuff: initialize a thread pool with nThreads workers.
  size_t nThreads=std::thread::hardware_concurrency(); //NV schleppi: 4
  thread_pool_t mtpool(nThreads);
  //
  // init parallel steppers with different initial conditions.
  std::vector<delay_stepper_t> vstepper(trials);
  for(size_t i=0;i<trials;++i)
    {
    bool bOk=vstepper[i].init(dim,delay,dt,J0/dim,J1/dim,0.0,0.0,Noise,Iext,0);
    if(bOk)
      ;
    else
      fprintf(stderr,"stepper %zu: FAILED!\n",i);
    }
  // init stimuli.
  std::vector<double> vnull(dim,0.0);
  std::vector<std::vector<double> > vvStimBase(npos); // #npos stimuli x nodes
  for(size_t i=0;i<npos;++i)                        //5 features => posStep=100/5 => i%20 => 0,20,40,60,80
    vvStimBase[i]=generate_stimulus(dim,i*posStep,JStim,SStim,MStim); // 5x100nodes: shifted amplitudes on ring //new Ampl & width
  fprintf(stderr,"generate stim end\n");
  // rate matrix.
  std::vector<std::vector<double> > vvRate(trials);// of all nodes: trials x nodes
  std::vector<std::vector<double> > vvStim(trials);// of all nodes: trials x nodes
  for(size_t i=0;i<trials;++i)
    {
    vvRate[i]=vstepper[i].get_rate();
    vvStim[i]=vnull;
    }
  // will save average firing rates to file for t=tis...nsteps.
  std::vector<std::vector<std::vector<double> > > vvvRateAvg;
  vvvRateAvg.reserve(nsteps-tis);
  size_t saveStep = 10; 
  //stimulus injection matrix init to zero => once write out as .npy, only every 10th value according to avFR
  std::vector<std::vector<std::vector<double> > > vvStimSave(npos);
  for(size_t np=0;np<npos;++np)
    {
    vvStimSave[np].resize((nsteps-tis)/saveStep); // time x nodes
    reset(vvStimSave[np],vnull);
    }
  //
  std::vector<std::vector<std::vector<double> > > vvvRateWin; //ring buffer for past rates: windowSize x trials x nodes
  std::vector<std::vector<std::vector<double> > > vvvStimWin; //ring buffer for stim:       windowSize x trials x nodes
  //
  // initialize struct to handle average and variance calculations.
  //
  double avgSpace=0.0,varSpace=0.0;
  std::vector<double> vZeroN(dim,0.0); //to reset
  node_avg_t avgFR;
  avgFR.init(trials);
  avgFR.reset(vZeroN);
  std::vector<double> vAvRt(vvRate[0].size(),0.0);
  std::vector<double> vAvRt2(vvRate[0].size(),0.0);
  //
  std::vector<std::vector<double> > vvzeroStimPast(trials); //(re-)set to zeros for ringbuffer for stim, full size now!
  for(size_t i=0;i<trials;++i)
    vvzeroStimPast[i].resize(dim,0.0);
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
  vvvRateWin.resize(windowSize);
  vvvStimWin.resize(windowSize);
  for(size_t i=0;i<vvvRateWin.size();++i)
    {
    vvvRateWin[i].resize(trials);
    vvvStimWin[i].resize(trials);
    for(size_t j=0;j<vvvRateWin[i].size();++j)
      {
      vvvRateWin[i][j].resize(dim,0.0);
      vvvStimWin[i][j].resize(dim,0.0);
      }
    }
  no_stats_task_t transientTask; // transient: integrate but collect no information at all.
  transientTask.tbeg=0;
  transientTask.tend=tis-windowSize;
  transientTask.pvs=&vstepper;
  std::vector<no_stats_task_t>  vTransientTask(trials,transientTask);
  stats_task_t runTask;         // transient: integrate the last windowSize steps, start filling the ring buffers.
  runTask.tbeg=tis-windowSize;
  runTask.tend=tis;
  runTask.pvs=&vstepper;
  runTask.pvvvRate=&vvvRateWin;
  runTask.pvvvStim=&vvvStimWin;
  std::vector<stats_task_t> vRunTask(trials,runTask);
  for(size_t i=0;i<vRunTask.size();++i)
    {
    vRunTask[i].index=i;
    vRunTask[i].taskgen.seed(rd());
    }
  stats1_task_t runTask1;       // real task: just integrate one time step 
  runTask1.end=trials;
  runTask1.step=nThreads;
  runTask1.t=tis;
  runTask1.pvs=&vstepper;
  runTask1.pvvRate=&vvRate;
  runTask1.pvvStim=&vvStim;
  std::vector<stats1_task_t> vRunTask1(nThreads,runTask1);
  for(size_t i=0;i<vRunTask1.size();++i)
    {
    vRunTask1[i].taskgen.seed(rd());
    vRunTask1[i].beg=i;
    }
  START(integrate); // timer macros to print out run time: START(name);do something; STOP(name);
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  // run transient task.
  //
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  START(transient);
  for(size_t i=0;i<vTransientTask.size();++i)
    { // transient task: until tis-windowSize, no update of ringbuffers required.
    vTransientTask[i].index=i;
    vTransientTask[i].taskgen.seed(rd());
    mtpool.run_task(std::bind(run_no_stats_task,&vTransientTask[i]));
    }
  mtpool.wait();
  for(size_t i=0;i<vRunTask.size();++i)
    {
    // run task: until tis, preparing ring buffer for firing rate. note that the most
    // current rate for ti=tis-1 is at the *back* of the ring buffer, and the oldest
    // rate at the front (so ringpos=0 indicates the oldest value!).
    mtpool.run_task(std::bind(run_stats_task,&vRunTask[i]));
    }
  mtpool.wait();
  STOPN(transient,tis);
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  // main integration loop.
  //
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  bool bUpdateStimRingBuffer=false,bStoppingStimRingBuffer=false;
  size_t ringpos=0,ringcur=0; // ringpos: oldest vvRate in ring buffer vvvRateWin, initially at zero -- also used for stimPast.
  size_t turnOffUpdateStimRingBuffer=0;
  for(size_t ti=tis;ti<=nsteps;++ti)
    {
    if( (ti>tis) && (0==(ti%saveStep)) )
      { // save FRs to write out data, happens every 10th time steps.	
      vvvRateAvg.push_back(get_ringbuffer_avg(vvvRateWin,ringcur,saveStep));  
      fprintf(stderr,"ti=%zu\n",ti);
      }
    // save stimulus injection matrix to write out just once.
    if( (ti>=tstimON) && (ti<tstimOFF) )  
      {
      if( (ti-tis)%saveStep == 0)
	{
	for(size_t np=0;np<npos;++np)
	  vvStimSave[np][(ti-tis)/saveStep]=vvvStimWin[ringcur][np];
	}
      }
    //
    if(ti==nsteps)
      break; // exit, properly handling average rates above at the end of the main loop.
    //
    // switch stimulus according to time intervals defined in vsi.
    //
    if( (nxtstim<vsi.size()) && (ti==vsi[nxtstim].beg) )      // beg => ON
      {
      //
      avgSpace=0.0,varSpace=0.0; // switching => evaluate space avg and var.
      avgFR.get_avg_var_space(avgSpace,varSpace);
      fprintf(stderr,"%zu: avgSpace=%f +/- %f (count %zu)\n",ti,avgSpace,sqrt(varSpace),avgFR.m_count);
      avgFR.reset(vZeroN);
      //
      // do stimulus ring buffer updates now. also remember when to stop the updates.
      bUpdateStimRingBuffer=true;
      bStoppingStimRingBuffer=false;
      turnOffUpdateStimRingBuffer=vsi[nxtstim].end+windowSize;
      fprintf(stderr,"bUpdateStimRingBuffer: %s\nbStoppingStimRingBuffer: %s turn off at %zu.\n",
	      bUpdateStimRingBuffer?"true":"false",bStoppingStimRingBuffer?"true":"false",
	      turnOffUpdateStimRingBuffer);
      }
    else if( (nxtstim<vsi.size()) && (ti==vsi[nxtstim].end) ) // end => OFF
      {
      avgSpace=0.0,varSpace=0.0; // switching => evaluate space avg and var.
      avgFR.get_avg_var_space(avgSpace,varSpace);
      fprintf(stderr,"%zu: avgSpace=%f +/- %f (count %zu)\n",ti,avgSpace,sqrt(varSpace),avgFR.m_count);
      avgFR.reset(vZeroN);
      // turn off stim in steppers
      for(size_t i=0;i<vstepper.size();++i)
	vstepper[i].viext = vnull;
      //ti0=ti+25;
      ++nxtstim;
      bUpdateStimRingBuffer=false;
      bStoppingStimRingBuffer=true;
      fprintf(stderr,"bUpdateStimRingBuffer: %s\nbStoppingStimRingBuffer: %s turn off at %zu.\n",
	      bUpdateStimRingBuffer?"true":"false",bStoppingStimRingBuffer?"true":"false",
	      turnOffUpdateStimRingBuffer);
      }
    if(bStoppingStimRingBuffer && (turnOffUpdateStimRingBuffer<=ti) )
      { // we need to update the stimulus ring buffer windowSize steps after turning off the stimulus...
      bStoppingStimRingBuffer=false;
      fprintf(stderr,"bUpdateStimRingBuffer: %s\nbStoppingStimRingBuffer: %s turn off at %zu.\n",
	      bUpdateStimRingBuffer?"true":"false",bStoppingStimRingBuffer?"true":"false",
	      turnOffUpdateStimRingBuffer);
      }
    //
    if(bUpdateStimRingBuffer) //to do if stim is ON
      {
      for(size_t i=0;i<vstepper.size();++i)
	{
	vstepper[i].viext = vvStimBase[i%npos];
	}
      }
    //
    // integrate one step, so now vvRate contains the firing rates for time ti.
    //    
    for(size_t i=0;i<vRunTask1.size();++i)
      {
      vRunTask1[i].t=ti;
      mtpool.run_task(std::bind(run_stats1_task,&vRunTask1[i]));
      }
    mtpool.wait();
    //
    // update the ring buffers...
    //
    vvvRateWin[ringpos]=vvRate;            // ringbuffer rates: replace oldest with current values.
    if(bUpdateStimRingBuffer)
      vvvStimWin[ringpos]=vvStim;          // ringbuffer stimulus: replace oldest with current values, if required.
    else if(bStoppingStimRingBuffer)       // ringbuffer stimulus: replace oldest with zeros, if required.
      {
      //fprintf(stderr,"setting stim ringbuffer to zero at %zu.\n",ringpos);
      vvvStimWin[ringpos]=vvzeroStimPast;
      }
    ringcur=ringpos;                       // index of current rates in ring buffer.
    ringpos = (ringpos+1)%windowSize;      // update index ringpos so that it again points to the oldest rates.
    //
    // ...then update averages and variances with the new rates...
    //
    avgFR.update(vvRate);
    }
  //
  // end of main integration loop, dump run time on stderr.
  //
  STOPN(integrate,(nsteps-tis)); // timer macros to print out run time: START(name);do something; STOP(name);  
  avgSpace=0.0,varSpace=0.0;     // end of integration => evaluate space avg and var.
  avgFR.get_avg_var_space(avgSpace,varSpace);
  fprintf(stderr,"%zu: avgSpace=%f +/- %f (count %zu)\n",nsteps,avgSpace,sqrt(varSpace),avgFR.m_count);
  avgFR.reset(vZeroN); // well...are we ever going to reuse it?
  
  // write out FR, stim  
  std::vector<std::vector<double> > vvRateR0(vvvRateAvg.size());
  std::vector<std::vector<double> > vvRateRout(vvvRateAvg.size());
  for(size_t i=0;i<vvvRateAvg.size();++i)
    {
    vvRateR0[i].resize(dim,0.0);
    vvRateRout[i].resize(dim,0.0);
    }
  for(size_t ti=0;ti<vvvRateAvg.size();++ti)
    for(size_t tr=0;tr<vvvRateAvg[ti].size();++tr)
      for(size_t n=0;n<vvvRateAvg[ti][tr].size();++n)
	{
	if(0==tr)
	  vvRateR0[ti][n] =  vvvRateAvg[ti][tr][n];
	vvRateRout[ti][n] += vvvRateAvg[ti][tr][n];
	}
  for(size_t ti=0;ti<vvRateRout.size();++ti)
    for(size_t n=0;n<vvRateRout[ti].size();++n)
      vvRateRout[ti][n] /= trials;
  // write out simple heatmap of FR averages.
  char bufFR[1024];
  snprintf(bufFR,1024,"%s/%s_FRavg_trials%zu_ws%zu.png",sBaseDir.c_str(),bIsSU?"su":"sb",trials,windowSize);
  write_color_png_file(bufFR,vvRateRout); // averages.
  snprintf(bufFR,1024,"%s/%s_FR0_ws%zu.png",sBaseDir.c_str(),bIsSU?"su":"sb",windowSize);
  write_color_png_file(bufFR,vvRateR0); // trial 0.
  // write out numpy files: FR and Stim: save FR all trials, all nodes, time averaged with deltaT=10:
  snprintf(bufFR,1024,"%s/%s_FRavg_trials%zu_ws%zu.npy",sBaseDir.c_str(),bIsSU?"su":"sb",trials,windowSize);
  npy_save(bufFR,vvvRateAvg); 
  
  char bufFRstim[1024];
  snprintf(bufFRstim,1024,"%s/stimTraceAvg10_trials%zu_ws%zu.npy",sBaseDir.c_str(),trials,windowSize);
  npy_save(bufFRstim,vvStimSave);
  
  fprintf(stderr,"saved FRs size= %zu %zu %zu\n",vvvRateAvg.size(),vvvRateAvg[0].size(),vvvRateAvg[0][0].size()); // time, trials, nodes
  fprintf(stderr,"saved Stims size= %zu %zu %zu\n",vvStimSave.size(),vvStimSave[0].size(),vvStimSave[0][0].size()); // npos, time, nodes
  
  return 0;
  }
