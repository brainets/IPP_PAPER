// -*- c++ -*-
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// Copyright 2019 Johannes Hausmann, Nicole Voges (drjoe@free.fr, nicole.voges@gmx.com)
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* modified: Sun Mar 26 16:47:08 CEST 2023 joe */
#include "util.h"
#include <thread>
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// stuff to handle creation of directories...
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::string append_slash(const std::string& crsDir)
  {
  std::string sDir;
  if(crsDir.empty())
    return sDir;
  sDir=crsDir;
  if(sDir.back()!='/')
    sDir+='/';
  return sDir;
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::string check_dir(const std::string& crsDir)
  {
  struct stat statBuf;
  std::string sDir;
  if(0==stat(crsDir.c_str(),&statBuf))
    {
    if(S_ISDIR(statBuf.st_mode)) 
      sDir=append_slash(crsDir);
    }
  return sDir;
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::string make_dir(const std::string& crsDir)
  {
  std::string sDir=check_dir(crsDir);
  if(sDir.empty())
    {
    int ret=mkdir(crsDir.c_str(),S_IRWXU | S_IRGRP | S_IXGRP);
    if(0==ret)
      sDir=check_dir(crsDir);
    }
  return sDir;
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// helper functions for averages.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// given a 3d vector crvvv: nSteps x nTrials x nNodes, representing a ring buffer of 2d vectors of
// size nTrials x nNodes, calculate the time average from wsz values up to and including the ring
// buffer's current index curr.
// note that the function returns a "regular" average across all nSteps if curr=0 and wsz=-1.
//
//    returns a 2d vector: nTrials x nNodes of time averaged values.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<std::vector<double> > get_ringbuffer_avg(const std::vector<std::vector<std::vector<double> > >& crvvv,
						     size_t curr, size_t wsz)
  {
  std::vector<std::vector<double> > vv;
  if(crvvv.empty())
    return vv;
  size_t nsteps=crvvv.size(),windowSize=nsteps,ntrials=0,nnodes=0;
  if(0<nsteps)
    ntrials=crvvv[0].size();
  if(0<ntrials)
    nnodes=crvvv[0][0].size();  
  if( (ntrials==0) || (nnodes==0) )
    return vv;
  if(nsteps>wsz)
    nsteps = wsz;
  vv.resize(ntrials);
  for(size_t i=0;i<vv.size();++i)
    vv[i].resize(nnodes,0.0);
  size_t index=curr;
  for(size_t i=0;i<nsteps;++i) //time
    {
    index = (curr+windowSize-i)%windowSize; //count backwards from currPos in ringbuffer to obtain last wsz current values with pbc
    for(size_t j=0;j<ntrials;++j)
      for(size_t k=0;k<nnodes;++k)
	vv[j][k] += crvvv[index][j][k];
    }
  for(size_t i=0;i<vv.size();++i)
    for(size_t j=0;j<vv[i].size();++j)
      vv[i][j] /= nsteps;

  return vv;
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// get indices of maximum values (i.e., index of maximum FR node per trial).
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<size_t> get_trial_max(const std::vector<std::vector<double> >& crvv)
  {
  std::vector<size_t> vmax;
  if(!crvv.empty() && !crvv[0].empty())
    {
    size_t ntrials=crvv.size(),nnodes=crvv[0].size();
    vmax.resize(ntrials,0);
    for(size_t i=0;i<ntrials;++i)  // trials
      {
      size_t maxn=0;
      double maxval=0.0;
      for(size_t n=0;n<nnodes;++n) // nodes
	if(maxval<crvv[i][n])
	  {
	  maxval=crvv[i][n];
	  maxn=n;
	  }
      vmax[i]=maxn;
      }
    }
  return vmax;
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// node_avg_t
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void node_avg_t::init(size_t ntrials)
  {
  std::vector<std::vector<double> >(ntrials).swap(m_vvavg);
  std::vector<std::vector<double> >(ntrials).swap(m_vvavg2);
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void node_avg_t::reset(const std::vector<double>& crvZero)
  {
  ::reset(m_vvavg,crvZero);
  ::reset(m_vvavg2,crvZero);
  m_count=0;
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void node_avg_t::update(const std::vector<std::vector<double> >& crvv)
  {
  if( m_vvavg.empty() || m_vvavg[0].empty() || (crvv.size()!=m_vvavg.size()) || (crvv[0].size()!=m_vvavg[0].size()) )
    return;
  size_t ntrials=m_vvavg.size(),nnodes=m_vvavg[0].size();
  for(size_t i=0;i<ntrials;++i) // trials
    {
    for(size_t n=0;n<nnodes;++n)       // nodes
      {
      m_vvavg[i][n]+=crvv[i][n];
      m_vvavg2[i][n]+=crvv[i][n]*crvv[i][n];
      }
    }
  ++m_count;
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// get straight average and variance of accumulated data, aka time average per trial and node:
//     rvvAvg[trial][node] = m_vvavg[trial][node]/m_count
//     rvvVar[trial][node] = m_vvavg2[trial][node]/m_count - rvvAvg[trial][node]*rvvAvg[trial][node]
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void node_avg_t::get_avg_var_time(std::vector<std::vector<double> >& rvvAvg, std::vector<std::vector<double> >& rvvVar)
  {
  rvvAvg=m_vvavg;
  rvvVar=m_vvavg2;
  if( (0<m_count) && !m_vvavg.empty() && !m_vvavg[0].empty() )
    {
    size_t ntrials=m_vvavg.size(),nnodes=m_vvavg[0].size();
    double nrm=1.0/m_count;
    for(size_t i=0;i<ntrials;++i)  // trials
      for(size_t n=0;n<nnodes;++n) // nodes
	{
	rvvAvg[i][n]*=nrm;
	rvvVar[i][n]*=m_vvavg2[i][n]*nrm-rvvAvg[i][n]*rvvAvg[i][n]; // var=<X²>-<X>²
	}
    }
  return;
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// get trial-averaged node average and variance of accumulated data, aka space average.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void node_avg_t::get_avg_var_space(double& rAvg, double& rVar)
  {
  double avg=0.0,var=0.0;
  if( (0<m_count) && !m_vvavg.empty() && !m_vvavg[0].empty() )
    {
    size_t ntrials=m_vvavg.size(),nnodes=m_vvavg[0].size();
    double nrm=1.0/(m_count*nnodes),nrm2=1.0/(m_count*m_count*nnodes);
    for(size_t i=0;i<ntrials;++i)  // trials
      {
      double trialavg=0.0,trialvar=0.0;
      for(size_t n=0;n<nnodes;++n) // nodes
	{
	trialavg+=m_vvavg[i][n];
	trialvar+=m_vvavg[i][n]*m_vvavg[i][n];
	}
      trialavg*=nrm;
      trialvar=nrm2*trialvar-trialavg*trialavg;
      avg+=trialavg;
      var+=trialvar;
      }
    avg/=ntrials;
    var/=ntrials;
    }
  rAvg=avg;
  rVar=var;
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// get trial-averaged node average and variance of accumulated data, aka space average.
// in addition, return the trial-average for nodes in the range pos-dlt...pos+dlt.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double node_avg_t::get_avg_var_space_max(double& rAvg, double& rVar, size_t pos, size_t dlt)
  {
  double rangeavg=0.0,avg=0.0,var=0.0;
  if( (0<m_count) && !m_vvavg.empty() && !m_vvavg[0].empty() )
    {
    size_t ntrials=m_vvavg.size(),nnodes=m_vvavg[0].size(),rangesz=1+2*dlt,rangeend=pos+2*dlt;
    double nrm=1.0/(m_count*nnodes),nrm2=1.0/(m_count*m_count*nnodes),rangenrm=1.0/(m_count*rangesz);
    for(size_t i=0;i<ntrials;++i)  // trials
      {
      double trialavg=0.0,trialvar=0.0,trialrangeavg=0.0;
      for(size_t n=0;n<nnodes;++n) // nodes
	{
	trialavg+=m_vvavg[i][n];
	trialvar+=m_vvavg[i][n]*m_vvavg[i][n];
	if( (pos<=n+dlt) && (n+dlt<=rangeend) )
	  trialrangeavg+=m_vvavg[i][n];
	}
      trialavg*=nrm;
      trialvar=nrm2*trialvar-trialavg*trialavg;
      avg+=trialavg;
      var+=trialvar;
      rangeavg+=rangenrm*trialrangeavg;
      }
    avg/=ntrials;
    var/=ntrials;
    rangeavg/=ntrials;
    }
  rAvg=avg;
  rVar=var;
  return rangeavg;
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void node_avg_t::get_node_avg(std::vector<double>& rvAvg)
  {
  std::vector<double> vavg;
  if( (0<m_count) && !m_vvavg.empty() && !m_vvavg[0].empty() )
    {
    size_t ntrials=m_vvavg.size(),nnodes=m_vvavg[0].size();
    double nrm=1.0/(m_count*ntrials);
    vavg.resize(nnodes,0.0);
    for(size_t i=0;i<ntrials;++i)  // trials
      for(size_t n=0;n<nnodes;++n) // nodes
	vavg[n]+=nrm*m_vvavg[i][n];
    }
  vavg.swap(rvAvg);
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// get indices of maximum accumulated values (i.e., index of maximum average FR node per trial).
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<size_t> node_avg_t::get_trial_max(void)
  {
  std::vector<size_t> vmax;
  if( (0<m_count) && !m_vvavg.empty() && !m_vvavg[0].empty() )
    {
    size_t ntrials=m_vvavg.size(),nnodes=m_vvavg[0].size();
    vmax.resize(ntrials,0);
    for(size_t i=0;i<ntrials;++i)  // trials
      {
      size_t maxn=0;
      double maxval=0.0;
      for(size_t n=0;n<nnodes;++n) // nodes
	if(maxval<m_vvavg[i][n])
	  {
	  maxval=m_vvavg[i][n];
	  maxn=n;
	  }
      vmax[i]=maxn;
      }
    }
  return vmax;
  }

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// heatmap color generator.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<uint32_t> heatmap()
{
  std::vector<uint32_t> vc;
  uint32_t a=0xFF,max=1+0x7FFF;
  double x=0.0,dx=M_PI/max;
  a<<=24;
  for(uint32_t n=0;n<max;++n,x+=dx)
    {
      uint32_t r=static_cast<uint32_t>(128.0-127.0*cos(x));
      uint32_t g=static_cast<uint32_t>(255.0*sin(x));
      uint32_t b=static_cast<uint32_t>(128.0+127.0*cos(x));
      //fprintf(stderr,"%u: %02x%02x%02x%02x",n,r,g,b,0xFF);
      r<<=16;
      g<<=8;
      vc.push_back(r|g|b|a);
      //fprintf(stderr," => %08x\n",vc.back());
    }
  return vc;
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// function to directly generate color png plot for a single ring.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double write_color_png_file(const char* pFileName, const std::vector<std::vector<double> >& crvv, size_t offset, bool bZero)
{
  std::vector<uint32_t> vc=heatmap();
  //int w=(offset<crvv.size()) ? crvv.size()-offset:0,h=(0<w) ? 1+crvv[0].size() : 0;
  int w=(offset<crvv.size()) ? crvv.size()-offset:0,h=(0<w) ? crvv[0].size() : 0; //08/11/2019
  int depth=8,color_type=PNG_COLOR_TYPE_RGB_ALPHA;
  size_t clim=(0<vc.size()) ? vc.size()-1:0;
  const uint32_t cblack=((uint32_t)0xFF)<<24;
  double dmax=0.0,dmin=10.0e10,scl=clim,realdmin=0.0;
  std::vector<uint32_t> vpix(w,0);
  if(0==h)
    return -1.0;
  for(size_t i=offset;i<crvv.size();++i)
    for(size_t j=0;j<crvv[i].size();++j)
      {
	if(dmax<crvv[i][j])
	  dmax=crvv[i][j];
	if(dmin>crvv[i][j])
	  dmin=crvv[i][j];
      }
  realdmin=dmin;
  if(bZero && (0.>dmin) )
    {
    fprintf(stderr,"write_color_png_file: dmin = %f < 0 in data, setting to zero!",dmin);
    dmin=0.0;
    }
  scl/=dmax-dmin;
  fprintf(stderr,"write_color_png_file: writing PNG image %dx%d depth %d (dmax %f dmin %f (real %f) => scl %f)\n",w,h,depth,dmax,dmin,realdmin,scl);
  //
  FILE* pFile=fopen(pFileName, "wb");
  if(NULL==pFile)
    return -1.0;
  //
  png_structp png_ptr=png_create_write_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);
  if (NULL==png_ptr)
    { fclose(pFile);return -1.0; }
  png_infop info_ptr=png_create_info_struct(png_ptr);
  if (NULL==info_ptr)
    { fclose(pFile);return -1.0; }
  if(setjmp(png_jmpbuf(png_ptr)))
    { fclose(pFile);return -1.0; }
  png_init_io(png_ptr, pFile);
  // write png header...
  if (setjmp(png_jmpbuf(png_ptr)))
    { fclose(pFile);return -1.0; }
  png_set_IHDR(png_ptr,info_ptr,w,h,depth,color_type,PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_BASE,PNG_FILTER_TYPE_BASE);
#if 1
#warning "Writing description text to PNG file"
  {
  png_text text_ptr[3];
  char key0[] = "Min";
  char text0[128];
  snprintf(text0,127,"%f",realdmin);
  char key1[] = "Max";
  char text1[128];
  snprintf(text1,127,"%f",dmax);
  char key2[] = "Scale";
  char text2[128];
  snprintf(text2,127,"%f",scl);
  //
  text_ptr[0].key = key0;
  text_ptr[0].text = text0;
  text_ptr[0].compression = PNG_TEXT_COMPRESSION_NONE;//PNG_TEXT_COMPRESSION_zTXt;
  text_ptr[0].itxt_length = 0;
  text_ptr[0].lang = NULL;
  text_ptr[0].lang_key = NULL;
  //
  text_ptr[1].key = key1;
  text_ptr[1].text = text1;
  text_ptr[1].compression = PNG_TEXT_COMPRESSION_NONE;//PNG_TEXT_COMPRESSION_zTXt;
  text_ptr[1].itxt_length = 0;
  text_ptr[1].lang = NULL;
  text_ptr[1].lang_key = NULL;
  //
  text_ptr[2].key = key2;
  text_ptr[2].text = text2;
  text_ptr[2].compression = PNG_TEXT_COMPRESSION_NONE;//PNG_TEXT_COMPRESSION_zTXt;
  text_ptr[2].itxt_length = 0;
  text_ptr[2].lang = NULL;
  text_ptr[2].lang_key = NULL;
  //
  png_set_text(png_ptr,info_ptr,text_ptr, 3);
  }
#endif
  png_write_info(png_ptr,info_ptr);
  // write bytes...
  if (setjmp(png_jmpbuf(png_ptr)))
    { fclose(pFile);return -1.0; }
  ///size_t ci=0;
  for(int j=0;j<h;++j)
    {
    for(size_t i=offset,k=0;i<crvv.size();++i,++k)
      {
      ///ci=static_cast<size_t>(scl*(dmax-crvv[i][j]));
      double val=scl*(crvv[i][j]-dmin);
      uint32_t col=(bZero && (0.>val) ) ? cblack : vc[clim-static_cast<size_t>(val)]; // NEW: invalid values < 0.0 are plotted in black!
      ///if(ci>=vc.size()) fprintf(stderr,"buggerSE...%zu >= %zu...\n",ci,vc.size());
      vpix[k]=col;
      }
    png_write_row(png_ptr,reinterpret_cast<png_bytep>(&vpix[0]));
    }  
  // end write...
  if (setjmp(png_jmpbuf(png_ptr)))
    { fclose(pFile);return -1.0; }
  png_write_end(png_ptr, NULL);
  fclose(pFile);
  return dmax-dmin;
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// function to directly generate color png plot for two rings.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*void write_2color_png_file(const char* pFileName, const std::vector<std::vector<double> >& crvv, size_t offset)
{
  std::vector<uint32_t> vc=heatmap();
  int w=(offset<crvv.size()) ? crvv.size()-offset:0,h=(0<w) ? 1+crvv[0].size() : 0;
  int depth=8,color_type=PNG_COLOR_TYPE_RGB_ALPHA;
  double dmax=0.0,dmin=10.0e10,scl=vc.size()-1;
  std::vector<uint32_t> vpix(w,0);
  if(0==h)
    return;
  for(size_t i=offset;i<crvv.size();++i)
    for(size_t j=0;j<crvv[i].size();++j)
      {
	if(dmax<crvv[i][j])
	  dmax=crvv[i][j];
	if(dmin>crvv[i][j])
	  dmin=crvv[i][j];
      }  
  scl/=dmax-dmin;
  fprintf(stderr,"writing PNG image %dx%d depth %d (dmax %f dmin %f => scl %f)\n",w,h,depth,dmax,dmin,scl);
  //
  FILE* pFile=fopen(pFileName, "wb");
  if(NULL==pFile)
    return;
  //
  png_structp png_ptr=png_create_write_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);
  if (NULL==png_ptr)
    { fclose(pFile);return; }
  png_infop info_ptr=png_create_info_struct(png_ptr);
  if (NULL==info_ptr)
    { fclose(pFile);return; }
  if(setjmp(png_jmpbuf(png_ptr)))
    { fclose(pFile);return; }
  png_init_io(png_ptr, pFile);
  // write png header...
  if (setjmp(png_jmpbuf(png_ptr)))
    { fclose(pFile);return; }
  png_set_IHDR(png_ptr,info_ptr,w,h,depth,color_type,PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_BASE,PNG_FILTER_TYPE_BASE);
  png_write_info(png_ptr,info_ptr);
  // write bytes...
  if (setjmp(png_jmpbuf(png_ptr)))
    { fclose(pFile);return; }
  size_t ci=0;
  for(int j=0;j<h;++j)
    {
      if(j<(h>>1))
	for(size_t i=offset,k=0;i<crvv.size();++i,++k)
	  {
	    ci=static_cast<size_t>(scl*(dmax-crvv[i][j]));
	    //if(ci>=vc.size()) fprintf(stderr,"buggerSE...%zu >= %zu...\n",ci,vc.size());
	    vpix[k]=vc[ci];
	  }
      else if(j==(h>>1))
	{
	  for(size_t i=offset,k=0;i<crvv.size();++i,++k) // draw a black line in the middle ;-)
	    vpix[k]=0x0;//FF00FFFF;
	}
      else
	for(size_t i=offset,k=0;i<crvv.size();++i,++k)
	  {
	    ci=static_cast<size_t>(scl*(dmax-crvv[i][j-1]));
	    //if(ci>=vc.size()) fprintf(stderr,"buggerWM...%zu >= %zu...\n",ci,vc.size());
	    vpix[k]=vc[ci];
	  }
      png_write_row(png_ptr,reinterpret_cast<png_bytep>(&vpix[0]));
    }  
  // end write...
  if (setjmp(png_jmpbuf(png_ptr)))
    { fclose(pFile);return; }
  png_write_end(png_ptr, NULL);
  fclose(pFile);
  return;
}*/
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// write numpy file stuff.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
char numpy_map_type(const std::type_info& t)
  {
  if(t == typeid(float) ) return 'f';
  if(t == typeid(double) ) return 'f';
  if(t == typeid(long double) ) return 'f';

  if(t == typeid(int) ) return 'i';
  if(t == typeid(char) ) return 'i';
  if(t == typeid(short) ) return 'i';
  if(t == typeid(long) ) return 'i';
  if(t == typeid(long long) ) return 'i';

  if(t == typeid(unsigned char) ) return 'u';
  if(t == typeid(unsigned short) ) return 'u';
  if(t == typeid(unsigned long) ) return 'u';
  if(t == typeid(unsigned long long) ) return 'u';
  if(t == typeid(unsigned int) ) return 'u';
  
  if(t == typeid(bool) ) return 'b';
  
  if(t == typeid(std::complex<float>) ) return 'c';
  if(t == typeid(std::complex<double>) ) return 'c';
  if(t == typeid(std::complex<long double>) ) return 'c';
  
  return '?';
}
