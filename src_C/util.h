// -*- c++ -*-
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// Copyright 2019 Johannes Hausmann, Nicole Voges (drjoe@free.fr, nicole.voges@gmx.com)
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* modified: Sat Apr  1 09:41:15 CEST 2023 joe */
#ifndef __UTIL_H
#define __UTIL_H

#include <sys/stat.h> // for directory creation stuff.
#include <unistd.h>
#include <dirent.h>

#include <cmath>
#include <cstring>
#include <chrono>
#include <complex>
#include <functional>
#include <typeinfo> // for numpy stuff...
#include <vector>
#include <png.h>    // for PNG image stuff.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// timing stuff.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#define START(s) std::chrono::steady_clock::time_point beg##s=std::chrono::steady_clock::now();
#define STOP(s) std::chrono::steady_clock::time_point end##s=std::chrono::steady_clock::now();fprintf(stderr,"timer %s took %.3fms.\n",#s,0.001*std::chrono::duration_cast<std::chrono::microseconds>(end##s-beg##s).count());
#define STOPN(s,n) std::chrono::steady_clock::time_point end##s=std::chrono::steady_clock::now();fprintf(stderr,"timer %s took %.3fms (%.3f ms/op for %zu ops).\n",#s,0.001*std::chrono::duration_cast<std::chrono::microseconds>(end##s-beg##s).count(),0.001*static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end##s-beg##s).count())/n,n);
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// stuff to handle creation of directories...
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::string append_slash(const std::string& crsDir);
std::string check_dir(const std::string& crsDir);
std::string make_dir(const std::string& crsDir);
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// helper template to reset vectors to zero, usage: reset(vvR,vzeroR); reset(vvvvRRpS,vvvzeroRRS);
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template<typename T>
void reset(std::vector<T>& v, const T& zero)
  {
  for(size_t i=0;i<v.size();++i) // loop nodes.
    v[i]=zero;
  }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// helper functions for averages.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<std::vector<double> > get_ringbuffer_avg(const std::vector<std::vector<std::vector<double> > >& crvvv,
						     size_t curr=0, size_t wsz=-1);
std::vector<size_t> get_trial_max(const std::vector<std::vector<double> >& crvv);
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
struct node_avg_t
{
  size_t m_count;
  std::vector<std::vector<double> > m_vvavg;
  std::vector<std::vector<double> > m_vvavg2;
  void init(size_t ntrials);
  void reset(const std::vector<double>& crvZero);
  void update(const std::vector<std::vector<double> >& crvvVal);
  void get_avg_var_time(std::vector<std::vector<double> >& rvvAvg, std::vector<std::vector<double> >& rvvVar);
  void get_avg_var_space(double& rAvg, double& rVar);
  double get_avg_var_space_max(double& rAvg, double& rVar, size_t pos, size_t dlt);
  void get_node_avg(std::vector<double>& rvAvg);
  std::vector<size_t> get_trial_max(void);
};
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// functions to directly generate png plots.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<uint32_t> heatmap();
double write_color_png_file(const char* pFileName, const std::vector<std::vector<double> >& crvv, size_t offset=0, bool bZero=true);
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// write numpy file...stolen from https://github.com/rogersce/cnpy/blob/master/cnpy.h
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
char numpy_map_type(const std::type_info& t);
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template<typename T> std::vector<char> create_npy_header(const std::vector<size_t>& crvShape)
  {
  std::string sDict;
  std::vector<char> vHeader;
  if(crvShape.empty())
    return vHeader;
  int bet=1;
  char endian=(((char *)&bet)[0]) ? '<' : '>';
  sDict += "{'descr': '";
  sDict += endian;
  sDict += numpy_map_type(typeid(T));
  sDict += std::to_string(sizeof(T));
  sDict += "', 'fortran_order': False, 'shape': (";
  sDict += std::to_string(crvShape[0]);
  for(size_t i=1;i<crvShape.size();++i)
    {
    sDict += ", ";
    sDict += std::to_string(crvShape[i]);
    }
  if(1==crvShape.size())
    sDict += ",";
  sDict += "), }";
  //pad with spaces so that preamble+dict is modulo 16 bytes. preamble is 10 bytes. dict needs to end with \n
  int remainder=16-((10+sDict.size())%16); 
  sDict.insert(sDict.end(),remainder,' ');
  sDict.back() = '\n';
  fprintf(stderr,"%s\n",sDict.c_str());
  if(UINT16_MAX<=sDict.size())
    return vHeader;
  uint16_t dictsz=sDict.size();
  vHeader.resize(10+dictsz);
  vHeader[0]=(char) 0x93;
  memcpy(&(vHeader[1]),"NUMPY",5);
  vHeader[6]=(char) 0x01; //major version of numpy format
  vHeader[7]=(char) 0x00; //minor version of numpy format
  memcpy(&(vHeader[8]),&dictsz,sizeof(uint16_t));
  memcpy(&(vHeader[10]),sDict.c_str(),dictsz);
  
  return vHeader;
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template<typename T>
bool npy_save(const std::string& crsFileName, const std::vector<T>& crvData)
  {
  std::vector<size_t> vShape(1,crvData.size());
  std::vector<char> vHeader=create_npy_header<T>(vShape);
  if(vHeader.empty())
    return false;
  FILE* pFile=fopen(crsFileName.c_str(),"wb");
  if(NULL==pFile)
    return false;
  //size_t nels = std::accumulate(shape.begin(),shape.end(),1,std::multiplies<size_t>());
  fwrite(&vHeader[0],sizeof(char),vHeader.size(),pFile);
  fwrite(&crvData[0],sizeof(T),crvData.size(),pFile);
  fclose(pFile);
  return true;
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template<typename T>
bool npy_save(const std::string& crsFileName, const std::vector<std::vector<T> >& crvvData)
  {
  std::vector<size_t> vShape(2,crvvData.size());
  if(!crvvData.empty())
    vShape[1]=crvvData[0].size();
  else
    return false;
  std::vector<char> vHeader=create_npy_header<T>(vShape);
  if(vHeader.empty())
    return false;
  FILE* pFile=fopen(crsFileName.c_str(),"wb");
  if(NULL==pFile)
    return false;
  fwrite(&vHeader[0],sizeof(char),vHeader.size(),pFile);
  for(size_t i=0;i<crvvData.size();++i)
    fwrite(&crvvData[i][0],sizeof(T),crvvData[i].size(),pFile);
  fclose(pFile);
  return true;
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template<typename T>
bool npy_save(const std::string& crsFileName, const std::vector<std::vector<std::vector<T> > >& crvvvData)
  {
  std::vector<size_t> vShape(3,crvvvData.size());
  if(!crvvvData.empty() && !crvvvData[0].empty() )
    {
    vShape[1]=crvvvData[0].size();
    vShape[2]=crvvvData[0][0].size();
    }
  else
    return false;
  std::vector<char> vHeader=create_npy_header<T>(vShape);
  if(vHeader.empty())
    return false;
  FILE* pFile=fopen(crsFileName.c_str(),"wb");
  if(NULL==pFile)
    return false;
  fwrite(&vHeader[0],sizeof(char),vHeader.size(),pFile);
  for(size_t i=0;i<crvvvData.size();++i)
    for(size_t j=0;j<crvvvData[i].size();++j)
      fwrite(&crvvvData[i][j][0],sizeof(T),crvvvData[i][j].size(),pFile);
  fclose(pFile);
  return true;
}
#endif
