#include <cstdio>
#include <chrono>
#include "threadpool.h"
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
thread_pool_t::thread_pool_t(size_t nThreads) : m_bShutdown(false),m_pending(0)
{
  m_vPool.reserve(nThreads);
  for(size_t i=0;i<nThreads;++i)
    m_vPool.emplace_back(std::bind(&thread_pool_t::thread_entry,this,i));
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
thread_pool_t::~thread_pool_t()
{
  { // Unblock any threads and tell them to stop     
    std::unique_lock<std::mutex> lock(m_mutex);      
    m_bShutdown=true;
    m_cond.notify_all();
  }
  // Wait for all threads to stop
  fprintf(stderr,"Joining threads for shutdown...");
  for(size_t i=0;i<m_vPool.size();++i)
    m_vPool[i].join();
  //fprintf(stderr,"done, bye.\n");
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void thread_pool_t::run_task(std::function<void(void)> fun)
{   
  std::unique_lock<std::mutex> lock(m_mutex);    
  m_queue.emplace(std::move(fun));
  m_pending++;
  m_cond.notify_one();
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void thread_pool_t::wait(void)
{
  std::unique_lock<std::mutex> lock(m_busymutex);
  int np=1;
  while(0<np)
    {    
    std::unique_lock<std::mutex> lock(m_mutex);    
    np=m_pending.load();
    if(0<np)
      {
      //fprintf(stderr,"Still have %d pending tasks and queue size %zu, waiting...\n",np,m_queue.size());
      m_busycond.wait(lock);
      }
    }
  //fprintf(stderr,"done, bye.\n");
}  
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void thread_pool_t::thread_entry(size_t i)
{
  std::function <void(void)> task;    
  while(true)
    {
      {
      std::unique_lock<std::mutex> lock(m_mutex);	  
      while(!m_bShutdown && m_queue.empty())
	m_cond.wait(lock);
      if(m_queue.empty())
	{ // No tasks left and we are shutting down.
	fprintf(stderr,"Thread %zu terminating (shutdown is %s).\n",i,m_bShutdown?"TRUE":"FALSE");
	return;
	}
      //fprintf(stderr,"Thread %zu running task.\n",i);
      task=std::move(m_queue.front());
      m_queue.pop();
      //if(--m_pending==0)
      //m_busycond.notify_one();
      }	
      task();
      //int np=-1;
      {
      std::unique_lock<std::mutex> lock(m_mutex);
      m_pending--;
      //np=m_pending.load();
      }
      int np=m_pending.load();
      if(0==np)
	{
	//fprintf(stderr,"We have %d pending tasks and queue size %zu, signaling...\n",np,m_queue.size());
	m_busycond.notify_all();
	}
    }
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// compile with: c++ -Wall -O6 -std=c++11 -pthread -D__TESTMAIN -o threadpool threadpool.cpp
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#if defined(__TESTMAIN)
#include <random>
static std::random_device rd;
static std::mt19937 gen(rd());

void testsimple(int n)
{
  fprintf(stderr,"Sleeping for %ds\n",n);
  std::this_thread::sleep_for(std::chrono::seconds(n));
}
struct task_t
{
  int n;
};
void testtask(task_t* pt)
{
  fprintf(stderr,"Sleeping for %dms\n",pt->n);
  std::this_thread::sleep_for(std::chrono::milliseconds(pt->n));
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int main()
{
  thread_pool_t p(2);
  /*p.run_task(std::bind(testsimple,1));
  p.run_task(std::bind(testsimple,2));
  p.run_task(std::bind(testsimple,3));
  p.run_task(std::bind(testsimple,4));*/

  size_t sz=100;
  std::uniform_int_distribution<uint8_t> unidist(0,255);
  std::vector<task_t> vtsk(sz);
  for(size_t i=0;i<sz;++i)
    vtsk[i].n=unidist(gen);
  for(size_t i=0;i<sz;++i)
    p.run_task(std::bind(testtask,&vtsk[i]));
}
#endif
