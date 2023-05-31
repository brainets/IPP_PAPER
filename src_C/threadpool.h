#ifndef __THREAD_POOL_T_H
#define __THREAD_POOL_T_H
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <queue>
#include <functional>

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class thread_pool_t
{
public:
  thread_pool_t(size_t nThreads);
  virtual ~thread_pool_t();
  void run_task(std::function<void(void)> fun);
  void wait(void);
  size_t get_waiting(void) const { return m_queue.size(); }
  int get_pending(void) const { return m_pending.load(); }
  size_t get_pool_size(void) const { return m_vPool.size(); }
protected:
  void thread_entry(size_t i);
private:
  bool m_bShutdown;
  std::atomic<int> m_pending;
  std::mutex m_mutex;
  std::condition_variable m_cond;
  std::mutex m_busymutex;
  std::condition_variable m_busycond;
  std::queue<std::function<void(void)> > m_queue;
  std::vector<std::thread> m_vPool;
};
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif
