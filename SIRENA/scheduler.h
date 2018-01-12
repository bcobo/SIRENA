
#ifndef SCHEDULER_H
#define SCHEDULER_H

#include <thread>
#include <future>
#include <vector>

#include "integraSIRENA.h"

#define detection_input ReconstructInitSIRENA
#define detection_output ReconstructInitSIRENA

#if 0
struct detection_input
{
  int pulse_length;//recontruct_init->pulse_length
  int eventsz;//record->trigger_size
  int tstart_pulse1;//reconstruct_init->tstartPulse1
  int tstart_pulse2;//reconstruct_init->tstartPulse2
  int tstart_pulse3;//reconstruct_init->tstartPulse3
  int intermediate;//reconstruct_init->intermediate
  
  
};

struct detection_output
{
};

struct energy_input
{
};

struct energy_output
{
};
#endif
class scheduler
{
 public:

  void push_detection(/*TODO: add parameters*/
                      const detection_input &input);
  void run_energy(/*TODO: add parameters*/);

  virtual ~scheduler();
  inline static scheduler* get()
  { 
    return instance ? instance : instance = new scheduler();
  }
 private:
  static scheduler* instance;
  scheduler();
  scheduler(const scheduler& copy){}
  scheduler& operator=(const scheduler&){}

  unsigned int num_cores;
  unsigned int max_detection_workers;
  unsigned int max_energy_workers;
  
  std::promise<bool>* detection_worker_status;
  std::promise<bool>* energy_worker_status;

  std::thread* detection_workers;
  std::thread* energy_workers;
};

#endif
