
#include "scheduler.h"

#include "log.h"

#include "threadsafe_queue.h"

std::mutex end_workers_mut;

threadsafe_queue<detection_input> detection_queue;
threadsafe_queue<detection_output> energy_queue;

scheduler* scheduler::instance = 0;

bool end_workers = false;

void run_detect()
{
#if 0
  if(tstart pulses > trigger size) error;
  
  create intermediate fits if required;
  calls createDetectFile();
  
  filderLibrary();
  
  store the input record in invector;
  calls loadRecord();

  convert I into R if energyMethod = I2R;
  calls convertI2R();

  procRecord();

  if intermediate fits write keywords;

  if lastrecord and pca;
  calls weightMatrix();
  calls eigenVV();

  Polyfitlinear();

#endif
}

void detection_worker()
{
  log_trace("Starting worker");
  while(1){
    detection_input data;
    if(detection_queue.wait_and_pop(data))
      log_trace("Data extracted from queue: %i", data.pulse_length);
    std::unique_lock<std::mutex> lk(end_workers_mut);
    if(end_workers){
      log_trace("Finishing worker");
      lk.unlock();
      break;
    }
    lk.unlock();
  }
}

void scheduler::push_detection(const detection_input &input)
{ 
  log_trace("pushing input");
  /*detection_input d1, d2, d3, d4;
  d1.pulse_length = 1;
  d2.pulse_length = 20;
  d3.pulse_length = 5;
  d4.pulse_length = 12;
  detection_queue.push(d1);
  std::this_thread::sleep_for(std::chrono::milliseconds(900));
  detection_queue.push(d2);
  std::this_thread::sleep_for(std::chrono::milliseconds(900));
  detection_queue.push(d3);
  std::this_thread::sleep_for(std::chrono::milliseconds(900));
  detection_queue.push(d4);*/
  log_trace("pushing param");
  detection_queue.push(input);
  log_trace("end");
}

void scheduler::run_energy()
{
  
}

scheduler::scheduler():
  num_cores(0),
  max_detection_workers(0),
  max_energy_workers(0)
{
  this->num_cores = std::thread::hardware_concurrency();
  this->max_detection_workers = this->num_cores - 1;
  this->detection_workers = new std::thread[this->max_detection_workers];
  log_trace("Num of cores %u", this->num_cores);
  for (unsigned int i = 0; i < this->max_detection_workers; ++i){
    this->detection_workers[i] = std::thread (detection_worker);
  }
}

scheduler::~scheduler()
{
  log_trace("scheduler destructor");
  instance = 0;
  std::unique_lock<std::mutex> lk(end_workers_mut);
  end_workers = true;
  lk.unlock();
  //TODO: clean the threads
  for(unsigned int i = 0; i < this->max_detection_workers; ++i){
    this->detection_workers[i].join();
  }
}

detection::detection()
{
  
}

detection::detection(const detection& other)
{
  
}

detection& detection::operator=(const detection& other)
{
  
}

detection::~detection()
{
  
}

energy::energy()
{
  
}

energy::energy(const energy& other)
{
  
}

energy& energy::operator=(const energy& other)
{
  
}

energy::~energy()
{
  
}
