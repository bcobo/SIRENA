
#include "scheduler.h"

#include "log.h"

#include "threadsafe_queue.h"
#include "tasksSIRENA.h"

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
    if(detection_queue.wait_and_pop(data)){
      log_trace("Data extracted from queue: %i", data.rec_init.pulse_length);
      //th_runDetect(data.rec,data.n_record,data.last_record,
      //             data.all_pulses,&data.rec_init,&data.record_pulses);
    }
    std::unique_lock<std::mutex> lk(end_workers_mut);
    if(end_workers){
      log_trace("Finishing worker");
      lk.unlock();
      break;
    }
    lk.unlock();
  }
}

void scheduler::push_detection(TesRecord* record, 
                               int nRecord, 
                               int lastRecord, 
                               PulsesCollection *pulsesAll, 
                               ReconstructInitSIRENA** reconstruct_init, 
                               PulsesCollection** pulsesInRecord)
{
  detection_input input;
  input.rec = record;
  input.n_record = nRecord;
  input.last_record = lastRecord;
  input.all_pulses = pulsesAll;
  input.record_pulses = *pulsesInRecord;
  input.rec_init = **reconstruct_init;
  this->push_detection(input);
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

void scheduler::end_detection()
{
  
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

detection::detection():
  n_record(0),
  last_record(0),
  all_pulses(0),
  record_pulses(0)
{
  
}

detection::detection(const detection& other):
  n_record(other.n_record),
  last_record(other.last_record),
  all_pulses(0),
  record_pulses(0),
  rec(other.rec),
  rec_init(other.rec_init)
{

}

detection& detection::operator=(const detection& other)
{
  if(this != &other){
    n_record = other.n_record;
    last_record = other.last_record;
    rec = other.rec;
    rec_init = other.rec_init;
  }
  return *this;
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
