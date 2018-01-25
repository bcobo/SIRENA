
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

/* Data structures implementation *********************************************/

phidlist::phidlist():
  phid_array(0),
  times(0),
  wait_list(0),
  n_elements(0),
  index(0),
  size(0)
{
  
}

phidlist::phidlist(const phidlist& other):
  phid_array(0),
  times(0),
  wait_list(other.wait_list),
  n_elements(other.n_elements),
  index(other.index),
  size(other.size)
{
  if (other.phid_array && other.size > 0){
    phid_array = new long[other.size];
    for (unsigned int i = 0; i < size; ++i){
      phid_array[i] = other.phid_array[i];
    }
  }
  if (other.times && other.size > 0){
    times = new double[other.size];
    for (unsigned int i = 0; i < size; ++i){
      times[i] = other.times[i];
    }
  }
}

phidlist::phidlist(PhIDList* other):
  phid_array(0),
  times(0),
  wait_list(other->wait_list),
  n_elements(other->n_elements),
  index(other->index),
  size(other->size)
{
  if (other->phid_array && size > 0){
    phid_array = new long[size];
    for (unsigned int i = 0; i < size; ++i){
      phid_array[i] = other->phid_array[i];
    }
  }
  if (other->times && size > 0){
    times = new double[size];
    for (unsigned int i = 0; i < size; ++i){
      times[i] = other->times[i];
    }
  }
  
}

phidlist& phidlist::operator=(const phidlist& other)
{
  if(this != &other){
    wait_list = other.wait_list;
    n_elements = other.n_elements;
    index = other.index;
    size = other.size;
    if(phid_array){
      delete [] phid_array; phid_array = 0;
    }
    if (times){
      delete [] times; times = 0;
    }
    if (other.phid_array && other.size > 0){
      phid_array = new long[other.size];
      for (unsigned int i = 0; i < size; ++i){
        phid_array[i] = other.phid_array[i];
      }
    }
    if (other.times && other.size > 0){
      times = new double[other.size];
      for (unsigned int i = 0; i < size; ++i){
        times[i] = other.times[i];
      }
    }
  }
  return *this;
}

phidlist::~phidlist()
{
  //TODO:
}

tesrecord::tesrecord():
  trigger_size(0),
  time(0),
  delta_t(0),
  adc_array(0),
  adc_double(0),
  pixid(0),
  phid_list(0)
{
  
}

tesrecord::tesrecord(TesRecord* other):
  trigger_size(other->trigger_size),
  time(other->time),
  delta_t(other->delta_t),
  adc_array(0),//trigger_size
  adc_double(0),//trigger_size
  pixid(other->pixid),
  phid_list(0)//MAXIMPACTNUMBER
{
  if(other->adc_array && trigger_size > 0){
    adc_array = new uint16_t[trigger_size];
    for (unsigned int i = 0; i < trigger_size; ++i){
      adc_array[i] = other->adc_array[i];
    }
  }
  if(other->adc_double && trigger_size > 0){
    adc_double = new double[trigger_size];
    for (unsigned int i = 0; i < trigger_size; ++i){
      adc_double[i] = other->adc_double[i];
    }
  }
  if(other->phid_list && other->phid_list->size > 0){
    phid_list = new phidlist(other->phid_list);
  }
}

tesrecord::tesrecord(const tesrecord& other):
  trigger_size(other.trigger_size),
  time(other.time),
  delta_t(other.delta_t),
  adc_array(0),
  adc_double(0),
  pixid(other.pixid),
  phid_list(0)
{
  if(other.adc_array && trigger_size > 0){
    adc_array = new uint16_t[trigger_size];
    for (unsigned int i = 0; i < trigger_size; ++i){
      adc_array[i] = other.adc_array[i];
    }
  }
  if(other.adc_double && trigger_size > 0){
    adc_double = new double[trigger_size];
    for (unsigned int i = 0; i < trigger_size; ++i){
      adc_double[i] = other.adc_double[i];
    }
  }
  if(other.phid_list && other.phid_list->size > 0){
    phid_list = new phidlist(*other.phid_list);
  }
}

tesrecord& tesrecord::operator=(const tesrecord& other)
{
  if(this != &other){
    trigger_size = other.trigger_size;
    time = other.time;
    delta_t = other.delta_t;
    pixid = other.pixid;
    
    if(adc_array){
      delete [] adc_array; adc_array = 0;
    }

    if(adc_double){
      delete [] adc_double; adc_double = 0;
    }

    if(phid_list){
      delete phid_list; phid_list = 0;
    }

    if(other.adc_array && trigger_size > 0){
      adc_array = new uint16_t[trigger_size];
      for (unsigned int i = 0; i < trigger_size; ++i){
        adc_array[i] = other.adc_array[i];
      }
    }
    if(other.adc_double && trigger_size > 0){
      adc_double = new double[trigger_size];
      for (unsigned int i = 0; i < trigger_size; ++i){
        adc_double[i] = other.adc_double[i];
      }
    }
    if(other.phid_list && other.phid_list->size > 0){
      phid_list = new phidlist(*other.phid_list);
    }
  }
  return *this;
}

tesrecord::~tesrecord()
{
  //TODO:
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
  //TODO:
}

detection& detection::operator=(const detection& other)
{
  if(this != &other){
    n_record = other.n_record;
    last_record = other.last_record;
    rec = other.rec;
    rec_init = other.rec_init;
    //TODO:
  }
  return *this;
}

detection::~detection()
{
  //TODO:
}

energy::energy()
{
  //TODO:
}

energy::energy(const energy& other)
{
  //TODO:
}

energy& energy::operator=(const energy& other)
{
  //TODO:
}

energy::~energy()
{
  //TODO:
}
