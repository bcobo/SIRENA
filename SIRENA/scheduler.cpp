
#include "scheduler.h"

#include "log.h"

#include "threadsafe_queue.h"
#include "tasksSIRENA.h"

std::mutex end_workers_mut;
std::mutex records_detected_mut;

threadsafe_queue<detection_input*> detection_queue;
threadsafe_queue<detection_input*> energy_queue;

scheduler* scheduler::instance = 0;

bool end_workers = false;

unsigned int records_detected = 0;

/* ****************************************************************************/
/* Workers ********************************************************************/
/* ****************************************************************************/
void detection_worker()
{
  log_trace("Starting detection worker");
  while(1){
    detection_input* data;
    if(detection_queue.wait_and_pop(data)){
      log_trace("Data extracted from queue: %i", data->rec_init->pulse_length);
      log_trace("Data record from queue: %f", data->rec->time);
      //ReconstructInitSIRENA* aux = &data.rec_init;
      th_runDetect(data->rec->get_TesRecord(),
                   data->n_record,data->last_record,
                   data->all_pulses,
                   &(data->rec_init),
                   &(data->record_pulses));
      energy_queue.push(data);
      std::unique_lock<std::mutex> lk(records_detected_mut);
      ++records_detected;
      lk.unlock();
    }
    std::unique_lock<std::mutex> lk(end_workers_mut);
    if(end_workers){
      log_trace("Finishing detection worker");
      lk.unlock();
      break;
    }
    lk.unlock();
  }
}
#if 0
void reconstruction_worker()
{
  log_trace("reconstruction worker");
  while(1){
    detection_input* data;
    if(energy_queue.wait_and_pop(data)){
      log_trace("Data extracted from energy queue: %i", 
                data->rec_init->pulse_length);
#if 0
      if(strcmp(data.rec_init.EnergyMethod, "PCA") != 0){
        ReconstructInitSIRENA* aux = &data.rec_init;
        runEnergy(data.rec, &aux, &data.record_pulses, &data.optimal_filter);
      }
      if(data.n_record == 1){
        
      }
#endif
    }
    std::unique_lock<std::mutex> lk(end_workers_mut);
    if(end_workers && energy_queue.empty()){
      log_trace("Finishing reconstruction worker");
      lk.unlock();
      break;
    }
    lk.unlock();
  }
}
#endif
void scheduler::push_detection(TesRecord* record, 
                               int nRecord, 
                               int lastRecord, 
                               PulsesCollection *pulsesAll, 
                               ReconstructInitSIRENA** reconstruct_init, 
                               PulsesCollection** pulsesInRecord,
                               OptimalFilterSIRENA** optimal,
                               TesEventList* event_list)
{
  detection_input* input = new detection_input;
  input->rec = new tesrecord(record);
  input->n_record = nRecord;
  input->last_record = lastRecord;
  input->all_pulses = pulsesAll;//TODO
  input->record_pulses = *pulsesInRecord;//TODO
  input->rec_init = *reconstruct_init;
  input->optimal_filter = *optimal;
  input->event_list = event_list;
  detection_queue.push(input);
  ++num_records;
  //this->push_detection(input);
}

void scheduler::push_detection(const detection_input &input)
{ 
#if 0
  log_trace("pushing input");
  //std::this_thread::sleep_for(std::chrono::milliseconds(900));
  log_trace("pushing param");
  detection_queue.push(input);
  log_trace("end");
#endif
}

void scheduler::finish_reconstruction(PulsesCollection** pulsesAll, 
                                      OptimalFilterSIRENA** optimalFilter)
{
  // Waits until all the records are detected
  // this works because this function should only be called
  // after all the records are queue
  log_trace("Waiting until all the workers finish");
  while(1){
    std::unique_lock<std::mutex> lk(records_detected_mut);
    if (records_detected == this->num_records){
      lk.unlock();
      break;
    }
    lk.unlock();
  }
  
  // Sortint the arrays by record number
  if ((*pulsesAll)->pulses_detected){
    delete [] (*pulsesAll)->pulses_detected;
  }
  (*pulsesAll)->size = this->num_records;
  (*pulsesAll)->ndetpulses = this->num_records;
  (*pulsesAll)->pulses_detected = new PulseDetected[this->num_records];
  detection_input** data_array = new detection_input*[this->num_records];
  while(!energy_queue.empty()){
    detection_input* data;
    if(energy_queue.wait_and_pop(data)){
      data_array[data->n_record] = data;
    }
  }
  
  for (unsigned int i = 0; i < this->num_records; ++i){
    log_trace("Pulses sorted %i", data_array[i]->n_record);
  }

  // Energy
#if 0
  while(!energy_queue.empty()){
    detection_input* data;
    if(energy_queue.wait_and_pop(data)){
      log_trace("Data extracted from energy queue: %i", 
                data->rec_init->pulse_length);
      //#if 0
      if(strcmp(data.rec_init.EnergyMethod, "PCA") != 0){
        ReconstructInitSIRENA* aux = &data.rec_init;
        runEnergy(data.rec, &aux, &data.record_pulses, &data.optimal_filter);
      }
      if(data.n_record == 1){
        
      }
      //#endif
    }
  }
#endif
}

scheduler::scheduler():
  threading(false),
  num_cores(0),
  max_detection_workers(0),
  max_energy_workers(0),
  num_records(0)
{
  if(threading){
    this->num_cores = std::thread::hardware_concurrency();
    this->max_detection_workers = this->num_cores - 2;
    this->detection_workers = new std::thread[this->max_detection_workers];
    log_trace("Num of cores %u", this->num_cores);
    for (unsigned int i = 0; i < this->max_detection_workers; ++i){
      this->detection_workers[i] = std::thread (detection_worker);
    }
    //this->reconstruct_worker = std::thread(reconstruction_worker);
  }
}

scheduler::~scheduler()
{
  if(threading){
    log_trace("scheduler destructor");
    instance = 0;
    std::unique_lock<std::mutex> lk(end_workers_mut);
    end_workers = true;
    lk.unlock();
    //TODO: clean the threads
    for(unsigned int i = 0; i < this->max_detection_workers; ++i){
      this->detection_workers[i].join();
    }
    //this->reconstruct_worker.join();
  }
}

/* ****************************************************************************/
/* Data structures implementation *********************************************/
/* ****************************************************************************/

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
  if (phid_array){
    delete [] phid_array; phid_array = 0;
  }
  if (times){
    delete [] times; times = 0;
  }
}

PhIDList* phidlist::get_PhIDList() const
{
  PhIDList* ret = new PhIDList;
  ret->wait_list = wait_list;
  ret->n_elements = n_elements;
  ret->index = index;
  ret->size = size;
  if (phid_array && size > 0){
    ret->phid_array = new long[size];
    for (unsigned int i = 0; i < size; ++i){
      ret->phid_array[i] = phid_array[i];
    }
  }
  if (times && size > 0){
    ret->times = new double[size];
    for (unsigned int i = 0; i < size; ++i){
      ret->times[i] = times[i];
    }
  }
  return ret;
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
  if (adc_array){
    delete [] adc_array; adc_array = 0;
  }
  if (adc_double){
    delete [] adc_double; adc_double = 0;
  }
  if (phid_list){
    delete phid_list; phid_list = 0;
  }
}

TesRecord* tesrecord::get_TesRecord() const
{
  TesRecord* ret = new TesRecord;
  ret->trigger_size = trigger_size;
  ret->time = time;
  ret->delta_t = delta_t;
  ret->pixid = pixid;

  if(adc_array && trigger_size > 0){
    ret->adc_array = new uint16_t[trigger_size];
    for (unsigned int i = 0; i < trigger_size; ++i){
      ret->adc_array[i] = adc_array[i];
    }
  }
  if(adc_double && trigger_size > 0){
    ret->adc_double = new double[trigger_size];
    for (unsigned int i = 0; i < trigger_size; ++i){
        ret->adc_double[i] = adc_double[i];
    }
  }
  if(phid_list && phid_list->size > 0){
    ret->phid_list = phid_list->get_PhIDList();
  }
  return ret;
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

energy::energy():
  record_pulses(0),
  optimal_filter(0)
{
  //TODO:
}

energy::energy(const energy& other):
  rec_init(other.rec_init),
  record_pulses(0),
  optimal_filter(0)
{
  //TODO:
}

energy& energy::operator=(const energy& other)
{
  if(this != &other){
    rec_init = other.rec_init;
    //TODO:
  }
  return *this;
}

energy::~energy()
{
  //TODO:
}
