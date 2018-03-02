
#include "scheduler.h"

#include "log.h"

#include "threadsafe_queue.h"
#include "tasksSIRENA.h"
//#include "teseventlist.h"
//#include "testriggerfile.h"

std::mutex end_workers_mut;
std::mutex records_detected_mut;
std::mutex records_energy_mut;

threadsafe_queue<detection_input*> detection_queue;
threadsafe_queue<detection_input*> detected_queue;
threadsafe_queue<detection_input*> energy_queue;
threadsafe_queue<detection_input*> end_queue;

scheduler* scheduler::instance = 0;

bool end_workers = false;

unsigned int records_detected = 0;
unsigned int records_energy = 0;

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
      th_runDetect(data->rec,
                   data->n_record,data->last_record,
                   data->all_pulses,
                   &(data->rec_init),
                   &(data->record_pulses));
      detected_queue.push(data);
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

void energy_worker()
{
  log_trace("Starting energy worker");
  while(1){
    detection_input* data;
    if(energy_queue.wait_and_pop(data)){
#if 0
      log_trace("run energy");
      log_trace("energy data:");
      log_trace("%f", data->rec_init->pulse_length);
      log_trace("%f", data->rec->time);
      log_trace("%i", data->record_pulses->ndetpulses);
      log_trace("%p",&data->optimal_filter);
      log_trace("%p",&data->record_pulses);
      log_trace("%p",&data->rec_init);
      log_trace("trigger %i",data->rec->trigger_size);
      log_trace("phdlist %i", data->rec->phid_list->size);
      //log_trace("%p", data->rec->get_TesRecord());
#endif
      
      th_runEnergy(data->rec, 
                   &(data->rec_init),
                   &(data->record_pulses),//copy
                   &(data->optimal_filter));
      
      log_trace("end run energy");
      end_queue.push(data);
      std::unique_lock<std::mutex> lk(records_energy_mut);
      ++records_energy;
      lk.unlock();
    }
    std::unique_lock<std::mutex> lk(end_workers_mut);
    if(end_workers){
      log_trace("Finishing energy worker");
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

/* ****************************************************************************/
/* Scheduler ******************************************************************/
/* ****************************************************************************/
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
  tesrecord* rec = new tesrecord(record);
  input->rec = rec->get_TesRecord();
  input->n_record = nRecord;
  input->last_record = lastRecord;
  input->all_pulses = new PulsesCollection;
  if (pulsesAll and pulsesAll->ndetpulses > 0){
    *input->all_pulses = *pulsesAll;
  }
  //input->all_pulses = pulsesAll;//TODO copy?
  input->record_pulses = *pulsesInRecord;
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

void scheduler::finish_reconstruction(ReconstructInitSIRENA* reconstruct_init,
                                      PulsesCollection** pulsesAll, 
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
  
  (*pulsesAll)->size = 10;//this->num_records;
  (*pulsesAll)->ndetpulses = 0;
  (*pulsesAll)->pulses_detected = new PulseDetected[10];//[this->num_records];
  detection_input** data_array = new detection_input*[this->num_records+1];
  while(!detected_queue.empty()){
    detection_input* data;
    if(detected_queue.wait_and_pop(data)){
      data_array[data->n_record] = data;
      //log_trace("trigger %i",data->rec->trigger_size);
      //log_trace("phdlist %i", data->rec->phid_list->size);
      //log_debug("allPulses %i", data->all_pulses->ndetpulses);
    }
  }

  std::unique_lock<std::mutex> lk(end_workers_mut);
  end_workers = true;
  lk.unlock();
  for(unsigned int i = 0; i < this->max_detection_workers; ++i){
    this->detection_workers[i].join();
  }

  //
  // Energy
  //
  if(this->is_running_energy){
    std::unique_lock<std::mutex> lk_end(end_workers_mut);
    end_workers = false;
    lk_end.unlock();
    this->energy_workers = new std::thread[this->max_detection_workers];
    for (unsigned int i = 0; i < 1; ++i){//this->max_detection_workers; ++i){//
      this->energy_workers[i] = std::thread (energy_worker);
    }
    
    for (unsigned int i = 1; i <= this->num_records; ++i){
      energy_queue.push(data_array[i]);
    }
  
    // Waits until all the energies are calculated
    log_trace("Waiting until all the energy workers finish");
    while(1){
      std::unique_lock<std::mutex> lk_energy(records_energy_mut);
      if (records_energy == this->num_records){
        lk_energy.unlock();
        break;
      }
      lk_energy.unlock();
    }
    
    std::unique_lock<std::mutex> lk_end2(end_workers_mut);
    end_workers = true;
    lk_end2.unlock();
    for(unsigned int i = 0; i < 1; ++i){//this->max_detection_workers; ++i){
      this->energy_workers[i].join();
    }
  }//end energy

  //
  // Reconstruction of the pulses array
  //
  for (unsigned int i = 1; i <= this->num_records; ++i){
    /*log_trace("Pulses sorted %i all pulses %i record pulses %i", 
      data_array[i]->n_record,
              data_array[i]->all_pulses->ndetpulses,
              data_array[i]->record_pulses->ndetpulses);*/
    PulsesCollection* all = data_array[i]->all_pulses;
    PulsesCollection* in_record = data_array[i]->record_pulses;
    PulsesCollection aux;
    //log_trace("Size: all %i",(*pulsesAll)->ndetpulses); 
    //log_trace("detected %i",  in_record->ndetpulses);
    if ((*pulsesAll)->size < 
        ((*pulsesAll)->ndetpulses + in_record->ndetpulses)){
      aux = *(*pulsesAll);
      delete *pulsesAll;
      *pulsesAll = new PulsesCollection;
      //**pulsesAll = aux;
      (*pulsesAll)->ndetpulses = aux.ndetpulses;
      (*pulsesAll)->size = (*pulsesAll)->ndetpulses + in_record->ndetpulses;
      (*pulsesAll)->pulses_detected =  new PulseDetected[(*pulsesAll)->size];
      for (unsigned int i = 0; i < aux.ndetpulses; ++i){
        (*pulsesAll)->pulses_detected[i] = aux.pulses_detected[i];
      }
    }
    for (unsigned int i = 0; i < in_record->ndetpulses; ++i){
      (*pulsesAll)->pulses_detected[i+(*pulsesAll)->ndetpulses] = 
        in_record->pulses_detected[i];
    }
    (*pulsesAll)->ndetpulses += in_record->ndetpulses;
  }// End reconstruction of the pulses array

  for (unsigned int i = 1; i < (*pulsesAll)->ndetpulses; ++i){
    log_trace("PulsesAll %d ", (*pulsesAll)->pulses_detected[i].quality);
  }// end reconstruction of the pulses array

  /*
  for (unsigned int i = 1; i <= this->num_records; ++i){
    log_trace("Pulses sorted %i all pulses %i record pulses %i", 
              data_array[i]->n_record,
              data_array[i]->all_pulses->ndetpulses,
              data_array[i]->record_pulses->ndetpulses);
  }
  */

  for(unsigned int i = 1; i <= this->num_records; ++i){
    TesEventList* event_list = data_array[i]->event_list;
    if (event_list->energies != 0) delete [] event_list->energies;
    if (event_list->avgs_4samplesDerivative != 0) delete [] event_list->avgs_4samplesDerivative;
    if (event_list->grades1 != 0) delete [] event_list->grades1;
    if (event_list->grades2 != 0) delete [] event_list->grades2;
    if (event_list->pulse_heights != 0) delete [] event_list->pulse_heights;
    if (event_list->ph_ids != 0) delete [] event_list->ph_ids;

    if (strcmp(data_array[i]->rec_init->EnergyMethod,"PCA") != 0){
      event_list->index = data_array[i]->record_pulses->ndetpulses;
      event_list->energies = new double[event_list->index];
      event_list->avgs_4samplesDerivative = new double[event_list->index];
      event_list->grades1  = new int[event_list->index];
      event_list->grades2  = new int[event_list->index];
      event_list->pulse_heights  = new double[event_list->index];
      event_list->ph_ids   = new long[event_list->index];

      for (int ip=0; ip<data_array[i]->record_pulses->ndetpulses; ip++) {
        event_list->event_indexes[ip] = 
          (data_array[i]->record_pulses->pulses_detected[ip].Tstart 
           - data_array[i]->rec->time)/data_array[i]->rec->delta_t;
        
        event_list->energies[ip] = data_array[i]->record_pulses->pulses_detected[ip].energy;
        
        event_list->avgs_4samplesDerivative[ip] = 
          data_array[i]->record_pulses->pulses_detected[ip].avg_4samplesDerivative;
        event_list->grades1[ip]  = data_array[i]->record_pulses->pulses_detected[ip].grade1;
        event_list->grades2[ip]  = data_array[i]->record_pulses->pulses_detected[ip].grade2;
        event_list->pulse_heights[ip]  = data_array[i]->record_pulses->pulses_detected[ip].pulse_height;
        event_list->ph_ids[ip]   = 0;
      }
      if (data_array[i]->last_record == 1) {       
        double numLagsUsed_mean;
        double numLagsUsed_sigma;
        gsl_vector *numLagsUsed_vector = gsl_vector_alloc((*pulsesAll)->ndetpulses);
        
        for (int ip=0; ip<(*pulsesAll)->ndetpulses; ip++) {
          gsl_vector_set(numLagsUsed_vector,ip,(*pulsesAll)->pulses_detected[ip].numLagsUsed);
        }
        if (findMeanSigma (numLagsUsed_vector, &numLagsUsed_mean, &numLagsUsed_sigma)) {
          EP_EXIT_ERROR("Cannot run findMeanSigma routine for calculating numLagsUsed statistics",EPFAIL);
        }
        gsl_vector_free(numLagsUsed_vector);
      }
    }else{
      if (data_array[i]->last_record == 1) {
        // Free & Fill TesEventList structure
        event_list->index = (*pulsesAll)->ndetpulses;
        event_list->event_indexes = new double[event_list->index];
        event_list->energies = new double[event_list->index];
        event_list->avgs_4samplesDerivative = new double[event_list->index];
        event_list->grades1  = new int[event_list->index];
        event_list->grades2  = new int[event_list->index];
        event_list->pulse_heights  = new double[event_list->index];
        event_list->ph_ids   = new long[event_list->index];
        
        for (int ip=0; ip<(*pulsesAll)->ndetpulses; ip++) {
          event_list->event_indexes[ip] = ((*pulsesAll)->pulses_detected[ip].Tstart 
                                           - data_array[i]->rec->time)/data_array[i]->rec->delta_t;
          
          event_list->energies[ip] = (*pulsesAll)->pulses_detected[ip].energy;
          
          event_list->avgs_4samplesDerivative[ip]  = (*pulsesAll)->pulses_detected[ip].avg_4samplesDerivative;
          event_list->grades1[ip]  = (*pulsesAll)->pulses_detected[ip].grade1;
          event_list->grades2[ip]  = (*pulsesAll)->pulses_detected[ip].grade2;
          event_list->pulse_heights[ip]  = (*pulsesAll)->pulses_detected[ip].pulse_height;
          event_list->ph_ids[ip]   = 0;    
        }
      }
    }
    int status=EXIT_SUCCESS;
    saveEventListToFile(this->outfile,
                        event_list,
                        data_array[i]->rec->time,
                        this->record_file_delta_t,
                        data_array[i]->rec->pixid,
                        &status);
    //CHECK_STATUS_BREAK(status);
    if (EXIT_SUCCESS!=status){
      EP_EXIT_ERROR("Something went wrong while saving the event_list",EPFAIL);
    }
  }// for event_list
#if 0
  // TODO: construct pulsesAll

  // TODO: here we're just saving the event_list info
  for(unsigned int i = 0; i < this->num_records; ++i){
    
    nRecord = data[i+1]->n_record;
    pullsesAll = data[i+1]->all_pulses;
    pulsesInRecord = data[i+1]->record_pulses;
    event_list = data[i+1]->event_list;
    
    if (nRecord == 1){
      (*pulsesAll)->ndetpulses = pulsesInRecord->ndetpulses;
      if((*pulsesAll)->pulses_detected != 0 && (*pulsesAll)->size < pulsesInRecord->ndetpulses){
        delete [] (*pulsesAll)->pulses_detected; (*pulsesAll)->pulses_detected = 0;
        (*pulsesAll)->size = resize_array((*pulsesAll)->size, (*pulsesAll)->ndetpulses);
        (*pulsesAll)->pulses_detected = new PulseDetected[(*pulsesAll)->size];
      }
      
    
#ifndef POOLS
      if((*pulsesAll)->pulses_detected == 0){
      (*pulsesAll)->pulses_detected = new PulseDetected[pulsesInRecord->ndetpulses];
      (*pulsesAll)->size = pulsesInRecord->ndetpulses;
    }
#endif
      for (int i=0;i<(*pulsesAll)->ndetpulses;i++){
        (*pulsesAll)->pulses_detected[i] = pulsesInRecord->pulses_detected[i];
      }            
    } else {
      if (event_list->energies != NULL) delete [] event_list->energies;
      if (event_list->avgs_4samplesDerivative != NULL) delete [] event_list->avgs_4samplesDerivative;
      if (event_list->grades1 != NULL) delete [] event_list->grades1;
      if (event_list->grades2 != NULL) delete [] event_list->grades2;
      if (event_list->pulse_heights != NULL) delete [] event_list->pulse_heights;
      if (event_list->ph_ids != NULL) delete [] event_list->ph_ids;
    
      pulsesAllAux->ndetpulses = (*pulsesAll)->ndetpulses;
      (*pulsesAll)->ndetpulses = (*pulsesAll)->ndetpulses + pulsesInRecord->ndetpulses;
      
               
      if ((*pulsesAll)->pulses_detected != NULL && (*pulsesAll)->size < (*pulsesAll)->ndetpulses){
        pulsesAllAux->pulses_detected = new PulseDetected[(*pulsesAll)->ndetpulses];
      
        for (int i=0;i<pulsesAllAux->ndetpulses;i++){
          pulsesAllAux->pulses_detected[i] = (*pulsesAll)->pulses_detected[i];
        }
      
        delete [] (*pulsesAll)->pulses_detected; (*pulsesAll)->pulses_detected = 0; 
        (*pulsesAll)->size = resize_array((*pulsesAll)->size, (*pulsesAll)->ndetpulses);     
        (*pulsesAll)->pulses_detected = new PulseDetected[(*pulsesAll)->size];
        
        for (int i=0;i<pulsesAllAux->ndetpulses;i++){
          (*pulsesAll)->pulses_detected[i] = pulsesAllAux->pulses_detected[i];
        }
        delete [] pulsesAllAux->pulses_detected; pulsesAllAux->pulses_detected = 0;
      }
                
#ifndef POOLS
      if((*pulsesAll)->pulses_detected == 0){
      (*pulsesAll)->pulses_detected = new PulseDetected[(*pulsesAll)->ndetpulses];
      (*pulsesAll)->size = (*pulsesAll)->ndetpulses;
    }
#endif
                
      // Save pulses detected in current record
      for (int i=0;i<pulsesInRecord->ndetpulses;i++) {
        (*pulsesAll)->pulses_detected[i+pulsesAllAux->ndetpulses] = pulsesInRecord->pulses_detected[i];
      }
  }
  
  event_list->index = pulsesInRecord->ndetpulses;
  event_list->energies = new double[event_list->index];
  event_list->avgs_4samplesDerivative = new double[event_list->index];
  event_list->grades1  = new int[event_list->index];
  event_list->grades2  = new int[event_list->index];
  event_list->pulse_heights  = new double[event_list->index];
  event_list->ph_ids   = new long[event_list->index];
  
  if (strcmp(reconstruct_init->EnergyMethod,"PCA") != 0)     // Different from PCA
    {
      for (int ip=0; ip<pulsesInRecord->ndetpulses; ip++) {            
        event_list->event_indexes[ip] = 
          (pulsesInRecord->pulses_detected[ip].Tstart - record->time)/record->delta_t;	
        
        if (reconstruct_init->mode == 1) {
          event_list->energies[ip] = pulsesInRecord->pulses_detected[ip].energy;
        }
        else if (reconstruct_init->mode == 0) {
          event_list->energies[ip] = 999.;
        }
        
        event_list->avgs_4samplesDerivative[ip] = pulsesInRecord->pulses_detected[ip].avg_4samplesDerivative;
        event_list->grades1[ip]  = pulsesInRecord->pulses_detected[ip].grade1;
        event_list->grades2[ip]  = pulsesInRecord->pulses_detected[ip].grade2;
        event_list->pulse_heights[ip]  = pulsesInRecord->pulses_detected[ip].pulse_height;
        event_list->ph_ids[ip]   = 0;
      }
      if (lastRecord == 1) {       
        double numLagsUsed_mean;
        double numLagsUsed_sigma;
        gsl_vector *numLagsUsed_vector = gsl_vector_alloc((*pulsesAll)->ndetpulses);
        
        for (int ip=0; ip<(*pulsesAll)->ndetpulses; ip++) {
          gsl_vector_set(numLagsUsed_vector,ip,(*pulsesAll)->pulses_detected[ip].numLagsUsed);
        }
        if (findMeanSigma (numLagsUsed_vector, &numLagsUsed_mean, &numLagsUsed_sigma)) {
          EP_EXIT_ERROR("Cannot run findMeanSigma routine for calculating numLagsUsed statistics",EPFAIL);
        }
        gsl_vector_free(numLagsUsed_vector);
      }
    } else {
    if (lastRecord == 1) {
      // Free & Fill TesEventList structure
      event_list->index = (*pulsesAll)->ndetpulses;
      event_list->event_indexes = new double[event_list->index];
      event_list->energies = new double[event_list->index];
      event_list->avgs_4samplesDerivative = new double[event_list->index];
      event_list->grades1  = new int[event_list->index];
      event_list->grades2  = new int[event_list->index];
      event_list->pulse_heights  = new double[event_list->index];
      event_list->ph_ids   = new long[event_list->index];
      
      for (int ip=0; ip<(*pulsesAll)->ndetpulses; ip++) {
        event_list->event_indexes[ip] = ((*pulsesAll)->pulses_detected[ip].Tstart - record->time)/record->delta_t;
        
        if (reconstruct_init->mode == 1) {
          event_list->energies[ip] = (*pulsesAll)->pulses_detected[ip].energy;
        }
        else if (reconstruct_init->mode == 0) {
          event_list->energies[ip] = 999.;
        }
        
        event_list->avgs_4samplesDerivative[ip]  = (*pulsesAll)->pulses_detected[ip].avg_4samplesDerivative;
        event_list->grades1[ip]  = (*pulsesAll)->pulses_detected[ip].grade1;
        event_list->grades2[ip]  = (*pulsesAll)->pulses_detected[ip].grade2;
        event_list->pulse_heights[ip]  = (*pulsesAll)->pulses_detected[ip].pulse_height;
        event_list->ph_ids[ip]   = 0;    
      }
    }
  }
  
  delete pulsesAllAux; pulsesAllAux = 0;
  delete [] pulsesInRecord->pulses_detected; pulsesInRecord->pulses_detected = 0;
  delete pulsesInRecord; pulsesInRecord = 0;
  }// records for
#endif
}

scheduler::scheduler():
  threading(true),
  num_cores(0),
  max_detection_workers(0),
  max_energy_workers(0),
  num_records(0),
  is_running_energy(false),
  outfile(0),
  record_file_delta_t(0.0f)
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
#if 0
    std::unique_lock<std::mutex> lk(end_workers_mut);
    end_workers = true;
    lk.unlock();
    //TODO: clean the threads
    for(unsigned int i = 0; i < this->max_detection_workers; ++i){
      this->detection_workers[i].join();
    }
#endif
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
