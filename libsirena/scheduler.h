/**********************************************************************
*   This software is part of the grant PID2021-122955OB-C41
*   and by 'ERDF A way of making Europe'.
*
***********************************************************************
*                      SCHEDULER
*
*  File:       scheduler.h
*  Developers: Beatriz Cobo
* 	           cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/

#ifndef SCHEDULER_H
#define SCHEDULER_H

#include <thread>
#include <future>
#include <vector>

#include "integraSIRENA.h"

#include <time.h>

struct phidlist{
  /** Array containing the phIDs */
  long* phid_array;
  
  /** Array containing the corresponding impact times (only relevant in case of wait list) */
  double* times;
  
  /** Boolean to state if this should be a wait list */
  unsigned char wait_list;
  
  /** Number of elements in the wait list */
  int n_elements;
  
  /** Current index in the list */
  int index;
  
  /** Size of the list */
  int size;

  phidlist();
  phidlist(const phidlist& other);
  phidlist(PhIDList* other);
  phidlist& operator=(const phidlist& other);
  ~phidlist();

  PhIDList* get_PhIDList() const;
};//PhIDList;

struct tesrecord{
  /** Number of ADC values in the record */
  unsigned long trigger_size;
  
  /** Start time of the record */
  double time;
  
  /** Time difference between two samples */
  double delta_t;

  /** Buffer to read a record of ADC values */
  uint16_t* adc_array;
  
  /** Double version of the record */
  double* adc_double;
  
  /** PIXID of the record */
  long pixid;
  
  /** Array of the PH_ID in the record */
  phidlist* phid_list;
  tesrecord();
  tesrecord(TesRecord* other);
  tesrecord(const tesrecord& other);
  tesrecord& operator=(const tesrecord& other);
  ~tesrecord();

  TesRecord* get_TesRecord() const;
};//TesRecord;

struct datax
{ 
  TesRecord* rec;
  int trig_reclength;
  ReconstructInitSIRENA* rec_init;
  int n_record;
  int last_record;
  PulsesCollection* all_pulses;
  PulsesCollection* record_pulses;
  OptimalFilterSIRENA* optimal_filter;
  TesEventListSIRENA* event_list;

  datax();
  datax(const datax& other);
  datax& operator=(const datax& other);
  ~datax();
};

#define sirena_data datax

class scheduler
{
 public:

  void push_detection(TesRecord* record, int trig_reclength,
                      int nRecord, 
                      int lastRecord, 
                      PulsesCollection *pulsesAll, 
                      ReconstructInitSIRENA** reconstruct_init, 
                      PulsesCollection** pulsesInRecord,
                      TesEventListSIRENA* event_list);
  void finish_reconstruction(PulsesCollection** pulsesAll);
  void finish_reconstruction_v2(PulsesCollection** pulsesAll);

  inline bool is_threading() const { return threading; }
  inline bool is_reentrant() const { return fits_is_reentrant(); }

  inline void set_is_running_energy(bool val){ this->is_running_energy = val; }

  inline bool has_records(){return(this->num_records > this->current_record);}

  void get_test_event(TesEventListSIRENA** test_event, TesRecord** record);

  virtual ~scheduler();
  inline static scheduler* get()
  { 
    return instance ? instance : instance = new scheduler();
  }
 private:
  static scheduler* instance;
  scheduler();
  //scheduler(const scheduler& copy){}
  scheduler(const scheduler& copy) = delete;
  //scheduler& operator=(const scheduler&){}
  scheduler& operator=(const scheduler&) = delete;

  void init();
  void init_v2();

  unsigned int num_cores;
  unsigned int max_detection_workers;
  unsigned int max_energy_workers;
  unsigned int num_records;
  unsigned int current_record;

  sirena_data** data_array;

  bool is_running_energy;
  
  bool threading;

  std::promise<bool>* detection_worker_status;
  std::promise<bool>* energy_worker_status;

  std::thread* detection_workers;
  std::thread* energy_workers;
  std::thread reconstruct_worker;
};

#endif
