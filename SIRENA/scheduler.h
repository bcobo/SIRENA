
#ifndef SCHEDULER_H
#define SCHEDULER_H

#include <thread>
#include <future>
#include <vector>

#include "integraSIRENA.h"

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

struct detection
{
#if 0
  int pulse_length;//recontruct_init->pulse_length
  int eventsz;//record->trigger_size
  int tstart_pulse1;//reconstruct_init->tstartPulse1
  int tstart_pulse2;//reconstruct_init->tstartPulse2
  int tstart_pulse3;//reconstruct_init->tstartPulse3
  int intermediate;//reconstruct_init->intermediate
#endif
  
  tesrecord* rec;
  ReconstructInitSIRENA* rec_init;
  int n_record;
  int last_record;
  PulsesCollection* all_pulses;
  PulsesCollection* record_pulses;
  OptimalFilterSIRENA* optimal_filter;
  TesEventList* event_list;

  detection();
  detection(const detection& other);
  detection& operator=(const detection& other);
  ~detection();
};

struct energy
{
  tesrecord* rec;
  ReconstructInitSIRENA* rec_init;
  int n_record;
  PulsesCollection* all_pulses;
  PulsesCollection** record_pulses;
  OptimalFilterSIRENA** optimal_filter;

  energy();
  energy(const energy& other);
  energy& operator=(const energy& other);
  ~energy();
};

#define detection_input detection
#define detection_output energy

class scheduler
{
 public:

  void push_detection(const detection_input &input);
  void push_detection(TesRecord* record, 
                      int nRecord, 
                      int lastRecord, 
                      PulsesCollection *pulsesAll, 
                      ReconstructInitSIRENA** reconstruct_init, 
                      PulsesCollection** pulsesInRecord,
                      OptimalFilterSIRENA** optimal,
                      TesEventList* event_list);
  void run_energy(/*TODO: add parameters*/);
  void finish_reconstruction(PulsesCollection** pulsesAll, 
                             OptimalFilterSIRENA** optimalFilter);

  inline bool is_threading() const { return threading; }
  inline bool is_reentrant() const { return fits_is_reentrant(); }

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
  unsigned int num_records;
  
  bool threading;

  std::promise<bool>* detection_worker_status;
  std::promise<bool>* energy_worker_status;

  std::thread* detection_workers;
  std::thread* energy_workers;
  std::thread reconstruct_worker;
};

#endif
