/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2025-12-09 17:07:45
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-16 15:40:30
*----------------------------------------------------------------------------*/
#ifndef MC_OBSERVABLES_H
#define MC_OBSERVABLES_H

#include <string>
#include <vector>
#include <stdexcept>
#include "../scheduler/task.h"
#include "../lattice/lattice.h"
#include "../model/model.h"
#include "./bandstruct.h"

namespace diag {

class ObservableSet 
{
public:
  ObservableSet();
  ~ObservableSet() {}
  std::stringstream& headstream(void) { return headstream_; }
  //void init(const input::Parameters& inputs, 
  //  void (&print_copyright)(std::ostream& os), const lattice::Lattice& lattice, 
  //  const model::Hamiltonian& model, const SysConfig& config);
  void init(const input::Parameters& inputs, const lattice::Lattice& lattice, 
    const kSpace& kspace, const Hamiltonian& ham, const std::string& prefix);
  void as_functions_of(const std::vector<std::string>& xvars=std::vector<std::string>());
  void as_functions_of(const std::string& xvar);
  void switch_off(void);
  void reset(void); 
  void reset_batch_limit(const int& sample_size);
  void reset_grand_data(void); 
  int compute(const lattice::Lattice& lattice, const kSpace& kspace, const Hamiltonian& ham);
  //int do_measurement(const lattice::Lattice& lattice, 
  //  const model::Hamiltonian& model, const SysConfig& config, const SiteDisorder& site_disorder);
  void save_results(void); 
  void avg_grand_data(void); 
  const BandStruct& band_struct(void) const { return band_struct_; }
  void finalize(void);
  void print_heading(void);
  void print_results(const std::vector<double>& xvals=std::vector<double>()); 
  void print_results(const double& xval); 
  void MPI_send_results(const mpi::mpi_communicator& mpi_comm, 
    const mpi::proc& proc, const int& msg_tag);
  void MPI_recv_results(const mpi::mpi_communicator& mpi_comm, 
    const mpi::proc& proc, const int& msg_tag); 
private:
  bool replace_mode_{true};
  std::stringstream headstream_;
  std::vector<std::string> xvars_;
  int num_xvars_{0};
  BandStruct band_struct_; 
};





} // end namespace vmc

#endif
