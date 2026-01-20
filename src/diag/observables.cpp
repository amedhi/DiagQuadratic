/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2025-12-09 17:14:46
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-20 10:26:47
*----------------------------------------------------------------------------*/
#include "observables.h"
#include <boost/algorithm/string.hpp>

namespace diag {

ObservableSet::ObservableSet() 
  : band_struct_("BandStruct")
  , wave_function_("WaveFunction")
{
}

void ObservableSet::init(const input::Parameters& inputs, 
  const lattice::Lattice& lattice, kSpace& kspace, const Hamiltonian& ham, 
  const std::string& prefix)
{
  // file open mode
  std::string mode = inputs.set_value("mode", "NEW");
  boost::to_upper(mode);
  if (mode=="APPEND") replace_mode_ = false;
  else replace_mode_ = true;
  num_xvars_ = 0; 

  // files
  band_struct_.set_ofstream(prefix);
  wave_function_.set_ofstream(prefix);

  // switch on required observables
  band_struct_.check_on(inputs,replace_mode_);
  wave_function_.check_on(inputs,replace_mode_);

  // set up observables
  if (band_struct_) band_struct_.setup(inputs,lattice,kspace);
  if (wave_function_) wave_function_.setup(inputs,lattice,kspace);
  print_heading();
}

void ObservableSet::reset(void)
{
  if (band_struct_) band_struct_.reset();
  if (wave_function_) wave_function_.reset();
}

void ObservableSet::reset_batch_limit(const int& sample_size)
{
}

void ObservableSet::reset_grand_data(void)
{
  if (band_struct_) band_struct_.reset_grand_data();
  if (wave_function_) wave_function_.reset_grand_data();
}

void ObservableSet::save_results(void)
{
  if (band_struct_) band_struct_.save_result();
  if (wave_function_) wave_function_.save_result();
}

void ObservableSet::avg_grand_data(void)
{
  if (band_struct_) band_struct_.avg_grand_data();
  if (wave_function_) wave_function_.avg_grand_data();
}


int ObservableSet::compute(const lattice::Lattice& lattice, const kSpace& kspace,
    const Hamiltonian& ham)
{
  if (band_struct_) band_struct_.compute(kspace, ham);
  if (wave_function_) wave_function_.compute(kspace, ham);
  return 0;
}

void ObservableSet::finalize(void)
{
}

void ObservableSet::as_functions_of(const std::vector<std::string>& xvars)
{
  xvars_ = xvars;
  num_xvars_ = xvars_.size();
}

void ObservableSet::as_functions_of(const std::string& xvar)
{
  xvars_ = {xvar};
  num_xvars_ = 1;
}

void ObservableSet::switch_off(void) {
  band_struct_.switch_off();
  wave_function_.switch_off();
}

void ObservableSet::print_heading(void)
{
  band_struct_.print_heading(headstream_.rdbuf()->str(),xvars_);
  wave_function_.print_heading(headstream_.rdbuf()->str(),xvars_);
}

void ObservableSet::print_results(const std::vector<double>& xvals) 
{
  if (num_xvars_ != xvals.size()) 
    throw std::invalid_argument("Observables::print_result: 'x-vars' size mismatch");
  if (band_struct_) {
    band_struct_.print_heading(headstream_.rdbuf()->str(),xvars_);
    band_struct_.print_result(xvals);
  }
  if (wave_function_) {
    wave_function_.print_heading(headstream_.rdbuf()->str(),xvars_);
    wave_function_.print_result(xvals);
  }
}

void ObservableSet::print_results(const double& xval) 
{
  if (num_xvars_ != 1) 
    throw std::invalid_argument("ObservableSet::print_result: 'x-vars' size mismatch");
  std::vector<double> xvals{xval};
  if (band_struct_) {
    band_struct_.print_heading(headstream_.rdbuf()->str(),xvars_);
    band_struct_.print_result(xvals);
  }
  if (wave_function_) {
    wave_function_.print_heading(headstream_.rdbuf()->str(),xvars_);
    wave_function_.print_result(xvals);
  }
}

void ObservableSet::MPI_send_results(const mpi::mpi_communicator& mpi_comm, 
  const mpi::proc& proc, const int& msg_tag)
{
  if (band_struct_) band_struct_.MPI_send_data(mpi_comm, proc, msg_tag);
  if (wave_function_) wave_function_.MPI_send_data(mpi_comm, proc, msg_tag);
}

void ObservableSet::MPI_recv_results(const mpi::mpi_communicator& mpi_comm, 
  const mpi::proc& proc, const int& msg_tag)
{
  if (band_struct_) band_struct_.MPI_add_data(mpi_comm, proc, msg_tag);
  if (wave_function_) wave_function_.MPI_add_data(mpi_comm, proc, msg_tag);
}

} // end namespace vmc

