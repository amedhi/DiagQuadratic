/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2025-12-09 13:23:17
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-16 11:22:15
*----------------------------------------------------------------------------*/
#include <iomanip>
#include <fstream>
#include <numeric>
#include <boost/filesystem.hpp>
#include "diag.h"

namespace diag {

Diag::Diag(const input::Parameters& inputs) 
  : lattice_(inputs) 
  , model_(inputs,lattice_) 
  , kspace_(lattice_)
  , ham_(model_,lattice_) 
{
  // output directory
  std::string outdir = inputs.set_value("prefix", "");
  prefix_ = "./"+outdir+"/";
  boost::filesystem::path prefix_dir(prefix_);
  boost::filesystem::create_directory(prefix_dir);

  // x-variables (assuming max 1)
  xvar_names_.resize(1);
  xvar_values_.resize(1);
  xvar_names_[0]=inputs.set_value("as_func_of","X");

  // observables
  make_info_str(inputs);
  observables_.headstream() << info_str_.str();
  observables_.init(inputs,lattice_,kspace_,ham_,prefix_);
  observables_.as_functions_of(xvar_names_);
}

int Diag::run(const input::Parameters& inputs) 
{
  observables_.compute(lattice_,kspace_,ham_);
  observables_.print_results(0.0);
  return 0;
}

int Diag::compute_chern_number(void) 
{
  return 0;
}

int Diag::compute_band_gap(void) 
{
  return 0;
}

void Diag::make_info_str(const input::Parameters& inputs)
{
  info_str_.clear();
  copyright_msg(info_str_);
  info_str_ << "# "<< inputs.job_id() <<"\n"; 
  info_str_ << model_.info_str(); 
  //info_str_ << ", min_interval = " << min_interval_;
  //info_str_ << ", max_interval = " << max_interval_ << "\n";
}

void Diag::copyright_msg(std::ostream& os)
{
  os << "#" << std::string(72,'-') << "\n";
  os << "#" << " Program: Diag Quadratic\n";
  os << "#" << "          (c) Amal Medhi <amedhi@iisertvm.ac.in>\n";
  os << "#" << std::string(72,'-') << "\n";
}

void Diag::print_copyright(std::ostream& os)
{
  std::cout << "Diag\n";
}


} // end namespace diag
