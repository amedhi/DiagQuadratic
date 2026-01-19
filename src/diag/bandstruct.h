/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2025-12-09 17:26:41
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-18 23:23:41
*----------------------------------------------------------------------------*/
#ifndef BANDSTURCT_H
#define BANDSTURCT_H

#include "../mcdata/mc_observable.h"
#include "./hamiltonian.h"
#include "./kspace.h"

namespace diag {

class BandStruct : public mcdata::MC_Observable
{
public:
  using MC_Observable::MC_Observable;
  void reset(void) override;
  void setup(const input::Parameters& inputs, const lattice::Lattice& lattice, 
    const kSpace& kspace, const Hamiltonian& ham);
  void compute(const kSpace& kspace, const Hamiltonian& ham);
  void print_heading(const std::string& header, 
    const std::vector<std::string>& xvars) override;
  void print_result(const std::vector<double>& xvals) override; 
private:
  enum spin {UP, DN, BOTH};
  enum kpath_type {SYMM, FULL};
  spin spin_sector_{UP};
  bool setup_done_{false};
  bool computation_done_{false};
  kpath_type kpath_{SYMM};
  int num_bands_{1};
  std::vector<kvector> symm_kpoints_;
  std::vector<int> symm_pidx_;
  std::vector<std::string> symm_pname_;
  std::vector<std::string> xvars_;
  //Eigen::MatrixXd E_kn_;
};


} // end namespace vmc

#endif