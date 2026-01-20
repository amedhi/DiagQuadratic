/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2026-01-18 23:06:00
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-20 09:41:30
*----------------------------------------------------------------------------*/
#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <map>
#include "../mcdata/mc_observable.h"
#include "./hamiltonian.h"
#include "./kspace.h"

namespace diag {

class WaveFunction : public mcdata::MC_Observable
{
public:
  using MC_Observable::MC_Observable;
  void reset(void) override;
  void setup(const input::Parameters& inputs, const lattice::Lattice& lattice,
    kSpace& kspace);
  void compute(const kSpace& kspace, const Hamiltonian& ham);
  void print_heading(const std::string& header, 
    const std::vector<std::string>& xvars) override;
  void print_result(const std::vector<double>& xvals) override; 
private:
  enum kpath_type {SYMM, FULL};
  struct qn_type {int k; int n; int s;};

  // quantum states
  int num_states_{0};
  std::map<int,std::set<int>> up_states_; // k-{n1,n2,...} (up)
  std::map<int,std::set<int>> dn_states_; // k-{n1,n2,...} (dn)

  bool setup_done_{false};
  bool computation_done_{false};
  kpath_type kpath_{SYMM};
  std::vector<std::string> xvars_;
  //Eigen::MatrixXd E_kn_;
};


} // end namespace vmc

#endif
