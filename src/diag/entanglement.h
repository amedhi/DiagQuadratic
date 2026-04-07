/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2026-04-06 10:45:18
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-04-06 15:16:42
*----------------------------------------------------------------------------*/
#ifndef ENTANGLEMENT_H
#define ENTANGLEMENT_H

#include "../mcdata/mc_observable.h"
#include "./hamiltonian.h"
#include "./kspace.h"

namespace diag {

class EntanglementEntropy : public mcdata::MC_Observable
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
  enum spin {UP, DN, BOTH};
  enum kpath_type {SYMM, FULL};
  bool setup_done_{false};
  bool computation_done_{false};
  int size_A_{0};
  int size_cAB_{0};
  std::vector<int> subsys_A_;
  std::vector<int> subsys_B_;
  std::vector<int> subsys_cAB_;
  double entropy_A_{0.0};
  double entropy_B_{0.0};
  double entropy_cAB_{0.0};
  std::vector<std::string> xvars_;
};


} // end namespace vmc

#endif

