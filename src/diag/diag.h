/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef DIAG_H
#define DIAG_H

#include <iostream>
#include <Eigen/Eigenvalues>
#include "../scheduler/worker.h"
#include "../lattice/lattice.h"
#include "../lattice/graph.h"
#include "../model/model.h"
#include "./mf_model.h"
#include "./blochbasis.h"

namespace diag {

class Diag : public scheduler::Worker
{
public:
  Diag(const input::Parameters& parms); 
  ~Diag() {}
  int start(const input::Parameters& parms) override { return 0; }
  int run(const input::Parameters& parms) override;
  void finish(void) override {} 
  void dostep(void) override {} 
  void halt(void) override {} 
  static void print_copyright(std::ostream& os);
private:
  lattice::LatticeGraph graph_;
  model::Hamiltonian model_;
  MF_Model mf_model_;
  basis::BlochBasis blochbasis_;
  unsigned num_kpoints_{1};
  unsigned kblock_dim_{1};
  mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> es_k_up;

  std::vector<Vector3d> symm_line_;
};


} // end namespace diag

#endif
