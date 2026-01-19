/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2025-12-06 08:45:51
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-12 17:31:50
*----------------------------------------------------------------------------*/
#ifndef DIAG_H
#define DIAG_H

#include <iostream>
#include <Eigen/Eigenvalues>
#include "../scheduler/worker.h"
#include "../scheduler/mpi_comm.h"
#include "../lattice/constants.h"
#include "../lattice/lattice.h"
#include "./kspace.h"
#include "./hamiltonian.h"
#include "./observables.h"
//#include "./mf_model.h"

namespace diag {

class Diag : public scheduler::Worker
{
public:
  Diag(const input::Parameters& parms); 
  ~Diag() {}
  int start(const mpi::mpi_communicator& mpi_comm) override { return 0; } 
  int run(const input::Parameters& parms) override;
  void finish(void) override {} 
  void dostep(void) override {} 
  void halt(void) override {} 
  static void copyright_msg(std::ostream& os);
  static void print_copyright(std::ostream& os);
private:
  lattice::Lattice lattice_;
  model::Model model_;
  kSpace kspace_;
  Hamiltonian ham_;
  int num_kpoints_{1};
  int kblock_dim_{1};
  mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> es_k_up_;
  std::vector<Vector3d> symm_line_;
  std::vector<int> symm_pidx_;
  std::vector<std::string> symm_pname_;
  // outputs
  bool need_chern_number_{false};
  bool need_ebands_full_{false};
  bool need_ebands_symm_{false};
  bool need_band_gap_{false};
  bool need_fermi_energy_{false};
  bool need_dos_{false};

  // observables
  ObservableSet observables_;
  std::string prefix_{"./"};
  std::vector<std::string> xvar_names_;
  std::vector<double> xvar_values_;

  mutable std::ostringstream info_str_;

  void make_info_str(const input::Parameters& inputs);
  int compute_chern_number(void);
  int compute_band_gap(void);
  int compute_fermi_energy(const double& hole_doping);
  int compute_dos(void);
};


} // end namespace diag

#endif
