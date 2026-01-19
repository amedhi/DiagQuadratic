/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2025-12-06 12:11:48
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-18 23:04:43
*----------------------------------------------------------------------------*/
#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <stdexcept>
#include "../scheduler/task.h"
#include "../lattice/lattice.h"
#include "../model/model.h"
//#include "./blochbasis.h"
#include "./matrix.h"

//constexpr std::complex<double> ii(void) { return std::complex<double>{0.0,static_cast<double>(1.0)}; }

namespace diag {

// map [site_idx, orbital_idx] pair to basis_idx
using basis_idx_map = std::vector<std::vector<int>>;

class UnitcellTerm
{
public:
  UnitcellTerm() {}
  ~UnitcellTerm() {}
  void build_bondterm(const model::ModelTerm& bterm, const lattice::Lattice& lattice,
      const basis_idx_map& idx_map);
  void build_siteterm(const model::ModelTerm& sterm, const lattice::Lattice& lattice,
      const basis_idx_map& idx_map);
  void eval_coupling_constant(const model::ModelParams& cvals, const model::ModelParams& pvals);
  const int& num_out_bonds(void) const { return num_out_bonds_; } 
  const Vector3d& bond_vector(const int& i) const { return bond_vectors_[i]; }
  const ComplexMatrix& coeff_matrix(const int& i=0) const { return coeff_matrices_[i]; }
  //const double& coupling(const unsigned& site_type) const; 
  const model::op::quantum_op& qn_operator(void) const { return op_; }
private:
  model::op::quantum_op op_;
  int num_out_bonds_;
  int num_basis_sites_;
  int dim_;
  std::vector<ComplexMatrix> coeff_matrices_;
  std::vector<strMatrix> expr_matrices_;
  std::vector<Vector3d> bond_vectors_;
};

class Hamiltonian : public model::Model
{
public:
  Hamiltonian() {}
  Hamiltonian(const model::Model& model, const lattice::Lattice& lattice);
  ~Hamiltonian() {}
  int init(const lattice::Lattice& lattice) override;
  int finalize(const lattice::Lattice& lattice);
  int update(const input::Parameters& inputs);
  void update_terms(void) override;
  int update_site_parameter(const std::string& pname, const double& pvalue);
  int construct_upspin_block(const Vector3d& kvec) const;
  int construct_dnspin_block(const Vector3d& kvec) const;
  int construct_kblock(const Vector3d& kvec) const;
  int construct_pairing_part(const Vector3d& kvec) const;
  const int& dimension(void) const { return dim_; }
  const ComplexMatrix& upspin_block(void) const { return quadratic_block_up_; }
  const ComplexMatrix& dnspin_block(void) const { return quadratic_block_dn_; }
  const ComplexMatrix& pairing_block(void) const { return pairing_block_; }
private:
  using Model = model::Model;
  // dimension of Ham matrices
  int dim_;
  int num_basis_sites_;
  // [site_idx, orbital_idx] pair to basis_idx mapping
  basis_idx_map site_orb_idx_;
  // Unitecell Ham matrices
  std::vector<UnitcellTerm> usite_terms_;
  std::vector<UnitcellTerm> ubond_terms_;
  // matrices in kspace representation
  mutable ComplexMatrix quadratic_block_up_;
  mutable ComplexMatrix quadratic_block_dn_;
  mutable ComplexMatrix pairing_block_;
  mutable ComplexMatrix work; //, work2;

  int init_ham(const lattice::Lattice& lattice);
  int build_unitcell_terms(const lattice::Lattice& lattice);
  int update_unitcell_terms(void);
};


} // end namespace diag

#endif