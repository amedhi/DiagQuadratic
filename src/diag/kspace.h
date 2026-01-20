/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2025-12-10 19:14:38
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-19 16:35:30
*----------------------------------------------------------------------------*/
#ifndef KSPACE_H
#define KSPACE_H
#include <iostream>
#include <vector>
#include <set>
#include <stdexcept>
#include "./matrix.h"
#include "../lattice/lattice.h"

namespace diag {

using kvector = Vector3d;

class special_kpoint
{
public:
  special_kpoint() {}
  special_kpoint(const bool& exists, const std::string& name, const kvector& vec) 
    { exists_=exists; name_=name; vec_=vec; }
 ~special_kpoint() {}
  void set(const std::string& name, const kvector& vec) 
    { exists_=true; name_=name; vec_=vec; }
  const bool exists(void) const  { return exists_; }
  const kvector& vec(void) const 
    { 
      if (!exists_) throw std::range_error("error: Special Kpoint '"+name_+"' not deined");
      return vec_;
    }
  const std::string& name(void) const 
    { 
      if (!exists_) throw std::range_error("error: Special Kpoint not deined");
      return name_;
    }
private:
  bool exists_{false};
  std::string name_{""};
  kvector vec_{kvector(0,0,0)};
};

class BrillouinZone
{
public:
  struct half_space {kvector G; kvector ncap; double d;};
  BrillouinZone() {}
  ~BrillouinZone() {}
  int construct(const lattice::Lattice& lattice, const Vector3d& b1, 
    const Vector3d& b2, const Vector3d& b3);
  const std::vector<half_space>& half_spaces(void) const { return half_spaces_; } 
private:
  /* 
   A half-space wrt to G-point is a tuple (G, ncap, d), where 
   'G' is a reciprocal lattice point, 'ncap' is a unit vector along G
   and 'd' is the distance from the origin to the perpendicular
   bisector of the line joining origin to G.
  */
  // The intersection of the half-spaces define the FBZ
  std::vector<half_space> half_spaces_;
}; 


class kSpace 
{
public:
  kSpace() {}
  kSpace(const lattice::Lattice& lattice);
  ~kSpace() {}
  int init(const lattice::Lattice& lattice);
  int construct_kmesh(const lattice::Lattice& lattice);
  int construct_kpath(const lattice::Lattice& lattice);
  const kvector& vector_b1(void) const { return b1_; }
  const kvector& vector_b2(void) const { return b2_; }
  const kvector& vector_b3(void) const { return b3_; }
  const std::vector<kvector>& kmesh(void) const { return kmesh_; }
  const std::vector<kvector>& kpath(void) const { return kpath_; }
  const std::vector<special_kpoint>& kpath_nodes(void) const { return kpath_nodes_; }
  const std::vector<int>& kpath_nodes_idx(void) const { return kpath_nodes_idx_; }
  const BrillouinZone& FBZ(void) const { return FBZ_; }
  const special_kpoint& Gamma(void) const { return Gamma_; }
  const special_kpoint& X(void) const { return X_; }
  const special_kpoint& M(void) const { return M_; }
  const special_kpoint& L(void) const { return L_; }
  const special_kpoint& K(void) const { return K_; }
  const special_kpoint& Kprime(void) const { return Kprime_; }
private:
  kvector b1_;
  kvector b2_;
  kvector b3_;
  BrillouinZone FBZ_;

  // full kmesh
  bool kmesh_constructed_{false};
  std::vector<kvector> kmesh_;
  int num_kpoints_;

  // symmetry kpath
  bool kpath_constructed_{false};
  std::vector<kvector> kpath_;
  std::vector<special_kpoint> kpath_nodes_;
  std::vector<int> kpath_nodes_idx_;

  // Special Kpoints
  special_kpoint Gamma_{false, "Gamma", kvector(0,0,0)};
  special_kpoint X_{false, "X", kvector(0,0,0)};
  special_kpoint M_{false, "M", kvector(0,0,0)};
  special_kpoint K_{false, "K", kvector(0,0,0)};
  special_kpoint Kprime_{false, "Kprime", kvector(0,0,0)};
  special_kpoint L_{false, "L", kvector(0,0,0)};

  // helper functions
  void make_kpoints(const lattice::Lattice& lattice);
};


} // end namespace diag

#endif
