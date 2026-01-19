/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2025-12-10 19:14:38
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-14 13:21:04
*----------------------------------------------------------------------------*/
#ifndef KSPACE_H
#define KSPACE_H
#include <iostream>
#include <vector>
#include <set>
#include <stdexcept>
#include "./matrix.h"
#include "./brillouin_zone.h"
#include "../lattice/lattice.h"

namespace diag {

class kSpace 
{
public:
  kSpace() {}
  kSpace(const lattice::Lattice& lattice);
  ~kSpace() {}
  int construct(const lattice::Lattice& lattice);
  const Vector3d& vector_b1(void) const { return b1_; }
  const Vector3d& vector_b2(void) const { return b2_; }
  const Vector3d& vector_b3(void) const { return b3_; }
  const int& num_kpoints(void) const { return num_kpoints_; }
  kvector kpoint(const int& k) const { return kpoints_[k]; }
  const std::vector<kvector>& kpoints(void) { return kpoints_; }
  const BrillouinZone& FBZ(void) const { return FBZ_; }
private:
  Vector3d b1_;
  Vector3d b2_;
  Vector3d b3_;
  BrillouinZone FBZ_;
  std::vector<kvector> kpoints_;
  int num_kpoints_;
  // helper functions
  void make_kpoints(const lattice::Lattice& lattice);
};


} // end namespace diag

#endif
