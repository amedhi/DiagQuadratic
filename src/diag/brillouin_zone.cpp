/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2026-01-12 20:40:43
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-18 23:24:12
*----------------------------------------------------------------------------*/
#include "brillouin_zone.h"

namespace diag {

int BrillouinZone::make_special_Kpoints(const lattice::Lattice& lattice)
{
	// List special Kpoints for each Bravais lattice
  special_Kpoints_.clear();
  if (lattice.brav()==lattice::brav_id::ZERO) {
    kvector Gamma = kvector(0,0,0);
    Gamma_.set("Gamma",Gamma);
  }

  else if (lattice.brav()==lattice::brav_id::CHAIN) {
    kvector Gamma = kvector(0,0,0);
    Gamma_.set("Gamma",Gamma);
    kvector X = 0.5*b1_;
    X_.set("X",X);
  }

	else if (lattice.brav()==lattice::brav_id::SQUARE) {
    kvector Gamma = kvector(0,0,0);
    Gamma_.set("Gamma",Gamma);
    kvector X = 0.5*b1_;
    X_.set("X",X);
    kvector M = X+0.5*b2_;
    M_.set("M",M);
	}

	else if (lattice.brav()==lattice::brav_id::HEXAGONAL) {
    kvector Gamma = kvector(0,0,0);
    Gamma_.set("Gamma",Gamma);
    kvector M = 0.5*b1_;
    M_.set("M",M);
  	kvector K = (2*b1_+b2_)/3;
    K_.set("K",K);
  	kvector Kprime = (2*b1_-b2_)/3;
    Kprime_.set("Kprime",Kprime);
	}

	else {
    throw std::range_error("error: BrillouinZone: special Kpoints not defined for this lattice");
	}

	return 0;
} 

int BrillouinZone::construct(const lattice::Lattice& lattice, 
  const Vector3d& b1, const Vector3d& b2, const Vector3d& b3)
{
	// primitive vectors of the reciprocal lattice
	b1_ = b1;
	b2_ = b2;
	b3_ = b3;

  /* 
   Take a few reciprocal lattice points (G-points) near the origin.
   A half-space wrt to G-point is a tuple (G, ncap, d), where 
   'G' is a reciprocal lattice point, 'ncap' is a unit vector along G
   and 'd' is the distance from the origin to the perpendicular
   bisector of the line joining origin to G.
  */
  half_spaces_.clear();
  // 1D lattices
  if (lattice.dimension() == 1) {
    std::list<int> idx = {1,-1};
    for (const auto& i : idx) {
      kvector G = i*b1;
      double norm = G.norm();
      kvector ncap = G/norm;
      double d = 0.5*norm;
      half_spaces_.push_back({G,ncap,d});
    }
  }

  // 2D lattices
  if (lattice.dimension() == 2) {
    std::list<int> idx = {2,1,0,-1,-2};
    for (const auto& i : idx) {
      for (const auto& j : idx) {
        if (i==0 && j==0) continue;
        kvector G = i*b1+j*b2;
        double norm = G.norm();
        kvector ncap = G/norm;
        double d = 0.5*norm;
        half_spaces_.push_back({G,ncap,d});
      }
    }
  }

  // 3D lattices
  if (lattice.dimension() == 3) {
    std::list<int> idx = {2,1,0,-1,-2};
    for (const auto& i : idx) {
      for (const auto& j : idx) {
        for (const auto& k : idx) {
          if (i==0 && j==0 && k==0) continue;
          kvector G = i*b1+j*b2+k*b3;
          double norm = G.norm();
          kvector ncap = G/norm;
          double d = 0.5*norm;
          half_spaces_.push_back({G,ncap,d});
        }
      }
    }
  }

  // Special K-points
  make_special_Kpoints(lattice);

  return 0;
}


} // end namespace diag
