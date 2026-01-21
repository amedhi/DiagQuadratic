/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2025-12-10 19:13:42
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-21 12:44:12
*----------------------------------------------------------------------------*/
//#include <iomanip>
#include "kspace.h"

namespace diag {

kSpace::kSpace(const lattice::Lattice& lattice) 
{
  init(lattice);
  kmesh_constructed_ = false;
  kpath_constructed_ = false;
}

int kSpace::init(const lattice::Lattice& lattice)
{
  Vector3d a1, a2, a3, n1, n2, n3;
  double v;
  using bc = lattice::boundary_type;

  // reciprocal lattice vectors 
  a1 = lattice.vector_a1();
  a2 = lattice.vector_a2();
  a3 = lattice.vector_a3();
  b1_ = Vector3d(0.0,0.0,0.0);
  b2_ = Vector3d(0.0,0.0,0.0);
  b3_ = Vector3d(0.0,0.0,0.0);
  /*
  std::cout << "a1= " << a1.transpose() << "\n";
  std::cout << "a2= " << a2.transpose() << "\n";
  std::cout << "a3= " << a3.transpose() << "\n";
  */

  int symmetry_type = 0;
  if (lattice.bc1() == bc::periodic) {
    b1_ = two_pi() * a1 / a1.dot(a1); 
    symmetry_type = symmetry_type + 1;
  } 

  if (lattice.bc2() == bc::periodic) {
    switch (symmetry_type) {
      case 0:
        b2_ = two_pi() * a2 / a2.dot(a2); 
        break;
      case 1:
        n3 = a1.cross(a2);
        v = a1.dot(a2.cross(n3));
        b1_ = two_pi() * a2.cross(n3) / v;
        b2_ = two_pi() * n3.cross(a1) / v;
        break;
      default: break;
    }
    symmetry_type = symmetry_type + 2;
  } 

  if (lattice.bc3() == bc::periodic) {
    switch (symmetry_type) {
      case 0:
        b3_ = two_pi() * a3 / a3.dot(a3); 
        break;
      case 1:
        n2 = a3.cross(a1);
        v = a1.dot(n2.cross(a3));
        b1_ = two_pi() * n2.cross(a3) / v;
        b3_ = two_pi() * a1.cross(n2) / v;
        break;
      case 2:
        n1 = a2.cross(a3);
        v = n1.dot(a2.cross(a3));
        b2_ = two_pi() * a3.cross(n1) / v;
        b3_ = two_pi() * n1.cross(a2) / v;
        break;
      case 3:
        v = a1.dot(a2.cross(a3));
        b1_ = two_pi() * a2.cross(a3) / v;
        b2_ = two_pi() * a3.cross(a1) / v;
        b3_ = two_pi() * a1.cross(a2) / v;
        break;
      default: break;
    }
  }
  /*
  std::cout << "b1= " << b1_.transpose() << "\n";
  std::cout << "b2= " << b2_.transpose() << "\n";
  std::cout << "b3= " << b3_.transpose() << "\n";
  */

  // Construct First Brillouin Zone
  FBZ_.construct(lattice,b1_,b2_,b3_);

  return 0;
}

int kSpace::construct_kmesh(const lattice::Lattice& lattice)
{
  if (kmesh_constructed_) return 0;

  // Shifts in k-points due to 'twists in BCs'
  kvector bc_shift(0.0,0.0,0.0);
  bc_shift(0) = lattice.bc1_twist()/(two_pi()*lattice.size1());
  bc_shift(1) = lattice.bc2_twist()/(two_pi()*lattice.size2());
  bc_shift(2) = lattice.bc3_twist()/(two_pi()*lattice.size3());
  //std::cout << "k-shift = " << twist_shift.transpose() << "\n";

  // Generate the k-points 
  double x1, x2, x3;
  Vector3i n = {0,0,0};
  Vector3i m = {-lattice.size1()/2, -lattice.size2()/2, -lattice.size3()/2};
  // Special for 'HONEYCOMB' lattice Ribbon Goemetry
  if (lattice.id()==lattice::lattice_id::HONEYCOMB ||
    lattice.id()==lattice::lattice_id::HONEYCOMB2) {
    if (lattice.dimension()==1) {
      m = {0,0,0};
    }
  }

  num_kpoints_ = lattice.num_unitcells();
  kmesh_.resize(num_kpoints_);
  for (int i=0; i<num_kpoints_; i++) {
    x1 = static_cast<double>(m(0)+n(0))/lattice.size1() + bc_shift(0);
    x2 = static_cast<double>(m(1)+n(1))/lattice.size2() + bc_shift(1);
    x3 = static_cast<double>(m(2)+n(2))/lattice.size3() + bc_shift(2);
    kmesh_[i] = x1*b1_ + x2*b2_ + x3*b3_;
    /*
    auto kvec = x1 * b1 + x2 * b2 + x3 * b3;
    std::cout << i << ": " << kvec.transpose() << "\n";
    getchar();
    */
    //translation_vectors.push_back(n);
    n = lattice.get_next_bravindex(n);
  }

  /*
   The k-points may not lie inside FBZ for certain Bravais lattices.
   For those lattices, construct the FBZ explicitly and translate the
   k-points appropriately so as to shift them into the FBZ.
  */
  if (lattice.brav()==lattice::brav_id::HEXAGONAL) {
    // translate of the k-points outsize FBZ
    for (auto& kvec : kmesh_) {
      /*
       Check if the point lies within the half-spaces. If not, shift
       it by the corresponding '-G' vector. 
      */
      for (const auto& hs : FBZ_.half_spaces()){
        // projection of 'kvec' along ncap
        double proj = kvec.dot(hs.ncap);
        // Test if the point lies outside BZ boundary within a Tol
        if ((proj-hs.d)>1.0E-12) { 
          // point lies outside the half-space: translate it
          kvec = kvec - hs.G;
        }
      }
    }
  }

  /*
  int i=0;
  for (const auto& kvec : kmesh_) {
    std::cout << i++ << "  " << kvec(0) << "  " << kvec(1) << "\n";
  }
  std::cout << "Exiting at kspace.cpp\n";
  std::exit(0);
  */

  kmesh_constructed_ = true;
  return 0;
}

//=====================================================================
int BrillouinZone::construct(const lattice::Lattice& lattice, 
  const Vector3d& b1, const Vector3d& b2, const Vector3d& b3)
{
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

  return 0;
}

} // end namespace diag
