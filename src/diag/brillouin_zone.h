/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2026-01-12 20:40:07
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-16 16:16:38
*----------------------------------------------------------------------------*/
#ifndef BRILLOUIN_ZONE_H
#define BRILLOUIN_ZONE_H

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include "../lattice/lattice.h"
#include "./matrix.h"

namespace diag {

using kvector = Vector3d;
class special_Kpoint
{
public:
  special_Kpoint() {}
  special_Kpoint(const bool& exists, const std::string& name, const kvector& vec) 
  	{ exists_=exists; name_=name; vec_=vec; }
 ~special_Kpoint() {}
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
  int make_special_Kpoints(const lattice::Lattice& lattice);
  const std::vector<half_space>& half_spaces(void) const { return half_spaces_; } 
  const std::vector<special_Kpoint>& special_points(void) const { return special_Kpoints_; }
  const special_Kpoint& Gamma(void) const { return Gamma_; }
  const special_Kpoint& X(void) const { return X_; }
  const special_Kpoint& M(void) const { return M_; }
  const special_Kpoint& L(void) const { return L_; }
  const special_Kpoint& K(void) const { return K_; }
  const special_Kpoint& Kprime(void) const { return Kprime_; }
private:
  kvector b1_;
  kvector b2_;
  kvector b3_;
  /* 
   A half-space wrt to G-point is a tuple (G, ncap, d), where 
   'G' is a reciprocal lattice point, 'ncap' is a unit vector along G
   and 'd' is the distance from the origin to the perpendicular
   bisector of the line joining origin to G.
  */
  // The intersection of the half-spaces define the FBZ
  std::vector<half_space> half_spaces_;

  // list of special Kpoints
  std::vector<special_Kpoint> special_Kpoints_;
  special_Kpoint Gamma_{false, "Gamma", kvector(0,0,0)};
  special_Kpoint X_{false, "X", kvector(0,0,0)};
  special_Kpoint M_{false, "M", kvector(0,0,0)};
  special_Kpoint K_{false, "K", kvector(0,0,0)};
  special_Kpoint Kprime_{false, "Kprime", kvector(0,0,0)};
  special_Kpoint L_{false, "L", kvector(0,0,0)};
}; 







} // end namespace diag

#endif
