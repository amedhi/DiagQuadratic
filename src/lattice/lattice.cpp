/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2025-12-06 11:43:46
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-03-06 15:00:01
*----------------------------------------------------------------------------*/
#include <iomanip>
#include <cassert>
#include "lattice.h"
#include <boost/algorithm/string.hpp>

namespace lattice {

/*---------------'site' class-----------------*/
int Site::num_site = 0;

Site::Site(const int& uid, const int& type, const int& num_orbitals, 
 const Vector3i& bravindex, const Vector3d& coord, 
 const Vector3d& cell_coord)
{
  assert(uid>=0);
  assert(type>=0);
  assert(num_orbitals>=1);
  id_ = num_site++;
  uid_ = uid;
  type_ = type;
  num_orbitals_ = num_orbitals;
  bravindex_ = bravindex;
  coord_ = coord;
  cell_coord_ = cell_coord;
}

void Site::translate_by(const int& id_offset, const Vector3i& bravindex_offset, const Vector3d& coord_offset)
{
  id_ += id_offset;
  bravindex_ += bravindex_offset; 
  coord_ += coord_offset; 
  cell_coord_ += coord_offset; 
}

// friends
std::ostream& operator<<(std::ostream& os, const Site& s) 
{
  os << std::fixed << std::showpoint;
  os << "site " << s.id_ << ": type = " << s.type_; 
  os << ", orbitals =" << s.num_orbitals_ << "\n";
  os << "bravindex = (" << s.bravindex_(0) << ", " << s.bravindex_(1) << ", " << s.bravindex_(2) << ")\n";
  os << "coord = (" << s.coord_(0) << ", " << s.coord_(1) << ", " << s.coord_(2) << ")\n";
  os << "cell_coord = (" << s.cell_coord_(0) << ", " << s.cell_coord_(1) << ", " << s.cell_coord_(2) << ")\n";
  return os;
}

/*---------------'bond' class-----------------*/
int Bond::num_bond = 0;

Bond::Bond() 
{
  id_ = num_bond++;
  type_ = 0;
  ngb_ = 1;
  bravindex_ = Vector3i(0,0,0);
  sign_ = 1;
  bc_state_[0] = 1;
  bc_state_[1] = 1;
  bc_state_[2] = 1;
}

Bond::Bond(const int& type, const int& ngb, const Site& src, const Vector3i& src_offset, 
  const Site& tgt, const Vector3i& tgt_offset)
  : std::pair<Site,Site>(src,tgt)
{
  assert(type>=0);
  assert(ngb>=1);
  this->first.reset_bravindex(src_offset);
  this->second.reset_bravindex(tgt_offset);
  id_ = num_bond++;
  type_ = type;
  ngb_ = ngb;
  bravindex_ = Vector3i(0,0,0);
  sign_ = 1;
  bc_state_[0] = 1; 
  bc_state_[1] = 1;
  bc_state_[2] = 1;
}

void Bond::connect(const Site& src, const Vector3i& src_offset, const Site& tgt, 
    const Vector3i& tgt_offset, const int& sign)
{
  this->first = src;
  this->second = tgt;
  // It seems, we need to SWITCH OFF the next 2 lines
  // 2026-01-16: The two lines are needed.
  this->first.reset_bravindex(src_offset);
  this->second.reset_bravindex(tgt_offset);
  sign_ = sign;
}

void Bond::connect(const Site& src, const Vector3i& src_offset, const Site& tgt, 
    const Vector3i& tgt_offset, const BC_state& bstate)
{
  this->first = src;
  this->second = tgt;
  // It seems, we need to SWITCH OFF the next 2 lines
  // 2026-01-16: The two lines are needed.
  this->first.reset_bravindex(src_offset);
  this->second.reset_bravindex(tgt_offset);
  bc_state_ = bstate;
}

// friends
std::ostream& operator<<(std::ostream& os, const Bond& b) 
{
  os << std::fixed << std::showpoint;
  os << b.src_id() << " - " << b.tgt_id();
  //os << "bond " << b.id_ << ": " << b.src_ << " (" << b.src_offset_(0) << "," <<  b.src_offset_(1) << ","
  //  << b.src_offset_(2) << ") --- " << b.tgt_ << " (" << b.tgt_offset_(0) << "," 
  //  <<  b.tgt_offset_(1) << "," << b.tgt_offset_(2) << ")" << std::endl;
  //os << "type = " << b.type_ << ", ngb = " << b.ngb_ << "\n";
  return os;
}

/*---------------Unitcell class-----------------*/
void Unitcell::clear(void) 
{ 
  max_site_type_val=0; 
  max_bond_type_val=0; 
  max_neighb_val=0; 
  sites_.clear(); 
  bonds_.clear(); 
  Site::reset_count(); 
  Bond::reset_count(); 
  num_orbitals_ = 0;
  a1 = Vector3d(0.0, 0.0, 0.0); 
  a2 = Vector3d(0.0, 0.0, 0.0); 
  a3 = Vector3d(0.0, 0.0, 0.0); 
  sitetypes_map_.clear();
  bondtypes_map_.clear();
}

void Unitcell::set_basis(const Vector3d& av1, const Vector3d& av2, const Vector3d& av3)
{ 
  a1 = av1; a2 = av2; a3 = av3;
}

int Unitcell::add_site(const int& type, const int& num_orbitals, const Vector3d& site_coord)
{
  int uid = sites_.size();
  sites_.push_back(Site(uid, type, num_orbitals, bravindex(), site_coord, coord())); 
  if (max_site_type_val < type) max_site_type_val = type;
  return sites_.back().id(); 
}

int Unitcell::add_site(const int& type, const Vector3d& site_coord)
{
  int uid = sites_.size();
  sites_.push_back(Site(uid, type, 1, bravindex(), site_coord, coord())); 
  if (max_site_type_val < type) max_site_type_val = type;
  return sites_.back().id(); 
}

int Unitcell::add_site(const Site& s, const Vector3i& bravindex, const Vector3d& cell_coord)
{
  sites_.push_back(s); 
  sites_.back().reset_bravindex(bravindex);
  sites_.back().reset_cell_coord(cell_coord);
  return sites_.back().id();
}

int Unitcell::add_bond(const int& type, const int& ngb, const int& src_id, 
  const Vector3i& src_offset, const int& tgt_id, const Vector3i& tgt_offset)
{
  if (src_id >= sites_.size()) throw std::range_error("*error: add_bond:: 'src' site does not exist");
  if (tgt_id >= sites_.size()) throw std::range_error("*error: add_bond:: 'tgt' site does not exist");
  //bonds.push_back(Bond(type, ngb, bravindex(), src_id, src_offset, tgt_id, tgt_offset, 1));
  bonds_.push_back(Bond(type, ngb, sites_[src_id], src_offset, sites_[tgt_id], tgt_offset));
  if (max_bond_type_val < type) max_bond_type_val = type;
  if (max_neighb_val < ngb) max_neighb_val = ngb;
  return bonds_.back().id();
}

/*
int Unitcell::add_bond(const int& type, const int& src_id, 
  const Vector3i& src_offset, const int& tgt_id, const Vector3i& tgt_offset)
{
  if (src_id >= sites_.size()) throw std::range_error("*error: add_bond:: 'src' site does not exist");
  if (tgt_id >= sites_.size()) throw std::range_error("*error: add_bond:: 'tgt' site does not exist");
  //bonds.push_back(Bond(type, ngb, bravindex(), src_id, src_offset, tgt_id, tgt_offset, 1));
  int ngb = 1;
  bonds_.push_back(Bond(type, ngb, sites_[src_id], src_offset, sites_[tgt_id], tgt_offset));
  if (max_bond_type_val < type) max_bond_type_val = type;
  if (max_neighb_val < ngb) max_neighb_val = ngb;
  return bonds_.back().id();
}
*/


/*int Unitcell::add_bond(const Bond& bond, const Vector3i& src_offset, const Vector3i& tgt_offset)
{
  //bonds.push_back(Bond(bond.type(), bond.ngb(), bond.src_id(), src_offset, tgt_id, tgt_offset, vector));
  bonds.push_back(bond);
  bonds.back().reset_src_offset(src_offset);
  bonds.back().reset_tgt_offset(tgt_offset);
  return bonds.back().id();
}*/

void Unitcell::finalize(void) 
{
  // Re-assign the sites 'type' values making them contiguous
  sitetypes_map_.clear();
  // store the sitetype map
  std::set<int> types;
  // type values sorted in the set
  for (const auto& s : sites_) types.insert({s.type()});
  // map to contiguous indices
  int i = 0;
  for (const auto& t : types) {
    auto status = sitetypes_map_.insert({t, i});
    if (status.second) i++;
  }
  // now reassign the values
  for (auto& s : sites_) {
    int i = sitetypes_map_[s.type()];
    s.reset_type(i);
  }
  // number of orbitals
  num_orbitals_ = 0;
  for (auto& s : sites_) num_orbitals_ += s.num_orbitals();

  // same for bonds
  bondtypes_map_.clear();
  types.clear();
  for (const auto& b : bonds_) types.insert({b.type()});
  i = 0;
  for (const auto& t : types) {
    auto status = bondtypes_map_.insert({t, i});
    if (status.second) i++;
  }
  for (auto& b : bonds_) {
    int i = bondtypes_map_[b.type()];
    b.reset_type(i);
  }
  // type values within range?
  if (sitetypes_map_.size() >= MAX_SITE_TYPES) 
    throw std::range_error("error: latticelibrary: number of 'site types' exceed limit");
  if (bondtypes_map_.size() >= MAX_BOND_TYPES) 
    throw std::range_error("error: latticelibrary: number of 'bond types' exceed limit");
}

void Unitcell::update_num_orbitals(void) 
{
  // number of orbitals
  num_orbitals_ = 0;
  for (auto& s : sites_) num_orbitals_ += s.num_orbitals();
}


/*void Unitcell::reset(const std::vector<Site>& new_sites, const std::vector<Bond>& new_bonds)
{
  sites.clear();
  bonds.clear();
  sites = new_sites;
  bonds = new_bonds;
  for (int i=0; i<sites.size(); ++i) {
    sites[i].reset_uid(i);
    sites[i].reset_bravindex(Vector3i(0,0,0));
    sites[i].reset_cell_coord(Vector3d(0.0,0.0,0.0));
  }
  for (int i=0; i<bonds.size(); ++i) bonds[i].reset_bravindex(Vector3i(0,0,0));
}*/

void Unitcell::translate_by(const Vector3i& bravindex_offset, const int& cell_id_offset) 
{
  Vector3d coord_offset = bravindex_offset(0) * a1 + bravindex_offset(1) * a2 
                        + bravindex_offset(2) * a3;
  int id_offset = cell_id_offset * sites_.size();
  for (auto& s : sites_) s.translate_by(id_offset, bravindex_offset, coord_offset);
  for (auto& b : bonds_) b.translate_by(bravindex_offset);
  bravindex_ += bravindex_offset;
  coord_ += coord_offset;
}

void Unitcell::rotate_by(const Eigen::Matrix3d& matrix)
{
  // rotate the basis vectors
  a1 = matrix * a1;
  a2 = matrix * a2;
  a3 = matrix * a3;
  // rotate site coordinates
  Vector3d rv;
  for (int i=0; i<sites_.size(); ++i) {
    rv = matrix * sites_[i].coord();
    sites_[i].reset_coord(rv);
  }
}

/*--------------- Implementation of 'Lattice' class-----------------*/
void Lattice::set_basis_vectors(const Vector3d& av1, const Vector3d& av2, const Vector3d& av3)
{
  Unitcell::set_basis(av1, av2, av3); 
}

int Lattice::add_basis_site(const int& type, const int& num_orbitals, const Vector3d& site_coord)
{
  return Unitcell::add_site(type, num_orbitals, site_coord);  
}

int Lattice::add_basis_site(const int& type, const Vector3d& site_coord)
{
  return Unitcell::add_site(type, site_coord);  
}

int Lattice::add_bond(const int& type, const int& ngb, const int& src_id, const Vector3i& src_offset,
    const int& tgt_id, const Vector3i& tgt_offset)
{
  return Unitcell::add_bond(type, ngb, src_id, src_offset, tgt_id, tgt_offset);
}

int Lattice::add_bond(const int& type, const int& src_id, const Vector3i& src_offset,
    const int& tgt_id, const Vector3i& tgt_offset)
{
  return Unitcell::add_bond(type, 1, src_id, src_offset, tgt_id, tgt_offset);
}

/*
int Lattice::add_bond(const int& type, const int& src_id, 
    const int& tgt_id, const Vector3i& tgt_offset, const int& ngb)
{
  return Unitcell::add_bond(type, src_id, tgt_id, tgt_offset, ngb);
}
*/

int Lattice::add_bond(const Bond& b)
{ 
  return Unitcell::add_bond(b);
}

Unitcell Lattice::get_translated_cell(const Vector3i& bravindex_offset) const
{
  Unitcell newcell(*this);
  int cell_id_offset = bravindex_offset[0] + bravindex_offset[1] * extent[dim1].size 
                     + bravindex_offset[2] * num_layer_cells_;
  newcell.translate_by(bravindex_offset, cell_id_offset);
  return newcell;
}


boundary_type Lattice::boundary_condition(std::string& bc) const
{
  boost::to_upper(bc);
  if (bc=="OPEN") {
    return boundary_type::open;
  } 
  else if (bc=="PERIODIC") {
    return boundary_type::periodic;
  }
  else if (bc=="ANTIPERIODIC") {
    return boundary_type::antiperiodic;
  }
  else {
    throw std::range_error("error: latticelibrary: invalid boundary condition");
  }
}

Vector3i Lattice::get_next_bravindex(const Vector3i& current_index) const
{
  /* Returns the next Bravais lattice index. 
  ! Index for first unit cell = (0,0,0)
  ! Index for last unit cell = (N1-1, N2-1, N3-1)
   */
  Vector3i next_index = current_index;
  if (++next_index[0] >= static_cast<int>(extent[dim1].size)) {
    next_index[0] = 0;
    if (++next_index[1] >= static_cast<int>(extent[dim2].size)) {
      next_index[1] = 0;
      if (++next_index[2] >= static_cast<int>(extent[dim3].size)) {
        next_index[2] = 0;
      }
    }
  }
  return next_index;
}


/*
bool Lattice::connect_bond(Bond& bond, const std::vector<Site>& sites) const
{
  // source site
  int sign1, sign2;
  Vector3i src_cell, src_offset;
  src_offset = bond.first.bravindex();
  //src_offset = bond.src_offset();
  src_cell = bond.bravindex() + src_offset;
  // if the source is inside the lattice, offset = 0
  for (int dim=dim1; dim<=dim3; ++dim) {
    if (src_cell[dim] >= 0 && src_cell[dim] < static_cast<int>(extent[dim].size)) src_offset[dim]=0;
  }
  // id of the source site
  boost::tie(src_cell, sign1) = boundary_wrap(src_cell);
  int src_id = mapped_site_id(bond.src_id(), src_cell);
  if (src_id < 0) return false;  // bond can't be connected due to open boundary

  // target site
  Vector3i tgt_cell, tgt_offset;
  tgt_offset = bond.second.bravindex();
  //tgt_offset = bond.tgt_offset();
  tgt_cell = bond.bravindex() + tgt_offset;
  // if the target is inside the lattice, offset = 0
  for (int dim=dim1; dim<=dim3; ++dim) {
    if (tgt_cell[dim] >= 0 && tgt_cell[dim] < static_cast<int>(extent[dim].size)) tgt_offset[dim]=0;
  }
  // id of the target site
  boost::tie(tgt_cell, sign2) = boundary_wrap(tgt_cell);
  int tgt_id = mapped_site_id(bond.tgt_id(), tgt_cell);
  if (tgt_id < 0) return false;  // bond can't be connected due to open boundary

  // connect the bond
  //int s = static_cast<int>(src_id);
  //int t = static_cast<int>(tgt_id);
  int sign = sign1 * sign2;
  bond.connect(sites[src_id], src_offset, sites[tgt_id], tgt_offset, sign);
  return true;
}
*/

std::pair<Vector3i, int> Lattice::boundary_wrap(const Vector3i& cell_idx) const
{
  Vector3i new_idx(cell_idx);
  int phase = 1;
  for (int dim=dim1; dim<=dim3; ++dim) {
    int size = static_cast<int>(extent[dim].size); 
    if (new_idx[dim]>=size) {
      if (extent[dim].bc == boundary_type::periodic) {
        while (new_idx[dim]>=size) {
          new_idx[dim] -= size;
          if (extent[dim].periodicity==boundary_type::antiperiodic) 
            phase=-phase;
        }
      }
      else new_idx[dim] = -1;
    }
    else if (new_idx[dim]<0) {
      if (extent[dim].bc == boundary_type::periodic) {
        while (new_idx[dim]<0) {
          new_idx[dim] += size;
          if (extent[dim].periodicity==boundary_type::antiperiodic) 
            phase=-phase;
        }
      }
      else new_idx[dim] = -1;
    }
  }
  return std::pair<Vector3i, int>(new_idx, phase);
}

bool Lattice::connect_bond2(Bond& bond, const std::vector<Site>& sites) const
{
  // source site
  BC_state bstate1, bstate2;
  Vector3i src_cell, src_offset;
  src_offset = bond.first.bravindex();
  //src_offset = bond.src_offset();
  src_cell = bond.bravindex() + src_offset;
  // if the source is inside the lattice, offset = 0
  for (int dim=dim1; dim<=dim3; ++dim) {
    if (src_cell[dim] >= 0 && src_cell[dim] < static_cast<int>(extent[dim].size)) src_offset[dim]=0;
  }
  // id of the source site
  boost::tie(src_cell, bstate1) = boundary_wrap2(src_cell);
  int src_id = mapped_site_id(bond.src_id(), src_cell);
  if (src_id < 0) return false;  // bond can't be connected due to open boundary

  // target site
  Vector3i tgt_cell, tgt_offset;
  tgt_offset = bond.second.bravindex();
  //tgt_offset = bond.tgt_offset();
  tgt_cell = bond.bravindex() + tgt_offset;
  // if the target is inside the lattice, offset = 0
  for (int dim=dim1; dim<=dim3; ++dim) {
    if (tgt_cell[dim] >= 0 && tgt_cell[dim] < static_cast<int>(extent[dim].size)) tgt_offset[dim]=0;
  }
  // id of the target site
  boost::tie(tgt_cell, bstate2) = boundary_wrap2(tgt_cell);
  int tgt_id = mapped_site_id(bond.tgt_id(), tgt_cell);
  if (tgt_id < 0) return false;  // bond can't be connected due to open boundary

  //std::cout << "connect_bond2: tgt_offset = " << tgt_offset.transpose() << "\n";
  //std::cout << "bond connected\n";
  //std::cout << "bond connected\n";
  //getchar();

  // connect the bond
  //int s = static_cast<int>(src_id);
  //int t = static_cast<int>(tgt_id);
  BC_state bstate;
  for (int dim=dim1; dim<=dim3; ++dim) {
    bstate[dim] = bstate1[dim] * bstate2[dim];
  }
  bond.connect(sites[src_id],src_offset,sites[tgt_id],tgt_offset,bstate);
  return true;
}

std::pair<Vector3i, BC_state> Lattice::boundary_wrap2(const Vector3i& cell_idx) const
{
  Vector3i new_idx(cell_idx);
  BC_state bstate;
  for (int dim=dim1; dim<=dim3; ++dim) {
    bstate[dim] = 1;
    int size = static_cast<int>(extent[dim].size); 
    if (new_idx[dim]>=size) {
      if (extent[dim].bc == boundary_type::periodic) {
        while (new_idx[dim]>=size) {
          new_idx[dim] -= size;
          bstate[dim] = -bstate[dim];
        }
      }
      else new_idx[dim] = -1;
    }
    else if (new_idx[dim]<0) {
      if (extent[dim].bc == boundary_type::periodic) {
        while (new_idx[dim]<0) {
          new_idx[dim] += size;
          bstate[dim] = -bstate[dim];
        }
      }
      else new_idx[dim] = -1;
    }
  }
  return std::pair<Vector3i, BC_state>(new_idx, bstate);
}


int Lattice::mapped_site_id(const int& local_id, const Vector3i& bravindex) const
{
  int cell_id = bravindex[0] + bravindex[1] * extent[dim1].size + bravindex[2] * num_layer_cells_;
  if (cell_id < 0) return -1;
  //return (static_cast<int>(local_id) + cell_id * num_basis_sites());
  // the above line caused a bug
  return (static_cast<int>(local_id) + cell_id * num_basis_sites_);
}

int Lattice::translation_mapped_site(const int& uid, 
  const Vector3i& bravindex, const Vector3i& translation_vec) const
{
  if (uid > Unitcell::num_sites()) return uid;
  Vector3i translated_cell = bravindex + translation_vec;
  int sign;
  boost::tie(translated_cell, sign) = boundary_wrap(translated_cell);
  int cell_id = translated_cell[0] + translated_cell[1] * extent[dim1].size 
                   + translated_cell[2] * num_layer_cells_;
  int mapped_id = uid + cell_id * num_basis_sites();
  if (mapped_id < this->num_sites()) return mapped_id;
  else return uid;
}

const Site& Lattice::translated_site(const Site& site, const Vector3i& translation_vec) const
{
  int id = translation_mapped_site(site.uid(), site.bravindex(), translation_vec); 
  return sites_[id];
}

Eigen::Matrix3d Lattice::rotation_matrix(const Vector3d& r, const Vector3d& rp)
{
  /* Calculates rotation matrix which would rotate vector 'r' 
  ! to align it to common origin vector 'rp'. Rotation axis is 
  ! along 'rp x r'. For the formula used to calculate the matrix,
  ! see H. Goldstein, Classical Mechanics, Eq. 4-92 through 4-96. 
  */

  double e1 = std::sqrt(r.dot(r));
  double e2 = std::sqrt(rp.dot(rp));
  double phi = std::acos(r.dot(rp)/(e1*e2));
  if (std::abs(phi) < dp_tol) return Eigen::Matrix<double,3,3>::Identity();

  // unit vector perpendicular to 'rp' and 'r' along \vec{r'} x \vec{r}
  Vector3d nhat = rp.cross(r);
  nhat = nhat/sqrt(nhat.dot(nhat));

  // the 'e' parameters
  double e0, e3;
  phi = phi*0.50;
  e0 = cos(phi); e1 = nhat(0)*sin(phi); e2 = nhat(1)*sin(phi); e3 = nhat(2)*sin(phi);

  // the rotation matrix
  Eigen::Matrix3d mat;
  mat(0,0) = e0*e0 + e1*e1 - e2*e2 - e3*e3;
  mat(0,1) = 2.0*(e1*e2 + e0*e3);
  mat(0,2) = 2.0*(e1*e3 - e0*e2);

  mat(1,0) = 2.0*(e1*e2 - e0*e3);
  mat(1,1) = e0*e0 - e1*e1 + e2*e2 - e3*e3;
  mat(1,2) = 2.0*(e2*e3 + e0*e1);

  mat(2,0) = 2.0*(e1*e3 + e0*e2);
  mat(2,1) = 2.0*(e2*e3 - e0*e1);
  mat(2,2) = e0*e0 - e1*e1 - e2*e2 + e3*e3;

  return mat;
}





} // end namespace lattice














