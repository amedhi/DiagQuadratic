/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2025-12-06 11:29:29
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-04-06 12:11:18
*----------------------------------------------------------------------------*/
#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <string>
#include <tuple>
#include <array>
#include <set>
#include <stdexcept>
#include "../scheduler/task.h"
#include "../diag/matrix.h"
#include <Eigen/Dense>
#include <boost/graph/adjacency_list.hpp>
#include "constants.h"

namespace lattice {

// Global Constants
const int MAX_SITE_TYPES = 100;
const int MAX_BOND_TYPES = 200;

/*---------------Bravais id-----------------*/
enum class brav_id {ZERO, CHAIN, SQUARE, RECTANGULAR, RECTANGULAR_C, HEXAGONAL, 
OBLIQUE, CUBIC, CUBIC_BCC, CUBIC_FCC
};

/*---------------lattice types-----------------*/
enum class lattice_id {
  UNDEFINED, SQUARE, SQUARE_NNN, SQUARE_2BAND, SQUARE_2SITE, SQUARE_4SITE, CHAIN, CHAIN_2SITE,
  HONEYCOMB, HONEYCOMB2, HONEYCOMB3, SW_GRAPHENE, SIMPLECUBIC, NICKELATE_2B, NICKELATE_2D, NICKELATE_2L,
  NICKELATE_4SITE, SQUARE_CDW4, SQUARE_STRIPE
};

/*---------------Lattice site class-----------------*/
using Vector3i = Eigen::Vector3i;
using Vector3d = Eigen::Vector3d;
using BC_state = std::array<int,3>; // boundary crossing state of a bond

//using BrvaisIdx = Eigen::Matrix<int, 3, 1>;

class Site 
{
public:
  // ctors
  Site() {}
  Site(const int& uid, const int& type, const int& num_orbitals, 
    const Vector3i& bravindex, const Vector3d& coord, const Vector3d& cell_coord);
  ~Site() {}
  // setter functions
  static void reset_count(void) { num_site=0; }
  void reset_id(const int& id) { id_= id; }
  void reset_type(const int& t) { type_=t; }
  void reset_uid(const int& uid) { uid_=uid; }
  void reset_bravindex(const Vector3i& idx) { bravindex_=idx; }
  void reset_coord(const Vector3d& v) { coord_=v; }
  void reset_cell_coord(const Vector3d& v) { cell_coord_=v; }
  void translate_by(const int& id_offset, const Vector3i& bravindex_offset, const Vector3d& coord_offset); 
  void clear_bonds(void) { out_bonds_.clear(); in_bonds_.clear(); }
  void add_out_bond(const int& bond_id) { out_bonds_.push_back(bond_id); }
  void add_in_bond(const int& bond_id) { in_bonds_.push_back(bond_id); }

  // getter functions
  const int& id(void) const {return id_;}
  const int& uid(void) const {return uid_;}
  const int& type(void) const {return type_;}
  const int& num_orbitals(void) const {return num_orbitals_;}
  const Vector3i& bravindex(void) const { return bravindex_; }
  const Vector3d& coord(void) const {return coord_;}
  const Vector3d& cell_coord(void) const {return cell_coord_;}
  const std::vector<int>& outbond_ids(void) const { return out_bonds_; }
  const std::vector<int>& inbond_ids(void) const { return out_bonds_; }
  // friends
  friend std::ostream& operator<<(std::ostream& os, const Site& site);
private:
  static int num_site;
  int id_{0};
  int uid_{0}; // local id within a unitcell
  int type_{0};
  int num_orbitals_{1};
  Vector3i bravindex_{Vector3i(0, 0, 0)};
  Vector3d coord_{Vector3d(0.0, 0.0, 0.0)};
  Vector3d cell_coord_{Vector3d(0.0, 0.0, 0.0)};
  std::vector<int> out_bonds_;
  std::vector<int> in_bonds_;
};

/*---------------Lattice bond class-----------------*/
class Bond : public std::pair<Site, Site> 
{
public:
  // ctors
  Bond(); // {}
  Bond(const int& type, const int& ngb, const Site& src, const Vector3i& bravindex,
    const Site& tgt, const Vector3i& tgt_offset);
  //Bond(const int& type, const int& ngb, const Vector3i& bravindex, const int& src_id, 
  //  const Vector3i& src_offset, const int& tgt_id, const Vector3i& tgt_offset, const int& sign);
  ~Bond() {}
  // setter functions
  static void reset_count(void) { num_bond=0; }
  void reset_id(const int& id) { id_=id; }
  void reset_type(const int& t) { type_=t; }
  void reset_sign(const int& s) { sign_=s; }
  void reset_phase(const std::complex<double>& ph) { phase_=ph; }
  //void reset_src_offset(const Vector3i& idx) { src_offset_=idx; }
  //void reset_tgt_offset(const Vector3i& idx) { tgt_offset_=idx; }
  void reset_bravindex(const Vector3i& idx) { bravindex_=idx; }
  void set_vector_id(const int& id) { vector_id_=id; }
  void set_vector(const Vector3d& R) { vector_=R; }
  //void shift_target_ids(const int& id_offset) { src_ += id_offset; tgt_ += id_offset; }
  void translate_by(const Vector3i& bravindex_offset) { bravindex_ += bravindex_offset; } 
  void connect(const Site& src, const Vector3i& src_offset, const Site& tgt, 
    const Vector3i& tgt_offset, const int& sign);
  void connect(const Site& src, const Vector3i& src_offset, const Site& tgt, 
    const Vector3i& tgt_offset, const BC_state& bstate);
  //void connect(const int& src_id, const Vector3i& src_offset, const int& tgt_id, 
  //  const Vector3i& tgt_offset, const int& sign);
  // getter functions
  const int& id(void) const { return id_; }
  const int& type(void) const {return type_;}
  const int& ngb(void) const {return ngb_;}
  int src_id(void) const { return first.id(); }
  int tgt_id(void) const { return second.id(); }
  const Site& src(void) const { return first; }
  const Site& tgt(void) const { return second; }
  const int& vector_id(void) const { return vector_id_; }
  const Vector3d& vector(void) const { return vector_; }
  const BC_state& bc_state(void) const { return bc_state_; }
  const int& sign(void) const { return sign_; }
  const std::complex<double>& phase(void) const { return phase_; }
  Vector3i bravindex(void) const { return bravindex_; }
  //Vector3i src_offset(void) const { return src_offset_; }
  //Vector3i tgt_offset(void) const { return tgt_offset_; }
  // friends
  friend std::ostream& operator<<(std::ostream& os, const Bond& bond);
private:
  static int num_bond;
  int id_{0};
  int type_{0};
  int ngb_{0};
  int vector_id_{0}; // integer id for the following vector
  Vector3d vector_{Vector3d(0,0,0)}; // coordinate of 'tgt cell' wrt 'src cell'
  //int src_ {0}; 
  //int tgt_ {0}; 
  int sign_{1}; // = -1 if across an antiperiodic boundary
  BC_state bc_state_; // bc_state_[dim]=1 for normal bond, 
                    // bc_state_[dim]=-1 for bond crossing 'dim' boundary
  std::complex<double> phase_; // bc phase
  Vector3i bravindex_ {Vector3i(0, 0, 0)};
  //Vector3i src_offset_ {Vector3i(0, 0, 0)};
  //Vector3i tgt_offset_ {Vector3i(0, 0, 0)};
};

/*---------------Unitcell class-----------------*/
class Unitcell 
{
public:
  // ctors
  Unitcell() {}
  ~Unitcell() {}
  // setter functions
  int add_site(const int& type, const Vector3d& site_coord); 
  int add_site(const int& type, const int& num_orbitals, const Vector3d& site_coord); 
  int add_site(const Site& s) { sites_.push_back(s); return sites_.back().id(); }
  int add_site(const Site& s, const Vector3i& bravindex, const Vector3d& cell_coord);
  int add_bond(const Bond& b) { bonds_.push_back(b); return bonds_.back().id(); }
  int add_bond(const int& type, const int& ngb, const int& src_id, const Vector3i& src_offset,
    const int& tgt_id, const Vector3i& tgt_offset); 
  int add_bond(const int& type, const int& src_id, const Vector3i& src_offset,
    const int& tgt_id, const Vector3i& tgt_offset); 
  //int add_bond(const int& type, const int& src_id, const int& tgt_id, const Vector3i& tgt_offset);
  void set_basis(const Vector3d& av1, const Vector3d& av2, const Vector3d& av3);
  void reset_a1(const Vector3d& av) { a1=av; }
  void reset_a2(const Vector3d& av) { a2=av; }
  void reset_a3(const Vector3d& av) { a3=av; }
  void clear(void); 
  void clear_sites(void) { sites_.clear(); }
  void clear_bonds(void) { bonds_.clear(); }
  void finalize(void);
  void update_num_orbitals(void);
  // getter functions
  int num_sites(void) const { return sites_.size(); }
  const int& num_orbitals(void) const { return num_orbitals_; }
  int num_bonds(void) const { return bonds_.size(); }
  const Vector3d& vector_a1(void) const { return a1; }
  const Vector3d& vector_a2(void) const { return a2; }
  const Vector3d& vector_a3(void) const { return a3; }
  const Vector3i& bravindex(void) const { return bravindex_; }
  const Vector3d& coord(void) const {return coord_;}
  const Site& site(const int& i) const { return sites_[i]; }
  const Bond& bond(const int& i) const { return bonds_[i]; }
  Bond& bond(const int& i) { return bonds_[i]; }
  int num_site_types(void) const { return sitetypes_map_.size(); }
  int num_bond_types(void) const { return bondtypes_map_.size(); }
  std::map<int,int> sitetypes_map(void) const { return sitetypes_map_; }
  std::map<int,int> bondtypes_map(void) const { return bondtypes_map_; }
  // transform
  void translate_by(const Vector3i& bravindex_offset, const int& cell_id_offset);
  void rotate_by(const Eigen::Matrix3d& matrix);
private:
  //int id {0};
  int max_site_type_val {0};
  int max_bond_type_val {0};
  int max_neighb_val {0};
  /*------------------------------------------------------- 
   * A 'site' means an 'orbital'. In single-orbital systems, 
   * an 'atom' is same as a 'site'. In multi-orbital systems,
   * the number of 'atom' is obviously not same as the 
   * number of sites in the unitcell.
   *-------------------------------------------------------*/
  int num_orbitals_{0};
  Vector3d a1 {Vector3d(0.0, 0.0, 0.0)};
  Vector3d a2 {Vector3d(0.0, 0.0, 0.0)};
  Vector3d a3 {Vector3d(0.0, 0.0, 0.0)};
  Vector3i bravindex_ {Vector3i(0, 0, 0)};
  Vector3d coord_ {Vector3d(0.0, 0.0, 0.0)};
  std::vector<Site> sites_;
  std::vector<Bond> bonds_;
  std::map<int,int> sitetypes_map_; // user set value to contiguous value 
  std::map<int,int> bondtypes_map_; // user set value to contiguous value 
};

/*---------------spatial dimension type-----------------*/
enum class boundary_type {open, periodic, antiperiodic};


class Lattice : public Unitcell
{
public:
  using Atom = std::vector<int>; 
  // ctors
  Lattice() {}
  Lattice(const input::Parameters& parms) { construct(parms); }
  ~Lattice() {}
  // setter functions
  int construct(const input::Parameters& parms);
  void set_basis_vectors(const Vector3d& av1, const Vector3d& av2, const Vector3d& av3);
  int add_basis_site(const int& type, const int& num_orbitals, const Vector3d& site_coord);
  int add_basis_site(const int& type, const Vector3d& site_coord);
  int add_bond(const int& type, const int& ngb, const int& src_id, const Vector3i& src_offset,
    const int& tgt_id, const Vector3i& tgt_offset); 
  int add_bond(const int& type, const int& src_id, const Vector3i& src_offset,
    const int& tgt_id, const Vector3i& tgt_offset); 
  //int add_bond(const int& type, const int& src_id, 
  //  const int& tgt_id, const Vector3i& tgt_offset, const int& ngb=1);
  int add_bond(const Bond& b); 
  /*
  void set_boundary_twist(const int& twist_id) 
  {
    for (int dim=dim1; dim<=dim3; ++dim) {
      extent[dim].bc_twist = twist_angles_(twist_id,dim);
      //std::cout << "twist = " << extent[dim].bc_twist << "\n";
    }
    //std::cout << "\n";
  }
  */
  int reset_boundary_twist(const int& twist_id);

  // getter functions
  std::string name(void) const { return lname; }
  lattice_id id(void) const { return lid; }
  brav_id brav(void) const { return ibrav; }
  int num_basis_sites(void) const { return Unitcell::num_sites(); }
  int num_basis_orbitals(void) const { return Unitcell::num_orbitals(); }
  int num_basis_bonds(void) const { return Unitcell::num_bonds(); }
  const int& dimension(void) const { return spatial_dim; }
  const int& num_sites(void) const { return num_total_sites_; }
  const int& num_unitcells(void) const { return num_total_cells_; }
  const int& num_bonds(void) const { return num_bonds_; }
  const int& num_atoms(void) const { return num_atoms_; }
  Vector3d basis_vector_a1(void) const { return Unitcell::vector_a1(); }
  Vector3d basis_vector_a2(void) const { return Unitcell::vector_a2(); }
  Vector3d basis_vector_a3(void) const { return Unitcell::vector_a3(); }
  int size(const int& idir) const {
    switch (idir) {
      case 0: return static_cast<int>(extent[dim1].size); break;
      case 1: return static_cast<int>(extent[dim2].size); break;
      case 3: return static_cast<int>(extent[dim3].size); break;
      default: throw std::range_error("error: Lattice::size: out-of-range argument");
    }
  }
  int size1(void) const { return static_cast<int>(extent[dim1].size); }
  int size2(void) const { return static_cast<int>(extent[dim2].size); }
  int size3(void) const { return static_cast<int>(extent[dim3].size); }
  int input_size1(void) const { return static_cast<int>(copy_extent[dim1].size); }
  int input_size2(void) const { return static_cast<int>(copy_extent[dim2].size); }
  int input_size3(void) const { return static_cast<int>(copy_extent[dim3].size); }
  const int& num_boundary_twists(void) const { return num_total_twists_; }
  const double& bc1_twist(void) const { return extent[dim1].bc_twist; }
  const double& bc2_twist(void) const { return extent[dim2].bc_twist; }
  const double& bc3_twist(void) const { return extent[dim3].bc_twist; }
  boundary_type bc1(void) const { return extent[dim1].bc; }
  boundary_type bc2(void) const { return extent[dim2].bc; }
  boundary_type bc3(void) const { return extent[dim3].bc; }
  //boundary_type bc2(void) const { return extent[dim2].bc; }
  //boundary_type bc3(void) const { return extent[dim3].bc; }
  boundary_type bc1_periodicity(void) const { return extent[dim1].periodicity; }
  boundary_type bc2_periodicity(void) const { return extent[dim2].periodicity; }
  boundary_type bc3_periodicity(void) const { return extent[dim3].periodicity; }
  const Site& site(const int& i) const { return sites_[i]; }
  const Bond& bond(const int& i) const { return bonds_[i]; }
  const std::vector<Site>& sites(void) const { return sites_; }
  const std::vector<Bond>& bonds(void) const { return bonds_; }

  // other methods 
  //Vector3i boundary_wrap(const Vector3i& cell_idx) const;
  std::pair<Vector3i,int> boundary_wrap(const Vector3i& cell_idx) const;
  std::pair<Vector3i,BC_state> boundary_wrap2(const Vector3i& cell_idx) const;
  Vector3i get_next_bravindex(const Vector3i& current_index) const;
  Unitcell get_translated_cell(const Vector3i& bravindex_offset) const;
  int mapped_site_id(const int& local_id, const Vector3i& bravindex) const;
  int translation_mapped_site(const int& uid, const Vector3i& bravindex,
    const Vector3i& translation_vec) const;
  const Site& translated_site(const Site& site, const Vector3i& translation_vec) const;
  //bool connect_bond(Bond& bond, const std::vector<Site>& sites) const;
  bool connect_bond2(Bond& bond, const std::vector<Site>& sites) const;
  const Site& basis_site(const int& i) const { return Unitcell::site(i); }
  const Bond& basis_bond(const int& i) const { return Unitcell::bond(i); }

  // for lattices with doped impurities
  //void add_new_bondtype(const int& type);
private:
  struct Extent {int size; boundary_type bc; boundary_type periodicity; double bc_twist;};
  enum Dimension {dim1, dim2, dim3};
  brav_id ibrav {brav_id::SQUARE};
  lattice_id lid {lattice_id::SQUARE};
  std::string lname {""};
  int spatial_dim {0};

  // lattice dimensions
  Extent extent[3] = {Extent{1, boundary_type::open, boundary_type::open, 0.0}, 
                      Extent{1, boundary_type::open, boundary_type::open, 0.0},
                      Extent{1, boundary_type::open, boundary_type::open, 0.0}
                     };

  // copy of user-set lattice dimensions
  Extent copy_extent[3] {Extent{1, boundary_type::open, boundary_type::open, 0.0}, 
                         Extent{1, boundary_type::open, boundary_type::open, 0.0},
                         Extent{1, boundary_type::open, boundary_type::open, 0.0}
                        };
  
  // number of unit cells in total and in one layer (for symmetrized lattice)
  int num_total_cells_{1};
  int num_layer_cells_{1};
  int num_total_sites_{0};
  int num_basis_sites_{0};
  int num_bonds_{0};
  int num_atoms_{0};

  int num_total_twists_{1};
  RealMatrix twist_angles_;

  // sites & bonds
  std::vector<Site> sites_;
  std::vector<Bond> bonds_;

  // set of type values
  //std::set<int> sitetype_set_;
  //std::set<int> bondtype_set_;

  // for lattices with impurities
  std::vector<Site> impurity_sites_;
  std::vector<Bond> impurity_bonds_;

  // helper functions
  int define_lattice(void); 
  int finalize_lattice(void); 
  int symmetrize_lattice(void);
  int construct_graph(void); 
  int print_lattice(void); 
  boundary_type boundary_condition(std::string& bc) const;
  Eigen::Matrix3d rotation_matrix(const Eigen::Vector3d& r, const Eigen::Vector3d& r_prime);
};


} // end namespace lattice

#endif







