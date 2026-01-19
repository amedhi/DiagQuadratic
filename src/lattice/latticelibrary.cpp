/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2016-01-17 21:32:15
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-17 21:40:51
*----------------------------------------------------------------------------*/
#include <stdexcept>
#include <iomanip>
#include <boost/algorithm/string.hpp>
#include "lattice.h"
//#include "graph.h"

namespace lattice {

// define the lattice
int Lattice::define_lattice(void) 
{
  using pos = Eigen::Vector3i;
  using vec = Eigen::Vector3d;
  int type, src, tgt, orbitals;
  vec a1, a2, a3, coord;
  pos offset, src_offset, tgt_offset, cell;

  /*------------- 'SQUARE' lattice--------------*/
  if (lname == "SQUARE") {
    // type
    lid = lattice_id::SQUARE;
    extent[dim3] = Extent{1, boundary_type::open, boundary_type::open, 0.0};
    // Bravais index
    ibrav = brav_id::SQUARE;
    // in case of 'open' BCs
    if (spatial_dim==0) ibrav = brav_id::ZERO;
    if (spatial_dim==1) ibrav = brav_id::CHAIN;

    // basis vectors
    set_basis_vectors(a1=vec(1,0,0), a2=vec(0,1,0), a3=vec(0,0,0));

    // add sites
    add_basis_site(type=0, coord=vec(0,0,0));

    // add bonds
    add_bond(type=0, src=0, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(1,0,0));
    add_bond(type=1, src=0, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(0,1,0));
  }

  else if (lname == "SQUARE_2BAND") {
    // type
    lid = lattice_id::SQUARE_2BAND;
    extent[dim3] = Extent{1, boundary_type::open, boundary_type::open, 0.0};
    // Bravais index
    ibrav = brav_id::SQUARE;
    // in case of 'open' BCs
    if (spatial_dim==0) ibrav = brav_id::ZERO;
    if (spatial_dim==1) ibrav = brav_id::CHAIN;

    // basis vectors
    set_basis_vectors(a1=vec(1,0,0), a2=vec(0,1,0), a3=vec(0,0,0));
    // add sites
    add_basis_site(type=0, orbitals=2, coord=vec(0,0,0));
    // add bonds
    add_bond(type=0,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,0,0));
    add_bond(type=1,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,1,0));
    add_bond(type=2,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,1,0));
    add_bond(type=3,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(-1,1,0));
  }

  else if (lname=="HONEYCOMB") {
    // type
    lid = lattice_id::HONEYCOMB;
    extent[dim3] = Extent{1, boundary_type::open, boundary_type::open, 0.0};
    ibrav = brav_id::HEXAGONAL;
    if (spatial_dim==0) ibrav = brav_id::ZERO;
    if (spatial_dim==1) ibrav = brav_id::CHAIN;

    // basis vectors
    double aa = 1.0; // A-B bond
    double x = std::sqrt(3.0)*aa;
    double y = 1.5*aa;
    set_basis_vectors(a1=vec(x,0,0), a2=vec(0.5*x,y,0), a3=vec(0,0,0));

    // add sites
    add_basis_site(type=0, orbitals=1, coord=vec(0,0,0));
    add_basis_site(type=1, orbitals=1, coord=vec(0.5*x,0.5*aa,0));

    // NN bonds
    add_bond(type=0,src=0,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,0,0));
    add_bond(type=1,src=1,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,0,0));
    add_bond(type=2,src=1,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,1,0));

    // NNN bonds
    add_bond(type=3,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,0,0));
    add_bond(type=4,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,1,0));
    add_bond(type=3,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(-1,1,0));
    add_bond(type=4,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(1,0,0));
    add_bond(type=3,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,1,0));
    add_bond(type=4,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(-1,1,0));
  }

  else if (lname=="HONEYCOMB2") {
    // type
    lid = lattice_id::HONEYCOMB;
    extent[dim3] = Extent{1, boundary_type::open, boundary_type::open, 0.0};
    ibrav = brav_id::HEXAGONAL;
    if (spatial_dim==0) ibrav = brav_id::ZERO;
    if (spatial_dim==1) ibrav = brav_id::CHAIN;

    // basis vectors
    double aa = 1.0; // A-B bond
    double x = std::sqrt(3.0)*aa;
    double y = 1.5*aa;
    set_basis_vectors(a1=vec(x,0,0), a2=vec(0.5*x,y,0), a3=vec(0,0,0));

    // add sites
    add_basis_site(type=0, orbitals=1, coord=vec(0,0,0));
    add_basis_site(type=1, orbitals=1, coord=vec(0,aa,0));

    // NN bonds
    add_bond(type=0,src=0,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,0,0));
    add_bond(type=1,src=0,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,-1,0));
    add_bond(type=2,src=0,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(1,-1,0));

    // NNN bonds
    add_bond(type=3,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,0,0));
    add_bond(type=4,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,1,0));
    add_bond(type=3,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(-1,1,0));
    add_bond(type=4,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(1,0,0));
    add_bond(type=3,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,1,0));
    add_bond(type=4,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(-1,1,0));
  }

  /*------------- undefined lattice--------------*/
  else {
    throw std::range_error("error: latticelibrary: undefined lattice");
  }
  return 0;
}

// read lattice parameters
int Lattice::construct(const input::Parameters& parms) 
{
  int info;
  // name
  lname = parms.set_value("lattice", "NULL");
  boost::to_upper(lname);

  // sizes
  for (int dim=dim1; dim<=dim3; ++dim) {
    std::string lsize = "lsize" + std::to_string(dim+1);
    extent[dim].size = parms.set_value(lsize, 1, info);
    if (extent[dim].size<1) throw std::range_error("error: latticelibrary: invalid lattice size");
  }

  // boundary conditions
  std::string bc; 
  for (int dim=dim1; dim<=dim3; ++dim) {
    std::string lbc = "bc" + std::to_string(dim+1);
    bc = parms.set_value(lbc, "open", info);
    extent[dim].periodicity = boundary_condition(bc);
    extent[dim].bc = extent[dim].periodicity;
    if (extent[dim].bc == boundary_type::antiperiodic) extent[dim].bc = boundary_type::periodic;
    if (extent[dim].periodicity == boundary_type::antiperiodic) {
      extent[dim].bc_twist = pi();
    }
    else {
      extent[dim].bc_twist = 0.0;
    }
  }

  // spatial dimension of the lattice
  spatial_dim = 0;
  for (int dim=dim1; dim<=dim3; ++dim) {
    if (extent[dim].bc==boundary_type::periodic) {
      spatial_dim += 1;
    }
  }

  //extent[0].bc_twist = 0.0;
  //extent[1].bc_twist += two_pi()/extent[dim2].size * 2;
  //extent[2].bc_twist = 0.0;
  //for (int dim=dim1; dim<=dim3; ++dim) {
  //  std::cout << "twist = " << extent[dim].bc_twist << "\n";
  //}
  //std::cout << "\n";

  
  // number of different boundary conditions
  int bc1_twists = parms.set_value("bc1_twists", 1, info);
  int bc2_twists = parms.set_value("bc2_twists", 1, info);
  int bc3_twists = parms.set_value("bc3_twists", 1, info);
  if (bc1_twists==0 || bc2_twists==0 || bc3_twists==0) {
    throw std::range_error("Lattice::construct: 'bc_twists' value must be > 0");
  }

  if (extent[dim1].size==1) bc1_twists = 1;
  if (extent[dim2].size==1) bc2_twists = 1;
  if (extent[dim3].size==1) bc3_twists = 1;
  num_total_twists_ = bc1_twists*bc2_twists*bc3_twists;
  twist_angles_.resize(num_total_twists_,3);
  // twist angles
  /*
  double dtheta1 = two_pi()/extent[dim1].size;
  double dtheta2 = two_pi()/extent[dim2].size;
  double dtheta3 = two_pi()/extent[dim3].size;
  */
  double dtheta1 = two_pi()/bc1_twists;
  double dtheta2 = two_pi()/bc2_twists;
  double dtheta3 = two_pi()/bc3_twists;


  int n = 0;
  for (int k=0; k<bc3_twists; ++k) {
    double step3 = k*dtheta3;
    for (int j=0; j<bc2_twists; ++j) {
      double step2 = j*dtheta2;
      for (int i=0; i<bc1_twists; ++i) {
        twist_angles_(n,0) = extent[dim1].bc_twist + i*dtheta1;
        twist_angles_(n,1) = extent[dim2].bc_twist + step2;
        twist_angles_(n,2) = extent[dim3].bc_twist + step3;
        n++;
      }
    }
  }
  // check
  /*
  for (int i=0; i<num_total_twists_; ++i) {
    std::cout << "twist["<<i<<"] = "<<twist_angles_.row(i)<<"\n";
  }
  */
  //getchar();

  // empty unitcell
  Unitcell::clear();

  // impurities
  //impurity_sites_.clear();
  //impurity_bonds_.clear();

  define_lattice();
  finalize_lattice();
  construct_graph();

  return 0;
}

int Lattice::finalize_lattice(void) 
{
  // Finalize the unit cell
  Unitcell::finalize();

  // copy the user set dimensions
  for (int dim=dim1; dim<=dim3; ++dim) copy_extent[dim] = extent[dim];

  // Is it necessary to construct 'symmetrized lattice'?
  for (int dim=dim1; dim<=dim3; ++dim) {
    if (extent[dim].bc==boundary_type::open && extent[dim].size>1) {
      symmetrize_lattice();
      break;
    }
  }

  // number of unit cells & sites
  num_layer_cells_ = extent[dim1].size * extent[dim2].size;
  num_total_cells_ = num_layer_cells_ * extent[dim3].size;
  num_basis_sites_ = Unitcell::num_sites();
  num_total_sites_ = num_total_cells_ * num_basis_sites_;


  // check
  /*
  std::cout << "------Sites-------\n";
  for (int i=0; i<Unitcell::num_sites(); ++i) {
    std::cout << Unitcell::site(i) << std::endl;
  }
  std::cout << "------Bonds-------\n";
  for (int i=0; i<Unitcell::num_bonds(); ++i) {
    std::cout << Unitcell::bond(i) << std::endl;
  }*/

  // 'vector' & 'vector_id' attributes of the bonds
  /*
  std::map<int,int> vecid_map;
  int id=0;
  for (int i=0; i<Unitcell::num_bonds(); ++i) {
    Vector3i ivec = Unitcell::bond(i).tgt().bravindex()-Unitcell::bond(i).src().bravindex();
    int key = ivec[0] + ivec[1]*extent[dim1].size + ivec[2]*num_layer_cells_;
    auto it = vecid_map.find(key);
    if (it != vecid_map.end()) Unitcell::bond(i).set_vector_id(it->second);
    else {
      vecid_map.insert({key, id});
      Unitcell::bond(i).set_vector_id(id);
      id++;
    }
    Unitcell::bond(i).set_vector((ivec[0]*vector_a1()+ivec[1]*vector_a2()+ivec[2]*vector_a3()));
    //std::cout << "bond " << i << ": vector_id = " << bond(i).vector_id() << "\n";
  }
  */
  
  // vector_id by making the cell-vector itself as a key
  using triplet = std::tuple<int,int,int>;
  std::map<triplet,int> vecid_map;
  int id=0;
  for (int i=0; i<Unitcell::num_bonds(); ++i) {
    Vector3i ivec = Unitcell::bond(i).tgt().bravindex()-Unitcell::bond(i).src().bravindex();
    //std::cout << "ivec=" << ivec.transpose() << "\n";
    triplet key = std::make_tuple(ivec[0], ivec[1], ivec[2]);
    auto it = vecid_map.find(key);
    if (it != vecid_map.end()) Unitcell::bond(i).set_vector_id(it->second);
    else {
      vecid_map.insert({key, id});
      Unitcell::bond(i).set_vector_id(id);
      id++;
    }
    Unitcell::bond(i).set_vector((ivec[0]*vector_a1()+ivec[1]*vector_a2()+ivec[2]*vector_a3()));
    //std::cout << "bond " << i << ": vector_id = " << Unitcell::bond(i).vector_id() << "\n";
  }
  //getchar();

  return 0;
}

int Lattice::symmetrize_lattice(void) 
{
  // initially, the 'dim' with periodic bc has size = 1
  int actual_spatial_dim = 0; // considering BCs
  Vector3d bvec; // will be used only if 'actual_spatial_dim=1'
  for (int dim=dim1; dim<=dim3; ++dim) {
    if (extent[dim].bc==boundary_type::periodic) {
      actual_spatial_dim++;
      // temporarily set size = 1 for dim with PBC
      extent[dim].size = 1;
      switch (dim) {
        case dim1: bvec = basis_vector_a1(); break;
        case dim2: bvec = basis_vector_a2(); break;
        case dim3: bvec = basis_vector_a3(); break;
      }
    }
  }
  // if 1 dimensional lattice, rotate the lattice to align 'bvec' along x-direction
  if (actual_spatial_dim == 1) {
    // rotation matrix to do that
    Eigen::Matrix3d matrix = rotation_matrix(bvec, Vector3d(1.0,0.0,0.0));
    // rotate the unitcell
    Unitcell::rotate_by(matrix);
  }

  // number of unit cells & sites
  num_layer_cells_ = extent[dim1].size * extent[dim2].size;
  num_total_cells_ = num_layer_cells_ * extent[dim3].size;
  num_basis_sites_ = Unitcell::num_sites();
  num_total_sites_ = num_total_cells_ * num_basis_sites_;

  // Add the sites & the bonds to the symmetrized unitcell
  std::vector<Site> sites;
  std::vector<Bond> bonds;
  Unitcell translated_cell;
  Vector3i bravindex(0,0,0);
  for (int i=0; i<num_total_cells_; ++i) {
    translated_cell = get_translated_cell(bravindex);
    // collect the sites
    for (int n=0; n<translated_cell.num_sites(); ++n) 
      sites.push_back(translated_cell.site(n));
    // collect the bonds
    for (int n=0; n<translated_cell.num_bonds(); ++n) 
      bonds.push_back(translated_cell.bond(n));
    bravindex = get_next_bravindex(bravindex);
  }

  // replace the old sites & bonds
  Unitcell::clear_sites();
  int i = 0;
  for (auto& s : sites) {
    s.reset_uid(i++); 
    s.reset_bravindex(Vector3i(0,0,0));
    s.reset_cell_coord(Vector3d(0.0,0.0,0.0));
    Unitcell::add_site(s);
  }
  Unitcell::update_num_orbitals();
  Unitcell::clear_bonds();
  for (auto& b : bonds) {
    if (connect_bond2(b, sites)) {
      b.reset_bravindex(Vector3i(0,0,0));
      Unitcell::add_bond(b);
    }
  }


  //std::cout << unitcell.vector_a1() << "\n";
  //std::cout << unitcell.vector_a2() << "\n";
  //std::cout << unitcell.vector_a3() << "\n";

  // extent & basis vectors of the symmetrized lattice
  for (int dim=dim1; dim<=dim3; ++dim) {
    extent[dim] = copy_extent[dim];
    if (extent[dim].bc != boundary_type::periodic) {
      extent[dim].size = 1;
      switch (dim) {
        case dim1: Unitcell::reset_a1(Vector3d(0.0,0.0,0.0)); break;
        case dim2: Unitcell::reset_a2(Vector3d(0.0,0.0,0.0)); break;
        case dim3: Unitcell::reset_a3(Vector3d(0.0,0.0,0.0)); break;
      }
    }
  }
  /*
  std::cout << unitcell.vector_a1() << "\n";
  std::cout << unitcell.vector_a2() << "\n";
  std::cout << unitcell.vector_a3() << "\n";
  getchar();
  */

  return 0;
}

int Lattice::construct_graph(void) 
{
  // all the sites and the bonds
  Unitcell translated_cell;
  Vector3i bravindex(0,0,0);
  std::vector<Site> allsites;
  std::vector<Bond> allbonds;
  for (int i=0; i<num_unitcells(); ++i) {
    translated_cell = get_translated_cell(bravindex);
    // collect the sites
    for (int n=0; n<translated_cell.num_sites(); ++n) {
      allsites.push_back(translated_cell.site(n));
    }
    // collect the bonds
    for (int n=0; n<translated_cell.num_bonds(); ++n) {
      allbonds.push_back(translated_cell.bond(n));
    }
    bravindex = get_next_bravindex(bravindex);
  }

  // save the sites
  sites_.clear();
  //sitetype_set_.clear(); // not used
  for (const auto& s : allsites) {
    sites_.push_back(s);
    sites_.back().clear_bonds(); // to be added later
    //sitetype_set_.insert(s.type());
    //std::cout << "uc = " << s.bravindex().transpose() << "\n";
    //std::cout << "id = " << s.id() << "\n";
    //std::cout << "uid = " << s.uid() << "\n";
    //std::cout << "aid = " << s.atomid() << "\n";
    //std::cout << "tp = " << s.type() << "\n\n";
    //getchar();
  }
  //getchar();

  /*
  int i=0;
  for (const auto& a : atoms_) {
    std::cout << "atom = "<<i<<": ";
    for (const auto& id: a) {
      std::cout << "  "<<id;
    }
    std::cout << "\n";
    i++;
    getchar();
  }
  */


  // Save the bonds.
  bonds_.clear();
  //bondtype_set_.clear(); // not used
  int id = 0;
  for (auto& b : allbonds) {
    // Connect the bond, discard if can't be connected 
    // (for bonds across open boundaries) 
    if (!connect_bond2(b, sites_)) continue;
    b.reset_id(id++);
    // Set boundary condition phase.
    // Also redefine meaning of 'sign' ('-'ve mean boundary bond)
    int sign = 1;
    std::complex<double> phase = 1.0;
    if (b.bc_state()[0]==-1) {
      sign = -1;
      phase *= std::exp(ii()*bc1_twist());
    }
    if (b.bc_state()[1]==-1) {
      sign = -1;
      phase *= std::exp(ii()*bc2_twist());
    } 
    if (b.bc_state()[2]==-1) {
      sign = -1;
      phase *= std::exp(ii()*bc3_twist());
    } 
    if (std::abs(phase.imag())<1.0E-15) phase.imag(0.0); 
    b.reset_sign(sign);
    b.reset_phase(phase);

    // site connections
    sites_[b.src_id()].add_out_bond(b.id());
    sites_[b.tgt_id()].add_in_bond(b.id());

    // save the bond
    //bondtype_set_.insert(b.type());
    bonds_.push_back(b);

    //std::cout << "b.vector = " << b.vector().transpose() << "\n";
    //getchar();
  }
  num_bonds_ = bonds_.size();

  // check consisency
  if (sites_.size() != num_total_sites_) {
    throw std::logic_error("Lattice::construct_graph: site count mismatch");
  }
  for (int i=0; i<sites_.size(); ++i) {
    if (i != sites_[i].id()) {
      throw std::logic_error("Lattice::construct_graph: site id mismatch");
    }
  }
  for (int i=0; i<bonds_.size(); ++i) {
    if (i != bonds_[i].id()) {
      throw std::logic_error("Lattice::construct_graph: bond id mismatch");
    }
  }

  return 0;
}

int Lattice::reset_boundary_twist(const int& twist_id)
{
  // twist angle 
  for (int dim=dim1; dim<=dim3; ++dim) {
    extent[dim].bc_twist = twist_angles_(twist_id,dim);
    //std::cout << "twist = " << extent[dim].bc_twist << "\n";
  }
  //std::cout << "twist = " << twist_angles_.row(twist_id) << "\n";
  //std::cout << "\n";

  // reset bond phase values
  for (Bond& b : bonds_) {
    if (b.sign() == -1) {
      std::complex<double> phase = 1.0;
      if (b.bc_state()[0]==-1) {
        phase *= std::exp(ii()*bc1_twist());
      }
      if (b.bc_state()[1]==-1) {
        phase *= std::exp(ii()*bc2_twist());
      } 
      if (b.bc_state()[2]==-1) {
        phase *= std::exp(ii()*bc3_twist());
      } 
      if (std::abs(phase.imag())<1.0E-15) phase.imag(0.0); 
      b.reset_phase(phase);
      //std::cout << "phase["<<b.id()<<"] = " << b.phase() << "\n";
    }
  }

  return 0;
}




} // end namespace lattice
