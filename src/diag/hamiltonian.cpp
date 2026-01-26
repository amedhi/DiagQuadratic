/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2025-12-06 12:15:35
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-22 11:32:04
*----------------------------------------------------------------------------*/
#include "./hamiltonian.h"
#include <boost/algorithm/string.hpp>
#include "../expression/complex_expression.h"

namespace diag {


/*--------------------------Hamiltonian------------------------*/
Hamiltonian::Hamiltonian(const model::Model& model, const lattice::Lattice& lattice)
: model::Model(model)
{
  init_ham(lattice);
}

int Hamiltonian::init(const lattice::Lattice& lattice)
{
  return Model::init(lattice);
}

int Hamiltonian::finalize(const lattice::Lattice& lattice)
{
  Model::finalize(lattice);
  init_ham(lattice);
  return 0;
}

int Hamiltonian::init_ham(const lattice::Lattice& lattice)
{
  // set Hamiltonian parameters
  num_basis_sites_ = lattice.num_basis_sites();
  dim_ = lattice.num_basis_orbitals();
  // indexing of the '[site, orbital] pairs'
  site_orb_idx_.clear();
  site_orb_idx_.resize(num_basis_sites_);
  int idx = 0;
  for (int i=0; i<num_basis_sites_; ++i) {
    for (int m=0; m<lattice.site(i).num_orbitals(); ++m) {
      site_orb_idx_[i].push_back(idx++);
    }
  }
  assert(dim_==idx);

  // labels and 'site' coordinates of the basis states
  basis_state_labels_.resize(dim_);
  basis_state_coords_.resize(dim_);
  idx = 0;
  for (int i=0; i<num_basis_sites_; ++i) {
    std::string str = std::to_string(i);
    for (int m=0; m<lattice.site(i).num_orbitals(); ++m) {
      basis_state_labels_[idx] = str+"-"+std::to_string(m); 
      basis_state_coords_[idx] = lattice.site(i).coord();
      idx++;
    }
  }

  // build unitcell terms
  build_unitcell_terms(lattice);

  // Hamiltonian matrices
  quadratic_block_up_.resize(dim_,dim_);
  quadratic_block_dn_.resize(dim_,dim_);
  pairing_block_.resize(dim_,dim_);
  work.resize(dim_,dim_);
  quadratic_block_up_.setZero();
  quadratic_block_dn_.setZero();
  pairing_block_.setZero();
  work.setZero();

  return 0;
}

int Hamiltonian::update(const input::Parameters& inputs)
{
  Model::update_parameters(inputs);
  update_terms();
  return 0;
}

void Hamiltonian::update_terms(void)
{
  update_unitcell_terms();
}

int Hamiltonian::update_site_parameter(const std::string& pname, const double& pvalue)
{
  Model::update_parameter(pname, pvalue);
  for (int i=0; i<usite_terms_.size(); ++i) 
    usite_terms_[i].eval_coupling_constant(Model::parameters(),Model::constants());
  return 0;
}

int Hamiltonian::build_unitcell_terms(const lattice::Lattice& lattice)
{
  // take only quadratic & pairing terms
  int num_siteterms = 0;
  int num_bondterms = 0;
  for (auto sterm=siteterms_begin(); sterm!=siteterms_end(); ++sterm) {
    if (sterm->qn_operator().is_quadratic() || sterm->qn_operator().is_pairing())
      num_siteterms++;
  }
  for (auto bterm=bondterms_begin(); bterm!=bondterms_end(); ++bterm) {
    if (bterm->qn_operator().is_quadratic() || bterm->qn_operator().is_pairing()) {
      //std::cout << bterm->qn_operator().name() << "\n"; getchar();
      num_bondterms++;
    }
  }
  usite_terms_.resize(num_siteterms);
  ubond_terms_.resize(num_bondterms);
  int i = 0;
  for (auto sterm=siteterms_begin(); sterm!=siteterms_end(); ++sterm) {
    if (sterm->qn_operator().is_quadratic() || sterm->qn_operator().is_pairing()) {
      usite_terms_[i].build_siteterm(*sterm, lattice, site_orb_idx_);
      i++;
    }
  }
  i = 0;
  for (auto bterm=bondterms_begin(); bterm!=bondterms_end(); ++bterm) {
    if (bterm->qn_operator().is_quadratic() || bterm->qn_operator().is_pairing()) {
      ubond_terms_[i].build_bondterm(*bterm, lattice, site_orb_idx_);
      i++;
    }
  }
  return 0;
}

int Hamiltonian::construct_upspin_block(const Vector3d& kvec) const
{
  work.setZero(); 
  for (const auto& term : ubond_terms_) {
    if (term.qn_operator().is_quadratic() && term.qn_operator().spin_up()) {
      //std::cout << term.num_out_bonds() << "\n";
      for (int i=0; i<term.num_out_bonds(); ++i) {
        Vector3d delta = term.bond_vector(i);
        work += term.coeff_matrix(i) * std::exp(ii()*kvec.dot(delta));
        //std::cout << "delta = " << delta.transpose() << "   " << 
        //kvec.dot(delta) << "\n"; getchar();
      }
    }
  }
  // add hermitian conjugate part
  quadratic_block_up_ = work + work.adjoint();
  // add site terms 
  for (const auto& term : usite_terms_) {
    if (term.qn_operator().spin_up()) {
      quadratic_block_up_ += term.coeff_matrix();
    }
  }
  return 0;
}

int Hamiltonian::construct_dnspin_block(const Vector3d& kvec) const
{
  work.setZero(); 
  for (const auto& term : ubond_terms_) {
    if (term.qn_operator().is_quadratic() && term.qn_operator().spin_dn()) {
      for (int i=0; i<term.num_out_bonds(); ++i) {
        Vector3d delta = term.bond_vector(i);
        work += term.coeff_matrix(i) * std::exp(ii()*kvec.dot(delta));
      }
    }
  }
  // add hermitian conjugate part
  quadratic_block_dn_ = work + work.adjoint();
  // site terms 
  for (const auto& term : usite_terms_) {
    if (term.qn_operator().spin_dn()) {
      quadratic_block_dn_ += term.coeff_matrix();
    }
  }
  return 0;
}

int Hamiltonian::construct_kblock(const Vector3d& kvec) const
{
  construct_upspin_block(kvec);
  construct_dnspin_block(kvec);
  return 0;
}

int Hamiltonian::construct_pairing_part(const Vector3d& kvec) const
{
  pairing_block_.setZero();
  for (const auto& term : ubond_terms_) {
    if (term.qn_operator().is_pairing()) {
      for (int i=0; i<term.num_out_bonds(); ++i) {
        Vector3d delta = term.bond_vector(i);
        pairing_block_ += term.coeff_matrix(i) * std::exp(ii()*kvec.dot(delta));
      }
    }
  }
  return 0;
}

/*
int Hamiltonian::construct_kblock(const Vector3d& kvec) const
{
  work.setZero(); 
  pairing_block_.setZero();
  //std::cout << "---------------------------------\n";
  //std::cout << "kvec = " << kvec.transpose() << "\n";
  //std::cout << "---------------------------------\n";
  // bond terms
  for (const auto& term : ubond_terms_) {
    if (term.qn_operator().is_quadratic() && term.qn_operator().spin_up()) {
      //std::cout << term.num_out_bonds() << "\n";
      for (int i=0; i<term.num_out_bonds(); ++i) {
        Vector3d delta = term.bond_vector(i);
        work += term.coeff_matrix(i) * std::exp(ii()*kvec.dot(delta));
        //std::cout << "delta = " << delta.transpose() << "   " << 
        //kvec.dot(delta) << "\n"; getchar();
      }
    }
    if (term.qn_operator().is_pairing()) {
      for (int i=0; i<term.num_out_bonds(); ++i) {
        Vector3d delta = term.bond_vector(i);
        pairing_block_ += term.coeff_matrix(i) * std::exp(ii()*kvec.dot(delta));
      }
    }
  }
  // add hermitian conjugate part
  quadratic_block_up_ = work + work.adjoint();
  // add site terms 
  for (const auto& term : usite_terms_) {
    if (term.qn_operator().spin_up()) {
      quadratic_block_up_ += term.coeff_matrix();
    }
    if (term.qn_operator().is_pairing()) {
      pairing_block_ += term.coeff_matrix();
    }
  }
  //std::cout << "pairing_block = " << pairing_block_ << "\n"; getchar();
  //quadratic_block_up_ += work1.adjoint();
  //pairing_block_ = work2;
  //pairing_block_ += work2.adjoint();
  // site terms
  //std::cout << "ek = " << quadratic_block_(0,0) << "\n";

  // quadratic 'dn'-spin block
  work.setZero(); 
  for (const auto& term : ubond_terms_) {
    if (term.qn_operator().is_quadratic() && term.qn_operator().spin_dn()) {
      for (int i=0; i<term.num_out_bonds(); ++i) {
        Vector3d delta = term.bond_vector(i);
        work += term.coeff_matrix(i) * std::exp(ii()*kvec.dot(delta));
      }
    }
  }
  // add hermitian conjugate part
  quadratic_block_dn_ = work + work.adjoint();
  // site terms 
  for (const auto& term : usite_terms_) {
    if (term.qn_operator().spin_dn()) {
      quadratic_block_dn_ += term.coeff_matrix();
    }
  }
  return 0;
}
*/


int Hamiltonian::update_unitcell_terms(void)
{
  for (int i=0; i<ubond_terms_.size(); ++i) 
    ubond_terms_[i].eval_coupling_constant(Model::parameters(),Model::constants());
  for (int i=0; i<usite_terms_.size(); ++i) 
    usite_terms_[i].eval_coupling_constant(Model::parameters(),Model::constants());
  return 0;
}


/* Write a bond term like,
 H = \sum_{Ia,Jb}c^{\dag}_{Ia} t_{Ia,Jb} c_{Jb}
 for lattices with multiple sites per unit cell as
 H = \sum_{I,delta} Psi^{\dag}_{I} M^{\delta} Psi_{I+delta}
 Assumption: 'site's in the Graph are numbered contigously. For
 sites in the unitcell, the 'sl number' is same as 'uid'.
*/

/*--------------------------UnitcellTerm------------------------*/
void UnitcellTerm::build_bondterm(const model::ModelTerm& bterm, 
  const lattice::Lattice& lattice, const basis_idx_map& siteorb_idx)
{
  num_basis_sites_ = lattice.num_basis_sites();
  dim_ = lattice.num_basis_orbitals();
  //std::cout << num_basis_sites_ << "  " << dim_ << "\n";

  // get number of unique 'cell bond vectors'
  num_out_bonds_ = 0;
  for (int i=0; i<num_basis_sites_; ++i) {
    std::vector<int> out_bond_ids = lattice.site(i).outbond_ids();
    for (const auto& id : out_bond_ids) {
      if (id > num_out_bonds_) num_out_bonds_ = id;
    }
  }
  num_out_bonds_++;
  //std::cout << "num_out_bonds_ = " << num_out_bonds_ << "\n";

  bond_vectors_.resize(num_out_bonds_);
  coeff_matrices_.resize(num_out_bonds_);
  for (auto& M : coeff_matrices_) {
    M.resize(dim_, dim_);
    M.setZero();
  }
  expr_matrices_.clear();
  expr_matrices_.resize(num_out_bonds_);
  for (auto& M : expr_matrices_) {
    M.resize(dim_,dim_);
  }
  // initialize expression matrices
  for (int id=0; id<num_out_bonds_; ++id) {
    for (int i=0; i<dim_; ++i) {
      for (int j=0; j<dim_; ++j) {
        expr_matrices_[id](i,j) = "0";
      }
    }
  }

  // operator
  op_ = bterm.qn_operator();
  // build the matrices (for each 'bond vector')
  for (int i=0; i<num_basis_sites_; ++i) {
    std::vector<int> out_bond_ids = lattice.site(i).outbond_ids();
    for (const auto& id : out_bond_ids) {
      int btype = lattice.bond(id).type();
      int j = lattice.bond(id).tgt().uid();
      strMatrix expr_mat = bterm.coupling_expr(btype);
      ComplexMatrix coeff_mat = bterm.coupling(btype);
      int rows = lattice.site(i).num_orbitals(); 
      int cols = lattice.site(j).num_orbitals(); 
      if ((expr_mat.rows()!=rows) && (expr_mat.cols()!=cols)) {
        throw std::runtime_error("Dimension mismatch in the 'bondterm'");
      }
      for (int m=0; m<rows; ++m) {
        for (int n=0; n<cols; ++n) {
          int p = siteorb_idx[i][m];  
          int q = siteorb_idx[j][n];  
          std::string cc_expr(expr_mat(m,n));
          if (cc_expr.size()>0) {
            expr_matrices_[id](p,q) += "+("+cc_expr+")";
          }
          // values
          coeff_matrices_[id](p,q) += coeff_mat(m,n);
        }
      }
      //std::cout << id << " " << coeff_matrices_[id](i,j) << "\n";
      bond_vectors_[id] = lattice.bond(id).vector();
      //std::cout << i << " " << j << "\n"; 
      //std::cout<<bond_vectors_[id].transpose()<<"\n"; getchar();
      //std::cout << "--------HERE--------\n";
    }
  }
  /*std::cout <<"------"<<op_.name()<<"-------\n"; getchar();
  for (int id=0; id<num_out_bonds_; ++id) {
    std::cout << "delta = " << id << "\n";
    for (int i=0; i<dim_; ++i) {
      for (int j=0; j<dim_; ++j) {
        std::cout<<i<<", "<<j<<"= "<<expr_matrices_[id][i][j]<<"\n";
      }
    }
    std::cout<<"\n";
  }
  getchar();
  */
}

void UnitcellTerm::build_siteterm(const model::ModelTerm& sterm, 
  const lattice::Lattice& lattice, const basis_idx_map& siteorb_idx)
{
  num_basis_sites_ = lattice.num_basis_sites();
  dim_ = lattice.num_basis_orbitals();

  num_out_bonds_ = 1; // dummy, no real contrib
  bond_vectors_.resize(1);
  coeff_matrices_.resize(1);
  coeff_matrices_[0].resize(dim_,dim_);
  coeff_matrices_[0].setZero();
  expr_matrices_.resize(1);
  expr_matrices_[0].resize(dim_,dim_);
  // initialize expression matrix
  for (int i=0; i<dim_; ++i) {
    for (int j=0; j<dim_; ++j) {
      expr_matrices_[0](i,j) = "0";
    }
  }

  // operator
  op_ = sterm.qn_operator();
  // build the matrix 
  for (int i=0; i<num_basis_sites_; ++i) {
    int stype = lattice.site(i).type();
    //coeff_matrices_[0](i,i) = model_term.coupling(stype);
    // expression
    strMatrix expr_mat = sterm.coupling_expr(stype);
    ComplexMatrix coeff_mat = sterm.coupling(stype);
    int rows = lattice.site(i).num_orbitals(); 
    int cols = rows;
    if (expr_mat.cols()!=cols) {
      throw std::runtime_error("Dimension mismatch in the 'siteterm'");
    }
    // diagonal site term
    if (expr_mat.rows()==1) {
      for (int m=0; m<cols; ++m) {
        int p = siteorb_idx[i][m];  
        std::string cc_expr(expr_mat(0,m));
        boost::trim(cc_expr);
        if (cc_expr.size()>0) {
          expr_matrices_[0](p,p) += "+("+cc_expr+")";
        }
        coeff_matrices_[0](p,p) = coeff_mat(0,m);
      }
    }
    else {
      // off-diagonal site term
      if (expr_mat.rows()!=rows) {
        throw std::runtime_error("Dimension mismatch in the 'siteterm'");
      }
      for (int m=0; m<rows; ++m) {
        for (int n=0; n<cols; ++n) {
          int p = siteorb_idx[i][m];  
          int q = siteorb_idx[i][n];  
          std::string cc_expr(expr_mat(m,n));
          if (cc_expr.size()>0) {
            expr_matrices_[0](p,q) += "+("+cc_expr+")";
          }
          // values
          coeff_matrices_[0](p,q) += coeff_mat(m,n);
        }
      }
    }
  }
  bond_vectors_[0] = Vector3d(0.0,0.0,0.0);
}


void UnitcellTerm::eval_coupling_constant(const model::ModelParams& pvals, const model::ModelParams& cvals)
{
  //std::cout << "------eval_coupling_constant---------\n"; getchar();
  expr::ComplexExpr expr;
  for (const auto& p : pvals) expr.add_var(p.first, p.second);
  for (const auto& c : cvals) expr.add_var(c.first, c.second);
  try { 
    for (int n=0; n<num_out_bonds_; ++n) {
      for (int i=0; i<dim_; ++i) {
        for (int j=0; j<dim_; ++j) {
          std::string cc_expr(expr_matrices_[n](i,j));
          if (cc_expr.size()>0) {
            expr.set_expr(cc_expr);
            coeff_matrices_[n](i,j) = expr.evaluate(); 
            //std::cout << "cc = " << cc_expr << "\n"; getchar();
          }
          else coeff_matrices_[n](i,j) = 0.0;
        }
      }
    }
  }
  catch (std::exception& e) 
  { 
    std::string msg = "UnitcellTerm::evaluate_coupling_constant:\n" + std::string(e.what());
    throw std::runtime_error(msg);
  }
}

/*void Hamiltonian::check_xml(void)
{
  std::cout << "Checking XML parser\n";
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("model.xml", pugi::parse_trim_pcdata);
  //std::cout << "Load result: " << result.description() << ", mesh name: " << "\n"; 
  pugi::xml_node model = doc.child("model");
  //std::cout << model.child("parameter").attribute("default").value() << std::endl;
  for (pugi::xml_node p = model.child("parameter"); p; p = p.next_sibling())
  {
    std::cout << "Parameter: ";
    for (pugi::xml_attribute attr = p.first_attribute(); attr; attr = attr.next_attribute())
    {
      std::cout << attr.name() << " = " << attr.as_double() << "\n";
    }
    std::cout << std::endl;
  }
  std::cout << "---------------------------------\n\n";
}*/


} // end namespace diag











