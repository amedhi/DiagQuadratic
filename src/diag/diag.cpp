/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-06-11 17:23:02
*----------------------------------------------------------------------------*/
#include <iomanip>
#include <fstream>
#include "diag.h"

namespace diag {

Diag::Diag(const input::Parameters& inputs) 
  : graph_(inputs) 
  , model_(inputs, graph_.lattice())
  , mf_model_(model_,graph_)
  , blochbasis_(graph_)
{
  num_kpoints_ = blochbasis_.num_kpoints();
  kblock_dim_ = blochbasis_.subspace_dimension();
  mf_model_.init(graph_.lattice());
  std::string name; 
  double defval;
  using namespace model;
  model::CouplingConstant cc;

  if (graph_.lattice().id()==lattice::lattice_id::KAGOME) {
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_parameter(name="t2", defval=1.0, inputs);
    mf_model_.add_parameter(name="lambda", defval=1.0, inputs);
    mf_model_.add_parameter(name="lambda2", defval=1.0, inputs);
    // upspin hop
    cc = CouplingConstant({0,"-t+i*lambda"},{1,"-t-i*lambda"},
      {2,"-t2+i*lambda2"},{3,"-t2-i*lambda2"});
    mf_model_.add_bondterm(name="hopping", cc, op::upspin_hop());
    // dnspin hop
    cc = CouplingConstant({0,"-t-i*lambda"},{1,"-t+i*lambda"},
      {2,"-t2-i*lambda2"},{3,"-t2+i*lambda2"});
    mf_model_.add_bondterm(name="hopping", cc, op::dnspin_hop());
    mf_model_.finalize(graph_);
  }
  else if (graph_.lattice().id()==lattice::lattice_id::PYROCHLORE) {
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_parameter(name="t2", defval=1.0, inputs);
    mf_model_.add_parameter(name="lambda", defval=1.0, inputs);
    mf_model_.add_parameter(name="lambda2", defval=1.0, inputs);
    mf_model_.add_parameter(name="th", defval=1.0, inputs);
    // upspin hop
    cc = CouplingConstant({0,"-t+i*lambda"},{1,"-t-i*lambda"},
      {2,"-t2+i*lambda2"},{3,"-t2-i*lambda2"},{4,"-th"});
    mf_model_.add_bondterm(name="hopping", cc, op::upspin_hop());
    // dnspin hop
    cc = CouplingConstant({0,"-t-i*lambda"},{1,"-t+i*lambda"},
      {2,"-t2-i*lambda2"},{3,"-t2+i*lambda2"},{4,"-th"});
    mf_model_.add_bondterm(name="hopping", cc, op::dnspin_hop());
    mf_model_.finalize(graph_);
  }

  std::cout << "b1 = " << blochbasis_.vector_b1().transpose() << "\n";
  std::cout << "b2 = " << blochbasis_.vector_b2().transpose() << "\n";
  std::cout << "b3 = " << blochbasis_.vector_b3().transpose() << "\n";
  Vector3d K_point = 0.5*(blochbasis_.vector_b1()-blochbasis_.vector_b2());
  Vector3d M_point = 0.5*blochbasis_.vector_b1();
  Vector3d G_point = Vector3d(0,0,0);
  int N = 100;
  Vector3d step = (K_point-G_point)/N;
  std::cout << "K = " << step.transpose() << "\n";
  for (int i=0; i<N; ++i) symm_line_.push_back(G_point+i*step);
  step = (M_point-K_point)/N;
  for (int i=0; i<N; ++i) symm_line_.push_back(K_point+i*step);
  step = (G_point-M_point)/N;
  for (int i=0; i<N; ++i) symm_line_.push_back(M_point+i*step);
  symm_line_.push_back(G_point);
}

int Diag::run(const input::Parameters& inputs) 
{
  //std::ios state(NULL);
  //state.copyfmt(std::cout);
  //std::cout << std::scientific << std::uppercase << std::setprecision(6) << std::right;
  std::ofstream of("energy.txt");
  of << std::scientific << std::uppercase << std::setprecision(6) << std::right;
  std::vector<double> ek;
  //for (unsigned k=0; k<num_kpoints_; ++k) {
  for (unsigned k=0; k<symm_line_.size(); ++k) {
    std::cout << k << " of " << symm_line_.size() << "\n";
    //Vector3d kvec = blochbasis_.kvector(k);
    Vector3d kvec = symm_line_[k];
    mf_model_.construct_kspace_block(kvec);
    es_k_up.compute(mf_model_.quadratic_spinup_block(), Eigen::EigenvaluesOnly);
    //ek.insert(ek.end(),es_k_up.eigenvalues().data(),
    //    es_k_up.eigenvalues().data()+kblock_dim_);
      //if (k%graph_.lattice().size1()==0) of << "\n";
    of << std::setw(6) << k; 
    of << std::setw(14) << kvec(0) << std::setw(14) << kvec(1); 
    of << std::setw(14) << es_k_up.eigenvalues().transpose() << "\n";
  }
  of << "\n\n"; 
  of.close();
  //std::sort(ek.begin(),ek.end());
  //std::cout.copyfmt(state);
  return 0;
}


void Diag::print_copyright(std::ostream& os)
{
  std::cout << "Diag\n";
}


} // end namespace diag
