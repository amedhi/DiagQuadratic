/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-06-17 13:17:42
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

  blochbasis_.gen_mesh_neighbors(graph_.lattice());
  compute_chern_number();


  //std::cout << "b1 = " << blochbasis_.vector_b1().transpose() << "\n";
  //std::cout << "b2 = " << blochbasis_.vector_b2().transpose() << "\n";
  //std::cout << "b3 = " << blochbasis_.vector_b3().transpose() << "\n";
  Vector3d K_point = 0.5*(blochbasis_.vector_b1()-blochbasis_.vector_b2());
  Vector3d M_point = 0.5*blochbasis_.vector_b1();
  Vector3d G_point = Vector3d(0,0,0);
  
  int N = 100;
  Vector3d step = (K_point-G_point)/N;
  //std::cout << "K = " << step.transpose() << "\n";
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
  for (unsigned k=0; k<num_kpoints_; ++k) {
  //for (unsigned k=0; k<symm_line_.size(); ++k) {
    //std::cout << k << " of " << symm_line_.size() << "\n";
    Vector3d kvec = blochbasis_.kvector(k);
    //Vector3d kvec = symm_line_[k];
    mf_model_.construct_kspace_block(kvec);
    //std::cout << mf_model_.quadratic_spinup_block() << "\n"; getchar();
    es_k_up.compute(mf_model_.quadratic_spinup_block(), Eigen::EigenvaluesOnly);
    //ek.insert(ek.end(),es_k_up.eigenvalues().data(),
    //    es_k_up.eigenvalues().data()+kblock_dim_);
    if (k%graph_.lattice().size1()==0) of << "\n";
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

int Diag::compute_chern_number(void) 
{
  //using xy_pair = std::pair<std::complex<double>,std::complex<double> >;
  //std::vector<xy_pair> BerryConnection(num_kpoints_);
  //ComplexVector psi_k1(kblock_dim_);
  //ComplexVector psi_k2(kblock_dim_);
  //int band_idx = 0;

  std::vector<ComplexVector> BerryConnection_dir1(num_kpoints_);
  std::vector<ComplexVector> BerryConnection_dir2(num_kpoints_);
  for (int k=0; k<num_kpoints_; ++k) {
    BerryConnection_dir1[k].resize(kblock_dim_);
    BerryConnection_dir2[k].resize(kblock_dim_);
  }
  ComplexMatrix psi_k(kblock_dim_,kblock_dim_);
  for (unsigned k=0; k<num_kpoints_; ++k) {
    Vector3d kv1 = blochbasis_.kvector(k);
    mf_model_.construct_kspace_block(kv1);
    es_k_up.compute(mf_model_.quadratic_spinup_block());
    //psi_k1 = es_k_up.eigenvectors().col(band_idx);
    psi_k = es_k_up.eigenvectors();

    int k_nn = blochbasis_.mesh_nn_xp(k);
    Vector3d kv2 = blochbasis_.kvector(k_nn);
    mf_model_.construct_kspace_block(kv2);
    es_k_up.compute(mf_model_.quadratic_spinup_block());
    /*
    psi_k2 = es_k_up.eigenvectors().col(band_idx);
    auto x = psi_k1.dot(psi_k2);
    BerryConnection[k].first  = x/std::abs(x);
    */
    // Berry Connection at 'k' along dir1 for all bands
    for (int n=0; n<kblock_dim_; ++n) {
      auto x = psi_k.col(n).dot(es_k_up.eigenvectors().col(n));
      BerryConnection_dir1[k][n] = x/std::abs(x);
    }


    k_nn = blochbasis_.mesh_nn_yp(k);
    kv2 = blochbasis_.kvector(k_nn);
    mf_model_.construct_kspace_block(kv2);
    es_k_up.compute(mf_model_.quadratic_spinup_block());
    /*
    psi_k2 = es_k_up.eigenvectors().col(band_idx);
    x = psi_k1.dot(psi_k2);
    BerryConnection[k].second = x/std::abs(x); // along dir-2
    */
    // Berry Connection at 'k' along dir2 for all bands
    for (int n=0; n<kblock_dim_; ++n) {
      auto x = psi_k.col(n).dot(es_k_up.eigenvectors().col(n));
      BerryConnection_dir2[k][n] = x/std::abs(x);
    }
  }

  std::complex<double> U_01, U_12, U_23, U_30;
  //std::complex<double> csum(0.0);
  ComplexVector BerrySum(kblock_dim_);
  BerrySum.setZero();
  for (unsigned k=0; k<num_kpoints_; ++k) {
    /*
    U_01 = BerryConnection[k].first;
    U_30 = std::conj(BerryConnection[k].second);
    int k_nn = blochbasis_.mesh_nn_xp(k);
    U_12 = BerryConnection[k_nn].second;
    k_nn = blochbasis_.mesh_nn_yp(k);
    U_23 = std::conj(BerryConnection[k_nn].first);
    csum += std::log(U_01 * U_12 * U_23 * U_30);
    */
    int k_nn1 = blochbasis_.mesh_nn_xp(k);
    int k_nn2 = blochbasis_.mesh_nn_yp(k);
    for (int n=0; n<kblock_dim_; ++n) {
      U_01 = BerryConnection_dir1[k][n];
      U_30 = std::conj(BerryConnection_dir2[k][n]);
      U_12 = BerryConnection_dir2[k_nn1][n];
      U_23 = std::conj(BerryConnection_dir1[k_nn2][n]);
      BerrySum[n] += std::log(U_01 * U_12 * U_23 * U_30);
    }
  }
  //std::cout << "csum = " << csum/two_pi() <<"\n";
  //int chern_number = std::nearbyint(std::imag(csum/two_pi()));
  std::vector<int> ChernNumber(kblock_dim_);
  for (int n=0; n<kblock_dim_; ++n) {
    ChernNumber[n] = std::nearbyint(std::imag(BerrySum[n]/two_pi()));
    std::cout << "Chern number of band-"<<n<<" = "<<ChernNumber[n]<<"\n";
  }

  return 0;
}

void Diag::print_copyright(std::ostream& os)
{
  std::cout << "Diag\n";
}


} // end namespace diag
