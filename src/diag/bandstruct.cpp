/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2025-12-09 17:26:56
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-18 23:10:55
*----------------------------------------------------------------------------*/
#include <boost/algorithm/string.hpp>
#include "./bandstruct.h"

namespace diag {

void BandStruct::setup(const input::Parameters& inputs, const lattice::Lattice& lattice, 
  const kSpace& kspace, const Hamiltonian& ham)
{
  MC_Observable::switch_on();
  if (setup_done_) return;

  // setup k-points along symmetry path
  num_bands_ = ham.dimension();
  symm_kpoints_.clear();
  symm_pidx_.clear();
  symm_pname_.clear();

  // spin sector
  int info;
  std::string spin_sector = inputs.set_value("spin_sector", "UP", info);
  boost::to_upper(spin_sector);
  if (spin_sector=="UP") {
    spin_sector_ = spin::UP;
  }
  else if (spin_sector=="DN") {
    spin_sector_ = spin::DN;
  }
  else if (spin_sector=="BOTH") {
    spin_sector_ = spin::BOTH;
  }
  else {
    throw std::range_error("error: BandStruct: invalid value for input 'spin_sector'");
  }


  // kpoints 
  std::string kpath = inputs.set_value("kpath", "FULL");
  boost::to_upper(kpath);
  if (kpath=="SYMM") {
    kpath_ = kpath_type::SYMM;
  }
  else if (kpath=="FULL") {
    kpath_ = kpath_type::FULL;
  }
  else {
    throw std::range_error("error: BandStruct: invalid value for input 'kpath'");
  }

  if (kpath_==kpath_type::SYMM) {
    int N = 100;
    int idx = 0;
    if (lattice.brav()==lattice::brav_id::CHAIN) {
      kvector Gamma = kspace.FBZ().Gamma().vec();
      kvector X = kspace.FBZ().X().vec();
      //---------------------------------
      symm_pidx_.push_back(idx);
      symm_pname_.push_back("Gamma");
      kvector step = (X-Gamma)/N;
      for (int i=0; i<N; ++i) symm_kpoints_.push_back(Gamma+i*step);
      //---------------------------------
      idx += N;
      symm_pidx_.push_back(idx);
      symm_pname_.push_back("X");
      step = (Gamma-X)/N;
      for (int i=0; i<N; ++i) symm_kpoints_.push_back(X+i*step);
      //---------------------------------
      idx += N;
      symm_pidx_.push_back(idx);
      symm_pname_.push_back("Gamma");
      symm_kpoints_.push_back(Gamma);
    }

    else if (lattice.brav()==lattice::brav_id::SQUARE) {
      kvector Gamma = kspace.FBZ().Gamma().vec();
      kvector X = kspace.FBZ().X().vec();
      kvector M = kspace.FBZ().M().vec();
      //---------------------------------
      symm_pidx_.push_back(idx);
      symm_pname_.push_back("Gamma");
      kvector step = (X-Gamma)/N;
      for (int i=0; i<N; ++i) symm_kpoints_.push_back(Gamma+i*step);
      //---------------------------------
      idx += N;
      symm_pidx_.push_back(idx);
      symm_pname_.push_back("X");
      step = (M-X)/N;
      for (int i=0; i<N; ++i) symm_kpoints_.push_back(X+i*step);
      //---------------------------------
      idx += N;
      symm_pidx_.push_back(idx);
      symm_pname_.push_back("M");
      step = (Gamma-M)/N;
      for (int i=0; i<N; ++i) symm_kpoints_.push_back(M+i*step);
      //---------------------------------
      idx += N;
      symm_pidx_.push_back(idx);
      symm_pname_.push_back("Gamma");
      symm_kpoints_.push_back(Gamma);
    }

    else if (lattice.brav()==lattice::brav_id::HEXAGONAL) {
      kvector Gamma = kspace.FBZ().Gamma().vec();
      kvector M = kspace.FBZ().M().vec();
      kvector K = kspace.FBZ().K().vec();
      //---------------------------------
      symm_pidx_.push_back(idx);
      symm_pname_.push_back("Gamma");
      kvector step = (M-Gamma)/N;
      for (int i=0; i<N; ++i) symm_kpoints_.push_back(Gamma+i*step);
      //---------------------------------
      idx += N;
      symm_pidx_.push_back(idx);
      symm_pname_.push_back("M");
      step = (K-M)/N;
      for (int i=0; i<N; ++i) symm_kpoints_.push_back(M+i*step);
      //---------------------------------
      idx += N;
      symm_pidx_.push_back(idx);
      symm_pname_.push_back("K");
      step = (Gamma-K)/N;
      for (int i=0; i<N; ++i) symm_kpoints_.push_back(K+i*step);
      //---------------------------------
      idx += N;
      symm_pidx_.push_back(idx);
      symm_pname_.push_back("Gamma");
      symm_kpoints_.push_back(Gamma);
    }

    else {
      throw std::range_error("error: BandStruct: not implemented for this lattice");
    }
  }

  //E_kn_.resize(symm_kpoints_.size(),num_bands_);
  this->resize(1); // not used, actually
  setup_done_ = true;
  computation_done_ = false;
}

void BandStruct::reset(void) 
{
  computation_done_ = false;
}

void BandStruct::compute(const kSpace& kspace, const Hamiltonian& ham)
{
  if (computation_done_) return;

  // header information
  if (!is_on()) {
    std::cout << " >>warning: BandStruct::compute:: observable not 'switch on'\n";
    return; 
  } 
  if (!is_open()) open_file();
  fs_ << std::right;
  fs_ << std::scientific << std::uppercase << std::setprecision(6);

  /*
  for (int i=0; i<xvars_.size(); ++i) {
    fs_ << "# ";
    fs_ << xvars_[i].substr(0,14)<<" =";
    fs_ << std::setw(14)<<xvals[i];
    fs_ << std::endl;
  }*/
  if (spin_sector_==spin::UP) {
    fs_ << "# bands: {spin=UP, n=[0:"<<(num_bands_-1)<<"], cols=[5:"<<(5+num_bands_-1)<<"]}\n";
  }
  if (spin_sector_==spin::DN) {
    fs_ << "# bands: {spin=DN, n=[0:"<<(num_bands_-1)<<"], cols=[5:"<<(5+num_bands_-1)<<"]}\n";
  }
  if (spin_sector_==spin::BOTH) {
    fs_ << "# bands: {spin=UP, n=[0:"<<(num_bands_-1)<<"], cols=[5:"<<(5+num_bands_-1)<<"]}; ";
    fs_ << "{spin=DN, n=["<<num_bands_<<":"<<(2*num_bands_-1)<<"], cols=["<<5+num_bands_<<":"<<(5+2*num_bands_-1)<<"]}\n";
  }

  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << "# ";
  fs_ << std::left;
  fs_ << std::setw(6)<<"k"<<std::setw(14)<<"kx"<<std::setw(14)<<"ky"<<std::setw(14)<<"kz";
  if (spin_sector_==spin::UP) {
    for (int n=0; n<num_bands_; ++n) {
      fs_<<std::setw(14)<<"ek_"+std::to_string(n)+"_0"; 
    }
  }
  if (spin_sector_==spin::DN) {
    for (int n=0; n<num_bands_; ++n) {
      fs_<<std::setw(14)<<"ek_"+std::to_string(n)+"_1"; 
    }
  }
  if (spin_sector_==spin::BOTH) {
    for (int n=0; n<num_bands_; ++n) {
      fs_<<std::setw(14)<<"ek_"+std::to_string(n)+"_0"; 
    }
    for (int n=0; n<num_bands_; ++n) {
      fs_<<std::setw(14)<<"ek_"+std::to_string(n)+"_1"; 
    }
  }
  fs_ << std::endl;
  fs_ << "#" << std::string(72, '-') << "\n";
  //---------------------------------------------------------------------


  // compute dispersion
  Eigen::SelfAdjointEigenSolver<ComplexMatrix> es;
  //---------------------------------------------------------------------
  // Energy dispersion along symmetry path
  if (kpath_==kpath_type::SYMM) {
    fs_ << std::right;
    // only UP-spin sector
    if (spin_sector_==spin::UP) {
      for (int k=0; k<symm_kpoints_.size(); ++k) {
        Vector3d kvec = symm_kpoints_[k];
        ham.construct_upspin_block(kvec);
        es.compute(ham.upspin_block(), Eigen::EigenvaluesOnly);
        auto e_kn = es.eigenvalues().transpose();
        fs_<<std::setw(6) << k; 
        fs_<<std::setw(14)<<kvec(0)<<std::setw(14)<<kvec(1)<<std::setw(14)<<kvec(2); 
        for (int n=0; n<num_bands_; ++n) {
          fs_<<std::setw(14)<<e_kn(n); 
        }
        fs_ << std::endl; 
      }
    }
    // only DN-spin sector
    if (spin_sector_==spin::DN) {
      for (int k=0; k<symm_kpoints_.size(); ++k) {
        Vector3d kvec = symm_kpoints_[k];
        ham.construct_dnspin_block(kvec);
        es.compute(ham.dnspin_block(), Eigen::EigenvaluesOnly);
        auto e_kn = es.eigenvalues().transpose();
        fs_<<std::setw(6) << k; 
        fs_<<std::setw(14)<<kvec(0)<<std::setw(14)<<kvec(1)<<std::setw(14)<<kvec(2); 
        for (int n=0; n<num_bands_; ++n) {
          fs_<<std::setw(14)<<e_kn(n); 
        }
        fs_ << std::endl; 
      }
    }
    // BOTH-spin sectors
    if (spin_sector_==spin::BOTH) {
      for (int k=0; k<symm_kpoints_.size(); ++k) {
        Vector3d kvec = symm_kpoints_[k];
        fs_<<std::setw(6) << k; 
        fs_<<std::setw(14)<<kvec(0)<<std::setw(14)<<kvec(1)<<std::setw(14)<<kvec(2); 
        ham.construct_kblock(kvec);
        es.compute(ham.upspin_block(), Eigen::EigenvaluesOnly);
        RealVector e_kn = es.eigenvalues().transpose();
        for (int n=0; n<num_bands_; ++n) {
          fs_<<std::setw(14)<<e_kn(n); 
        }
        es.compute(ham.dnspin_block(), Eigen::EigenvaluesOnly);
        e_kn = es.eigenvalues().transpose();
        for (int n=0; n<num_bands_; ++n) {
          fs_<<std::setw(14)<<e_kn(n); 
        }
        fs_ << std::endl; 
      }
    }
    fs_ << std::endl; 
  }
  //---------------------------------------------------------------------

  //---------------------------------------------------------------------
  // Energy dispersion in the full FBZ
  if (kpath_==kpath_type::FULL) {
    fs_ << std::right;
    // only UP-spin sector
    if (spin_sector_==spin::UP) {
      for (int k=0; k<kspace.num_kpoints(); ++k) {
        Vector3d kvec = kspace.kpoint(k);
        ham.construct_upspin_block(kvec);
        es.compute(ham.upspin_block(), Eigen::EigenvaluesOnly);
        auto e_kn = es.eigenvalues().transpose();
        fs_<<std::setw(6) << k; 
        fs_<<std::setw(14)<<kvec(0)<<std::setw(14)<<kvec(1)<<std::setw(14)<<kvec(2); 
        for (int n=0; n<num_bands_; ++n) {
          fs_<<std::setw(14)<<e_kn(n); 
        }
        fs_ << std::endl; 
      }
    }
    // only DN-spin sector
    if (spin_sector_==spin::DN) {
      for (int k=0; k<kspace.num_kpoints(); ++k) {
        Vector3d kvec = kspace.kpoint(k);
        ham.construct_dnspin_block(kvec);
        es.compute(ham.dnspin_block(), Eigen::EigenvaluesOnly);
        auto e_kn = es.eigenvalues().transpose();
        fs_<<std::setw(6) << k; 
        fs_<<std::setw(14)<<kvec(0)<<std::setw(14)<<kvec(1)<<std::setw(14)<<kvec(2); 
        for (int n=0; n<num_bands_; ++n) {
          fs_<<std::setw(14)<<e_kn(n); 
        }
        fs_ << std::endl; 
      }
    }
    // BOTH-spin sectors
    if (spin_sector_==spin::BOTH) {
      for (int k=0; k<kspace.num_kpoints(); ++k) {
        Vector3d kvec = kspace.kpoint(k);
        fs_<<std::setw(6) << k; 
        fs_<<std::setw(14)<<kvec(0)<<std::setw(14)<<kvec(1)<<std::setw(14)<<kvec(2); 
        ham.construct_kblock(kvec);
        es.compute(ham.upspin_block(), Eigen::EigenvaluesOnly);
        RealVector e_kn = es.eigenvalues().transpose();
        for (int n=0; n<num_bands_; ++n) {
          fs_<<std::setw(14)<<e_kn(n); 
        }
        es.compute(ham.dnspin_block(), Eigen::EigenvaluesOnly);
        e_kn = es.eigenvalues().transpose();
        for (int n=0; n<num_bands_; ++n) {
          fs_<<std::setw(14)<<e_kn(n); 
        }
        fs_ << std::endl; 
      }
    }
    fs_ << std::endl; 
  }
  //---------------------------------------------------------------------

  computation_done_ = true;
}

void BandStruct::print_heading(const std::string& header,
  const std::vector<std::string>& xvars) 
{
  if (!is_on()) return;
  if (heading_printed_) return;
  if (!replace_mode_) return;
  if (!is_open()) open_file();
  xvars_ = xvars;
  fs_ << header;
  fs_ << "# Results: " << name() << "\n";
  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << "# Special k-points:\n";
  fs_ << "# ";
  fs_ << std::right;
  for (const auto& name : symm_pname_) fs_<<std::setw(8)<<name; 
  fs_ << std::endl;
  fs_ << "# ";
  for (const auto& idx : symm_pidx_) fs_<<std::setw(8)<<idx; 
  fs_ << std::endl;
  fs_ << "#" << std::string(72, '-') << "\n";

  fs_ << std::flush;
  heading_printed_ = true;
  close_file();
}

void BandStruct::print_result(const std::vector<double>& xvals) 
{
  // Already done, nothing to be done
  /*
  if (!is_on()) return;
  if (!is_open()) open_file();
  fs_ << std::right;
  fs_ << std::scientific << std::uppercase << std::setprecision(6);

  for (int i=0; i<xvars_.size(); ++i) {
    fs_ << "# ";
    fs_ << xvars_[i].substr(0,14)<<" =";
    fs_ << std::setw(14)<<xvals[i];
    fs_ << std::endl;
  }

  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << "# ";
  fs_ << std::left;
  //for (const auto& p : xvars) fs_ << std::setw(14)<<p.substr(0,14);
  fs_ << std::setw(6)<<"k"<<std::setw(14)<<"kx"<<std::setw(14)<<"ky"<<std::setw(14)<<"kz";
  for (int n=0; n<num_bands_; ++n) {
    fs_<<std::setw(14)<<"Ek_"+std::to_string(n); 
  }
  fs_ << std::endl;
  fs_ << "#" << std::string(72, '-') << "\n";

  // Energy dispersion along symmetry path
  if (kpath_=kpath_type::SYMM) {
    fs_ << std::right;
    for (int k=0; k<symm_kpoints_.size(); ++k) {
      Vector3d kvec = symm_kpoints_[k];
      fs_<<std::setw(6) << k; 
      fs_<<std::setw(14)<<kvec(0)<<std::setw(14)<<kvec(1)<<std::setw(14)<<kvec(2); 
      for (int n=0; n<num_bands_; ++n) {
        fs_<<std::setw(14)<<E_kn_(k,n); 
      }
      fs_ << std::endl; 
    }
    fs_ << std::endl; 
  }

  // Energy dispersion in the full FBZ
  if (kpath_=kpath_type::FULL) {
    for (int k=0; k<kspace.num_kpoints(); ++k) {
      Vector3d kvec = kspace.kpoint(k);
      fs_<<std::setw(6) << k; 
      fs_<<std::setw(14)<<kvec(0)<<std::setw(14)<<kvec(1)<<std::setw(14)<<kvec(2); 
      for (int n=0; n<num_bands_; ++n) {
        fs_<<std::setw(14)<<E_kn_(k,n); 
      }
      fs_ << std::endl; 
    }
    fs_ << std::endl; 
  }

  fs_ << std::flush;
  close_file();
  */
}


} // end namespace vmc