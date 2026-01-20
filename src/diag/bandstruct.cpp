/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2025-12-09 17:26:56
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-19 22:15:22
*----------------------------------------------------------------------------*/
#include <boost/algorithm/string.hpp>
#include "./bandstruct.h"

namespace diag {

void BandStruct::setup(const input::Parameters& inputs, const lattice::Lattice& lattice,
  kSpace& kspace)
{
  MC_Observable::switch_on();
  if (setup_done_) return;

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

  // kpoints list 
  std::string kpath = inputs.set_value("kpath", "FULL");
  boost::to_upper(kpath);
  if (kpath=="SYMM") {
    kspace.construct_kpath(lattice);
    kpath_ = kpath_type::SYMM;
    kpath_nodes_ = kspace.kpath_nodes();
    kpath_nodes_idx_ = kspace.kpath_nodes_idx();
  }
  else if (kpath=="FULL") {
    kspace.construct_kmesh(lattice);
    kpath_ = kpath_type::FULL;
  }
  else {
    throw std::range_error("error: BandStruct: invalid value for input 'kpath'");
  }

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
  int num_bands=ham.dimension();

  // header information
  if (!is_on()) {
    std::cout << " >>warning: BandStruct::compute:: observable not 'switched on'\n";
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
    fs_ << "# bands: {spin=UP, n=[0:"<<(num_bands-1)<<"], cols=[5:"<<(5+num_bands-1)<<"]}\n";
  }
  if (spin_sector_==spin::DN) {
    fs_ << "# bands: {spin=DN, n=[0:"<<(num_bands-1)<<"], cols=[5:"<<(5+num_bands-1)<<"]}\n";
  }
  if (spin_sector_==spin::BOTH) {
    fs_ << "# bands: {spin=UP, n=[0:"<<(num_bands-1)<<"], cols=[5:"<<(5+num_bands-1)<<"]}; ";
    fs_ << "{spin=DN, n=[0:"<<(num_bands-1)<<"], cols=["<<5+num_bands<<":"<<(5+2*num_bands-1)<<"]}\n";
  }

  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << "# ";
  fs_ << std::left;
  fs_ << std::setw(6)<<"k"<<std::setw(14)<<"kx"<<std::setw(14)<<"ky"<<std::setw(14)<<"kz";
  if (spin_sector_==spin::UP) {
    for (int n=0; n<num_bands; ++n) {
      fs_<<std::setw(14)<<"ek_"+std::to_string(n)+"_0"; 
    }
  }
  if (spin_sector_==spin::DN) {
    for (int n=0; n<num_bands; ++n) {
      fs_<<std::setw(14)<<"ek_"+std::to_string(n)+"_1"; 
    }
  }
  if (spin_sector_==spin::BOTH) {
    for (int n=0; n<num_bands; ++n) {
      fs_<<std::setw(14)<<"ek_"+std::to_string(n)+"_0"; 
    }
    for (int n=0; n<num_bands; ++n) {
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
      int kidx = 0;
      for (const auto& kvec : kspace.kpath()) {
        ham.construct_upspin_block(kvec);
        es.compute(ham.upspin_block(), Eigen::EigenvaluesOnly);
        auto e_kn = es.eigenvalues().transpose();
        fs_<<std::setw(6) << kidx++; 
        fs_<<std::setw(14)<<kvec(0)<<std::setw(14)<<kvec(1)<<std::setw(14)<<kvec(2); 
        for (int n=0; n<num_bands; ++n) {
          fs_<<std::setw(14)<<e_kn(n); 
        }
        fs_ << std::endl; 
      }
    }
    // only DN-spin sector
    if (spin_sector_==spin::DN) {
      int kidx = 0;
      for (const auto& kvec : kspace.kpath()) {
        ham.construct_dnspin_block(kvec);
        es.compute(ham.dnspin_block(), Eigen::EigenvaluesOnly);
        auto e_kn = es.eigenvalues().transpose();
        fs_<<std::setw(6) << kidx++; 
        fs_<<std::setw(14)<<kvec(0)<<std::setw(14)<<kvec(1)<<std::setw(14)<<kvec(2); 
        for (int n=0; n<num_bands; ++n) {
          fs_<<std::setw(14)<<e_kn(n); 
        }
        fs_ << std::endl; 
      }
    }
    // BOTH-spin sectors
    if (spin_sector_==spin::BOTH) {
      int kidx = 0;
      for (const auto& kvec : kspace.kpath()) {
        fs_<<std::setw(6) << kidx++; 
        fs_<<std::setw(14)<<kvec(0)<<std::setw(14)<<kvec(1)<<std::setw(14)<<kvec(2); 
        ham.construct_kblock(kvec);
        es.compute(ham.upspin_block(), Eigen::EigenvaluesOnly);
        RealVector e_kn = es.eigenvalues().transpose();
        for (int n=0; n<num_bands; ++n) {
          fs_<<std::setw(14)<<e_kn(n); 
        }
        es.compute(ham.dnspin_block(), Eigen::EigenvaluesOnly);
        e_kn = es.eigenvalues().transpose();
        for (int n=0; n<num_bands; ++n) {
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
      int kidx = 0;
      for (const auto& kvec : kspace.kmesh()) {
        ham.construct_upspin_block(kvec);
        es.compute(ham.upspin_block(), Eigen::EigenvaluesOnly);
        auto e_kn = es.eigenvalues().transpose();
        fs_<<std::setw(6) << kidx++; 
        fs_<<std::setw(14)<<kvec(0)<<std::setw(14)<<kvec(1)<<std::setw(14)<<kvec(2); 
        for (int n=0; n<num_bands; ++n) {
          fs_<<std::setw(14)<<e_kn(n); 
        }
        fs_ << std::endl; 
      }
    }
    // only DN-spin sector
    if (spin_sector_==spin::DN) {
      int kidx = 0;
      for (const auto& kvec : kspace.kmesh()) {
        ham.construct_dnspin_block(kvec);
        es.compute(ham.dnspin_block(), Eigen::EigenvaluesOnly);
        auto e_kn = es.eigenvalues().transpose();
        fs_<<std::setw(6) << kidx++; 
        fs_<<std::setw(14)<<kvec(0)<<std::setw(14)<<kvec(1)<<std::setw(14)<<kvec(2); 
        for (int n=0; n<num_bands; ++n) {
          fs_<<std::setw(14)<<e_kn(n); 
        }
        fs_ << std::endl; 
      }
    }
    // BOTH-spin sectors
    if (spin_sector_==spin::BOTH) {
      int kidx = 0;
      for (const auto& kvec : kspace.kmesh()) {
        fs_<<std::setw(6) << kidx++; 
        fs_<<std::setw(14)<<kvec(0)<<std::setw(14)<<kvec(1)<<std::setw(14)<<kvec(2); 
        ham.construct_kblock(kvec);
        es.compute(ham.upspin_block(), Eigen::EigenvaluesOnly);
        RealVector e_kn = es.eigenvalues().transpose();
        for (int n=0; n<num_bands; ++n) {
          fs_<<std::setw(14)<<e_kn(n); 
        }
        es.compute(ham.dnspin_block(), Eigen::EigenvaluesOnly);
        e_kn = es.eigenvalues().transpose();
        for (int n=0; n<num_bands; ++n) {
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
  if (kpath_==kpath_type::SYMM) {
    fs_ << "# Special k-points:\n";
    fs_ << "# ";
    fs_ << std::right;
    for (const auto& G : kpath_nodes_) fs_<<std::setw(8)<<G.name(); 
    fs_ << std::endl;
    fs_ << "# ";
    for (const auto& idx : kpath_nodes_idx_) fs_<<std::setw(8)<<idx; 
    fs_ << std::endl;
    fs_ << "#" << std::string(72, '-') << "\n";
  }

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