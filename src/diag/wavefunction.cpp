/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2026-01-18 23:08:54
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-18 23:22:24
*----------------------------------------------------------------------------*/
#include <boost/algorithm/string.hpp>
#include "./bandstruct.h"

namespace diag {

void WaveFunction::setup(const input::Parameters& inputs, const lattice::Lattice& lattice, 
  const kSpace& kspace, const Hamiltonian& ham)
{
  MC_Observable::switch_on();
  if (setup_done_) return;

  // setup k-points along symmetry path
  num_bands_ = ham.dimension();
  symm_kpoints_.clear();
  symm_pidx_.clear();
  symm_pname_.clear();


  this->resize(1); // not used, actually
  setup_done_ = true;
  computation_done_ = false;
}

void WaveFunction::reset(void) 
{
  computation_done_ = false;
}

void WaveFunction::compute(const kSpace& kspace, const Hamiltonian& ham)
{
  if (computation_done_) return;

  // header information
  if (!is_on()) {
    std::cout << " >>warning: WaveFunction::compute:: observable not 'switch on'\n";
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

void WaveFunction::print_heading(const std::string& header,
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

void WaveFunction::print_result(const std::vector<double>& xvals) 
{
  // Already done, nothing to be done
}


} // end namespace diag
