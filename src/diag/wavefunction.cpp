/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2026-01-18 23:08:54
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-22 08:25:55
*----------------------------------------------------------------------------*/
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include "./wavefunction.h"

namespace diag {

void WaveFunction::setup(const input::Parameters& inputs, const lattice::Lattice& lattice,
  kSpace& kspace)
{
  MC_Observable::switch_on();
  if (setup_done_) return;

  // kpath type 
  std::string kpath = inputs.set_value("kpath", "FULL");
  boost::to_upper(kpath);
  if (kpath=="SYMM") {
    kspace.construct_kpath(lattice);
    kpath_ = kpath_type::SYMM;
  }
  else if (kpath=="FULL") {
    kspace.construct_kmesh(lattice);
    kpath_ = kpath_type::FULL;
  }
  else {
    throw std::range_error("*error: WaveFunction: invalid value for input 'kpath'");
  }

  // Set of QNs for which to compute WFs
  std::string qn_specs = inputs.set_value("qn_list", "");
  // parse the string to read QN values
  std::string::size_type pos1, pos2;
  boost::char_separator<char> comma(",");
  boost::tokenizer<boost::char_separator<char> >::iterator it;

  boost::trim(qn_specs);
  if (qn_specs.front()!='[' || qn_specs.back()!=']') {
    throw std::range_error("*error: WaveFunction::setup: invalid format for input 'qn_list'");
  }
  qn_specs.erase(0,1);

  // read inner lists
  up_states_.clear();
  dn_states_.clear();
  while (true) {
    boost::trim(qn_specs);
    pos1 = qn_specs.find_first_of("(");
    pos2 = qn_specs.find_first_of(")");
    if (pos1==std::string::npos || pos2==std::string::npos) {
      throw std::range_error("*error: WaveFunction::setup: invalid format for input 'qn_list'");
    }
    pos1++;
    int len = pos2-pos1;
    std::string qn_str = qn_specs.substr(pos1,len);
    boost::tokenizer<boost::char_separator<char> > qn_set(qn_str, comma);
    std::vector<int> qn_list;
    for (it=qn_set.begin(); it!=qn_set.end();++it) {
      qn_list.push_back(std::stoi(*it));
      //std::cout << "qn = " <<  *it << "\n";
    }
    if (qn_list.size()!=3) {
      throw std::range_error("*error: WaveFunction::setup: invalid quantum no in 'qn_list'");
    }
    int k = qn_list[0];
    int n = qn_list[1];
    int s = qn_list[2];
    if (s==0) {
      if (up_states_.find(k)!=up_states_.end()) up_states_[k].insert(n);
      else up_states_.insert({k,{n}});
    }
    else if (s==1) {
      if (dn_states_.find(k)!=dn_states_.end()) dn_states_[k].insert(n);
      else dn_states_.insert({k,{n}});
    }
    else {
      throw std::range_error("*error: WaveFunction::setup: invalid spin quantum no in 'qn_list'");
    }

    // parse next
    pos2++;
    qn_specs = qn_specs.substr(pos2);
    boost::trim(qn_specs);
    if (qn_specs.front()==',') {
      qn_specs.erase(0,1);
      continue;
    }
    if (qn_specs.front()==']') break;

    // should not reach here
    throw std::range_error("*error: WaveFunction::setup: invalid format for input 'qn_list'");
  }

  // total number of states read
  num_states_ = 0;
  for (auto it=up_states_.begin(); it!=up_states_.end(); ++it) {
    auto n_set = it->second;
    for (auto it2=n_set.begin(); it2!=n_set.end(); ++it2) {
      num_states_++;
    }
  }
  for (auto it=dn_states_.begin(); it!=dn_states_.end(); ++it) {
    auto n_set = it->second;
    for (auto it2=n_set.begin(); it2!=n_set.end(); ++it2) {
      num_states_++;
    }
  }

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
  int num_bands=ham.dimension();

  // header information
  if (!is_on()) {
    std::cout << " >>warning: WaveFunction::compute:: observable not 'switched on'\n";
    return; 
  } 
  if (!is_open()) open_file();
  fs_ << std::right;
  fs_ << std::scientific << std::uppercase << std::setprecision(6);

  fs_ << std::left;
  int idx = 0;
  for (auto it=up_states_.begin(); it!=up_states_.end();++it) {
    int k = it->first;
    auto n_set = it->second;
    for (auto it2=n_set.begin(); it2!=n_set.end();++it2) {
      int n = *it2;
      fs_ << "# ";
      fs_ << std::setw(8)<<"wf_"+std::to_string(idx)+":";
      fs_ << std::setw(20) << "psi("+std::to_string(k)+","+std::to_string(n)+",0)";
      fs_ << std::setw(10) << "cols: ("+std::to_string(6+2*idx)+","+std::to_string(6+2*idx+1)+")";
      fs_ << std::endl;
      idx++;
    }
  }
  for (auto it=dn_states_.begin(); it!=dn_states_.end();++it) {
    int k = it->first;
    auto n_set = it->second;
    for (auto it2=n_set.begin(); it2!=n_set.end();++it2) {
      int n = *it2;
      fs_ << "# ";
      fs_ << std::setw(8)<<"wf_"+std::to_string(idx)+":";
      fs_ << std::setw(20) << "psi("+std::to_string(k)+","+std::to_string(n)+",1)";
      fs_ << std::setw(10) << "cols: ("+std::to_string(1+idx)+","+std::to_string(2+idx)+")";
      fs_ << std::endl;
      idx++;
    }
  }

  /*
  for (int i=0; i<xvars_.size(); ++i) {
    fs_ << "# ";
    fs_ << xvars_[i].substr(0,14)<<" =";
    fs_ << std::setw(14)<<xvals[i];
    fs_ << std::endl;
  }*/
  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << "# ";
  fs_ << std::left;
  fs_ << std::setw(8)<<"basis";
  fs_ << std::setw(8)<<"label";
  fs_ << std::setw(14)<<"xcoord";
  fs_ << std::setw(14)<<"ycoord";
  fs_ << std::setw(14)<<"zcoord";
  idx = 1;
  for (auto it=up_states_.begin(); it!=up_states_.end();++it) {
    auto n_set = it->second;
    for (auto it2=n_set.begin(); it2!=n_set.end();++it2) {
      fs_ << std::setw(14) << "Re[wf_"+std::to_string(idx)+"]";
      fs_ << std::setw(14) << "Im[wf_"+std::to_string(idx)+"]";
      idx++;
    }
  }
  for (auto it=dn_states_.begin(); it!=dn_states_.end();++it) {
    auto n_set = it->second;
    for (auto it2=n_set.begin(); it2!=n_set.end();++it2) {
      fs_ << std::setw(14) << "Re[wf_"+std::to_string(idx)+"]";
      fs_ << std::setw(14) << "Im[wf_"+std::to_string(idx)+"]";
      idx++;
    }
  }
  fs_ << std::endl;
  fs_ << "#" << std::string(72, '-') << "\n";
  //---------------------------------------------------------------------

  // compute wavefunctions
  Matrix Psi(num_bands,num_states_);
  Eigen::SelfAdjointEigenSolver<ComplexMatrix> es;
  int istate = 0;
  //---------------------------------------------------------------------
  // Energy dispersion along symmetry path
  if (kpath_==kpath_type::SYMM) {
    // UP-spin states
    for (auto it=up_states_.begin(); it!=up_states_.end();++it) {
      int k = it->first;
      if (k>=kspace.kpath().size()) {
        throw std::range_error("*error: WaveFunction::compute: quantum no 'k' out of range");
      }
      kvector kvec = kspace.kpath()[k];
      // compute the eigensystem 
      ham.construct_upspin_block(kvec);
      es.compute(ham.upspin_block());
      // print the eigen functions for the set of 'n' values
      auto n_set = it->second;
      for (auto it2=n_set.begin(); it2!=n_set.end();++it2) {
        int n = *it2;
        if (n>=num_bands) {
          throw std::range_error("*error: WaveFunction::compute: quantum no 'n' out of range");
        }
        Psi.col(istate) = es.eigenvectors().col(n);
        istate++;
      }
    }
    // DN-spin states
    for (auto it=dn_states_.begin(); it!=dn_states_.end();++it) {
      int k = it->first;
      if (k>=kspace.kpath().size()) {
        throw std::range_error("*error: WaveFunction::compute: quantum no 'k' out of range");
      }
      kvector kvec = kspace.kpath()[k];
      // compute the eigensystem 
      ham.construct_dnspin_block(kvec);
      es.compute(ham.dnspin_block());
      // print the eigen functions for the set of 'n' values
      auto n_set = it->second;
      for (auto it2=n_set.begin(); it2!=n_set.end();++it2) {
        int n = *it2;
        if (n>=num_bands) {
          throw std::range_error("*error: WaveFunction::compute: quantum no 'n' out of range");
        }
        Psi.col(istate) = es.eigenvectors().col(n);
        istate++;
      }
    }
  }
  //---------------------------------------------------------------------

  //---------------------------------------------------------------------
  // Energy dispersion in the full FBZ
  if (kpath_==kpath_type::FULL) {
    // UP-spin states
    for (auto it=up_states_.begin(); it!=up_states_.end();++it) {
      int k = it->first;
      if (k>=kspace.kmesh().size()) {
        throw std::range_error("*error: WaveFunction::compute: quantum no 'k' out of range");
      }
      kvector kvec = kspace.kmesh()[k];
      // compute the eigensystem 
      ham.construct_upspin_block(kvec);
      es.compute(ham.upspin_block());
      // print the eigen functions for the set of 'n' values
      auto n_set = it->second;
      for (auto it2=n_set.begin(); it2!=n_set.end();++it2) {
        int n = *it2;
        if (n>=num_bands) {
          throw std::range_error("*error: WaveFunction::compute: quantum no 'n' out of range");
        }
        Psi.col(istate) = es.eigenvectors().col(n);
        istate++;
      }
    }
    // DN-spin states
    for (auto it=dn_states_.begin(); it!=dn_states_.end();++it) {
      int k = it->first;
      if (k>=kspace.kmesh().size()) {
        throw std::range_error("*error: WaveFunction::compute: quantum no 'k' out of range");
      }
      kvector kvec = kspace.kmesh()[k];
      // compute the eigensystem 
      ham.construct_dnspin_block(kvec);
      es.compute(ham.dnspin_block());
      // print the eigen functions for the set of 'n' values
      auto n_set = it->second;
      for (auto it2=n_set.begin(); it2!=n_set.end();++it2) {
        int n = *it2;
        if (n>=num_bands) {
          throw std::range_error("*error: WaveFunction::compute: quantum no 'n' out of range");
        }
        Psi.col(istate) = es.eigenvectors().col(n);
        istate++;
      }
    }
  }
  //---------------------------------------------------------------------

  // Print the wave functions
  fs_ << std::right;
  for (int i=0; i<num_bands; ++i) {
    fs_<<std::setw(8) << i; 
    fs_<<std::setw(8) << ham.basis_state_labels()[i]; 
    auto R = ham.basis_state_coords()[i];
    fs_<<std::setw(14)<<R(0)<<std::setw(14)<<R(1)<<std::setw(14)<<R(2); 
    for (int j=0; j<num_states_; ++j) {
      fs_<<std::setw(14)<<std::real(Psi(i,j)); 
      fs_<<std::setw(14)<<std::imag(Psi(i,j)); 
    }
    fs_ << std::endl; 
  }
  fs_ << std::endl; 

  close_file();
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
  fs_ << std::flush;
  heading_printed_ = true;
  close_file();
}

void WaveFunction::print_result(const std::vector<double>& xvals) 
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


} // end namespace diag
