/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2025-12-06 11:34:31
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-19 11:29:49
*----------------------------------------------------------------------------*/
#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <sstream>
#include <string>
#include <complex>
#include <vector>
#include <map>
#include <stdexcept>
#include "strmatrix.h"
#include "../lattice/lattice.h"
#include "./model_term.h"

namespace model {

enum class model_id {
  UNDEFINED, HUBBARD, KM, MKM, TJ, TBI_HUBBARD, 
  PYROCHLORE, NICKELATE
};

class Model 
{
public:
  using siteterm_iterator = std::vector<ModelTerm>::const_iterator; 
  using bondterm_iterator = std::vector<ModelTerm>::const_iterator; 
  Model() {}
  Model(const input::Parameters& inputs, const lattice::Lattice& lattice)
  { construct(inputs, lattice); }
  ~Model() {}
  int construct(const input::Parameters& inputs, const lattice::Lattice& lattice);
  virtual int init(const lattice::Lattice& lattice);
  int define_model(const input::Parameters& inputs, const lattice::Lattice& lattice);
  int finalize(const lattice::Lattice& lattice);

  //unsigned add_sitebasis(SiteBasis& sitebasis);
  //unsigned add_sitebasis(const unsigned& type, SiteBasis& sitebasis);
  int add_parameter(const std::string& pname, const double& defval, 
    const input::Parameters& inputs)
    { parms_[pname] = inputs.set_value(pname, defval); return parms_.size(); }
  int add_parameter(const std::string& pname, const double& defval, 
    const input::Parameters& inputs, int& info)
    { parms_[pname] = inputs.set_value(pname, defval, info); return parms_.size(); }
  int add_parameter(const std::string& pname, const double& val) 
    { parms_[pname] = val; return parms_.size(); }
  void update_parameters(const input::Parameters& inputs);
  void update_parameter(const std::string& pname, const double& val); 
  virtual void update_terms(void);
  //void change_parameter_value(const std::string& pname, const double& pval);
  double get_parameter_value(const std::string& pname) const;
  int add_constant(const std::string& cname, const double& val) 
    { constants_.insert({cname, val}); return constants_.size(); }
  int add_siteterm(const std::string& name, const CouplingConstant& cc, const op::quantum_op& op);
  int add_bondterm(const std::string& name, const CouplingConstant& cc, const op::quantum_op& op);
  void set_no_dbloccupancy(void) { double_occupancy_=false; }

  //const BasisDescriptor& basis(void) const { return basis_; }
  //const SiteBasis& site_basis(const unsigned& site_type) const { return basis_.at(site_type); }
  //unsigned sitebasis_dimension(const unsigned& site_type) const
  //{ return basis_.dimension(site_type); }
  const ModelParams& parameters(void) const { return parms_; }
  const ModelParams& constants(void) const { return constants_; }
  const bool& double_occupancy(void) const { return double_occupancy_; }
  const bool& have_disorder_term(void) const { return have_disorder_term_; }
  const bool& have_siteterm(void) const { return have_siteterm_; }
  const bool& have_bondterm(void) const { return have_bondterm_; }
  const siteterm_iterator& siteterms_begin(void) const { return st_begin_; }
  const siteterm_iterator& siteterms_end(void) const { return st_end_; }
  const bondterm_iterator& bondterms_begin(void) const { return bt_begin_; }
  const bondterm_iterator& bondterms_end(void) const { return bt_end_; }
  const siteterm_iterator& disorder_term_begin(void) const { return dterm_begin_; }
  const siteterm_iterator& disorder_term_end(void) const { return dterm_end_; }
  std::pair<siteterm_iterator, siteterm_iterator> site_terms(void) const 
    { return std::make_pair(site_terms_.cbegin(), site_terms_.cend()); }
  std::pair<bondterm_iterator, bondterm_iterator> bond_terms(void) const 
    { return std::make_pair(bond_terms_.cbegin(), bond_terms_.cend()); }
  std::pair<siteterm_iterator, siteterm_iterator> disorder_terms(void) const 
    { return std::make_pair(disorder_terms_.cbegin(), disorder_terms_.cend()); }
  int num_siteterms(void) const { return site_terms_.size(); }
  int num_bondterms(void) const { return bond_terms_.size(); }
  int num_disorder_terms(void) const { return disorder_terms_.size(); }
  int num_terms(void) const 
    { return site_terms_.size()+bond_terms_.size()+disorder_terms_.size(); }
  void get_term_names(std::vector<std::string>& term_names) const;
  const std::string& signature_str(void) const { return signature_str_; }
  const std::string& info_str(void) const { return info_str_; }
  std::ostream& print_info(std::ostream& os) const { return os << info_str_; }
private:
  model_id mid {model_id::UNDEFINED};
  std::string model_name;
  //BasisDescriptor basis_;
  std::map<int,int> sitetypes_map_;
  std::map<int,int> bondtypes_map_;
  //std::map<int,int> type_dim_map_;
  //BondTerm::BondSiteMap bond_sites_map_;  
  std::vector<ModelTerm> bond_terms_;
  std::vector<ModelTerm> site_terms_;

  bool double_occupancy_{true};
  bool have_siteterm_{false};
  bool have_bondterm_{false};
  siteterm_iterator st_begin_;
  siteterm_iterator st_end_;
  bondterm_iterator bt_begin_;
  bondterm_iterator bt_end_;
  //std::vector<SiteTermDescriptor> siteterms_;
  //std::vector<BondTermDescriptor> bondterms_;
  ModelParams parms_;
  ModelParams constants_;

  // disorder term
  bool have_disorder_term_{false};
  //SiteDisorder site_disorder_;
  std::vector<ModelTerm> disorder_terms_;
  siteterm_iterator dterm_begin_;
  siteterm_iterator dterm_end_;

  //std::ostringstream info_str_;
  //std::ostringstream signature_str_;
  std::string info_str_;
  std::string signature_str_;
  void set_info_string(const lattice::Lattice& lattice); 
};


} // end namespace model

#endif