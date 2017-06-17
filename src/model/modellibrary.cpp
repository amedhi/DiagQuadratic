/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-11 13:02:35
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-06-17 12:02:29
*----------------------------------------------------------------------------*/
#include <cmath>
#include "model.h"
#include <boost/algorithm/string.hpp>

namespace model {

int Hamiltonian::define_model(const input::Parameters& inputs, 
  const lattice::Lattice& lattice)
{
  //int info;
  //unsigned ntypes;
  //std::vector<MatrixElement> matrix_elem(20);
  double defval;
  //unsigned sitetype, change, type, src_type, tgt_type;
  std::string name; //, matrixelem, op, qn, site, src, tgt, fact;
  //SiteBasis site_basis;
  //BasisDescriptor basis;
  //QuantumNumber::value_type min, max, step;
  CouplingConstant cc;

  // define the models 
  model_name = inputs.set_value("model", "HUBBARD");
  boost::to_upper(model_name);

  if (model_name == "HUBBARD") {
    mid = model_id::HUBBARD;
    // model parameters
    add_parameter(name="t", defval=1.0, inputs);
    add_parameter(name="U", defval=0.0, inputs);
    // bond operator terms
    add_bondterm(name="hopping", cc="-t", op::spin_hop());
    add_siteterm(name="hubbard", cc="U", op::hubbard_int());
  }

  else if (model_name == "TJ") {
    mid = model_id::TJ;
    int nowarn;
    if (inputs.set_value("no_double_occupancy",true,nowarn))
      set_no_dbloccupancy();
    // model parameters
    add_parameter(name="t", defval=1.0, inputs);
    add_parameter(name="J", defval=0.0, inputs);
    // bond operator terms
    add_bondterm(name="hopping", cc="-t", op::spin_hop());
    add_bondterm(name="exchange", cc="J", op::sisj_plus());
  }

  else if (model_name == "TBI_HUBBARD") {
    mid = model_id::TBI_HUBBARD;
    switch (lattice.id()) {
      case lattice::lattice_id::SQUARE_2BAND:
        add_parameter(name="e0", defval=1.0, inputs);
        add_parameter(name="t", defval=1.0, inputs);
        add_parameter(name="tsp", defval=1.0, inputs);
        // site operators
        cc = CouplingConstant({0,"-(e0-4*t)"}, {1,"(e0-4*t)"});
        add_siteterm(name="hopping", cc, op::ni_up());
        add_siteterm(name="hopping", cc, op::ni_dn());

        // upspin hop
        cc = CouplingConstant({0,"-t"}, {1,"-i*tsp"},
          {2,"-i*tsp"}, {3,"t"}, {4,"-tsp"},{5,"tsp"});
        add_bondterm(name="hopping", cc, op::upspin_hop());
        // dnspin hop
        cc = CouplingConstant({0,"-t"}, {1,"i*tsp"},
          {2,"i*tsp"}, {3,"t"}, {4,"-tsp"},{5,"tsp"});
        add_bondterm(name="hopping", cc, op::dnspin_hop());
        // Hubbard interaction
        add_parameter(name="U", defval=0.0, inputs);
        add_siteterm(name="hubbard", cc="U", op::hubbard_int());
        break;

      case lattice::lattice_id::HONEYCOMB:
        add_parameter(name="t", defval=1.0, inputs);
        add_parameter(name="lambda", defval=1.0, inputs);
        // upspin hop
        cc = CouplingConstant({0,"-t"}, {1,"i*lambda"},
          {2,"-i*lambda"}, {3,"-i*lambda"}, {4,"i*lambda"});
        add_bondterm(name="hopping", cc, op::upspin_hop());
        // dnspin hop
        cc = CouplingConstant({0,"-t"}, {1,"-i*lambda"},
          {2,"i*lambda"}, {3,"i*lambda"}, {4,"-i*lambda"});
        add_bondterm(name="hopping", cc, op::dnspin_hop());
        // Hubbard interaction
        add_parameter(name="U", defval=0.0, inputs);
        add_siteterm(name="hubbard", cc="U", op::hubbard_int());
        break;

      case lattice::lattice_id::KAGOME:
        add_parameter(name="t", defval=1.0, inputs);
        add_parameter(name="t2", defval=1.0, inputs);
        add_parameter(name="lambda", defval=1.0, inputs);
        add_parameter(name="lambda2", defval=1.0, inputs);
        // upspin hop
        cc = CouplingConstant({0,"-t+i*lambda"},{1,"-t-i*lambda"},
          {2,"-t2+i*lambda2"},{3,"-t2-i*lambda2"});
        add_bondterm(name="hopping", cc, op::upspin_hop());
        // dnspin hop
        cc = CouplingConstant({0,"-t-i*lambda"},{1,"-t+i*lambda"},
          {2,"-t2-i*lambda2"},{3,"-t2+i*lambda2"});
        add_bondterm(name="hopping", cc, op::dnspin_hop());
        // Hubbard interaction
        add_parameter(name="U", defval=0.0, inputs);
        add_siteterm(name="hubbard", cc="U", op::hubbard_int());
        break;

      case lattice::lattice_id::PYROCHLORE:
        add_parameter(name="t", defval=1.0, inputs);
        add_parameter(name="t2", defval=1.0, inputs);
        add_parameter(name="lambda", defval=1.0, inputs);
        add_parameter(name="lambda2", defval=1.0, inputs);
        add_parameter(name="th", defval=1.0, inputs);
        // upspin hop
        cc = CouplingConstant({0,"-t+i*lambda"},{1,"-t-i*lambda"},
          {2,"-t2+i*lambda2"},{3,"-t2-i*lambda2"},{4,"-th"});
        add_bondterm(name="hopping", cc, op::upspin_hop());
        // dnspin hop
        cc = CouplingConstant({0,"-t-i*lambda"},{1,"-t+i*lambda"},
          {2,"-t2-i*lambda2"},{3,"-t2+i*lambda2"},{4,"-th"});
        add_bondterm(name="hopping", cc, op::dnspin_hop());
        // Hubbard interaction
        add_parameter(name="U", defval=0.0, inputs);
        add_siteterm(name="hubbard", cc="U", op::hubbard_int());
        break;

      default:
        throw std::range_error("*error: modellibrary: model not defined for the lattice"); 
    }
  }


  /*------------- undefined model--------------*/
  else {
    throw std::range_error("*error: modellibrary: undefined model");
  }

  // if the model has site disorder
  /*
  if (site_disorder) {
    add_disorder_term(name="disorder", op::ni_sigma());
  }*/
  
  return 0;
}

int Hamiltonian::construct(const input::Parameters& inputs, 
  const lattice::Lattice& lattice)
{
  init(lattice);
  define_model(inputs, lattice);
  finalize(lattice);
  return 0;
}


} // end namespace model
