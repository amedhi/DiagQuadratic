/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2025-12-06 11:31:10
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-17 21:43:36
*----------------------------------------------------------------------------*/
#include <cmath>
#include "model.h"
#include <boost/algorithm/string.hpp>

namespace model {

int Model::define_model(const input::Parameters& inputs, const lattice::Lattice& lattice)
{
  //int info;
  //std::vector<MatrixElement> matrix_elem(20);
  double defval;
  std::string name, path; //, matrixelem, op, qn, site, src, tgt, fact;
  CouplingConstant cc;

  // define the models 
  model_name = inputs.set_value("model", "HUBBARD");
  boost::to_upper(model_name);

  strMatrix expr_mat;
  strMatrix::row_t expr_vec;

  if (model_name == "HUBBARD") {
    mid = model_id::HUBBARD;
    if (lattice.id()==lattice::lattice_id::SQUARE) {
      // model parameters
      add_parameter(name="t", defval=1.0, inputs);
      // bond operator terms
      add_bondterm(name="hopping", cc="-t", op::spin_hop());
    }
    else {
      throw std::range_error("*error: modellibrary: model not defined for this lattice");
    }
  }

  else if (model_name == "KM") {
    mid = model_id::KM;
    if (lattice.id()==lattice::lattice_id::HONEYCOMB) {
      // model parameters
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="lambda", defval=1.0, inputs);

      // bond operator terms
      cc.create(5);
      cc.add_type(0,"-t");
      cc.add_type(1,"-t");
      cc.add_type(2,"-t");
      cc.add_type(3,"-i*lambda");
      cc.add_type(4,"i*lambda");
      add_bondterm(name="hopping", cc, op::upspin_hop());

      cc.create(5);
      cc.add_type(0,"-t");
      cc.add_type(1,"-t");
      cc.add_type(2,"-t");
      cc.add_type(3,"i*lambda");
      cc.add_type(4,"-i*lambda");
      add_bondterm(name="hopping", cc, op::dnspin_hop());
    }
    else {
      throw std::range_error("*error: modellibrary: model not defined for this lattice");
    }
  }

  /*------------- undefined model--------------*/
  else {
    throw std::range_error("*error: modellibrary: undefined model");
  }

  return 0;
}

int Model::construct(const input::Parameters& inputs, 
  const lattice::Lattice& lattice)
{
  init(lattice);
  define_model(inputs, lattice);
  finalize(lattice);
  return 0;
}


} // end namespace model
