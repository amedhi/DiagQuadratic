/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2025-12-06 11:31:10
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-03-22 20:03:58
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

  else if (model_name == "MKM") {
    mid = model_id::MKM;
    if (lattice.id()==lattice::lattice_id::HONEYCOMB2) {
      // model parameters
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="eta", defval=1.0, inputs);
      add_parameter(name="lambda", defval=1.0, inputs);
      add_parameter(name="lambda_R", defval=1.0, inputs);

      // constants
      double aa = 1.0;
      double dx = 0.5*std::sqrt(3.0)*aa;
      double dy = 0.5*aa;
      add_constant(name="dx", dx);
      add_constant(name="dy", dy);

      // bond operator terms
      cc.create(5);
      /*
      // first NN-bond
      expr_mat.resize(2,2);
      expr_mat(0,0) = "-t";
      expr_mat(0,1) = "lambda_R*(-dx+i*dy)";
      expr_mat(1,0) = "lambda_R*(dx+i*dy)";
      expr_mat(1,1) = "-t";
      cc.add_type(0, expr_mat); 

      // second NN-bond
      expr_mat.resize(2,2);
      expr_mat(0,0) = "-t*eta";
      expr_mat(0,1) = "lambda_R*(-dx-i*dy)";
      expr_mat(1,0) = "lambda_R*(dx-i*dy)";
      expr_mat(1,1) = "-t*eta";
      cc.add_type(1, expr_mat); 

      // third NN-bond
      expr_mat.resize(2,2);
      expr_mat(0,0) = "-t*eta";
      expr_mat(0,1) = "i*lambda_R";
      expr_mat(1,0) = "i*lambda_R";
      expr_mat(1,1) = "-t*eta";
      cc.add_type(2, expr_mat); 
      */
      // vertical bond is taken as the "strong" bond
      // first NN-bond
      expr_mat.resize(2,2);
      expr_mat(0,0) = "-t*eta";
      expr_mat(0,1) = "lambda_R*(-dx+i*dy)";
      expr_mat(1,0) = "lambda_R*(dx+i*dy)";
      expr_mat(1,1) = "-t*eta";
      cc.add_type(0, expr_mat); 

      // second NN-bond
      expr_mat.resize(2,2);
      expr_mat(0,0) = "-t*eta";
      expr_mat(0,1) = "lambda_R*(-dx-i*dy)";
      expr_mat(1,0) = "lambda_R*(dx-i*dy)";
      expr_mat(1,1) = "-t*eta";
      cc.add_type(1, expr_mat); 

      // third NN-bond
      expr_mat.resize(2,2);
      expr_mat(0,0) = "-t";
      expr_mat(0,1) = "i*lambda_R";
      expr_mat(1,0) = "i*lambda_R";
      expr_mat(1,1) = "-t";
      cc.add_type(2, expr_mat); 

      // NNN-bonds
      expr_mat.resize(2,2);
      expr_mat(0,0) = "-i*lambda";
      expr_mat(0,1) = "0.0";
      expr_mat(1,0) = "0.0";
      expr_mat(1,1) = "i*lambda";
      cc.add_type(3, expr_mat); 

      expr_mat.resize(2,2);
      expr_mat(0,0) = "i*lambda";
      expr_mat(0,1) = "0.0";
      expr_mat(1,0) = "0.0";
      expr_mat(1,1) = "-i*lambda";
      cc.add_type(4, expr_mat); 
      add_bondterm(name="hopping", cc, op::upspin_hop());
    }

    else if (lattice.id()==lattice::lattice_id::HONEYCOMB3) {
      // model parameters
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="eta", defval=1.0, inputs);
      add_parameter(name="lambda", defval=1.0, inputs);
      add_parameter(name="lambda_R", defval=1.0, inputs);

      // bond operator terms
      cc.create(5);
      cc.add_type(0,"-t");
      cc.add_type(1,"-t*eta");
      cc.add_type(2,"-t*eta");
      cc.add_type(3,"-i*lambda");
      cc.add_type(4,"i*lambda");
      add_bondterm(name="hopping", cc, op::upspin_hop());

      /*
      // constants
      double aa = 1.0;
      double dx = 0.5*std::sqrt(3.0)*aa;
      double dy = 0.5*aa;
      add_constant(name="dx", dx);
      add_constant(name="dy", dy);

      // bond operator terms
      cc.create(5);
      // first NN-bond
      expr_mat.resize(2,2);
      expr_mat(0,0) = "-t";
      expr_mat(0,1) = "lambda_R*(-dx+i*dy)";
      expr_mat(1,0) = "lambda_R*(dx+i*dy)";
      expr_mat(1,1) = "-t";
      cc.add_type(0, expr_mat); 

      // second NN-bond
      expr_mat.resize(2,2);
      expr_mat(0,0) = "-t*eta";
      expr_mat(0,1) = "lambda_R*(-dx-i*dy)";
      expr_mat(1,0) = "lambda_R*(dx-i*dy)";
      expr_mat(1,1) = "-t*eta";
      cc.add_type(1, expr_mat); 

      // third NN-bond
      expr_mat.resize(2,2);
      expr_mat(0,0) = "-t*eta";
      expr_mat(0,1) = "i*lambda_R";
      expr_mat(1,0) = "i*lambda_R";
      expr_mat(1,1) = "-t*eta";
      cc.add_type(2, expr_mat); 

      // NNN-bonds
      expr_mat.resize(2,2);
      expr_mat(0,0) = "-i*lambda";
      expr_mat(0,1) = "0.0";
      expr_mat(1,0) = "0.0";
      expr_mat(1,1) = "i*lambda";
      cc.add_type(3, expr_mat); 

      expr_mat.resize(2,2);
      expr_mat(0,0) = "i*lambda";
      expr_mat(0,1) = "0.0";
      expr_mat(1,0) = "0.0";
      expr_mat(1,1) = "-i*lambda";
      cc.add_type(4, expr_mat); 
      add_bondterm(name="hopping", cc, op::upspin_hop());
      */
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
