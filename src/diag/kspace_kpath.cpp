/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2026-01-19 14:43:15
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-01-21 16:42:00
*----------------------------------------------------------------------------*/
//#include <iomanip>
#include "kspace.h"

namespace diag {

int kSpace::construct_kpath(const lattice::Lattice& lattice)
{
  if (kpath_constructed_) return 0;

	// List special Kpoints for each Bravais lattice
	kpath_.clear();
	kpath_nodes_.clear();
	kpath_nodes_idx_.clear();

  int N = 100;
  int idx = 0;
  if (lattice.brav()==lattice::brav_id::ZERO) {
    kvector Gamma = kvector(0,0,0);
    Gamma_.set("Gamma",Gamma);
  }

  else if (lattice.brav()==lattice::brav_id::CHAIN) {
    kvector Gamma = kvector(0,0,0);
    Gamma_.set("Gamma",Gamma);
    kvector X = 0.5*b1_;
    X_.set("X",X);
    //---------------------------------
		kpath_nodes_.push_back(Gamma_);
		kpath_nodes_idx_.push_back(idx);
    kvector step = (X-Gamma)/N;
    for (int i=0; i<N; ++i) kpath_.push_back(Gamma+i*step);
    //---------------------------------
    idx += N;
		kpath_nodes_.push_back(X_);
		kpath_nodes_idx_.push_back(idx);
    step = (Gamma-X)/N;
    for (int i=0; i<N; ++i) kpath_.push_back(X+i*step);
    //---------------------------------
    idx += N;
		kpath_nodes_.push_back(Gamma_);
		kpath_nodes_idx_.push_back(idx);
    kpath_.push_back(Gamma);
    //---------------------------------
  }

	else if (lattice.brav()==lattice::brav_id::SQUARE) {
    kvector Gamma = kvector(0,0,0);
    Gamma_.set("Gamma",Gamma);
    kvector X = 0.5*b1_;
    X_.set("X",X);
    kvector M = X+0.5*b2_;
    M_.set("M",M);
    //---------------------------------
		kpath_nodes_.push_back(Gamma_);
		kpath_nodes_idx_.push_back(idx);
    kvector step = (X-Gamma)/N;
    for (int i=0; i<N; ++i) kpath_.push_back(Gamma+i*step);
    //---------------------------------
    idx += N;
		kpath_nodes_.push_back(X_);
		kpath_nodes_idx_.push_back(idx);
    step = (M-X)/N;
    for (int i=0; i<N; ++i) kpath_.push_back(X+i*step);
    //---------------------------------
    idx += N;
		kpath_nodes_.push_back(M_);
		kpath_nodes_idx_.push_back(idx);
    step = (Gamma-M)/N;
    for (int i=0; i<N; ++i) kpath_.push_back(M+i*step);
    //---------------------------------
    idx += N;
		kpath_nodes_.push_back(Gamma_);
		kpath_nodes_idx_.push_back(idx);
    kpath_.push_back(Gamma);
	}

	else if (lattice.brav()==lattice::brav_id::HEXAGONAL) {
    kvector Gamma = kvector(0,0,0);
    Gamma_.set("Gamma",Gamma);
    //kvector M = 0.5*(b1_+b2_);
    kvector M = 0.5*b2_;
    M_.set("M",M);
  	//kvector K = (b1_-b2_)/3;
    kvector K = (2*b1_+b2_)/3;
    K_.set("K",K);
  	//kvector Kprime = (2*b1_+b2_)/3;
    kvector Kprime = (2*b2_+b1_)/3;
    Kprime_.set("Kprime",Kprime);
    //---------------------------------
		kpath_nodes_.push_back(Gamma_);
		kpath_nodes_idx_.push_back(idx);
    kvector step = (M-Gamma)/N;
    for (int i=0; i<N; ++i) kpath_.push_back(Gamma+i*step);
    //---------------------------------
    idx += N;
		kpath_nodes_.push_back(M_);
		kpath_nodes_idx_.push_back(idx);
    step = (Kprime-M)/N;
    for (int i=0; i<N; ++i) kpath_.push_back(M+i*step);
    //---------------------------------
    idx += N;
		kpath_nodes_.push_back(Kprime_);
		kpath_nodes_idx_.push_back(idx);
    step = (Gamma-Kprime)/N;
    for (int i=0; i<N; ++i) kpath_.push_back(Kprime+i*step);
    //---------------------------------
    idx += N;
		kpath_nodes_.push_back(Gamma_);
		kpath_nodes_idx_.push_back(idx);
    step = (K-Gamma)/N;
    for (int i=0; i<N; ++i) kpath_.push_back(Gamma+i*step);
    //---------------------------------
    idx += N;
    kpath_nodes_.push_back(K_);
    kpath_nodes_idx_.push_back(idx);
    kpath_.push_back(K);
	}

	else {
    throw std::range_error("error: kSpace::construct_kpath: not defined for this lattice");
	}

  kpath_constructed_ = true;
	return 0;
} 


} // end namespace diag