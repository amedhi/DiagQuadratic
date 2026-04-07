/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2026-04-06 10:53:30
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-04-07 13:10:26
*----------------------------------------------------------------------------*/
#include <boost/algorithm/string.hpp>
#include "./entanglement.h"

namespace diag {

void EntanglementEntropy::setup(const input::Parameters& inputs, const lattice::Lattice& lattice,
  kSpace& kspace)
{
  MC_Observable::switch_on();
  if (setup_done_) return;


	if (lattice.id()!=lattice::lattice_id::HONEYCOMB3) {
    throw std::range_error("error: EntanglementEntropy: not implemented for the lattice");
	}

	if ((lattice.input_size1()!=10)&&(lattice.input_size2()!=10)) {
    throw std::range_error("error: EntanglementEntropy: implemented only for a 10x10 lattice");
	}

	// List the sites belonging to subsystems A & B
	subsys_A_.clear();
	for (int i=12; i<=18; ++i) subsys_A_.push_back(i);
	for (int i=31; i<=38; ++i) subsys_A_.push_back(i);
	for (int i=51; i<=58; ++i) subsys_A_.push_back(i);
	for (int i=71; i<=78; ++i) subsys_A_.push_back(i);
	subsys_B_.clear();
	for (int i=119; i<=126; ++i) subsys_B_.push_back(i);
	for (int i=139; i<=146; ++i) subsys_B_.push_back(i);
	for (int i=159; i<=166; ++i) subsys_B_.push_back(i);
	for (int i=179; i<=185; ++i) subsys_B_.push_back(i);
	assert(subsys_A_.size()==subsys_B_.size());
	size_A_ = subsys_A_.size();

	// subsystem: Complement A+B
  subsys_cAB_.clear();
	for (int n=0; n<lattice.num_sites(); ++n) {
		bool in_A_union_B = false;
		for (int i=0; i<size_A_; ++i) {
			if (n==subsys_A_[i]) {
				in_A_union_B = true;
				break;
			}
		}
		for (int i=0; i<size_A_; ++i) {
			if (n==subsys_B_[i]) {
				in_A_union_B = true;
				break;
			}
		}
		if (!in_A_union_B) subsys_cAB_.push_back(n);
	}
	size_cAB_ = subsys_cAB_.size();
	assert(2*size_A_+size_cAB_==lattice.num_sites());


  this->resize(1); // not used, actually
  setup_done_ = true;
  computation_done_ = false;
}

void EntanglementEntropy::reset(void) 
{
  computation_done_ = false;
}

void EntanglementEntropy::compute(const kSpace& kspace, const Hamiltonian& ham)
{
	if (kspace.kmesh().size()!=1) {
    throw std::range_error("error: EntanglementEntropy:compute: expecting kmesh.size=1");
	}

	// construct the hamiltonian (UP-SPIN block)
  Eigen::SelfAdjointEigenSolver<ComplexMatrix> es;
	kvector kvec(0.0,0.0,0.0);
  ham.construct_upspin_block(kvec);
  // diagonalize
  es.compute(ham.upspin_block());

  int N = ham.dimension()/2;
  ComplexVector U(N-1);
  ComplexVector V(N-1);
  ComplexVector U_degen(2);
  ComplexVector V_degen(2);

  // Correlation Matrix A
  ComplexMatrix CorrMatrix_A(size_A_,size_A_);
  for (int i=0; i<size_A_; ++i) {
  	int r_alpha = subsys_A_[i];
  	U = es.eigenvectors().row(r_alpha)(Eigen::seqN(0,N-1));
    U_degen = es.eigenvectors().row(r_alpha)(Eigen::seq(N-1,N));

    //std::cout << i << "\n"; 
    //std::cout << U_degen.size() << "\n"; getchar();

  	for (int j=0; j<size_A_; ++j) {
  		int rprime_beta = subsys_A_[j];
  		V = es.eigenvectors().row(rprime_beta)(Eigen::seqN(0,N-1));
      V_degen = es.eigenvectors().row(rprime_beta)(Eigen::seq(N-1,N));
  		CorrMatrix_A(i,j) = U.dot(V) + 0.5*U_degen.dot(V_degen);
  	}
  }

  // Correlation Matrix B
  ComplexMatrix CorrMatrix_B(size_A_,size_A_);
  for (int i=0; i<size_A_; ++i) {
  	int r_alpha = subsys_B_[i];
  	U = es.eigenvectors().row(r_alpha)(Eigen::seqN(0,N-1));
    U_degen = es.eigenvectors().row(r_alpha)(Eigen::seq(N-1,N));
  	for (int j=0; j<size_A_; ++j) {
  		int rprime_beta = subsys_B_[j];
  		V = es.eigenvectors().row(rprime_beta)(Eigen::seqN(0,N-1));
      V_degen = es.eigenvectors().row(rprime_beta)(Eigen::seq(N-1,N));
  		CorrMatrix_B(i,j) = U.dot(V) + 0.5*U_degen.dot(V_degen);
  	}
  }

  // Correlation Matrix cAB
  ComplexMatrix CorrMatrix_cAB(size_cAB_,size_cAB_);
  for (int i=0; i<size_cAB_; ++i) {
  	int r_alpha = subsys_cAB_[i];
  	U = es.eigenvectors().row(r_alpha)(Eigen::seqN(0,N-1));
    U_degen = es.eigenvectors().row(r_alpha)(Eigen::seq(N-1,N));
  	for (int j=0; j<size_cAB_; ++j) {
  		int rprime_beta = subsys_cAB_[j];
  		V = es.eigenvectors().row(rprime_beta)(Eigen::seqN(0,N-1));
      V_degen = es.eigenvectors().row(rprime_beta)(Eigen::seq(N-1,N));
  		CorrMatrix_cAB(i,j) = U.dot(V) + 0.5*U_degen.dot(V_degen);
  	}
  }

  // Entropy of subsystem-A
  es.compute(CorrMatrix_A, Eigen::EigenvaluesOnly);
  entropy_A_ = 0.0;
  for (int i=0; i<size_A_; ++i) {
  	double lambda = es.eigenvalues()(i);
  	//if (std::abs(lambda)<1.0E-15) lambda=0.0;
  	double lambda1 = 1.0-lambda;
    if (std::abs(lambda)<1.0E-12) {
      entropy_A_ += -lambda1*std::log(lambda1);
    }
    else if (std::abs(lambda1)<1.0E-12) {
      entropy_A_ += -lambda*std::log(lambda);
    }
    else {
      entropy_A_ += -lambda*std::log(lambda)-lambda1*std::log(lambda1);
    }
  }
  std::cout<<"S_A = "<<entropy_A_<<"\n"; 

  // Entropy of subsystem-A
  es.compute(CorrMatrix_B, Eigen::EigenvaluesOnly);
  entropy_B_ = 0.0;
  for (int i=0; i<size_A_; ++i) {
  	double lambda = es.eigenvalues()(i);
  	double lambda1 = 1.0-lambda;
    if (std::abs(lambda)<1.0E-12) {
      entropy_B_ += -lambda1*std::log(lambda1);
    }
    else if (std::abs(lambda1)<1.0E-12) {
      entropy_B_ += -lambda*std::log(lambda);
    }
    else {
      entropy_B_ += -lambda*std::log(lambda)-lambda1*std::log(lambda1);
    }
  }
  std::cout<<"S_B = "<<entropy_B_<<"\n"; 

  // Entropy of subsystem-cAB
  es.compute(CorrMatrix_cAB, Eigen::EigenvaluesOnly);
  entropy_cAB_ = 0.0;
  for (int i=0; i<size_cAB_; ++i) {
  	double lambda = es.eigenvalues()(i);
  	double lambda1 = 1.0-lambda;
  	if (std::abs(lambda)<1.0E-12) {
  		entropy_cAB_ += -lambda1*std::log(lambda1);
  	}
  	else if (std::abs(lambda1)<1.0E-12) {
  		entropy_cAB_ += -lambda*std::log(lambda);
  	}
    else {
      entropy_cAB_ += -lambda*std::log(lambda)-lambda1*std::log(lambda1);
      //std::cout << lambda << " " << lambda1 << "\n";
      //std::cout << -lambda*std::log(lambda) << "\n";
      //std::cout << -lambda1*std::log(lambda1) << "\n";
      //getchar();
    }
  }
  std::cout<<"S_AB = "<<entropy_cAB_<<"\n"; 
  std::cout<<"Invariant = "<<entropy_A_+entropy_B_-entropy_cAB_<<"\n"; 
}

void EntanglementEntropy::print_heading(const std::string& header,
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

void EntanglementEntropy::print_result(const std::vector<double>& xvals) 
{
  // Already done, nothing to be done
}



} // end namespace diag
