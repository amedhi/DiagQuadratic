/*---------------------------------------------------------------------------
* Copyright (C) 2025 by Amal Medhi <amedhi@iisertvm.ac.in>
* All rights reserved.
* Date:   2026-04-06 10:53:30
* Last Modified by:   Amal Medhi
* Last Modified time: 2026-04-06 13:13:30
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
  ComplexVector U(N);
  ComplexVector V(N);

  // --------------S_A-------------------
  // Correlation Matrix A
  CorrMatrix_.resize(size_A_,size_A_);
  for (int i=0; i<size_A_; ++i) {
  	int r_alpha = subsys_A_[i];
  	U = es.eigenvectors().row(r_alpha)(Eigen::seqN(0,N));
  	for (int j=0; j<size_A_; ++j) {
  		int rprime_beta = subsys_A_[j];
  		V = es.eigenvectors().row(rprime_beta)(Eigen::seqN(0,N));
  		CorrMatrix_(i,j) = U.dot(V);
  	}
  }
  // Eigenvalues of CorrMatrix_A
  es.compute(CorrMatrix_, Eigen::EigenvaluesOnly);
  for (int i=0; i<N; ++i) {
  	std::cout << es.eigenvalues()(i) << "\n"; 
  }
  std::cout << "\n"; 


  // --------------S_B-------------------
  // Correlation Matrix B
  for (int i=0; i<size_A_; ++i) {
  	int r_alpha = subsys_B_[i];
  	U = es.eigenvectors().row(r_alpha)(Eigen::seqN(0,N));
  	for (int j=0; j<size_A_; ++j) {
  		int rprime_beta = subsys_B_[j];
  		V = es.eigenvectors().row(rprime_beta)(Eigen::seqN(0,N));
  		CorrMatrix_(i,j) = U.dot(V);
  	}
  }

  // --------------S_cAB-------------------

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
