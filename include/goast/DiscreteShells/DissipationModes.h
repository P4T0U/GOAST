// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Dissipation modes as eigenmodes of the Hessian of an elastic shell energy.
 * \author Heeren
 *
 */
 #ifndef __DISSIPATIONMODES_H
#define __DISSIPATIONMODES_H

#include <goast/Core.h>

//===============================================================================================================================
//===============================================================================================================================

//! \brief Abstract base class for matrix operator
//! \author Heeren
template<typename ConfiguratorType>
class AbstractMatrixOperator : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType>{    
    
   typedef typename ConfiguratorType::VectorType VectorType;       
   typedef typename ConfiguratorType::SparseMatrixType MatrixType;
   typedef typename ConfiguratorType::TripletType TripletType;
   typedef std::vector<TripletType> TripletListType;
  
public:
  AbstractMatrixOperator() { }
  virtual ~AbstractMatrixOperator(){}    
  
  virtual void apply( const VectorType& Arg, MatrixType& Matrix ) const {

    Matrix.resize( Arg.size(), Arg.size() );
    TripletListType tripletList;
    pushTriplets( Arg, tripletList );
    Matrix.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
  
  virtual void pushTriplets( const VectorType&, TripletListType& ) const = 0;
  
  virtual int getDimOfKernel() const = 0;
};
  
//! \brief Matrix operator representing the Hessian of an elastic shell energy
//! \author Heeren
template<typename ConfiguratorType>
class ShellHessianOperator : public AbstractMatrixOperator<ConfiguratorType> {    

   typedef typename ConfiguratorType::VectorType VectorType;       
   typedef typename ConfiguratorType::TripletType TripletType;
   typedef std::vector<TripletType> TripletListType;

  const DeformationBase<ConfiguratorType>& _W;

public:
  ShellHessianOperator( const DeformationBase<ConfiguratorType>& W ) : _W(W) {}
    
  void pushTriplets( const VectorType& Arg, TripletListType& triplets ) const {
      _W.pushTripletsDefHessian ( Arg, Arg, triplets, 0, 0, 1.0 );
  }
    
  int getDimOfKernel() const {
      return 6;
  }
};

  
//===============================================================================================================================
//===============================================================================================================================

/**
 * \brief Abstract base class for a metric
 * \author Heeren
 *
 * The metric is given by \f$ g(v,w) := v^TMw \f$ for some matrix representation \f$M\f$.
 * The pure virtual apply function maps \f$v\f$ to \f$Mv\f$ and has to be provided in all derived classes.
 */
template<typename ConfiguratorType>
class AbstractMetricClass {    
   typedef typename ConfiguratorType::RealType RealType;
   typedef typename ConfiguratorType::VectorType VectorType; 
   
public:
  AbstractMetricClass( ){ }                       
  virtual ~AbstractMetricClass(){}  
  
  virtual void apply( const VectorType&, VectorType& ) const {
    throw BasicException ( "AbstractMetricClass(): abstract base class!" );
  }    
  
  virtual RealType evaluate( const VectorType& V, const VectorType& W ) const {
    VectorType temp;
    apply( V, temp );
    return temp.dot(W);
  }  

  virtual void normalize( VectorType& W ) const {
    W /= std::sqrt( evaluate( W, W ) );
  } 
};
  
//! \brief Euclidean metric, i.e. M = Id
//! \author Heeren
template<typename ConfiguratorType>
class EuclideanMetric : public AbstractMetricClass<ConfiguratorType> {    
  typedef typename ConfiguratorType::VectorType VectorType;   
public:    
  void apply( const VectorType& V, VectorType& W ) const { W = V; }
};
  
//! \brief L2 metric, i.e. M is (lumped) mass matrix
//! \author Heeren
template<typename ConfiguratorType>
class L2Metric : public AbstractMetricClass<ConfiguratorType> {    
    
   typedef typename ConfiguratorType::RealType RealType;
   typedef typename ConfiguratorType::VectorType VectorType;
   typedef typename ConfiguratorType::SparseMatrixType MatrixType;
   
   MatrixType _lumpedMassMatrix;
   std::vector<int> _bdryMask;
   
public:
  L2Metric( const MeshTopologySaver& Topology, const VectorType& Geometry ) {
      computeLumpedMassMatrix<ConfiguratorType>(Topology, Geometry, _lumpedMassMatrix);
  }

  void apply( const VectorType& V, VectorType& W ) const {
    W = _lumpedMassMatrix * V;
    // take care of Dirichlet nodes!?
    if( _bdryMask.size() > 0 )
      applyMaskToVector( _bdryMask, W );
  }
};

//===============================================================================================================================
//===============================================================================================================================

/**
 * \brief Operator to compute eigenmodes wrt. a particular metric
 * \author Heeren
 *
 * For a given metric \f$ g(v,w) = v^TMw \f$ and some matrix (operator) \f$ A \f$,
 * compute sequence \f$ (v_i, \lambda_i)_{i \geq 0}\f$ such that \f$ A v_i = \lambda_i M v_i \f$.
 *
 *  Here we have (approximately) \f$ 0 \leq \lambda_\min = \lambda_0 \leq \lambda_1 \leq \ldots \f$,
 *  i.e. we compute the smallest eigenvalues of \f$ A \f$ by means of an inverse power iteration.
 */
template<typename ConfiguratorType>
class EigenmodesOp {

protected:
   typedef typename ConfiguratorType::RealType RealType;
   typedef typename ConfiguratorType::VectorType VectorType;
   typedef typename ConfiguratorType::SparseMatrixType MatrixType;

  const MeshTopologySaver& _topology;
  const AbstractMatrixOperator<ConfiguratorType>& _matrixOp;
  const AbstractMetricClass<ConfiguratorType>& _metric;
  std::vector<int> _bdryMask;
  bool _quiet;
  int _maxIterations;

public:
  EigenmodesOp( const MeshTopologySaver& Topology, 
                const AbstractMatrixOperator<ConfiguratorType>& MatrixOp, 
                const AbstractMetricClass<ConfiguratorType>& Metric,
                bool quiet = true ) :
      _topology( Topology ),
      _matrixOp( MatrixOp ),
      _metric( Metric ),    
      _quiet( quiet ),           
      _maxIterations(1000){ }
  
  void setBoundaryMask( const std::vector<int>& Mask ) {
    _bdryMask.resize( Mask.size() );
    _bdryMask = Mask;    
  }

  //
  void execute( const VectorType& Geometry, int numEigenvalues, VectorType& Eigenvalues, std::vector<VectorType>& Eigenvectors ) const {
    
    // initialization of eigenvalues and eigenvectors
    Eigenvectors.resize( numEigenvalues );
    Eigenvalues.resize( numEigenvalues );
    
    int dim = 3 * _topology.getNumVertices();
    if( dim != Geometry.size() )
        throw BasicException("EigenmodesOp::execute(): geometry has wrong dimsension!");
    
    int numConstraints = _matrixOp.getDimOfKernel();
    int totalNumDofs = _bdryMask.size() > 0 ? dim : dim + numConstraints;
    MatrixType matrix( totalNumDofs, totalNumDofs );     
    
    if( _bdryMask.size() == 0 ){
        if(!_quiet) std::cerr << "Use Lagrange setup to account for rigid body motions..." << std::endl;
        if(!_quiet) std::cerr << "Number of constraints is = " << numConstraints << std::endl;
        RigidBodyMotionsConstraintHandler<ConfiguratorType> constHandler( _topology, Geometry, 1 );
        LagrangeHessian<ConfiguratorType, AbstractMatrixOperator<ConfiguratorType>, RigidBodyMotionsConstraintHandler<ConfiguratorType> > LagrangeHessian( _matrixOp, constHandler, dim );
        VectorType Arg( Geometry );
        Arg.conservativeResize( totalNumDofs );   
        LagrangeHessian.apply( Arg, matrix );
    }
    else{
      if(!_quiet) std::cerr << "Consider bdry mask with " << _bdryMask.size()/3 << " boundary nodes."  << std::endl;
      _matrixOp.apply( Geometry, matrix );
      applyMaskToSymmetricMatrix( _bdryMask, matrix );      
    }

    // calculate eigenvalues and modes
    if(!_quiet) std::cerr<<"================================="<< std::endl;
    if(!_quiet) std::cerr<<"Start inverse vector iteration. "<< std::endl;
    if(!_quiet) std::cerr<<"================================="<< std::endl;
    
    auto t_start = std::chrono::high_resolution_clock::now();  
    performInversePowerIteration( matrix, Eigenvectors, Eigenvalues );
    auto t_end = std::chrono::high_resolution_clock::now();
    std::cerr << std::fixed << "Inverse vector iteration done in " << std::chrono::duration<double, std::ratio<1> >(t_end - t_start).count() << " seconds." << std::endl;   
      
    // orthonormality check
    if(!_quiet) checkOrthonormality( Eigenvectors ); 
  }
  
  // Perform inverse power iteration to get smallest (in absolute value) eigenvalues a_i and eigenvectors b_i of some matrix A,
  // i.e. we calculate iteratively (A - \mu Id) b = b^k and b^{k+1} := b / |b| with \mu = 0 (as we are looking for the smallest Eigenvalue)
  void performInversePowerIteration( const MatrixType& Matrix, std::vector< VectorType >& eigenvectors, VectorType& eigenvalues ) const {
                                         
        int numEigenvalues = eigenvalues.size();
        int dim = 3 * _topology.getNumVertices();
        int numConstraints = _bdryMask.size() > 0 ? 0 : _matrixOp.getDimOfKernel();
	
        // initialize solver   
        if(!_quiet) std::cerr << "Factorize matrix..." << std::endl;
        LinearSolver<ConfiguratorType> linSolver;
        auto t_start = std::chrono::high_resolution_clock::now(); 
        linSolver.prepareSolver( Matrix );
        auto t_end = std::chrono::high_resolution_clock::now();
        std::cerr << std::fixed << "done in " << std::chrono::duration<double, std::ratio<1> >(t_end - t_start).count() << " seconds." << std::endl;  

        // compute eigenvalues via inverse power iteration
        for ( int i = 0; i < numEigenvalues; ++i ){
	  
	    if(!_quiet) std::cerr << "Start to compute " << i+1 << "th of " << numEigenvalues << " eigenvectors..." << std::endl;

	    if( eigenvectors[i].size() == 0 ){
              eigenvectors[i] = VectorType::Constant( dim, 1. );
	      _metric.normalize( eigenvectors[i] );
	    }

            RealType err = 1.;
	    int iter = 0;
            eigenvalues[i] = (i>0) ? eigenvalues[i-1] : 0;            
            VectorType rhs;
            
	    // start inverse vector iteration
            for (; iter < _maxIterations && err > 1e-10; ++iter){	        
            		
                // rhs = Mx^k, where M is matrix that represents metric
                _metric.apply( eigenvectors[i], rhs );
                rhs.conservativeResize( dim + numConstraints );
                if( _bdryMask.size() > 0. )
                    applyMaskToVector( _bdryMask, rhs );

                // solve A x^{k+1} = Mx^k
                VectorType newEV;
                linSolver.backSubstitute( rhs, newEV );
                newEV.conservativeResize( dim );

                // projection
                VectorType factors(i);
                for (int j = 0; j < i; ++j)
                    factors[j] = _metric.evaluate( newEV, eigenvectors[j] );
                for (int j = 0; j < i; ++j)
                  newEV -= factors[j] * eigenvectors[j];
                
                // approximation by Rayleigh coefficients             
                RealType newEigenValue = 1. / _metric.evaluate( newEV, eigenvectors[i] );                 
                                
                // normalize
		_metric.normalize( newEV );
		
		// update error
                err = std::abs( (newEigenValue - eigenvalues[i]) / eigenvalues[i] );
                
                // update
                eigenvalues[i] = newEigenValue;
                eigenvectors[i] = newEV;
            }
            
            // final console output
            if(!_quiet) std::cerr << iter << " iterations: lambda[" << i << "] = " << eigenvalues[i] << ", error = " << err << std::endl;
	    if(!_quiet) std::cerr << "-----------------------------------------------------------" << std::endl;
        }
  }
   
  //
  void checkOrthonormality( const std::vector<VectorType>& eigenvectors ) const {
    std::cerr << std::endl << "Othonormality check." << std::endl;
    for( int i = 0; i < eigenvectors.size(); i++ ){
      for( int j = 0; j < eigenvectors.size(); j++){
          if( j == i )
            std::cerr << _metric.evaluate( eigenvectors[i], eigenvectors[j] ) << " ";
          else
            std::cerr << (_metric.evaluate( eigenvectors[i], eigenvectors[j]) > 1e-10 ? "x " : "0 ");
      }
      std::cerr << std::endl;
    }
    std::cerr << std::endl;
  }

};


#endif