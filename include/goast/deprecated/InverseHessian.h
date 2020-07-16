// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
//!==========================================================================================================
//! \brief Operator that maps x to D2F(x_0)^{-1}x where D2F is the Hessian matrix of an energy functional F
//! \author Heeren
//! This class can be usefull e.g. for a warm-start BFGS method.
template<typename ConfiguratorType, typename HessianType>
class InverseHessianOperator : public BaseOp<typename ConfiguratorType::VectorType> {
    
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  
  const MeshTopologySaver& _topology;
  const HessianType& _D2E;
  const VectorType& _point;
  const VectorType* _refShape;
  const std::vector<int>* _mask;
  LinearSolver<ConfiguratorType> _solver;
  int _numDofs;
  
public:
  InverseHessianOperator( const MeshTopologySaver& Topology, const HessianType& D2E, const VectorType& Point ) 
    : _topology(Topology), _D2E(D2E), _point(Point), _refShape(NULL), _mask(NULL), _numDofs(Point.size()){ }
  
  InverseHessianOperator( const MeshTopologySaver& Topology, const HessianType& D2E, const VectorType& Point, const VectorType& refShape ) 
    : _topology(Topology), _D2E(D2E), _point(Point), _refShape(&refShape), _mask(NULL), _numDofs(Point.size()){ updateSolver(Point); }  
  
  InverseHessianOperator( const MeshTopologySaver& Topology, const HessianType& D2E, const VectorType& Point, const std::vector<int>& Mask ) 
    : _topology(Topology), _D2E(D2E), _point(Point), _refShape(NULL), _mask(Mask), _numDofs(Point.size()){ updateSolver(Point); }
  
  void apply( const VectorType& Arg, VectorType& Dest ) const {      
      
      if( Arg.size() != _numDofs )
          throw BasicException("InverseHessianOperator::apply(): argument has wrong size!");
      
      if( _mask )
        _solver.backSubstitute( Arg, Dest );
      
      if( _refShape ){
          VectorType augmArg( Arg );          
          int numConstraints = 6;
          int numConstrainedShapes = _numDofs / (3 * _topology.getNumVertices());
          augmArg.conservativeResize( _numDofs + numConstrainedShapes * numConstraints  );
          _solver.backSubstitute( augmArg, Dest );
          Dest.conservativeResize( _numDofs );
      }
  }

  void updateSolver( const VectorType& Arg ) {          
      
      MatrixType Hessian;
      
      if( _mask ){
        _D2E.apply( Arg, Hessian );
        applyMaskToSymmetricMatrix( *_mask, Hessian );
      }
      
      if( _refShape ){
          int numConstrainedShapes = _numDofs / (3 * _topology.getNumVertices());
          int numConstraints = 6;
          RigidBodyMotionsConstraintHandler<ConfiguratorType> constHandler( _topology, *_refShape, numConstrainedShapes );          
          LagrangeHessian<ConfiguratorType, HessianType, RigidBodyMotionsConstraintHandler<ConfiguratorType> > D2L(  _D2E, constHandler, _numDofs );
          VectorType augmArg( Arg );
          augmArg.conservativeResize( _numDofs + numConstrainedShapes * numConstraints );
          D2L.apply( augmArg, Hessian );
      }
      
      _solver.prepareSolver( Hessian );      
  }
  
  void setBoundaryMask( const std::vector<int>& Mask ) {
    _mask = &Mask;    
    updateSolver(_point); 
  }
    
};

