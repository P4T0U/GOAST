// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef CONSTRAINEDOPTIMIZATION_H
#define CONSTRAINEDOPTIMIZATION_H

#include <vector>

#include <goast/Optimization/LagrangeInterface.h>

#include "DeformationInterface.h"


//===========================================================================================================================
// CONSTRAINT HANDLER FOR RIGID BODY MOTIONS (RBM)
//===========================================================================================================================

/**
 * \brief Handles rigid body motion constraints for a single shell.
 * \author Heeren
 *
 * There are usually 2 x Dim = 6 constraints, as we have to fix translations (one in each dimension) and rotations (one in each dimension).
 */
template <typename ConfiguratorType>
class RigidBodyMotionsConstraintHandler {
  
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef  typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;
  typedef typename ConfiguratorType::VecType    VecType;
  
  static const int DIM = 3;
  int _numV, _numObjects, _numOfConstraints;
  const VectorType& _identity;
  VectorType _allOnes, _massOpAllOnes;
  std::vector<VectorType> _massOpIdentity;
  MatrixType _massOp, _lumpedMassOp; 

public:
  RigidBodyMotionsConstraintHandler( const MeshTopologySaver& Topology, const VectorType& geometry, int numOfObjects ) :
    _numV( Topology.getNumVertices() ),
    _numObjects( numOfObjects ),
    _numOfConstraints(6),
    _identity( geometry ),
    _allOnes( _numV ),
    _massOpAllOnes( _numV ),
    _massOpIdentity( DIM ),
    _massOp( _numV, _numV ), 
    _lumpedMassOp( _numV, _numV )
    {
      _allOnes = VectorType::Constant( _numV, 1. );
      assembleMassOps( Topology, geometry );
    }
 
  int getNumOfTotalConstraints( ) const { return _numObjects * _numOfConstraints; }
 
  int getNumOfSingleConstraints( ) const { return _numOfConstraints; }
 
  //
  void addConstraintEnergy( const VectorType& geomArg, RealType & Energy, const VectorType& LagrangeMultiplier ) const {
      
    int numDofsSingleObject = DIM*_numV;
    if( geomArg.size() != _numObjects*numDofsSingleObject )
      throw BasicException("RigidBodyMotionsConstraintHandler::addConstraintEnergy: argument has wronf size!");  
      
    if( LagrangeMultiplier.size() != _numObjects * _numOfConstraints )
      throw BasicException("RigidBodyMotionsConstraintHandler::addConstraintEnergy: wrong number of Lagrange multipliers!");      
    
    // run over all objects
    for( int i = 0; i < _numObjects; i++ ){
        
        //std::cerr << "i = " << i << std::endl;
        //std::cerr << "off = " << i*_numOfConstraints << std::endl;
      VectorType evalConstraints( _numOfConstraints );
      Eigen::Ref<const VectorType> geomArgRef = geomArg.segment(i*numDofsSingleObject, numDofsSingleObject );
      Eigen::Ref<const VectorType> multRef    = LagrangeMultiplier.segment(i*_numOfConstraints, _numOfConstraints );
      evaluateConstraints( geomArgRef, evalConstraints );
      // add all constraints
      for( int k = 0; k < _numOfConstraints; k++ ){
        Energy += multRef[k] * evalConstraints[k];
        //std::cerr << "k = " << k << ", lambda_k = " << multRef[k] << ", e_k = " << evalConstraints[k] << std::endl;
      }
      //std::cerr << "----------------------" << std::endl;
    }

  } 
  
  // derivative of the i-th constraint
  void addConstraintGrad( const VectorType& geomArg, VectorType& Grad, const VectorType& LagrangeMultiplier  ) const {
      
    int numDofsSingleObject = DIM*_numV;  
    if( geomArg.size() != _numObjects*numDofsSingleObject )
      throw BasicException("RigidBodyMotionsConstraintHandler::addConstraintGrad: argument has wrong size!");   
      
    if( LagrangeMultiplier.size() != _numObjects * _numOfConstraints )
      throw BasicException("RigidBodyMotionsConstraintHandler::addConstraintGrad: wrong number of Lagrange multipliers!");
      
    if( Grad.size() != _numObjects * (numDofsSingleObject + _numOfConstraints) )
      throw BasicException("RigidBodyMotionsConstraintHandler::addConstraintGrad: gradient has wrong size!");

      
    // run over all objects
    for( int i = 0; i < _numObjects; i++ ){
        //std::cerr << "i = " << i << std::endl;
        
       // bring into more convenient form 
      std::vector<Eigen::Ref<VectorType> > gradRefs;
      gradRefs.reserve(DIM);
      for( int k = 0; k < DIM; k++ ){
          //std::cerr << "k = " << k << "th grad push starts at " << i*numDofsSingleObject + k*_numV << std::endl;
        gradRefs.push_back( Grad.segment(i*numDofsSingleObject + k*_numV, _numV) );
      }
      

      VectorType evalConstraints( _numOfConstraints );
      //std::cerr << "geom push starts at " << i*numDofsSingleObject << std::endl;
      Eigen::Ref<const VectorType> geomArgRef = geomArg.segment(i*numDofsSingleObject, numDofsSingleObject );
      //std::cerr << "mult push starts at " << i*_numOfConstraints << std::endl;
      Eigen::Ref<const VectorType> multRef    = LagrangeMultiplier.segment(i*_numOfConstraints, _numOfConstraints );
      evaluateConstraints( geomArgRef, evalConstraints );
      
      // geometry blocks of gradient
      for( int k = 0; k < DIM; k++ ){
        gradRefs[k] += multRef[k] * _massOpAllOnes;
        gradRefs[k] += multRef[DIM+k] * _massOpIdentity[(k+1)%DIM];
        gradRefs[k] -= multRef[DIM+((k+2)%DIM)] * _massOpIdentity[(k+2)%DIM];
      }
    
      // translation and rotation blocks
      for( int k = 0; k < _numOfConstraints; k++ ){
          //std::cerr << "grad mult part " << _numObjects*numDofsSingleObject + i*_numOfConstraints + k << std::endl;
        Grad[_numObjects*numDofsSingleObject + i*_numOfConstraints + k ] += evalConstraints[k];
      }
      //std::cerr << "----------------------" << std::endl;
    }

  }

  // derivative of the i-th constraint and its position in Hessian
  void addConstraintHessian( const VectorType& /*geomArg*/, TripletListType& tripletList, const VectorType& /*LagrangeMultiplier*/  ) const {
      
    int numDofsSingleObject = DIM*_numV;  
    
    // run over all objects
    for( int i = 0; i < _numObjects; i++ ){
      int colOffset = _numObjects * numDofsSingleObject + i * _numOfConstraints;  
      int rowOffset = i*numDofsSingleObject;
        
      //translations
      for( int j = 0; j < DIM; j++ ){
        for( int k = 0; k < _numV; k++ ){        
          tripletList.push_back( TripletType(  rowOffset + j*_numV + k, colOffset + j, _massOpAllOnes[k] ) );
          // transposed block
          tripletList.push_back( TripletType(  colOffset + j, rowOffset + j*_numV + k, _massOpAllOnes[k] ) );
        }
      }
    
      colOffset += DIM;
    
      //rotations
      for( int j = 0; j < DIM; j++ ){
        for( int k = 0; k < _numV; k++ ){        
          tripletList.push_back( TripletType(  rowOffset + j*_numV + k, colOffset + j, _massOpIdentity[(j+1)%DIM][k] ) );
          // transposed block
          tripletList.push_back( TripletType(  colOffset + j, rowOffset + j*_numV + k, _massOpIdentity[(j+1)%DIM][k] ) );
        
          tripletList.push_back( TripletType(  rowOffset + ((j+1)%DIM)*_numV + k, colOffset + j, -1. * _massOpIdentity[j][k] ) );
          // transposed block
          tripletList.push_back( TripletType(  colOffset + j, rowOffset + ((j+1)%DIM)*_numV + k, -1. * _massOpIdentity[j][k] ) );
        }
      }
    }
  }
  
  //
  void addConstraintHessian( TripletListType& tripletList ) const {
    VectorType geomArgDummy, LagrangeMultiplierDummy;
    addConstraintHessian( geomArgDummy, tripletList, LagrangeMultiplierDummy );
  }
  
protected:
  void assembleMassOps( const MeshTopologySaver& Topology, const VectorType& geometry ) {
           
      // run over all faces and compute area
      VectorType areas( _numV );
      areas.setZero(); 
      VecType Pi, Pj, Pk, temp;
      for( int k = 0; k < Topology.getNumFaces(); k++ ){
        getXYZCoord<VectorType, VecType>( geometry, Pi, Topology.getNodeOfTriangle(k,0) );
        getXYZCoord<VectorType, VecType>( geometry, Pj, Topology.getNodeOfTriangle(k,1) );
        getXYZCoord<VectorType, VecType>( geometry, Pk, Topology.getNodeOfTriangle(k,2) );
        temp.makeCrossProduct( Pk-Pj, Pi-Pk );
        RealType area = temp.norm() / 2.;
        for( int j = 0; j < 3; j++ )
          areas[Topology.getNodeOfTriangle(k,j)] += area;
      }
      
      // run over all vertices and assign areas
      TripletListType tripletList;
      for( int k = 0; k < _numV; k++ )
          tripletList.push_back( TripletType( k, k, areas[k] ) );
      _lumpedMassOp.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
          
          
      _massOpAllOnes = _lumpedMassOp * _allOnes;
      for( int k = 0; k < DIM; k++ )
          _massOpIdentity[k] = _lumpedMassOp * geometry.segment(k*_numV, _numV);
  }
  
  //
  void evaluateConstraints( const VectorType& geomArg, VectorType& Evals ) const {
    Evals.resize( _numOfConstraints );    
       
    // bring into more convenient form
    VectorType displacement( geomArg ); 
    std::vector< Eigen::Ref<VectorType> > dispRefs;
    dispRefs.reserve(DIM);    
    for( int k = 0; k < DIM; k++ ){
      dispRefs.push_back( displacement.segment(k*_numV, _numV) );
      dispRefs[k] -= _identity.segment(k*_numV, _numV);
    }
    
    // translation: F_i[\lambda_i,X] = M (X^i - R^i) * 1, where Arg = (X^1, X^2, X^3) and Ref = (R^1, R^2, R^3) is reference Object
    // Furthermore: X^i, R^i, 1 \in \R^n, and M is the nxn - mass matrix.
    for( int k = 0; k < DIM; k++ )
      Evals[k] =  _massOpAllOnes.dot( dispRefs[k] );
    
    // rotation
    for( int k = 0; k < DIM; k++ )
       Evals[DIM+k] = _massOpIdentity[ (k + 1)%DIM ].dot( dispRefs[k] ) - _massOpIdentity[k].dot( dispRefs[ (k + 1)%DIM ] );
 
  }

};


#endif