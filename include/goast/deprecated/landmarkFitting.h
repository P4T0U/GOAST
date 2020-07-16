// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by sassen on 22.06.18.
//

#ifndef NRIC_LANDMARKFITTING_H
#define NRIC_LANDMARKFITTING_H

#include "Subspace.h"

//! \brief Mismatch of discrete shell with given landmarks considering distance to LΘ-subspace as prior
//! \author Sassen
template<typename ConfiguratorType>
class LThetaSubspaceFittingFunctional
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;

  const MeshTopologySaver &_Topology;
  const std::vector<int> &_markerIndices;
  const std::vector<VecType> &_markerPositions;
  const RealType _gamma;

  const int _numLThetaDofs, _numLocalDofs, _numInputData;

  const SubspaceDistanceFunctional<ConfiguratorType> _distF;


public:
  LThetaSubspaceFittingFunctional( const MeshTopologySaver &Topology,
                                   const std::vector<int> &markerIndices,
                                   const std::vector<VecType> &markerPositions,
                                   const VectorType &refMeshLengthsAngles,
                                   const VectorType &IntegrationWeights,
                                   const RealType edgeLengthWeight,
                                   const RealType bendingWeight,
                                   const FullMatrixType &subspaceBasis,
                                   RealType priorWeight = 10. )
          : _Topology( Topology ),
            _markerIndices( markerIndices ), _markerPositions( markerPositions ),
            _gamma( priorWeight / 2. ),
            _numLThetaDofs( refMeshLengthsAngles.size()),
            _numInputData( subspaceBasis.cols()), _numLocalDofs( 3 * Topology.getNumVertices()),
            _distF( IntegrationWeights, refMeshLengthsAngles, subspaceBasis )
  {
    if ( markerIndices.size() != markerPositions.size())
      throw BasicException( "LThetaSubspaceFittingFunctional: number of marker indices and positions have to agree!" );


  }

  void apply( const VectorType &Arg, RealType &Dest ) const override {
    if ( Arg.size() != _numLocalDofs )
      throw BasicException( "LThetaLandmarkFittingFunctional::apply: wrong size of dofs!" );

    Dest.setZero();

    // Energy part
    VectorType lengthsAngles;
    computeLengthsAngles<ConfiguratorType>(_Topology, Arg, lengthsAngles);

    _distF.apply( lengthsAngles, Dest );

//    std::cerr << std::endl << "LSFF >> Prior contribution: " << _gamma * Dest[0] << " (" << Dest[0] << ")" << std::endl;

    Dest *= _gamma;

    // Model fitting part
    VecType vertexPosition;
    RealType differenceNorm;
    RealType fitContrib = 0;
//#ifdef _OPENMP
//#pragma omp parallel for private(vertexPosition, differenceNorm)
//#endif
    for ( int m = 0; m < _markerIndices.size(); m++ ) {
      getXYZCoord( Arg.segment( 0, _numLocalDofs ), vertexPosition, _markerIndices[m] );

      differenceNorm = (vertexPosition - _markerPositions[m]).squaredNorm();

#pragma omp critical
      {
        fitContrib += differenceNorm;
      }
    }

    Dest[0] += fitContrib;
//    std::cerr << "LSFF >> Fitting contribution: " << fitContrib << std::endl;
  }
};

//! \brief Derivative of LandmarkFittingFunctional w.r.t. s
//! \author Sassen
template<typename ConfiguratorType>
class LThetaSubspaceFittingGradient
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;

  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver &_Topology;
  const std::vector<int> &_markerIndices;
  const std::vector<VecType> &_markerPositions;
  const RealType _gamma;

  const int _numLThetaDofs, _numLocalDofs, _numInputData, _numVertices;

  const SubspaceDistanceGradient<ConfiguratorType> _distDF;
public:
  LThetaSubspaceFittingGradient( const MeshTopologySaver &Topology,
                                 const std::vector<int> &markerIndices,
                                 const std::vector<VecType> &markerPositions,
                                 const VectorType &refMeshLengthsAngles,
                                 const VectorType &IntegrationWeights,
                                 const RealType edgeLengthWeight,
                                 const RealType bendingWeight,
                                 const FullMatrixType &subspaceBasis,
                                 RealType priorWeight = 10. )
          : _Topology( Topology ),
            _markerIndices( markerIndices ), _markerPositions( markerPositions ),
            _gamma( priorWeight / 2. ),
            _numLThetaDofs( refMeshLengthsAngles.size()),
            _numInputData( subspaceBasis.cols()), _numLocalDofs( 3 * Topology.getNumVertices()),
            _numVertices( Topology.getNumVertices() ),
            _distDF( IntegrationWeights, refMeshLengthsAngles, subspaceBasis )
  {
    if ( markerIndices.size() != markerPositions.size())
      throw BasicException( "LThetaSubspaceFittingGradient: number of marker indices and positions have to agree!" );

  }

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    if ( Arg.size() != _numLocalDofs )
      throw BasicException( "LThetaSubspaceFittingGradient::apply: wrong size of dofs!" );

    Dest.resize(Arg.size());
    Dest.setZero();

    // Energy contribution
    VectorType lengthsAngles;
    computeLengthsAngles<ConfiguratorType>(_Topology, Arg, lengthsAngles);

    MatrixType LThetaDeriv;
    LThetaDeriv.resize( 3 * _Topology.getNumVertices(), 2 * _Topology.getNumEdges());
    LThetaDeriv.setZero();

    // set up triplet list
    TripletListType tripletList;
    pushLThetaTriplets( Arg, tripletList, true );
    // fill matrix from triplets
    LThetaDeriv.setFromTriplets( tripletList.cbegin(), tripletList.cend());

//    std::cerr << std::endl << "LSFG >> LThetaDeriv: " << LThetaDeriv.norm() << std::endl;

    VectorType distanceGrad;
    _distDF.apply( lengthsAngles, distanceGrad );
    Dest = LThetaDeriv * distanceGrad;

//    std::cerr << std::endl << "LSFG >> Prior contribution: " << (_gamma * Dest).norm() << " (" << _gamma << " * " <<  Dest.norm() << "); dG= " << distanceGrad.norm() << std::endl;

    Dest *=  _gamma;

    // fitting contribution
    VecType vertexPosition;
    RealType differenceNorm;
    VectorType fitContrib = VectorType::Zero(Arg.size());
//#ifdef _OPENMP
//#pragma omp parallel for private(vertexPosition, differenceNorm)
//#endif
    for ( int m = 0; m < _markerIndices.size(); m++ ) {
      getXYZCoord( Arg.segment(0, _numLocalDofs), vertexPosition, _markerIndices[m] );
      vertexPosition -=  _markerPositions[m];

#pragma omp critical
      {
        fitContrib[_markerIndices[m]] += 2 * vertexPosition[0];
        fitContrib[_numVertices + _markerIndices[m]] += 2 * vertexPosition[1];
        fitContrib[2 * _numVertices + _markerIndices[m]] += 2 * vertexPosition[2];
      }
    }

//    std::cerr << "LSFG >> Fitting contribution: " << fitContrib.norm() << std::endl;

//    std::cerr << "LSFG >> Angle: " << std::acos(fitContrib.dot(Dest) / (fitContrib.norm() * Dest.norm())) << std::endl;

    Dest += fitContrib;

//    std::cerr << "LSFG >> Whole norm: " << Dest.norm() << std::endl;
  }

  void pushLThetaTriplets( const VectorType &Arg, TripletListType &tripletList, bool transposed = false ) const {

    tripletList.clear();
    tripletList.reserve( 3 * (2 + 4) * _Topology.getNumEdges());
    int colOffset = _Topology.getNumVertices();
    int rowOffset = _Topology.getNumEdges();

    double maxEdge = 0, maxTheta = 0;
    int absurdTheta = 0;

    // run over all edges
    for ( int edgeIdx = 0; edgeIdx < _Topology.getNumEdges(); ++edgeIdx ) {

      int pi( _Topology.getAdjacentNodeOfEdge( edgeIdx, 0 )),
              pj( _Topology.getAdjacentNodeOfEdge( edgeIdx, 1 )),
              pk( _Topology.getOppositeNodeOfEdge( edgeIdx, 0 )),
              pl( _Topology.getOppositeNodeOfEdge( edgeIdx, 1 ));

      // set up vertices
      VecType Pi, Pj;
      getXYZCoord<VectorType, VecType>( Arg, Pi, pi );
      getXYZCoord<VectorType, VecType>( Arg, Pj, pj );
      VecType edge( Pi - Pj );
      RealType edgeLength = edge.norm();

      // assemble in global matrix
      for ( int i = 0; i < 3; i++ ) {
        pushTriplet( tripletList, edgeIdx, i * colOffset + pi,
                      edge[i] / edgeLength, transposed );
        if (std::abs(edge[i] / edgeLength) > maxEdge) maxEdge = std::abs(edge[i] / edgeLength);

        pushTriplet( tripletList, edgeIdx, i * colOffset + pj,
                     -1. * edge[i] / edgeLength, transposed );
      }

      // no bending at boundary edges
      if ( std::min( pl, pk ) < 0 )
        continue;

      // add bending contribution at all?
      VecType Pk, Pl;
      getXYZCoord<VectorType, VecType>( Arg, Pk, pk );
      getXYZCoord<VectorType, VecType>( Arg, Pl, pl );

      // compute first derivatives of dihedral angle
      VecType thetak, thetal, thetai, thetaj;
      getThetaGradK( Pi, Pj, Pk, thetak );
      getThetaGradK( Pj, Pi, Pl, thetal );
      getThetaGradI( Pi, Pj, Pk, Pl, thetai );
      getThetaGradJ( Pi, Pj, Pk, Pl, thetaj );

      if (thetai.norm() > maxTheta) maxTheta = std::abs(thetai.norm());
      if (thetaj.norm() > maxTheta) maxTheta = std::abs(thetaj.norm());
      if (thetak.norm() > maxTheta) maxTheta = std::abs(thetak.norm());
      if (thetal.norm() > maxTheta) maxTheta = std::abs(thetal.norm());

      if (thetai.norm() > 100000 || thetaj.norm() > 100000 || thetak.norm() > 100000 || thetal.norm() > 100000) {
        std::cout << "absurd: " << edgeIdx << " (" << pi << "," << pj << "," << pk << "," << pl << ")" << std::endl;
        absurdTheta++;
//        continue;
      }


      // assemble in global matrix
      for ( int i = 0; i < 3; i++ ) {
        pushTriplet( tripletList, rowOffset + edgeIdx, i * colOffset + pi, thetai[i], transposed );
        pushTriplet( tripletList, rowOffset + edgeIdx, i * colOffset + pj, thetaj[i], transposed );
        pushTriplet( tripletList, rowOffset + edgeIdx, i * colOffset + pk, thetak[i], transposed );
        pushTriplet( tripletList, rowOffset + edgeIdx, i * colOffset + pl, thetal[i], transposed );
      }

    }
//    std::cout << "pushLThetaTriplets: " << maxEdge << " / " << maxTheta << " / " << absurdTheta << std::endl;

  }

protected:
  void pushTriplet( TripletListType &tripletList, int row, int col, RealType value, bool transposed ) const {
    if ( transposed )
      tripletList.push_back( TripletType( col, row, value ));
    else
      tripletList.push_back( TripletType( row, col, value ));
  }
};

template<typename ConfiguratorType>
class LThetaSubspaceFittingHessian
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver &_Topology;
  const std::vector<int> &_markerIndices;
  const std::vector<VecType> &_markerPositions;
  const RealType _gamma;

  const int _numLThetaDofs, _numLocalDofs, _numInputData, _numVertices;

  const SubspaceDistanceHessian<ConfiguratorType> _distD2F;
  const SubspaceDistanceGradient<ConfiguratorType> _distDF;


public:
  LThetaSubspaceFittingHessian( const MeshTopologySaver &Topology,
                                const std::vector<int> &markerIndices,
                                const std::vector<VecType> &markerPositions,
                                const VectorType &refMeshLengthsAngles,
                                const VectorType &IntegrationWeights,
                                const RealType edgeLengthWeight,
                                const RealType bendingWeight,
                                const FullMatrixType &subspaceBasis,
                                RealType priorWeight = 10. )
          : _Topology( Topology ),
            _markerIndices( markerIndices ), _markerPositions( markerPositions ),
            _gamma( priorWeight ),
            _numLThetaDofs( refMeshLengthsAngles.size()),
            _numInputData( subspaceBasis.cols()), _numLocalDofs( 3 * Topology.getNumVertices()),
            _numVertices( Topology.getNumVertices() ),
            _distD2F( IntegrationWeights, refMeshLengthsAngles, subspaceBasis ),
            _distDF( IntegrationWeights, refMeshLengthsAngles, subspaceBasis )
  {
    if ( markerIndices.size() != markerPositions.size())
      throw BasicException( "LThetaSubspaceFittingHessian: number of marker indices and positions have to agree!" );
  }

  void apply( const VectorType &Arg, MatrixType &Dest ) const override {
    if ( Arg.size() != _numLocalDofs )
      throw BasicException( "LThetaSubspaceFittingHessian::apply: wrong size of dofs!" );

    if ((Dest.cols() != Arg.size()) || (Dest.rows() != Arg.size()))
      Dest.resize( Arg.size(), Arg.size());

    Dest.setZero();

    // fill triplet lists for marker fitting hessian
    TripletListType tripletList;
    tripletList.reserve(_markerPositions.size() * 3);
    pushFittingTriplets( Arg, tripletList );
    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend());

    // Energy contribution
    MatrixType temp;
    temp.resize( Arg.size(), Arg.size());
    _distD2F.apply( Arg, temp );

    Dest += _gamma * temp;
  }

  void pushFittingTriplets( const VectorType &Arg, TripletListType &tripletList ) const {
    if ( Arg.size() != _numLocalDofs + _numInputData )
      throw BasicException( "LThetaLandmarkFittingHessian::pushFittingTriplets: wrong size of dofs!" );

    VecType vertexPosition;
    RealType differenceNorm;
//#ifdef _OPENMP
//#pragma omp parallel for private(vertexPosition, differenceNorm)
//#endif
    for ( int m = 0; m < _markerIndices.size(); m++ ) {
      TripletListType localTripletList;
      localTripletList.reserve( 3 );

      localTripletList.push_back(TripletType(_markerIndices[m], _markerIndices[m], 2));
      localTripletList.push_back(TripletType(_numVertices + _markerIndices[m], _numVertices + _markerIndices[m], 2));
      localTripletList.push_back(TripletType(2 * _numVertices + _markerIndices[m], 2 * _numVertices + _markerIndices[m], 2));

#ifdef _OPENMP
#pragma omp critical
#endif
      tripletList.insert( tripletList.end(), localTripletList.begin(), localTripletList.end());
    }

  }

};

//! \brief Mismatch of discrete shell with given landmarks considering distance to a point in an LΘ-subspace as prior
//! \author Sassen
//!
//! Given landmark positions x_0,...,x_l, corresponding vertex positions X_0(s),..,X_l(s) and weight \gamma this
//! implements
//! F_x[s] = \sum_{l=1}^L ||X_l(s) - x_l||^2 + \gamma \sum_{e \in E} (\alpha_e (l_e[s] - (\bar{l}_e + \sum_j \omega_j u_e^j) )^2 + \eta\, \beta_e (\theta_e[s] - (\bar{\theta}_e + \sum_j \omega_j u_\theta^j) )^2)
//! for a shell s.
template<typename ConfiguratorType>
class LThetaLandmarkFittingFunctional
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;

  const DeformationBase<ConfiguratorType> &_W;
  const MeshTopologySaver &_Topology;
  const std::vector<int> &_markerIndices;
  const std::vector<VecType> &_markerPositions;
  const RealType _gamma;

  const int _numLThetaDofs, _numLocalDofs, _numInputData;

  LinearBlendingResidual<ConfiguratorType> _LThetaResidualEnergy;
  SquaredFunctional<ConfiguratorType> _squaredResidualEnergy;

public:
  LThetaLandmarkFittingFunctional( const MeshTopologySaver &Topology,
                                   const DeformationBase<ConfiguratorType> &W,
                                   const std::vector<int> &markerIndices,
                                   const std::vector<VecType> &markerPositions,
                                   const VectorType &refMeshLengthsAngles,
                                   const VectorType &IntegrationWeights,
                                   const RealType edgeLengthWeight,
                                   const RealType bendingWeight,
                                   const std::vector<VectorType> &inputDataLengthsAngles,
                                   RealType priorWeight = 10. )
          : _Topology( Topology ), _W( W ), _markerIndices( markerIndices ), _markerPositions( markerPositions ),
            _gamma( priorWeight ),
            _numLThetaDofs( refMeshLengthsAngles.size()),
            _numInputData( inputDataLengthsAngles.size()), _numLocalDofs( 3 * Topology.getNumVertices()),
            _LThetaResidualEnergy( Topology, refMeshLengthsAngles, IntegrationWeights, edgeLengthWeight, bendingWeight,
                                   inputDataLengthsAngles ),
            _squaredResidualEnergy( _LThetaResidualEnergy ) {
    if ( markerIndices.size() != markerPositions.size())
      throw std::length_error( "LThetaLandmarkFittingFunctional: number of marker indices and positions have to "
                               "agree!" );


  }

  void apply( const VectorType &Arg, RealType &Dest ) const override {
    if ( Arg.size() != _numInputData + _numLocalDofs )
      throw std::length_error( "LThetaLandmarkFittingFunctional::apply: wrong size of dofs!" );

    Dest.setZero();

    // Energy part
    _squaredResidualEnergy.apply( Arg, Dest );
    Dest *= _gamma;

    // Model fitting part
    VecType vertexPosition;
    RealType differenceNorm;
    RealType fitContrib = 0;

    for ( int m = 0; m < _markerIndices.size(); m++ ) {
      getXYZCoord<typename VectorType::ConstSegmentReturnType, VecType>( Arg.segment( 0, _numLocalDofs ),
                                                                         vertexPosition, _markerIndices[m] );

      differenceNorm = (vertexPosition - _markerPositions[m]).squaredNorm();

#pragma omp critical
      {
        fitContrib += differenceNorm;
      }
    }

    Dest[0] += fitContrib;
  }
};

//! \brief Derivative of LThetaLandmarkFittingFunctional w.r.t. s and \omega
//! \author Sassen
template<typename ConfiguratorType>
class LThetaLandmarkFittingGradient
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;


  const DeformationBase<ConfiguratorType> &_W;
  const MeshTopologySaver &_Topology;
  const std::vector<int> &_markerIndices;
  const std::vector<VecType> &_markerPositions;
  const RealType _gamma;

  const int _numLThetaDofs, _numLocalDofs, _numVertices;
  const unsigned long _numInputData;

  LinearBlendingResidual<ConfiguratorType> _LThetaResidualEnergy;
  LinearBlendingResidualGradient<ConfiguratorType> _LThetaResidualGradient;
  SquaredDerivative<ConfiguratorType> _squaredResidualGradient;
public:
  LThetaLandmarkFittingGradient( const MeshTopologySaver &Topology,
                                 const DeformationBase<ConfiguratorType> &W,
                                 const std::vector<int> &markerIndices,
                                 const std::vector<VecType> &markerPositions,
                                 const VectorType &refMeshLengthsAngles,
                                 const VectorType &IntegrationWeights,
                                 const RealType edgeLengthWeight,
                                 const RealType bendingWeight,
                                 const std::vector<VectorType> &inputDataLengthsAngles,
                                 RealType priorWeight = 10.)
          : _Topology( Topology ), _W( W ), _markerIndices( markerIndices ), _markerPositions( markerPositions ),
            _gamma( priorWeight ),
            _numLThetaDofs( refMeshLengthsAngles.size()), _numVertices( Topology.getNumVertices()),
            _numInputData( inputDataLengthsAngles.size()), _numLocalDofs( 3 * Topology.getNumVertices()),
            _LThetaResidualEnergy( Topology, refMeshLengthsAngles, IntegrationWeights, edgeLengthWeight, bendingWeight,
                                   inputDataLengthsAngles ),
            _LThetaResidualGradient( Topology, refMeshLengthsAngles, IntegrationWeights, edgeLengthWeight,
                                     bendingWeight, inputDataLengthsAngles, false ),
            _squaredResidualGradient(_LThetaResidualEnergy, _LThetaResidualGradient)
  {
    if ( markerIndices.size() != markerPositions.size())
      throw std::length_error( "LThetaLandmarkFittingGradient: number of marker indices and positions have to agree!" );

  }

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    if ( Arg.size() != _numLocalDofs + _numInputData )
      throw std::length_error( "LandmarkFittingGradient::apply: wrong size of dofs!" );

    Dest.resize(Arg.size());
    Dest.setZero();

    // Energy contribution
    _squaredResidualGradient.apply( Arg, Dest );
    Dest *= _gamma;

//    for (int i = 0; i < Arg.size(); i++)
//      if (Dest[i] > 1e3)
//        Dest[i] = 0;

//    int i;
//    std::cout << "Absurd: " << (Dest.array() > 1e6).sum() << " " << Dest.maxCoeff(&i) << std::endl;
//    std::cout << "Max entry: " << i << std::endl;

    // Model fitting contribution
    VecType vertexPosition;
    RealType differenceNorm;
    VectorType fitContrib = VectorType::Zero(Arg.size());
//#ifdef _OPENMP
//#pragma omp parallel for private(vertexPosition, differenceNorm)
//#endif
    for ( int m = 0; m < _markerIndices.size(); m++ ) {
      getXYZCoord<typename VectorType::ConstSegmentReturnType, VecType>( Arg.segment( 0, _numLocalDofs ),
                                                                         vertexPosition, _markerIndices[m] );
      vertexPosition -=  _markerPositions[m];

#pragma omp critical
      {
        fitContrib[_markerIndices[m]] += 2 * vertexPosition[0];
        fitContrib[_numVertices + _markerIndices[m]] += 2 * vertexPosition[1];
        fitContrib[2 * _numVertices + _markerIndices[m]] += 2 * vertexPosition[2];
      }
    }

    Dest += fitContrib;

//    std::cout << "Absurd: " << (Dest.array() > 1e6).sum() << " " << Dest.maxCoeff(&i) << std::endl;
//    std::cout << "Max entry: " << i << std::endl;
  }
};

//! \brief Hessian of LThetaLandmarkFittingFunctional w.r.t. s and \omega
//! \author Sassen
template<typename ConfiguratorType>
class LThetaLandmarkFittingHessian
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef std::vector<TripletType> TripletListType;

  const DeformationBase<ConfiguratorType> &_W;
  const MeshTopologySaver &_Topology;
  const std::vector<int> &_markerIndices;
  const std::vector<VecType> &_markerPositions;
  const RealType _gamma;

  const int _numLThetaDofs, _numLocalDofs, _numVertices;
  const unsigned long _numInputData;

  LinearBlendingResidual<ConfiguratorType> _LThetaResidualEnergy;
  LinearBlendingResidualGradient<ConfiguratorType> _LThetaResidualGradient;
  ReducedSquaredHessian<ConfiguratorType> _squaredResidualHessian;


public:
  LThetaLandmarkFittingHessian( const MeshTopologySaver &Topology,
                                const DeformationBase<ConfiguratorType> &W,
                                const std::vector<int> &markerIndices,
                                const std::vector<VecType> &markerPositions,
                                const VectorType &refMeshLengthsAngles,
                                const VectorType &IntegrationWeights,
                                const RealType edgeLengthWeight,
                                const RealType bendingWeight,
                                const std::vector<VectorType> &inputDataLengthsAngles,
                                RealType priorWeight = 10. )
          : _Topology( Topology ), _W( W ), _markerIndices( markerIndices ), _markerPositions( markerPositions ),
            _gamma( priorWeight ),
            _numLThetaDofs( refMeshLengthsAngles.size()), _numVertices( Topology.getNumVertices()),
            _numInputData( inputDataLengthsAngles.size()), _numLocalDofs( 3 * Topology.getNumVertices()),
            _LThetaResidualEnergy( Topology, refMeshLengthsAngles, IntegrationWeights, edgeLengthWeight, bendingWeight,
                                   inputDataLengthsAngles ),
            _LThetaResidualGradient( Topology, refMeshLengthsAngles, IntegrationWeights, edgeLengthWeight,
                                     bendingWeight, inputDataLengthsAngles, false ),
            _squaredResidualHessian(_LThetaResidualEnergy, _LThetaResidualGradient)
  {
    if ( markerIndices.size() != markerPositions.size())
      throw std::length_error( "LThetaLandmarkFittingHessian: number of marker indices and positions have to agree!" );
  }

  void apply( const VectorType &Arg, MatrixType &Dest ) const override {
    if ( Arg.size() != _numLocalDofs + _numInputData )
      throw std::length_error( "LThetaLandmarkFittingHessian::apply: wrong size of dofs!" );

    if ((Dest.cols() != Arg.size()) || (Dest.rows() != Arg.size()))
      Dest.resize( Arg.size(), Arg.size());

    Dest.setZero();

    // Model fitting contribution
    // fill triplet lists for marker fitting hessian
    TripletListType tripletList;
    tripletList.reserve(_markerPositions.size() * 3);
    pushFittingTriplets( Arg, tripletList );
    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend());

    // Prior contribution
    MatrixType temp;
    temp.resize( Arg.size(), Arg.size());
    _squaredResidualHessian.apply( Arg, temp );

    Dest += _gamma * temp;
  }

  void pushFittingTriplets( const VectorType &Arg, TripletListType &tripletList ) const {
    if ( Arg.size() != _numLocalDofs + _numInputData )
      throw std::length_error( "LThetaLandmarkFittingHessian::pushFittingTriplets: wrong size of dofs!" );

    for ( int m = 0; m < _markerIndices.size(); m++ ) {
      tripletList.push_back( TripletType( _markerIndices[m], _markerIndices[m], 2 ));
      tripletList.push_back( TripletType( _numVertices + _markerIndices[m], _numVertices + _markerIndices[m], 2 ));
      tripletList.push_back( TripletType( 2 * _numVertices + _markerIndices[m], 2 * _numVertices + _markerIndices[m], 2 ));
    }
  }

};

//! \brief Mismatch of discrete shell with given landmarks considering distance to a set of lengths and angles as prior
//! \author Sassen
//!
//! Given landmark positions x_0,...,x_l, corresponding vertex positions X_0(s),..,X_l(s) this implements
//! F_x[s] = \sum_{l=1}^L ||X_l(s) - x_l||^2
//!          + \gamma \sum_{e \in E} (\alpha_e (l_e[s] - l_e^* )^2 + \eta\, \beta_e (\theta_e[s] - \theta_e^* )^2)
//! for a shell s.
template<typename ConfiguratorType>
class SoftReconstructionFunctional
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;

  const DeformationBase<ConfiguratorType> &_W;
  const MeshTopologySaver &_Topology;
  const std::vector<int> &_markerIndices;
  const std::vector<VecType> &_markerPositions;
  const RealType _gamma;

  const int _numLThetaDofs, _numLocalDofs;

  LinearReconstructionFunctional<ConfiguratorType> _LThetaResidualEnergy;
  SquaredFunctional<ConfiguratorType> _squaredResidualEnergy;

public:
  SoftReconstructionFunctional( const MeshTopologySaver &Topology,
                                const DeformationBase<ConfiguratorType> &W,
                                const std::vector<int> &markerIndices,
                                const std::vector<VecType> &markerPositions,
                                const VectorType &refMeshLengthsAngles,
                                const VectorType &IntegrationWeights,
                                const RealType edgeLengthWeight,
                                const RealType bendingWeight,
                                RealType priorWeight = 10. )
          : _Topology( Topology ), _W( W ), _markerIndices( markerIndices ), _markerPositions( markerPositions ),
            _gamma( priorWeight ), _numLThetaDofs( refMeshLengthsAngles.size()),
            _numLocalDofs( 3 * Topology.getNumVertices()),
            _LThetaResidualEnergy( Topology, refMeshLengthsAngles, IntegrationWeights, edgeLengthWeight, bendingWeight  ),
            _squaredResidualEnergy( _LThetaResidualEnergy ){
    if ( markerIndices.size() != markerPositions.size())
      throw std::length_error( "SoftReconstructionFunctional: number of marker indices and positions have to "
                               "agree!" );


  }

  void apply( const VectorType &Arg, RealType &Dest ) const override {
    if ( Arg.size() != _numLocalDofs )
      throw std::length_error( "SoftReconstructionFunctional::apply: wrong size of dofs!" );

    Dest.setZero();

    // Energy part
    _squaredResidualEnergy.apply( Arg, Dest );
    Dest *= _gamma;

    // Model fitting part
    VecType vertexPosition;
    RealType differenceNorm;
    RealType fitContrib = 0;

    for ( int m = 0; m < _markerIndices.size(); m++ ) {
      getXYZCoord<VectorType, VecType>( Arg, vertexPosition, _markerIndices[m] );
      differenceNorm = (vertexPosition - _markerPositions[m]).squaredNorm();

#pragma omp critical
      {
        fitContrib += differenceNorm;
      }
    }

//    std::cout << "Energy: " << Dest[0] << " " << fitContrib << std::endl;

    Dest[0] += fitContrib;
  }
};

//! \brief Derivative of SoftReconstructionFunctional w.r.t. s
//! \author Sassen
template<typename ConfiguratorType>
class SoftReconstructionGradient
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;


  const DeformationBase<ConfiguratorType> &_W;
  const MeshTopologySaver &_Topology;
  const std::vector<int> &_markerIndices;
  const std::vector<VecType> &_markerPositions;
  const RealType _gamma;

  const int _numLThetaDofs, _numLocalDofs, _numVertices;

  LinearReconstructionFunctional<ConfiguratorType> _LThetaResidualEnergy;
  LinearReconstructionDerivative<ConfiguratorType> _LThetaResidualGradient;
  SquaredDerivative<ConfiguratorType> _squaredResidualGradient;
public:
  SoftReconstructionGradient( const MeshTopologySaver &Topology,
                              const DeformationBase<ConfiguratorType> &W,
                              const std::vector<int> &markerIndices,
                              const std::vector<VecType> &markerPositions,
                              const VectorType &refMeshLengthsAngles,
                              const VectorType &IntegrationWeights,
                              const RealType edgeLengthWeight,
                              const RealType bendingWeight,
                              RealType priorWeight = 10. )
          : _Topology( Topology ), _W( W ), _markerIndices( markerIndices ), _markerPositions( markerPositions ),
            _gamma( priorWeight ),
            _numLThetaDofs( refMeshLengthsAngles.size()), _numVertices( Topology.getNumVertices()),
            _numLocalDofs( 3 * Topology.getNumVertices()),
            _LThetaResidualEnergy( Topology, refMeshLengthsAngles, IntegrationWeights, edgeLengthWeight,
                                   bendingWeight ),
            _LThetaResidualGradient( Topology, refMeshLengthsAngles, IntegrationWeights, edgeLengthWeight,
                                     bendingWeight ),
            _squaredResidualGradient( _LThetaResidualEnergy, _LThetaResidualGradient ) {
    if ( markerIndices.size() != markerPositions.size())
      throw std::length_error( "SoftReconstructionGradient: number of marker indices and positions have to agree!" );

  }

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    if ( Arg.size() != _numLocalDofs )
      throw std::length_error( "SoftReconstructionGradient::apply: wrong size of dofs!" );

    Dest.resize( _numLocalDofs );
    Dest.setZero();

    // Energy contribution
    _squaredResidualGradient.apply( Arg, Dest );
    Dest *= _gamma;
//
//    for (int i = 0; i < Arg.size(); i++)
//      if (Dest[i] > 1e3)
//        Dest[i] = 0;

//    int i;
//    std::cout << "Absurd: " << (Dest.array() > 1e6).sum() << " " << Dest.maxCoeff(&i) << std::endl;
//    std::cout << "Max entry: " << i << std::endl;


    // Model fitting contribution
    VecType vertexPosition;
    RealType differenceNorm;
    VectorType fitContrib = VectorType::Zero( _numLocalDofs );
//#ifdef _OPENMP
//#pragma omp parallel for private(vertexPosition, differenceNorm)
//#endif
    for ( int m = 0; m < _markerIndices.size(); m++ ) {
      getXYZCoord<VectorType, VecType>( Arg, vertexPosition, _markerIndices[m] );
      vertexPosition -=  _markerPositions[m];

#pragma omp critical
      {
        fitContrib[_markerIndices[m]] += 2 * vertexPosition[0];
        fitContrib[_numVertices + _markerIndices[m]] += 2 * vertexPosition[1];
        fitContrib[2 * _numVertices + _markerIndices[m]] += 2 * vertexPosition[2];
      }
    }

//    std::cout << "Gradient: " << Dest.norm() << " " << fitContrib.norm() << std::endl;

    Dest += fitContrib;
//    std::cout << "Absurd: " << (Dest.array() > 1e6).sum() << " " << Dest.maxCoeff(&i) << std::endl;
//    std::cout << "Max entry: " << i << std::endl;
  }
};

//! \brief Hessian of SoftReconstructionFunctional w.r.t. s
//! \author Sassen
template<typename ConfiguratorType>
class SoftReconstructionHessian
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef std::vector<TripletType> TripletListType;

  const DeformationBase<ConfiguratorType> &_W;
  const MeshTopologySaver &_Topology;
  const std::vector<int> &_markerIndices;
  const std::vector<VecType> &_markerPositions;
  const RealType _gamma;

  const int _numLThetaDofs, _numLocalDofs, _numVertices;

  LinearReconstructionFunctional<ConfiguratorType> _LThetaResidualEnergy;
  LinearReconstructionDerivative<ConfiguratorType> _LThetaResidualGradient;
  LinearReconstructionHessian<ConfiguratorType> _LThetaResidualHessian;
  ReducedSquaredHessian<ConfiguratorType> _squaredResidualHessian;

public:
  SoftReconstructionHessian( const MeshTopologySaver &Topology,
                             const DeformationBase<ConfiguratorType> &W,
                             const std::vector<int> &markerIndices,
                             const std::vector<VecType> &markerPositions,
                             const VectorType &refMeshLengthsAngles,
                             const VectorType &IntegrationWeights,
                             const RealType edgeLengthWeight,
                             const RealType bendingWeight,
                             RealType priorWeight = 10. )
          : _Topology( Topology ), _W( W ), _markerIndices( markerIndices ), _markerPositions( markerPositions ),
            _gamma( priorWeight ),
            _numLThetaDofs( refMeshLengthsAngles.size()), _numVertices( Topology.getNumVertices()),
            _numLocalDofs( 3 * Topology.getNumVertices()),
            _LThetaResidualEnergy( Topology, refMeshLengthsAngles, IntegrationWeights, edgeLengthWeight,
                                   bendingWeight ),
            _LThetaResidualGradient( Topology, refMeshLengthsAngles, IntegrationWeights, edgeLengthWeight,
                                     bendingWeight ),
            _LThetaResidualHessian( Topology, refMeshLengthsAngles, IntegrationWeights, edgeLengthWeight,
                                     bendingWeight ),
            _squaredResidualHessian( _LThetaResidualEnergy, _LThetaResidualGradient ) {
    if ( markerIndices.size() != markerPositions.size())
      throw std::length_error( "SoftReconstructionHessian: number of marker indices and positions have to agree!" );
  }

  void apply( const VectorType &Arg, MatrixType &Dest ) const override {
    if ( Arg.size() != _numLocalDofs )
      throw std::length_error( "SoftReconstructionHessian::apply: wrong size of dofs!" );

    if ((Dest.cols() != _numLocalDofs) || (Dest.rows() != _numLocalDofs))
      Dest.resize( _numLocalDofs, _numLocalDofs );

    Dest.setZero();

    // Model fitting contribution
    // fill triplet lists for marker fitting hessian
    TripletListType tripletList;
    tripletList.reserve( _markerPositions.size() * 3 );
    pushFittingTriplets( Arg, tripletList );
    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend());

    // Prior contribution
    MatrixType temp;
    temp.resize( _numLocalDofs, _numLocalDofs );
    _squaredResidualHessian.apply( Arg, temp );

    Dest += _gamma * temp;
  }

  void pushFittingTriplets( const VectorType &Arg, TripletListType &tripletList ) const {
    if ( Arg.size() != _numLocalDofs )
      throw std::length_error( "SoftReconstructionHessian::pushFittingTriplets: wrong size of dofs!" );

    for ( int m = 0; m < _markerIndices.size(); m++ ) {
      tripletList.push_back( TripletType( _markerIndices[m], _markerIndices[m], 2 ));
      tripletList.push_back( TripletType( _numVertices + _markerIndices[m], _numVertices + _markerIndices[m], 2 ));
      tripletList.push_back(
              TripletType( 2 * _numVertices + _markerIndices[m], 2 * _numVertices + _markerIndices[m], 2 ));
    }
  }

};

#endif //NRIC_LANDMARKFITTING_H

