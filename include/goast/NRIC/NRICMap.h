// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Operators mapping vertex positions to lengths and angles
 * \author Sassen
 */

#ifndef NRIC_MAP_H
#define NRIC_MAP_H

#include <goast/Core/Auxiliary.h>
#include <goast/Core/LocalMeshGeometry.h>
#include <goast/Core/Topology.h>

#include <utility>

/**
 * \brief Operator mapping vertex positions to the corresponding edge lengths and dihedral angles
 * \tparam ConfiguratorType Container with data types
 */
template<typename ConfiguratorType>
class NRICMap : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;

  const MeshTopologySaver &_topology;
  const std::vector<int> _consideredIndices;
  const int _numEdges;

public:
  /**
   * \param topology The topology of the mesh.
   */
  explicit NRICMap( const MeshTopologySaver &topology ) : _topology( topology ), _numEdges( _topology.getNumEdges()),
                                                          _consideredIndices( 0 ) {}

  /**
   * \param topology The topology of the mesh.
   * \param consideredIndices A subset of the edges to which the output will be restricted.
   */
  NRICMap( const MeshTopologySaver &topology, std::vector<int> consideredIndices ) :
          _topology( topology ), _consideredIndices( std::move( consideredIndices )),
          _numEdges( _topology.getNumEdges()) {}

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    if ( Arg.size() != 3 * _topology.getNumVertices())
      throw std::length_error( "NRICMap::apply(): Arg has wrong size!" );

    if ( _consideredIndices.size() > 0 ) {
      Dest.resize( _consideredIndices.size());
      Dest.setZero();
    }
    else {
      Dest.resize( 2 * _topology.getNumEdges());
      Dest.setZero();
    }
    if ( _consideredIndices.size() > 0 ) {
      for ( int j = 0; j < _consideredIndices.size(); j++ ) {
        int edgeIdx = _consideredIndices[j] % _numEdges;


        int pi( _topology.getAdjacentNodeOfEdge( edgeIdx, 0 )),
                pj( _topology.getAdjacentNodeOfEdge( edgeIdx, 1 )),
                pk( _topology.getOppositeNodeOfEdge( edgeIdx, 0 )),
                pl( _topology.getOppositeNodeOfEdge( edgeIdx, 1 ));

        // set up vertices
        VecType Pi, Pj;
        getXYZCoord<VectorType, VecType>( Arg, Pi, pi );
        getXYZCoord<VectorType, VecType>( Arg, Pj, pj );

        if ( _consideredIndices[j] >= _numEdges ) {
          // compute dihedral angle  (no bending at boundary edges)
          if ( std::min( pl, pk ) == -1 )
            continue;

          VecType Pk, Pl;
          getXYZCoord<VectorType, VecType>( Arg, Pk, pk );
          getXYZCoord<VectorType, VecType>( Arg, Pl, pl );
          Dest[j] = getDihedralAngle( Pi, Pj, Pk, Pl );
        }
        else {
          Dest[j] = std::sqrt( dotProduct( Pj - Pi, Pj - Pi ));
        }
      }
    }
    else {
      for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ) {

        int pi( _topology.getAdjacentNodeOfEdge( edgeIdx, 0 )),
                pj( _topology.getAdjacentNodeOfEdge( edgeIdx, 1 )),
                pk( _topology.getOppositeNodeOfEdge( edgeIdx, 0 )),
                pl( _topology.getOppositeNodeOfEdge( edgeIdx, 1 ));

        // set up vertices
        VecType Pi, Pj;
        getXYZCoord<VectorType, VecType>( Arg, Pi, pi );
        getXYZCoord<VectorType, VecType>( Arg, Pj, pj );

        // compute length of edge
        Dest[edgeIdx] = std::sqrt( dotProduct( Pj - Pi, Pj - Pi ));

        // compute dihedral angle  (no bending at boundary edges)
        if ( std::min( pl, pk ) == -1 )
          continue;

        VecType Pk, Pl;
        getXYZCoord<VectorType, VecType>( Arg, Pk, pk );
        getXYZCoord<VectorType, VecType>( Arg, Pl, pl );
        Dest[_topology.getNumEdges() + edgeIdx] = getDihedralAngle( Pi, Pj, Pk, Pl );
      }
    }
  }

  int getTargetDimension() const override {
    if ( _consideredIndices.size() > 0 ) {
      return _consideredIndices.size();
    }
    else {
      return 2 * _topology.getNumEdges();
    }
  }
};

/**
 * \brief Derivative of operator mapping vertex positions to edge lengths and dihedral angles
 * \tparam ConfiguratorType Container with data types
 * \see NRICMap
 */
template<typename ConfiguratorType>
class NRICMapGradient
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver &_topology;
  const std::vector<int> _consideredIndices;
  const int _numEdges;

public:
  /**
   * \param topology The topology of the mesh.
   */
  explicit NRICMapGradient( const MeshTopologySaver &topology ) : _topology( topology ),
                                                                  _numEdges( _topology.getNumEdges()) {}

  /**
   * \param topology The topology of the mesh.
   * \param consideredIndices A subset of the edges to which the output will be restricted.
   */
  NRICMapGradient( const MeshTopologySaver &topology, std::vector<int> consideredIndices ) :
          _topology( topology ), _consideredIndices( std::move( consideredIndices )),
          _numEdges( _topology.getNumEdges()) {}

  void apply( const VectorType &Arg, SparseMatrixType &Dest ) const override {
    if ( Arg.size() != 3 * _topology.getNumVertices())
      throw std::length_error( "NRICMapGradient::apply(): Arg has wrong size!" );

    if ( _consideredIndices.size() > 0 ) {
      Dest.resize( _consideredIndices.size(), 3 * _topology.getNumVertices());
      Dest.setZero();
    }
    else {
      Dest.resize( 2 * _topology.getNumEdges(), 3 * _topology.getNumVertices());
      Dest.setZero();
    }
    TripletListType tripletList;
    pushTriplets( Arg, tripletList );

    Dest.setFromTriplets( tripletList.begin(), tripletList.end());
  }

  void pushTriplets( const VectorType &Arg, TripletListType &tripletList, bool transposed = false ) const {
    if ( _consideredIndices.size() > 0 ) {
      tripletList.reserve( 3 * ( 2 + 4 ) * _topology.getNumEdges());

      const int numVertices = _topology.getNumVertices(); // colOffset
      const int numEdges = _topology.getNumEdges(); // rowOffset

      for ( int j = 0; j < _consideredIndices.size(); j++ ) {
        int edgeIdx = _consideredIndices[j] % _numEdges;

        int pi( _topology.getAdjacentNodeOfEdge( edgeIdx, 0 )),
                pj( _topology.getAdjacentNodeOfEdge( edgeIdx, 1 )),
                pk( _topology.getOppositeNodeOfEdge( edgeIdx, 0 )),
                pl( _topology.getOppositeNodeOfEdge( edgeIdx, 1 ));

        // set up vertices
        VecType Pi, Pj;
        getXYZCoord<VectorType, VecType>( Arg, Pi, pi );
        getXYZCoord<VectorType, VecType>( Arg, Pj, pj );
        VecType edge( Pi - Pj );
        RealType edgeLength = edge.norm();

        if ( _consideredIndices[j] < _numEdges ) {
          // assemble in global matrix
          for ( int i = 0; i < 3; i++ ) {
            pushTriplet( tripletList, j, i * numVertices + pi, edge[i] / edgeLength, transposed );
            pushTriplet( tripletList, j, i * numVertices + pj, -1. * edge[i] / edgeLength, transposed );
          }
        }
        else {
          // no dihedral angle at boundary edges
          if ( std::min( pl, pk ) < 0 )
            continue;

          VecType Pk, Pl;
          getXYZCoord<VectorType, VecType>( Arg, Pk, pk );
          getXYZCoord<VectorType, VecType>( Arg, Pl, pl );

          // compute first derivatives of dihedral angle
          VecType thetak, thetal, thetai, thetaj;
          getThetaGradK( Pi, Pj, Pk, thetak );
          getThetaGradK( Pj, Pi, Pl, thetal );
          getThetaGradI( Pi, Pj, Pk, Pl, thetai );
          getThetaGradJ( Pi, Pj, Pk, Pl, thetaj );

          // assemble in global matrix
          for ( int i = 0; i < 3; i++ ) {
            pushTriplet( tripletList, j, i * numVertices + pi, thetai[i], transposed );
            pushTriplet( tripletList, j, i * numVertices + pj, thetaj[i], transposed );
            pushTriplet( tripletList, j, i * numVertices + pk, thetak[i], transposed );
            pushTriplet( tripletList, j, i * numVertices + pl, thetal[i], transposed );
          }
        }
      }
    }
    else {
      tripletList.reserve( 3 * ( 2 + 4 ) * _topology.getNumEdges());

      const int numVertices = _topology.getNumVertices(); // colOffset
      const int numEdges = _topology.getNumEdges(); // rowOffset

      for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ) {
        int pi( _topology.getAdjacentNodeOfEdge( edgeIdx, 0 )),
                pj( _topology.getAdjacentNodeOfEdge( edgeIdx, 1 )),
                pk( _topology.getOppositeNodeOfEdge( edgeIdx, 0 )),
                pl( _topology.getOppositeNodeOfEdge( edgeIdx, 1 ));

        // set up vertices
        VecType Pi, Pj;
        getXYZCoord<VectorType, VecType>( Arg, Pi, pi );
        getXYZCoord<VectorType, VecType>( Arg, Pj, pj );
        VecType edge( Pi - Pj );
        RealType edgeLength = edge.norm();

        // assemble in global matrix
        for ( int i = 0; i < 3; i++ ) {
          pushTriplet( tripletList, edgeIdx, i * numVertices + pi, edge[i] / edgeLength, transposed );
          pushTriplet( tripletList, edgeIdx, i * numVertices + pj, -1. * edge[i] / edgeLength, transposed );
        }

        // no dihedral angle at boundary edges
        if ( std::min( pl, pk ) < 0 )
          continue;

        VecType Pk, Pl;
        getXYZCoord<VectorType, VecType>( Arg, Pk, pk );
        getXYZCoord<VectorType, VecType>( Arg, Pl, pl );

        // compute first derivatives of dihedral angle
        VecType thetak, thetal, thetai, thetaj;
        getThetaGradK( Pi, Pj, Pk, thetak );
        getThetaGradK( Pj, Pi, Pl, thetal );
        getThetaGradI( Pi, Pj, Pk, Pl, thetai );
        getThetaGradJ( Pi, Pj, Pk, Pl, thetaj );

        // assemble in global matrix
        for ( int i = 0; i < 3; i++ ) {
          pushTriplet( tripletList, numEdges + edgeIdx, i * numVertices + pi, thetai[i], transposed );
          pushTriplet( tripletList, numEdges + edgeIdx, i * numVertices + pj, thetaj[i], transposed );
          pushTriplet( tripletList, numEdges + edgeIdx, i * numVertices + pk, thetak[i], transposed );
          pushTriplet( tripletList, numEdges + edgeIdx, i * numVertices + pl, thetal[i], transposed );
        }
      }
    }


//    tripletList.clear();

  }

  int getTargetDimension() const override {
    if ( _consideredIndices.size() > 0 ) {
      return _consideredIndices.size();
    }
    else {
      return 2 * _topology.getNumEdges();
    }
  }

  int getNNZ() const override {
    return 3 * 4 * _consideredIndices.size();
  }

protected:
  void pushTriplet( TripletListType &tripletList, int row, int col, RealType value, bool transposed ) const {
    if ( transposed )
      tripletList.push_back( TripletType( col, row, value ));
    else
      tripletList.push_back( TripletType( row, col, value ));
  }
};

/**
 * \brief Hessian of operator mapping vertex positions to edge lengths and dihedral angles
 * \tparam ConfiguratorType Container with data types
 * \see NRICMap
 */
template<typename ConfiguratorType>
class NRICMapHessian : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::TensorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;
  typedef typename ConfiguratorType::TensorType TensorType;

  const MeshTopologySaver &_topology;
  const std::vector<int> _consideredIndices;
  const int _numEdges;

public:
  /**
   * \param topology The topology of the mesh.
   */
  explicit NRICMapHessian( const MeshTopologySaver &topology ) : _topology( topology ),
                                                                 _numEdges( _topology.getNumEdges()) {}

  /**
   * \param topology The topology of the mesh.
   * \param consideredIndices A subset of the edges to which the output will be restricted.
   */
  NRICMapHessian( const MeshTopologySaver &topology, std::vector<int> consideredIndices ) :
          _topology( topology ), _consideredIndices( std::move( consideredIndices )),
          _numEdges( _topology.getNumEdges()) {}

  void apply( const VectorType &Arg, TensorType &Dest ) const override {
    if ( Arg.size() != 3 * _topology.getNumVertices())
      throw std::length_error( "NRICMapHessian::apply(): Arg has wrong size!" );

    if ( _consideredIndices.size() > 0 ) {
      Dest.resize( _consideredIndices.size(), 3 * _topology.getNumVertices(), 3 * _topology.getNumVertices());
      Dest.setZero();
    }
    else {
      Dest.resize( 2 * _topology.getNumEdges(), 3 * _topology.getNumVertices(), 3 * _topology.getNumVertices());
      Dest.setZero();
    }


    std::vector<TripletListType> tripletLists;
    pushTriplets( Arg, tripletLists );

    Dest.setFromTriplets( tripletLists );
  }

  void apply( const VectorType &Arg, std::vector<SparseMatrixType> &Dest ) const {
    if ( Arg.size() != 3 * _topology.getNumVertices())
      throw std::length_error( "NRICMapHessian::apply(): Arg has wrong size!" );

    std::vector<TripletListType> tripletLists;
    if ( _consideredIndices.size() > 0 ) {
      Dest.resize( _consideredIndices.size());
      tripletLists.resize( _consideredIndices.size());
    }
    else {
      Dest.resize( 2 * _topology.getNumEdges());
      tripletLists.resize( 2 * _topology.getNumEdges());
    }


    pushTriplets( Arg, tripletLists );

#ifdef GOAST_WITH_OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < tripletLists.size(); i++ ) {

      Dest[i].resize( 3 * _topology.getNumVertices(), 3 * _topology.getNumVertices());

      Dest[i].reserve( tripletLists[i].size());
      Dest[i].setFromTriplets( tripletLists[i].begin(), tripletLists[i].end());
    }
  }

  void pushTriplets( const VectorType &Arg, std::vector<TripletListType> &tripletLists ) const {
    if ( _consideredIndices.size() > 0 ) {
      tripletLists.resize( _consideredIndices.size());

      const int numVertices = _topology.getNumVertices(); // colOffset
      const int numEdges = _topology.getNumEdges(); // rowOffset

      for ( int j = 0; j < _consideredIndices.size(); j++ ) {
        int edgeIdx = _consideredIndices[j] % _numEdges;

        int pi( _topology.getAdjacentNodeOfEdge( edgeIdx, 0 )),
                pj( _topology.getAdjacentNodeOfEdge( edgeIdx, 1 )),
                pk( _topology.getOppositeNodeOfEdge( edgeIdx, 0 )),
                pl( _topology.getOppositeNodeOfEdge( edgeIdx, 1 ));

        // set up vertices
        VecType Pi, Pj;
        getXYZCoord<VectorType, VecType>( Arg, Pi, pi );
        getXYZCoord<VectorType, VecType>( Arg, Pj, pj );
        VecType edge( Pi - Pj );
        RealType edgeLength = edge.norm();

        if ( _consideredIndices[j] < _numEdges ) {
          MatType tensorProduct;
          tensorProduct.makeTensorProduct( edge, edge );
          tensorProduct *= -1. / ( edgeLength * edgeLength * edgeLength );
          tensorProduct.addToDiagonal( 1 / edgeLength );

          localToGlobal( tripletLists, j, pi, pi, tensorProduct );
          localToGlobal( tripletLists, j, pj, pj, tensorProduct );
          tensorProduct *= 1;
          localToGlobal( tripletLists, j, pi, pj, tensorProduct );
        }
        else {
          // no dihedral angle at boundary edges
          if ( std::min( pl, pk ) == -1 )
            continue;

          VecType Pk, Pl;
          getXYZCoord<VectorType, VecType>( Arg, Pk, pk );
          getXYZCoord<VectorType, VecType>( Arg, Pl, pl );

          // compute first derivatives of dihedral angle
          VecType thetak, thetal, thetai, thetaj;
          getThetaGradK( Pi, Pj, Pk, thetak );
          getThetaGradK( Pj, Pi, Pl, thetal );
          getThetaGradI( Pi, Pj, Pk, Pl, thetai );
          getThetaGradJ( Pi, Pj, Pk, Pl, thetaj );

          RealType delThetaDouble = 1.;

          // now compute second derivatives of dihedral angle
          MatType H;


          //kk
          getHessThetaKK( Pi, Pj, Pk, H );
          localToGlobal( tripletLists, j, pk, pk, H );

          //ik & ki (Hki = Hik)
          getHessThetaIK( Pi, Pj, Pk, H );
          localToGlobal( tripletLists, j, pi, pk, H );

          //jk & kj (Hkj = Hjk)
          getHessThetaJK( Pi, Pj, Pk, H );
          localToGlobal( tripletLists, j, pj, pk, H );

          //ll
          getHessThetaKK( Pj, Pi, Pl, H );
          localToGlobal( tripletLists, j, pl, pl, H );

          //il & li (Hli = Hil)
          getHessThetaJK( Pj, Pi, Pl, H );
          localToGlobal( tripletLists, j, pi, pl, H );

          //jl & lj (Hlj = Hjl)
          getHessThetaIK( Pj, Pi, Pl, H );
          localToGlobal( tripletLists, j, pj, pl, H );


          //ii
          getHessThetaII( Pi, Pj, Pk, Pl, H );
          localToGlobal( tripletLists, j, pi, pi, H );

          //jj
          getHessThetaII( Pj, Pi, Pl, Pk, H );
          localToGlobal( tripletLists, j, pj, pj, H );

          //ij & ji (Hij = Hji)
          getHessThetaJI( Pi, Pj, Pk, Pl, H );
          localToGlobal( tripletLists, j, pi, pj, H );
        }

      }
    }
    else {
      tripletLists.resize( 2 * _topology.getNumEdges());

      const int numVertices = _topology.getNumVertices(); // colOffset
      const int numEdges = _topology.getNumEdges(); // rowOffset

      for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ) {
        int pi( _topology.getAdjacentNodeOfEdge( edgeIdx, 0 )),
                pj( _topology.getAdjacentNodeOfEdge( edgeIdx, 1 )),
                pk( _topology.getOppositeNodeOfEdge( edgeIdx, 0 )),
                pl( _topology.getOppositeNodeOfEdge( edgeIdx, 1 ));

        // set up vertices
        VecType Pi, Pj;
        getXYZCoord<VectorType, VecType>( Arg, Pi, pi );
        getXYZCoord<VectorType, VecType>( Arg, Pj, pj );
        VecType edge( Pi - Pj );
        RealType edgeLength = edge.norm();

        MatType tensorProduct;
        tensorProduct.makeTensorProduct( edge, edge );
        tensorProduct *= -1. / ( edgeLength * edgeLength * edgeLength );
        tensorProduct.addToDiagonal( 1 / edgeLength );

        localToGlobal( tripletLists, edgeIdx, pi, pi, tensorProduct );
        localToGlobal( tripletLists, edgeIdx, pj, pj, tensorProduct );
        tensorProduct *= 1;
        localToGlobal( tripletLists, edgeIdx, pi, pj, tensorProduct );

        // no dihedral angle at boundary edges
        if ( std::min( pl, pk ) == -1 )
          continue;

        VecType Pk, Pl;
        getXYZCoord<VectorType, VecType>( Arg, Pk, pk );
        getXYZCoord<VectorType, VecType>( Arg, Pl, pl );

        // compute first derivatives of dihedral angle
        VecType thetak, thetal, thetai, thetaj;
        getThetaGradK( Pi, Pj, Pk, thetak );
        getThetaGradK( Pj, Pi, Pl, thetal );
        getThetaGradI( Pi, Pj, Pk, Pl, thetai );
        getThetaGradJ( Pi, Pj, Pk, Pl, thetaj );

        RealType delThetaDouble = 1.;

        // now compute second derivatives of dihedral angle
        MatType H;


        //kk
        getHessThetaKK( Pi, Pj, Pk, H );
        localToGlobal( tripletLists, numEdges + edgeIdx, pk, pk, H );

        //ik & ki (Hki = Hik)
        getHessThetaIK( Pi, Pj, Pk, H );
        localToGlobal( tripletLists, numEdges + edgeIdx, pi, pk, H );

        //jk & kj (Hkj = Hjk)
        getHessThetaJK( Pi, Pj, Pk, H );
        localToGlobal( tripletLists, numEdges + edgeIdx, pj, pk, H );

        //ll
        getHessThetaKK( Pj, Pi, Pl, H );
        localToGlobal( tripletLists, numEdges + edgeIdx, pl, pl, H );

        //il & li (Hli = Hil)
        getHessThetaJK( Pj, Pi, Pl, H );
        localToGlobal( tripletLists, numEdges + edgeIdx, pi, pl, H );

        //jl & lj (Hlj = Hjl)
        getHessThetaIK( Pj, Pi, Pl, H );
        localToGlobal( tripletLists, numEdges + edgeIdx, pj, pl, H );


        //ii
        getHessThetaII( Pi, Pj, Pk, Pl, H );
        localToGlobal( tripletLists, numEdges + edgeIdx, pi, pi, H );

        //jj
        getHessThetaII( Pj, Pi, Pl, Pk, H );
        localToGlobal( tripletLists, numEdges + edgeIdx, pj, pj, H );

        //ij & ji (Hij = Hji)
        getHessThetaJI( Pi, Pj, Pk, Pl, H );
        localToGlobal( tripletLists, numEdges + edgeIdx, pi, pj, H );
      }
    }

  }

  void pushTriplets( const VectorType &Arg, TripletListType tripletList, const VectorType &Lambda ) const {
    throw std::runtime_error( "NRICMapHessian::pushTriplets: Not implemented!" );
  }


  int getTargetDimension() const override {
    if ( _consideredIndices.size() > 0 ) {
      return _consideredIndices.size();
    }
    else {
      return 2 * _topology.getNumEdges();
    }
  }

protected:
  void pushTriplet( std::vector<TripletListType> &tripletLists, int i, int row, int col, RealType value,
                    bool transposed = false ) const {
    if ( transposed )
      tripletLists[i].emplace_back( col, row, value );
    else
      tripletLists[i].emplace_back( row, col, value );
  }

  void localToGlobal( std::vector<TripletListType> &tripletLists, int h, int k, int l,
                      const MatType &localMatrix ) const {
    const int numVertices = _topology.getNumVertices(); // colOffset

    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ )
        pushTriplet( tripletLists, h, i * numVertices + k, j * numVertices + l, localMatrix( i, j ));

    if ( k != l ) {
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          pushTriplet( tripletLists, h, i * numVertices + l, j * numVertices + k, localMatrix( j, i ));
    }
  }
};

#endif //NRIC_MAP_H
