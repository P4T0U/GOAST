// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Quadratic deformation energies for the space of lengths and angles
 * \author Sassen
 */

#ifndef NRIC_QUADRATICDEFORMATIONS_H
#define NRIC_QUADRATICDEFORMATIONS_H

#include <goast/Core/Topology.h>
#include <goast/Core/DeformationInterface.h>

/**
 * \brief Quadratic deformation energy on lengths and angles with constant weights
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 * \note This is nothing but a weighted Euclidean product
 */
template<typename ConfiguratorType>
class QuadraticDeformation : public DeformationBase<ConfiguratorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;

  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver &_topology;
  RealType _bendWeight, _memWeight;
  VectorType _weights;
  const int _numEdges;

public:
  /**
   * \brief Construct energy with membrane weight 1
   * \param Topology class containing information on the topology (i.e. connectivity) of the mesh.
   * \param Weights the weights as vector
   * \param bendWeight weight of the bending contribution, i.e. factor by which the second half of the weights are multiplied
   */
  QuadraticDeformation( const MeshTopologySaver &Topology, const VectorType &Weights, RealType bendWeight )
          : _topology( Topology ), _bendWeight( bendWeight ), _memWeight( 1. ), _weights( Weights ),
            _numEdges( Topology.getNumEdges()) {
    if ( Weights.size() != 2 * _numEdges )
      throw std::length_error( "QuadraticDeformation: Number of weights is not twice the number of edges!" );

    _weights.segment( 0, _numEdges ) *= _memWeight;
    _weights.segment( _numEdges, _numEdges ) *= _bendWeight;
  }

  /**
   * \brief Construct energy
   * \param Topology class containing information on the topology (i.e. connectivity) of the mesh.
   * \param Weights the weights as vector
   * \param memWeight weight of the membrane contribution, i.e. factor by which the first half of the weights are multiplied
   * \param bendWeight weight of the bending contribution, i.e. factor by which the second half of the weights are multiplied
   */
  QuadraticDeformation( const MeshTopologySaver &Topology, const VectorType &Weights, RealType memWeight,
                        RealType bendWeight ) : _topology( Topology ), _bendWeight( bendWeight ),
                                                _memWeight( memWeight ), _weights( Weights ),
                                                _numEdges( Topology.getNumEdges()) {
    if ( Weights.size() != 2 * _numEdges )
      throw std::length_error( "QuadraticDeformation: Number of weights is not twice the number of edges!" );

    _weights.segment( 0, _numEdges ) *= _memWeight;
    _weights.segment( _numEdges, _numEdges ) *= _bendWeight;
  }

  /**
   * \brief Evaluate the energy
   * \param[in] UndeformedGeom lengths and angles of the undeformed geometry
   * \param[in] DeformedGeom lengths and angles of the deformed geometry
   * \param Dest energy value as scalar
   */
  void applyEnergy( const VectorType &UndeformedGeom, const VectorType &DeformedGeom, RealType &Dest ) const {
    if ( UndeformedGeom.size() != 2 * _numEdges || DeformedGeom.size() != 2 * _numEdges )
      throw std::length_error( "QuadraticDeformation::applyEnergy: Number of variables don't match topology!" );

    VectorType diff = UndeformedGeom - DeformedGeom;

    Dest = diff.transpose() * _weights.asDiagonal() * diff;
  }

  /**
   * \brief This constructs the gradient w.r.t. the undeformed geometry, i.e. D_1 W[.,.]
   * \param[in] UndeformedGeom lengths and angles of the undeformed geometry
   * \param[in] DeformedGeom lengths and angles of the deformed geometry
   * \param[out] Dest vector which will contain the gradient afterwards
   */
  void applyUndefGradient( const VectorType &UndeformedGeom, const VectorType &DeformedGeom, VectorType &Dest ) const {
    if ( UndeformedGeom.size() != 2 * _numEdges || DeformedGeom.size() != 2 * _numEdges )
      throw std::length_error( "QuadraticDeformation::applyEnergy: Number of variables don't match topology!" );

    VectorType diff = UndeformedGeom - DeformedGeom;

    Dest.resize( 2 * _numEdges );
    Dest = 2. * _weights.cwiseProduct( diff );
  }

  /**
   * \brief This constructs the gradient w.r.t. the deformed geometry, i.e. D_2 W[.,.]
   * \param[in] UndeformedGeom lengths and angles of the undeformed geometry
   * \param[in] DeformedGeom lengths and angles of the deformed geometry
   * \param[out] Dest vector which will contain the gradient afterwards
   */
  void applyDefGradient( const VectorType &UndeformedGeom, const VectorType &DeformedGeom, VectorType &Dest ) const {
    if ( UndeformedGeom.size() != 2 * _numEdges || DeformedGeom.size() != 2 * _numEdges )
      throw std::length_error( "QuadraticDeformation::applyEnergy: Number of variables don't match topology!" );

    VectorType diff = UndeformedGeom - DeformedGeom;

    Dest.resize( 2 * _numEdges );
    Dest = -2. * _weights.cwiseProduct( diff );
  }

  /**
   * \brief This assembles the Hessian w.r.t. the deformed geometry, i.e. D_2 D_2 W[.,.], as triplet list
   * \param[in] UndeformedGeom lengths and angles of the undeformed geometry
   * \param[in] DeformedGeom lengths and angles of the deformed geometry
   * \param[out] triplets list of triplets to which the Hessian entries will be added
   * \param[in] rowOffset offset to add to the row indices
   * \param[in] colOffset offset to add to hte column indices
   * \param[in] factor the scaling factor of the Hessian
   */
  void pushTripletsDefHessian( const VectorType &UndeformedGeom, const VectorType &DeformedGeom,
                               TripletListType &triplets, int rowOffset, int colOffset, RealType factor = 1.0 ) const {
    for ( int i = 0; i < 2 * _numEdges; i++ )
      triplets.emplace_back( i + rowOffset, i + colOffset, 2. * factor * _weights[i] );
  }

  /**
   * \brief This assembles the Hessian w.r.t. the undeformed geometry, i.e. D_1 D_1 W[.,.], as triplet list
   * \param[in] UndeformedGeom lengths and angles of the undeformed geometry
   * \param[in] DeformedGeom lengths and angles of the deformed geometry
   * \param[out] triplets list of triplets to which the Hessian entries will be added
   * \param[in] rowOffset offset to add to the row indices
   * \param[in] colOffset offset to add to hte column indices
   * \param[in] factor the scaling factor of the Hessian
   */
  void pushTripletsUndefHessian( const VectorType &UndeformedGeom, const VectorType &DeformedGeom,
                                 TripletListType &triplets, int rowOffset, int colOffset,
                                 RealType factor = 1.0 ) const {
    for ( int i = 0; i < 2 * _numEdges; i++ )
      triplets.emplace_back( i + rowOffset, i + colOffset, 2. * factor * _weights[i] );
  }

  /**
   * \brief This assembles a mixed Hessian (i.e. D_1 D_2 W[.,.] or D_2 D_1 W[.,.]) as triplet list
   * \param[in] UndeformedGeom lengths and angles of the undeformed geometry
   * \param[in] DeformedGeom lengths and angles of the deformed geometry
   * \param[out] triplets list of triplets to which the Hessian entries will be added
   * \param[in] rowOffset offset to add to the row indices
   * \param[in] colOffset offset to add to hte column indices
   * \param[in] FirstDerivWRTDef whether to assemble D_1 D_2 W[.,.] (true) or D_2 D_1 W[.,.] (false)
   * \param[in] factor the scaling factor of the Hessian
   */
  void pushTripletsMixedHessian( const VectorType &UndeformedGeom, const VectorType &DeformedGeom,
                                 TripletListType &triplets, int rowOffset, int colOffset, const bool FirstDerivWRTDef,
                                 RealType factor = 1.0 ) const {
    for ( int i = 0; i < 2 * _numEdges; i++ )
      triplets.emplace_back( i + rowOffset, i + colOffset, -2. * factor * _weights[i] );
  }

  /**
   * \return the number of nonzero entries of the sparse hessian (twice the number of edges)
   */
  int numOfNonZeroHessianEntries() const {
    return 2 * _numEdges;
  }

  /**
   * \return the bending weight of the energy
   */
  RealType getBendingWeight() const {
    return _bendWeight;
  }
};


#endif //NRIC_QUADRATICDEFORMATIONS_H
