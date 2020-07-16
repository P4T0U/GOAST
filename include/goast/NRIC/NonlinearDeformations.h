// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Nonlinear deformation energies for the space of lengths and angles
 * \author Sassen
 *
 * See also
 * Sassen, J., Heeren, B., Hildebrandt, K., & Rumpf, M. (2020). Geometric optimization using nonlinear
 * rotation-invariant coordinates. Computer Aided Geometric Design, 77, 101829.
 * for a detailed formulation.
 */

#ifndef NRIC_NONLINEARDEFORMATIONS_H
#define NRIC_NONLINEARDEFORMATIONS_H

#include <goast/Core/Topology.h>
#include <goast/Core/DeformationInterface.h>
#include "TrigonometryOperators.h"


/**
 * \brief Bending deformation energy on lengths and angles with varying weights
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 *
 * This energy is taken from  Grinspun et al. \cite GrHiDeSc03
 *
 * Let \f$ z, \tilde z \f$ be the NRIC of two meshes that are in dense correspondence, let \f $ E \f$ be the set of edges.
 * Then this class realizes the energy
 * \f[ E[z, \tilde z] = |sum_{e \in E}  \frac{ (\theta_e[x] - \theta_e[\tilde x])^2 }{d_e[x]} l_e^2[x]\, , \f]
 * where \f$ \theta_e \f$ is the dihedral angle at the edge, \f$ l_e \f$ the length of the edge
 * and \f$ d_e = \frac13 (a_1 + a_2) \f$, if \f$ a_1 \f$ and \f$ a_2 \f$ denote the face area of the two adjacent
 * triangles, respectively, which are computed from edge lengths via Heron's formula.
 */
template<typename ConfiguratorType>
class NRICBendingDeformation : public DeformationBase<ConfiguratorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;

  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver &_topology;
  RealType _bendWeight;
  const int _numEdges;
  const int _numFaces;

  // Topological information
  std::vector<std::array<int, 4>> auxiliaryEdges;
  std::vector<int> interiorEdges;

  // Trigonometric operators
  TriangleAreaOp<ConfiguratorType> _A;

public:

  /**
   * \brief Construct energy
   * \param Topology class containing information on the topology (i.e. connectivity) of the mesh.
   * \param bendWeight weight of the bending contribution
   */
  NRICBendingDeformation( const MeshTopologySaver &Topology, RealType bendWeight )
          : _topology( Topology ), _bendWeight( bendWeight ), _numEdges( Topology.getNumEdges()),
            _numFaces( Topology.getNumFaces()), _A( Topology ) {

    auxiliaryEdges.resize( _numEdges, std::array<int, 4>{{ -1, -1, -1, -1 }} );
    // For each edge, store the indices of edges of neighboring triangles
    for ( int edgeIdx = 0; edgeIdx < _numEdges; edgeIdx++ ) {
      int f_r = _topology.getAdjacentTriangleOfEdge( edgeIdx, 0 );
      int f_l = _topology.getAdjacentTriangleOfEdge( edgeIdx, 1 );

      if ( f_r == -1 || f_l == -1 )
        continue;


      interiorEdges.emplace_back( edgeIdx );

      std::array<int, 3> e_r = { _topology.getEdgeOfTriangle( f_r, 0 ),
                                 _topology.getEdgeOfTriangle( f_r, 1 ),
                                 _topology.getEdgeOfTriangle( f_r, 2 ) };
      std::array<int, 3> e_l = { _topology.getEdgeOfTriangle( f_l, 0 ),
                                 _topology.getEdgeOfTriangle( f_l, 1 ),
                                 _topology.getEdgeOfTriangle( f_l, 2 ) };

      std::array<int, 4> auxEdges( { -1, -1, -1, -1 } ); // r_1, r_2, l_1, l_2
      int i = 0;
      for ( const auto &e : e_r ) {
        if ( e == edgeIdx ) continue;
        auxEdges[i] = e;
        i++;
      }
      for ( const auto &e : e_l ) {
        if ( e == edgeIdx ) continue;
        auxEdges[i] = e;
        i++;
      }

      auxiliaryEdges[edgeIdx] = auxEdges;
    }
  }

  /**
   * \brief Evaluate the energy
   * \param[in] UndeformedGeom lengths and angles of the undeformed geometry
   * \param[in] DeformedGeom lengths and angles of the deformed geometry
   * \param Dest energy value as scalar
   */
  void applyEnergy( const VectorType &UndeformedGeom, const VectorType &DeformedGeom, RealType &Dest ) const {
    if ( UndeformedGeom.size() != 2 * _numEdges || DeformedGeom.size() != 2 * _numEdges )
      throw std::length_error( "NRICBendingDeformation::applyEnergy: Number of variables don't match topology!" );

    Dest = 0.;

    VectorType triangleAreas( _numFaces );
    _A.apply( UndeformedGeom, triangleAreas );

    for ( const int &edgeIdx : interiorEdges ) {
      int f_r = _topology.getAdjacentTriangleOfEdge( edgeIdx, 0 );
      int f_l = _topology.getAdjacentTriangleOfEdge( edgeIdx, 1 );

      const RealType &theta = UndeformedGeom[_numEdges + edgeIdx];
      const RealType &theta_def = DeformedGeom[_numEdges + edgeIdx];
      const RealType &l_e = UndeformedGeom[edgeIdx];
      const RealType &r_1 = UndeformedGeom[auxiliaryEdges[edgeIdx][0]];
      const RealType &r_2 = UndeformedGeom[auxiliaryEdges[edgeIdx][1]];
      const RealType &l_1 = UndeformedGeom[auxiliaryEdges[edgeIdx][2]];
      const RealType &l_2 = UndeformedGeom[auxiliaryEdges[edgeIdx][3]];

      const RealType &a_r = triangleAreas[f_r];
      const RealType &a_l = triangleAreas[f_l];

      const RealType d_e = ( a_r + a_l ) / 3.;

      const RealType thetaDiff = theta - theta_def;

      Dest += ( thetaDiff * thetaDiff ) * ( l_e * l_e ) / d_e;
    }

    Dest *= _bendWeight;
  }

  /**
   * \brief This constructs the gradient w.r.t. the undeformed geometry, i.e. D_1 W[.,.]
   * \param[in] UndeformedGeom lengths and angles of the undeformed geometry
   * \param[in] DeformedGeom lengths and angles of the deformed geometry
   * \param[out] Dest vector which will contain the gradient afterwards
   */
  void applyUndefGradient( const VectorType &UndeformedGeom, const VectorType &DeformedGeom, VectorType &Dest ) const {
    if ( UndeformedGeom.size() != 2 * _numEdges || DeformedGeom.size() != 2 * _numEdges )
      throw std::length_error( "NRICBendingDeformation::applyEnergy: Number of variables don't match topology!" );

    Dest.resize( 2 * _numEdges );
    Dest.setZero();

    VectorType triangleAreas( _numFaces );
    _A.apply( UndeformedGeom, triangleAreas );

    for ( const int &edgeIdx : interiorEdges ) {
      int f_r = _topology.getAdjacentTriangleOfEdge( edgeIdx, 0 );
      int f_l = _topology.getAdjacentTriangleOfEdge( edgeIdx, 1 );

      const RealType &theta = UndeformedGeom[_numEdges + edgeIdx];
      const RealType &theta_def = DeformedGeom[_numEdges + edgeIdx];
      const RealType &l_e = UndeformedGeom[edgeIdx];
      const RealType &r_1 = UndeformedGeom[auxiliaryEdges[edgeIdx][0]];
      const RealType &r_2 = UndeformedGeom[auxiliaryEdges[edgeIdx][1]];
      const RealType &l_1 = UndeformedGeom[auxiliaryEdges[edgeIdx][2]];
      const RealType &l_2 = UndeformedGeom[auxiliaryEdges[edgeIdx][3]];

      const RealType &a_r = triangleAreas[f_r];
      const RealType &a_l = triangleAreas[f_l];

      const RealType l_eSqr = l_e * l_e;
      const RealType r_1Sqr = r_1 * r_1;
      const RealType r_2Sqr = r_2 * r_2;
      const RealType l_1Sqr = l_1 * l_1;
      const RealType l_2Sqr = l_2 * l_2;

      const RealType d_e = ( a_r + a_l ) / 3.;
      const RealType d_eSqr = d_e * d_e;

      const RealType thetaDiff = theta - theta_def;
      const RealType thetaDiffSqr = thetaDiff * thetaDiff;

      Dest[_numEdges + edgeIdx] += 2 * l_eSqr * thetaDiff / d_e;
      Dest[edgeIdx] += l_e * thetaDiffSqr / ( 6 * d_eSqr ) *
                       ( -l_eSqr *
                         (( -l_eSqr + l_1Sqr + l_2Sqr ) / ( 4 * a_l ) + ( -l_eSqr + r_1Sqr + r_2Sqr ) / ( 4 * a_r )) +
                         12 * d_e );
      Dest[auxiliaryEdges[edgeIdx][0]] += -( l_eSqr * r_1 * thetaDiffSqr * ( l_eSqr - r_1Sqr + r_2Sqr )) /
                                          ( 24 * d_eSqr * a_r );
      Dest[auxiliaryEdges[edgeIdx][1]] += -( l_eSqr * r_2 * thetaDiffSqr * ( l_eSqr + r_1Sqr - r_2Sqr )) /
                                          ( 24 * d_eSqr * a_r );
      Dest[auxiliaryEdges[edgeIdx][2]] += -( l_eSqr * l_1 * thetaDiffSqr * ( l_eSqr - l_1Sqr + l_2Sqr )) /
                                          ( 24 * d_eSqr * a_l );
      Dest[auxiliaryEdges[edgeIdx][3]] += -( l_eSqr * l_2 * thetaDiffSqr * ( l_eSqr + l_1Sqr - l_2Sqr )) /
                                          ( 24 * d_eSqr * a_l );
    }

    Dest.array() *= _bendWeight;

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

    Dest.resize( 2 * _numEdges );
    Dest.setZero();

    VectorType triangleAreas( _numFaces );
    _A.apply( UndeformedGeom, triangleAreas );

    for ( const int &edgeIdx : interiorEdges ) {
      int f_r = _topology.getAdjacentTriangleOfEdge( edgeIdx, 0 );
      int f_l = _topology.getAdjacentTriangleOfEdge( edgeIdx, 1 );

      const RealType &theta = UndeformedGeom[_numEdges + edgeIdx];
      const RealType &theta_def = DeformedGeom[_numEdges + edgeIdx];
      const RealType &l_e = UndeformedGeom[edgeIdx];
      const RealType &r_1 = UndeformedGeom[auxiliaryEdges[edgeIdx][0]];
      const RealType &r_2 = UndeformedGeom[auxiliaryEdges[edgeIdx][1]];
      const RealType &l_1 = UndeformedGeom[auxiliaryEdges[edgeIdx][2]];
      const RealType &l_2 = UndeformedGeom[auxiliaryEdges[edgeIdx][3]];

      const RealType &a_r = triangleAreas[f_r];
      const RealType &a_l = triangleAreas[f_l];

      const RealType d_e = ( a_r + a_l ) / 3.;

      const RealType thetaDiff = theta - theta_def;

      Dest[_numEdges + edgeIdx] = -2. * ( l_e * l_e ) * thetaDiff / d_e;
    }

    Dest.array() *= _bendWeight;
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
    VectorType triangleAreas( _numFaces );
    _A.apply( UndeformedGeom, triangleAreas );

    for ( const int &edgeIdx : interiorEdges ) {
      int f_r = _topology.getAdjacentTriangleOfEdge( edgeIdx, 0 );
      int f_l = _topology.getAdjacentTriangleOfEdge( edgeIdx, 1 );

      const RealType &theta = UndeformedGeom[_numEdges + edgeIdx];
      const RealType &theta_def = DeformedGeom[_numEdges + edgeIdx];
      const RealType &l_e = UndeformedGeom[edgeIdx];
      const RealType &r_1 = UndeformedGeom[auxiliaryEdges[edgeIdx][0]];
      const RealType &r_2 = UndeformedGeom[auxiliaryEdges[edgeIdx][1]];
      const RealType &l_1 = UndeformedGeom[auxiliaryEdges[edgeIdx][2]];
      const RealType &l_2 = UndeformedGeom[auxiliaryEdges[edgeIdx][3]];

      const RealType &a_r = triangleAreas[f_r];
      const RealType &a_l = triangleAreas[f_l];

      const RealType d_e = ( a_r + a_l ) / 3.;

      triplets.emplace_back( _numEdges + edgeIdx + rowOffset, _numEdges + edgeIdx + colOffset,
                             _bendWeight * factor * 2. * ( l_e * l_e ) / d_e );
    }
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
    VectorType triangleAreas( _numFaces );
    _A.apply( UndeformedGeom, triangleAreas );

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for ( int ivIdx = 0; ivIdx < interiorEdges.size(); ivIdx++ ) {
      const int &edgeIdx = interiorEdges[ivIdx];
      TripletListType localTriplets;

      int f_r = _topology.getAdjacentTriangleOfEdge( edgeIdx, 0 );
      int f_l = _topology.getAdjacentTriangleOfEdge( edgeIdx, 1 );

      const RealType &theta = UndeformedGeom[_numEdges + edgeIdx];
      const RealType &theta_def = DeformedGeom[_numEdges + edgeIdx];
      const RealType &l_e = UndeformedGeom[edgeIdx];
      const RealType &r_1 = UndeformedGeom[auxiliaryEdges[edgeIdx][0]];
      const RealType &r_2 = UndeformedGeom[auxiliaryEdges[edgeIdx][1]];
      const RealType &l_1 = UndeformedGeom[auxiliaryEdges[edgeIdx][2]];
      const RealType &l_2 = UndeformedGeom[auxiliaryEdges[edgeIdx][3]];

      const RealType &a_r = triangleAreas[f_r];
      const RealType &a_l = triangleAreas[f_l];

      // Often used powers
      const RealType l_eSqr = l_e * l_e;
      const RealType l_ePow4 = l_eSqr * l_eSqr;
      const RealType l_ePow6 = l_ePow4 * l_eSqr;
      const RealType l_ePow8 = l_ePow6 * l_eSqr;
      const RealType l_ePow10 = l_ePow8 * l_eSqr;

      const RealType r_1Sqr = r_1 * r_1;
      const RealType r_1Pow4 = r_1Sqr * r_1Sqr;

      const RealType r_2Sqr = r_2 * r_2;
      const RealType r_2Pow4 = r_2Sqr * r_2Sqr;

      const RealType l_1Sqr = l_1 * l_1;
      const RealType l_1Pow4 = l_1Sqr * l_1Sqr;

      const RealType l_2Sqr = l_2 * l_2;
      const RealType l_2Pow4 = l_2Sqr * l_2Sqr;

      const RealType a_lCub = a_l * a_l * a_l;
      const RealType a_rCub = a_r * a_r * a_r;

      // Edge-associated area
      const RealType d_e = ( a_r + a_l ) / 3.;
      const RealType d_eSqr = d_e * d_e;
      const RealType d_eCub12 = 12. * 12. * 12. * d_e * d_e * d_e;

      // Angle difference
      const RealType thetaDiff = theta - theta_def;
      const RealType thetaDiffSqr = thetaDiff * thetaDiff;

      // Helper Polynomials
      const RealType X_1 = l_1Sqr - l_2Sqr;
      const RealType X_1Sqr = X_1 * X_1;
      const RealType X_2 = l_1Sqr + l_2Sqr;
      const RealType X_3 = r_1Sqr - r_2Sqr;
      const RealType X_3Sqr = X_3 * X_3;
      const RealType X_4 = r_1Sqr + r_2Sqr;
      const RealType X_5 = r_1Sqr - 3 * r_2Sqr;
      const RealType X_6 = r_2Sqr - 3 * r_1Sqr;

      const RealType Y_1 = l_eSqr + X_1;
      const RealType Y_2 = l_eSqr + X_3;
      const RealType Y_3 = l_eSqr - X_1;
      const RealType Y_4 = l_eSqr - X_3;
      const RealType Y_5 = -l_eSqr + X_4;
      const RealType Y_6 = -l_eSqr + X_2;

      const RealType Z_1 = 16. * a_l * a_r;
      const RealType Z_2 = 4. * a_l;
      const RealType Z_3 = 4. * a_r;

      const RealType Z_4 = ( Y_6 / Z_2 + Y_5 / Z_3 );

      // temp variables
      RealType prefactor, term1, term2, term3, term4, term5, term6;

      // D_x D_theta
      localTriplets.emplace_back( _numEdges + edgeIdx + rowOffset, _numEdges + edgeIdx + colOffset,
                                  _bendWeight * factor * 2. * l_eSqr / d_e );

      localTriplets.emplace_back( edgeIdx + rowOffset, _numEdges + edgeIdx + colOffset,
                                  _bendWeight * factor * l_e * thetaDiff / ( 3. * d_eSqr ) *
                                  ( -l_eSqr * Z_4 + 12. * d_e ));
      localTriplets.emplace_back( _numEdges + edgeIdx + rowOffset, edgeIdx + colOffset,
                                  _bendWeight * factor * l_e * thetaDiff / ( 3. * d_eSqr ) *
                                  ( -l_eSqr * Z_4 + 12. * d_e ));

      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][0] + rowOffset, _numEdges + edgeIdx + colOffset,
                                  _bendWeight * factor * -l_eSqr * r_1 * thetaDiff * Y_4 / ( 12 * a_r * d_eSqr ));
      localTriplets.emplace_back( _numEdges + edgeIdx + rowOffset, auxiliaryEdges[edgeIdx][0] + colOffset,
                                  _bendWeight * factor * -l_eSqr * r_1 * thetaDiff * Y_4 / ( 12 * a_r * d_eSqr ));

      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][1] + rowOffset, _numEdges + edgeIdx + colOffset,
                                  _bendWeight * factor * -l_eSqr * r_2 * thetaDiff * Y_2 / ( 12 * a_r * d_eSqr ));
      localTriplets.emplace_back( _numEdges + edgeIdx + rowOffset, auxiliaryEdges[edgeIdx][1] + colOffset,
                                  _bendWeight * factor * -l_eSqr * r_2 * thetaDiff * Y_2 / ( 12 * a_r * d_eSqr ));

      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][2] + rowOffset, _numEdges + edgeIdx + colOffset,
                                  _bendWeight * factor * -l_eSqr * l_1 * thetaDiff * Y_3 / ( 12 * a_l * d_eSqr ));
      localTriplets.emplace_back( _numEdges + edgeIdx + rowOffset, auxiliaryEdges[edgeIdx][2] + colOffset,
                                  _bendWeight * factor * -l_eSqr * l_1 * thetaDiff * Y_3 / ( 12 * a_l * d_eSqr ));

      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][3] + rowOffset, _numEdges + edgeIdx + colOffset,
                                  _bendWeight * factor * -l_eSqr * l_2 * thetaDiff * Y_1 / ( 12 * a_l * d_eSqr ));
      localTriplets.emplace_back( _numEdges + edgeIdx + rowOffset, auxiliaryEdges[edgeIdx][3] + colOffset,
                                  _bendWeight * factor * -l_eSqr * l_2 * thetaDiff * Y_1 / ( 12 * a_l * d_eSqr ));


      // D_l D_l
      prefactor = 24 * thetaDiffSqr / d_eCub12;
      term1 = -1. * 48. * l_eSqr * d_e * Z_4;
      term2 = -1. * 12. * d_e * l_eSqr * ( -2. * l_eSqr * Y_6 * Y_6 / ( 64. * a_lCub ) + ( -3. * l_eSqr + X_2 ) / Z_2
                                           - 2. * l_eSqr * Y_5 * Y_5 / ( 64. * a_rCub ) +
                                           ( -3. * l_eSqr + X_4 ) / Z_3 );
      term3 = 4. * l_ePow4 * Z_4 * Z_4;
      term4 = 12. * 12 * d_eSqr;

      localTriplets.emplace_back( edgeIdx + rowOffset, edgeIdx + colOffset,
                                  _bendWeight * factor * prefactor * ( term1 + term2 + term3 + term4 ));

      // D_r1 D_l
      prefactor = ( -48. * l_e * r_1 * thetaDiffSqr ) / ( 4. * a_l * 64. * a_rCub * d_eCub12 );
      term1 = -1. * l_ePow10;
      term2 = 1. * l_ePow8 * ( 3 * r_1Sqr - r_2Sqr );
      term3 = l_ePow6 * ( Z_1 - 3 * r_1Pow4 + 3 * r_2Pow4 + 4 * r_2Sqr * X_2 + X_1Sqr );

      term4 = -1. * l_ePow4 * ( r_1Sqr * ( Z_1 - 3 * r_2Pow4 - 4 * r_2Sqr * X_2 + 3 * X_1Sqr )
                                + r_2Sqr * ( -3. * Z_1 + r_2Pow4 + 4 * r_2Sqr * X_2 + 3 * X_1Sqr )
                                - r_1Pow4 * r_1Sqr + 3 * r_1Pow4 * r_2Sqr );
      term5 = -1. * l_eSqr * X_3 * ( r_1Sqr * ( Z_1 - 3 * X_1Sqr ) + r_2Sqr * ( X_1Sqr - 3 * Z_1 ));
      term6 = X_3Sqr * X_3 * ( Z_1 - X_1Sqr );

      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][0] + rowOffset, edgeIdx + colOffset,
                                  _bendWeight * factor * prefactor * ( term1 + term2 + term3 + term4 + term5 + term6 ));
      localTriplets.emplace_back( edgeIdx + rowOffset, auxiliaryEdges[edgeIdx][0] + colOffset,
                                  _bendWeight * factor * prefactor * ( term1 + term2 + term3 + term4 + term5 + term6 ));

      // D_r2 D_l
      prefactor = ( -48. * l_e * r_2 * thetaDiffSqr ) / ( 4. * a_l * 64. * a_rCub * d_eCub12 );
      term1 = -1. * l_ePow10;
      term2 = 1. * l_ePow8 * ( 3 * r_2Sqr - r_1Sqr );
      term3 = l_ePow6 * ( Z_1 - 3 * r_2Pow4 + 3 * r_1Pow4 + 4 * r_1Sqr * X_2 + X_1Sqr );

      term4 = l_ePow4 * ( r_2Sqr * ( -Z_1 + r_2Pow4 - 3 * X_1Sqr )
                          + r_1Sqr * ( 3 * Z_1 - 3 * r_2Pow4 + 4 * r_2Sqr * X_2 - 3 * X_1Sqr )
                          - r_1Pow4 * r_1Sqr + 3 * r_1Pow4 * r_2Sqr - 4 * r_1Pow4 * X_2 );
      term5 = -1. * l_eSqr * X_3 * ( r_2Sqr * ( 3 * X_1Sqr - Z_1 ) + r_1Sqr * ( 3 * Z_1 - X_1Sqr ));
      term6 = -1. * X_3Sqr * X_3 * ( Z_1 - X_1Sqr );

      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][1] + rowOffset, edgeIdx + colOffset,
                                  _bendWeight * factor * prefactor * ( term1 + term2 + term3 + term4 + term5 + term6 ));
      localTriplets.emplace_back( edgeIdx + rowOffset, auxiliaryEdges[edgeIdx][1] + colOffset,
                                  _bendWeight * factor * prefactor * ( term1 + term2 + term3 + term4 + term5 + term6 ));

      // D_l1 D_l
      prefactor = ( -48. * l_e * l_1 * thetaDiffSqr ) / ( 4. * a_r * 64. * a_lCub * d_eCub12 );
      term1 = -1. * l_ePow10;
      term2 = 1. * l_ePow8 * ( 3 * l_1Sqr - l_2Sqr );
      term3 = l_ePow6 * ( Z_1 - 3 * l_1Pow4 + 3 * l_2Pow4 + 4 * l_2Sqr * X_4 + X_3Sqr );

      term4 = -1. * l_ePow4 * ( l_1Sqr * ( Z_1 - 3 * l_2Pow4 - 4 * l_2Sqr * X_4 + 3 * X_3Sqr )
                                + l_2Sqr * ( -3. * Z_1 + l_2Pow4 + 4 * l_2Sqr * X_4 + 3 * X_3Sqr )
                                - l_1Pow4 * l_1Sqr + 3 * l_1Pow4 * l_2Sqr );
      term5 = -1. * l_eSqr * X_1 * ( l_1Sqr * ( Z_1 - 3 * X_3Sqr ) + l_2Sqr * ( X_3Sqr - 3 * Z_1 ));
      term6 = X_1Sqr * X_1 * ( Z_1 - X_3Sqr );

      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][2] + rowOffset, edgeIdx + colOffset,
                                  _bendWeight * factor * prefactor * ( term1 + term2 + term3 + term4 + term5 + term6 ));
      localTriplets.emplace_back( edgeIdx + rowOffset, auxiliaryEdges[edgeIdx][2] + colOffset,
                                  _bendWeight * factor * prefactor * ( term1 + term2 + term3 + term4 + term5 + term6 ));

      // D_l2 D_l
      prefactor = ( -48. * l_e * l_2 * thetaDiffSqr ) / ( 4. * a_r * 64. * a_lCub * d_eCub12 );
      term1 = -1. * l_ePow10;
      term2 = 1. * l_ePow8 * ( 3 * l_2Sqr - l_1Sqr );
      term3 = l_ePow6 * ( Z_1 - 3 * l_2Pow4 + 3 * l_1Pow4 + 4 * l_1Sqr * X_4 + X_3Sqr );

      term4 = l_ePow4 * ( l_2Sqr * ( -Z_1 + l_2Pow4 - 3 * X_3Sqr )
                          + l_1Sqr * ( 3 * Z_1 - 3 * l_2Pow4 + 4 * l_2Sqr * X_4 - 3 * X_3Sqr )
                          - l_1Pow4 * l_1Sqr + 3 * l_1Pow4 * l_2Sqr - 4 * l_1Pow4 * X_4 );
      term5 = -1. * l_eSqr * X_1 * ( l_2Sqr * ( 3 * X_3Sqr - Z_1 ) + l_1Sqr * ( 3 * Z_1 - X_3Sqr ));
      term6 = -1. * X_1Sqr * X_1 * ( Z_1 - X_3Sqr );

      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][3] + rowOffset, edgeIdx + colOffset,
                                  _bendWeight * factor * prefactor * ( term1 + term2 + term3 + term4 + term5 + term6 ));
      localTriplets.emplace_back( edgeIdx + rowOffset, auxiliaryEdges[edgeIdx][3] + colOffset,
                                  _bendWeight * factor * prefactor * ( term1 + term2 + term3 + term4 + term5 + term6 ));

      // D_r1 D_r1
      prefactor = ( 6. * l_eSqr * thetaDiffSqr ) / d_eCub12;
      term1 = -1. * ( 12. * d_e * ( l_eSqr + X_6 )) / a_r;
      term2 = 1. * r_1Sqr * Y_4 * Y_4 / ( a_r * a_r );
      term3 = 1.5 * ( d_e * r_1Sqr * Y_4 * Y_4 ) / ( a_rCub );
      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][0] + rowOffset, auxiliaryEdges[edgeIdx][0] + colOffset,
                                  _bendWeight * factor * prefactor * ( term1 + term2 + term3 ));

      // D_r2 D_r1
      prefactor = ( -3. * l_eSqr * r_1 * r_2 * thetaDiffSqr ) / ( 2. * a_rCub * d_eCub12 );
      term1 = -12. * l_ePow4 * d_e - 4. * l_ePow4 * a_r + 12. * l_eSqr * X_4 * d_e + 4 * X_3Sqr * a_r;
      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][1] + rowOffset, auxiliaryEdges[edgeIdx][0] + colOffset,
                                  _bendWeight * factor * prefactor * term1 );
      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][0] + rowOffset, auxiliaryEdges[edgeIdx][1] + colOffset,
                                  _bendWeight * factor * prefactor * term1 );

      // D_l1 D_r1
      prefactor = ( 96. * l_eSqr * l_1 * r_1 * thetaDiffSqr ) / d_eCub12;
      term1 = Y_3 * Y_4 / Z_1;
      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][2] + rowOffset, auxiliaryEdges[edgeIdx][0] + colOffset,
                                  _bendWeight * factor * prefactor * term1 );
      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][0] + rowOffset, auxiliaryEdges[edgeIdx][2] + colOffset,
                                  _bendWeight * factor * prefactor * term1 );

      // D_l2 D_r1
      prefactor = ( 96. * l_eSqr * l_2 * r_1 * thetaDiffSqr ) / d_eCub12;
      term1 = Y_1 * Y_4 / Z_1;
      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][3] + rowOffset, auxiliaryEdges[edgeIdx][0] + colOffset,
                                  _bendWeight * factor * prefactor * term1 );
      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][0] + rowOffset, auxiliaryEdges[edgeIdx][3] + colOffset,
                                  _bendWeight * factor * prefactor * term1 );

      // D_r2 D_r2
      prefactor = ( 6. * l_eSqr * thetaDiffSqr ) / d_eCub12;
      term1 = -1. * ( 12. * d_e * ( l_eSqr + X_5 )) / a_r;
      term2 = 1. * r_2Sqr * Y_2 * Y_2 / ( a_r * a_r );
      term3 = 1.5 * ( d_e * r_2Sqr * Y_2 * Y_2 ) / ( a_rCub );
      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][1] + rowOffset, auxiliaryEdges[edgeIdx][1] + colOffset,
                                  _bendWeight * factor * prefactor * ( term1 + term2 + term3 ));

      // D_l1 D_r2
      prefactor = ( 96. * l_eSqr * l_1 * r_2 * thetaDiffSqr ) / d_eCub12;
      term1 = Y_3 * Y_2 / Z_1;
      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][2] + rowOffset, auxiliaryEdges[edgeIdx][1] + colOffset,
                                  _bendWeight * factor * prefactor * term1 );
      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][1] + rowOffset, auxiliaryEdges[edgeIdx][2] + colOffset,
                                  _bendWeight * factor * prefactor * term1 );

      // D_l2 D_r2
      prefactor = ( 96. * l_eSqr * l_2 * r_2 * thetaDiffSqr ) / d_eCub12;
      term1 = Y_1 * Y_2 / Z_1;
      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][3] + rowOffset, auxiliaryEdges[edgeIdx][1] + colOffset,
                                  _bendWeight * factor * prefactor * term1 );
      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][1] + rowOffset, auxiliaryEdges[edgeIdx][3] + colOffset,
                                  _bendWeight * factor * prefactor * term1 );

      // D_l1 D_l1
      prefactor = ( 6. * l_eSqr * thetaDiffSqr ) / d_eCub12;
      term1 = -1. * ( 12. * d_e * ( l_eSqr + l_2Sqr - 3. * l_1Sqr )) / a_l;
      term2 = 1. * l_1Sqr * Y_3 * Y_3 / ( a_l * a_l );
      term3 = 1.5 * ( d_e * l_1Sqr * Y_3 * Y_3 ) / ( a_lCub );
      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][2] + rowOffset, auxiliaryEdges[edgeIdx][2] + colOffset,
                                  _bendWeight * factor * prefactor * ( term1 + term2 + term3 ));

      // D_l2 D_l1
      prefactor = ( -3. * l_eSqr * l_1 * l_2 * thetaDiffSqr ) / ( 2. * a_lCub * d_eCub12 );
      term1 = -12. * l_ePow4 * d_e - 4. * l_ePow4 * a_l + 12. * l_eSqr * X_2 * d_e + 4 * X_1Sqr * a_l;
      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][3] + rowOffset, auxiliaryEdges[edgeIdx][2] + colOffset,
                                  _bendWeight * factor * prefactor * term1 );
      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][2] + rowOffset, auxiliaryEdges[edgeIdx][3] + colOffset,
                                  _bendWeight * factor * prefactor * term1 );

      // D_l2 D_l2
      prefactor = ( 6. * l_eSqr * thetaDiffSqr ) / d_eCub12;
      term1 = -1. * ( 12. * d_e * ( l_eSqr + l_1Sqr - 3. * l_2Sqr )) / a_l;
      term2 = 1. * l_2Sqr * Y_1 * Y_1 / ( a_l * a_l );
      term3 = 1.5 * ( d_e * l_2Sqr * Y_1 * Y_1 ) / ( a_lCub );
      localTriplets.emplace_back( auxiliaryEdges[edgeIdx][3] + rowOffset, auxiliaryEdges[edgeIdx][3] + colOffset,
                                  _bendWeight * factor * prefactor * ( term1 + term2 + term3 ));

#ifdef _OPENMP
#pragma omp critical
#endif
      triplets.insert( triplets.end(), localTriplets.begin(), localTriplets.end());
    }

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
    VectorType triangleAreas( _numFaces );
    _A.apply( UndeformedGeom, triangleAreas );

    for ( int ivIdx = 0; ivIdx < interiorEdges.size(); ivIdx++ ) {
      const int &edgeIdx = interiorEdges[ivIdx];

      int f_r = _topology.getAdjacentTriangleOfEdge( edgeIdx, 0 );
      int f_l = _topology.getAdjacentTriangleOfEdge( edgeIdx, 1 );

      const RealType &theta = UndeformedGeom[_numEdges + edgeIdx];
      const RealType &theta_def = DeformedGeom[_numEdges + edgeIdx];
      const RealType &l_e = UndeformedGeom[edgeIdx];
      const RealType &r_1 = UndeformedGeom[auxiliaryEdges[edgeIdx][0]];
      const RealType &r_2 = UndeformedGeom[auxiliaryEdges[edgeIdx][1]];
      const RealType &l_1 = UndeformedGeom[auxiliaryEdges[edgeIdx][2]];
      const RealType &l_2 = UndeformedGeom[auxiliaryEdges[edgeIdx][3]];

      const RealType &a_r = triangleAreas[f_r];
      const RealType &a_l = triangleAreas[f_l];

      // Often used powers
      const RealType l_eSqr = l_e * l_e;
      const RealType r_1Sqr = r_1 * r_1;
      const RealType r_2Sqr = r_2 * r_2;
      const RealType l_1Sqr = l_1 * l_1;
      const RealType l_2Sqr = l_2 * l_2;

      const RealType a_lCub = a_l * a_l * a_l;
      const RealType a_rCub = a_r * a_r * a_r;

      // Edge-associated area
      const RealType d_e = ( a_r + a_l ) / 3.;
      const RealType d_eSqr = d_e * d_e;

      // Angle difference
      const RealType thetaDiff = theta - theta_def;

      // Helper Polynomials
      const RealType X_1 = l_1Sqr - l_2Sqr;
      const RealType X_1Sqr = X_1 * X_1;
      const RealType X_2 = l_1Sqr + l_2Sqr;
      const RealType X_3 = r_1Sqr - r_2Sqr;
      const RealType X_3Sqr = X_3 * X_3;
      const RealType X_4 = r_1Sqr + r_2Sqr;

      const RealType Y_1 = l_eSqr + X_1;
      const RealType Y_2 = l_eSqr + X_3;
      const RealType Y_3 = l_eSqr - X_1;
      const RealType Y_4 = l_eSqr - X_3;
      const RealType Y_5 = -l_eSqr + X_4;
      const RealType Y_6 = -l_eSqr + X_2;

      const RealType Z_2 = 4. * a_l;
      const RealType Z_3 = 4. * a_r;

      const RealType Z_4 = ( Y_6 / Z_2 + Y_5 / Z_3 );

      // temp variables
      RealType prefactor, term1, term2, term3, term4, term5, term6;

      // D_x D_theta_tilde
      if ( FirstDerivWRTDef ) {
        triplets.emplace_back( _numEdges + edgeIdx + rowOffset, _numEdges + edgeIdx + colOffset,
                               _bendWeight * factor * -2. * l_eSqr / d_e );
        triplets.emplace_back( _numEdges + edgeIdx + rowOffset, edgeIdx + colOffset,
                               _bendWeight * factor * -l_e * thetaDiff / ( 3. * d_eSqr ) *
                               ( -l_eSqr * Z_4 + 12. * d_e ));
        triplets.emplace_back( _numEdges + edgeIdx + rowOffset, auxiliaryEdges[edgeIdx][0] + colOffset,
                               _bendWeight * factor * l_eSqr * r_1 * thetaDiff * Y_4 / ( 12 * a_r * d_eSqr ));
        triplets.emplace_back( _numEdges + edgeIdx + rowOffset, auxiliaryEdges[edgeIdx][1] + colOffset,
                               _bendWeight * factor * l_eSqr * r_2 * thetaDiff * Y_2 / ( 12 * a_r * d_eSqr ));
        triplets.emplace_back( _numEdges + edgeIdx + rowOffset, auxiliaryEdges[edgeIdx][2] + colOffset,
                               _bendWeight * factor * l_eSqr * l_1 * thetaDiff * Y_3 / ( 12 * a_l * d_eSqr ));
        triplets.emplace_back( _numEdges + edgeIdx + rowOffset, auxiliaryEdges[edgeIdx][3] + colOffset,
                               _bendWeight * factor * l_eSqr * l_2 * thetaDiff * Y_1 / ( 12 * a_l * d_eSqr ));
      }
      else {
        triplets.emplace_back( _numEdges + edgeIdx + rowOffset, _numEdges + edgeIdx + colOffset,
                               _bendWeight * factor * -2. * l_eSqr / d_e );
        triplets.emplace_back( edgeIdx + rowOffset, _numEdges + edgeIdx + colOffset,
                               _bendWeight * factor * -l_e * thetaDiff / ( 3. * d_eSqr ) *
                               ( -l_eSqr * Z_4 + 12. * d_e ));
        triplets.emplace_back( auxiliaryEdges[edgeIdx][0] + rowOffset, _numEdges + edgeIdx + colOffset,
                               _bendWeight * factor * l_eSqr * r_1 * thetaDiff * Y_4 / ( 12 * a_r * d_eSqr ));
        triplets.emplace_back( auxiliaryEdges[edgeIdx][1] + rowOffset, _numEdges + edgeIdx + colOffset,
                               _bendWeight * factor * l_eSqr * r_2 * thetaDiff * Y_2 / ( 12 * a_r * d_eSqr ));
        triplets.emplace_back( auxiliaryEdges[edgeIdx][2] + rowOffset, _numEdges + edgeIdx + colOffset,
                               _bendWeight * factor * l_eSqr * l_1 * thetaDiff * Y_3 / ( 12 * a_l * d_eSqr ));
        triplets.emplace_back( auxiliaryEdges[edgeIdx][3] + rowOffset, _numEdges + edgeIdx + colOffset,
                               _bendWeight * factor * l_eSqr * l_2 * thetaDiff * Y_1 / ( 12 * a_l * d_eSqr ));
      }
    }
  }

  /**
   * \return the number of nonzero entries of the sparse hessian (twice the number of edges)
   */
  int numOfNonZeroHessianEntries() const {
    return 49 * _numEdges; // for each edge we have 7 DOFs, 5 undeformed lengths, 1 undeformed angle, 1 deformed angle
  }

  /**
   * \return the bending weight of the energy
   */
  RealType getBendingWeight() const {
    return _bendWeight;
  }
};

/**
 * \brief Nonlinear membrane deformation energy on lengths and angles
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 *
 * See eq. (8) in \cite HeRuWa12.
 *
 * Let \f$ z, \tilde z \f$ be the NRIC of two meshes that are in dense correspondence, let \f $ E \f$ be the set of edges.
 * For a face \f$ f \in F \f$ we consider discrete first fundamental forms \f$ g_f \f$ and \f$ \tilde g_f \f$, respectively.
 * Furthermore, for \f$ f \in F \f$ we define the geometric distortion tensor as \f$ G_f = g_f^{-1}\tilde g_f \f$. Both
 * can be formulated completely in NRIC, cf. \cite SaHeHi20
 *
 * Then the hyperelastic energy is given as
 * \f[ E[z, \tilde z] = \sum_{f \in F} a_f W( \tr G_f[z], \det G_f[z]) \f]
 * where \f$ a_f \f$ is the face area and the hyperelastic energy density
 * \f[ W(a,d) = \frac\mu2 a + \frac\lambda4 d - (\frac\mu2 + \frac\lambda4) \log d - \mu -\frac\lambda4 \f]
 * for physical parameters \f$ \mu, \lambda > 0 \f$.
 */
template<typename ConfiguratorType>
class NRICMembraneDeformation : public DeformationBase<ConfiguratorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;

  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver &_topology;
  RealType _memWeight;
  const int _numEdges;
  const int _numFaces;

  const RealType mu = 1.;
  const RealType lam = 1.;

  const RealType muHalf;
  const RealType lambdaQuarter;

  // Trigonometric operators
  TriangleAreaOp<ConfiguratorType> _A;

public:

  /**
   * \brief Construct energy
   * \param Topology class containing information on the topology (i.e. connectivity) of the mesh.
   * \param bendWeight weight of the bending contribution
   */
  NRICMembraneDeformation( const MeshTopologySaver &Topology, RealType memWeight )
          : _topology( Topology ), _memWeight( memWeight ), _numEdges( Topology.getNumEdges()),
            _numFaces( Topology.getNumFaces()), _A( Topology ), muHalf( mu / 2. ), lambdaQuarter( lam / 4. ) {

  }

  /**
   * \brief Construct energy
   * \param Topology class containing information on the topology (i.e. connectivity) of the mesh.
   * \param bendWeight weight of the bending contribution
   * \param mu Material parameter mu
   * \param lambda Material paramter lambda
   */
  NRICMembraneDeformation( const MeshTopologySaver &Topology, RealType memWeight, RealType mu, RealType lambda )
          : _topology( Topology ), _memWeight( memWeight ), _numEdges( Topology.getNumEdges()),
            _numFaces( Topology.getNumFaces()), _A( Topology ), mu( mu ), lam( lambda ), muHalf( mu / 2. ),
            lambdaQuarter( lam / 4. ) {

  }

  /**
   * \brief Evaluate the energy
   * \param[in] UndeformedGeom lengths and angles of the undeformed geometry
   * \param[in] DeformedGeom lengths and angles of the deformed geometry
   * \param Dest energy value as scalar
   */
  void applyEnergy( const VectorType &UndeformedGeom, const VectorType &DeformedGeom, RealType &Dest ) const {
    if ( UndeformedGeom.size() != 2 * _numEdges || DeformedGeom.size() != 2 * _numEdges )
      throw std::length_error( "NRICMembraneDeformation::applyEnergy: Number of variables don't match topology!" );

    Dest = 0.;

    VectorType undefTriangleAreas( _numFaces );
    _A.apply( UndeformedGeom, undefTriangleAreas );

    VectorType defTriangleAreas( _numFaces );
    _A.apply( DeformedGeom, defTriangleAreas );

    for ( int faceIdx = 0; faceIdx < _numFaces; faceIdx++ ) {
      int e0 = _topology.getEdgeOfTriangle( faceIdx, 0 );
      int e1 = _topology.getEdgeOfTriangle( faceIdx, 1 );
      int e2 = _topology.getEdgeOfTriangle( faceIdx, 2 );

      const RealType &a = UndeformedGeom[e0];
      const RealType &b = UndeformedGeom[e1];
      const RealType &c = UndeformedGeom[e2];
      const RealType &a_t = DeformedGeom[e0];
      const RealType &b_t = DeformedGeom[e1];
      const RealType &c_t = DeformedGeom[e2];

      const RealType &a_u = undefTriangleAreas[faceIdx];
      const RealType &a_d = defTriangleAreas[faceIdx];

      const RealType cgamma_u = ( a * a + b * b - c * c ) / ( 2 * a * b );
      const RealType cgamma_d = ( a_t * a_t + b_t * b_t - c_t * c_t ) / ( 2 * a_t * b_t );

      const RealType detTerm = a_d * a_d / ( a_u * a_u );
      const RealType traceTerm = 1 / ( 4 * a_u * a_u ) * ( a * a * b_t * b_t -
                                                           2 * a * b * cgamma_u * a_t * b_t * cgamma_d +
                                                           a_t * a_t * b * b );

      Dest += a_u *
              (( mu / 2. ) * traceTerm + ( lam / 4. ) * detTerm - ( mu / 2. + lam / 4. ) * std::log( detTerm ) - mu -
               ( lam / 4. ));

    }
    Dest *= _memWeight;
  }

  /**
   * \brief This constructs the gradient w.r.t. the undeformed geometry, i.e. D_1 W[.,.]
   * \param[in] UndeformedGeom lengths and angles of the undeformed geometry
   * \param[in] DeformedGeom lengths and angles of the deformed geometry
   * \param[out] Dest vector which will contain the gradient afterwards
   */
  void applyUndefGradient( const VectorType &UndeformedGeom, const VectorType &DeformedGeom, VectorType &Dest ) const {
    if ( UndeformedGeom.size() != 2 * _numEdges || DeformedGeom.size() != 2 * _numEdges )
      throw std::length_error( "NRICMembraneDeformation::applyUndefGradient: "
                               "Number of variables doesn't match topology!" );

    Dest.resize( 2 * _numEdges );
    Dest.setZero();

    VectorType undefTriangleAreas( _numFaces );
    _A.apply( UndeformedGeom, undefTriangleAreas );

    VectorType defTriangleAreas( _numFaces );
    _A.apply( DeformedGeom, defTriangleAreas );

    for ( int faceIdx = 0; faceIdx < _numFaces; faceIdx++ ) {
      int e0 = _topology.getEdgeOfTriangle( faceIdx, 0 );
      int e1 = _topology.getEdgeOfTriangle( faceIdx, 1 );
      int e2 = _topology.getEdgeOfTriangle( faceIdx, 2 );

      // Undeformed Lengths
      const RealType &a = UndeformedGeom[e0];
      const RealType &b = UndeformedGeom[e1];
      const RealType &c = UndeformedGeom[e2];

      // Deformed Lengths
      const RealType &a_t = DeformedGeom[e0];
      const RealType &b_t = DeformedGeom[e1];
      const RealType &c_t = DeformedGeom[e2];

      // Squared Lengths
      const RealType aSqr = a * a;
      const RealType bSqr = b * b;
      const RealType cSqr = c * c;
      const RealType a_tSqr = a_t * a_t;
      const RealType b_tSqr = b_t * b_t;
      const RealType c_tSqr = c_t * c_t;

      // Areas
      const RealType &a_u = undefTriangleAreas[faceIdx];
      const RealType &a_d = defTriangleAreas[faceIdx];

      // Relative change of area = determinant term
      const RealType detTerm = a_d * a_d / ( a_u * a_u );
      const RealType logDetTerm = std::log( detTerm );

      // Helper polynomials / rational functions
      const RealType X_3 = ( -aSqr - bSqr + cSqr );
      const RealType X_2 = ( -a_tSqr - b_tSqr + c_tSqr );
      const RealType X_1 = aSqr * bSqr - X_3 * X_3 / 4;
      const RealType X_4 = 2 * a * bSqr + a * X_3;
      const RealType X_5 = 2 * aSqr * b + b * X_3;

      const RealType Z_1 = ( a_tSqr * bSqr - 0.5 * X_3 * X_2 + aSqr * b_tSqr ) / X_1;

      const RealType Y_1 = lambdaQuarter * detTerm - lambdaQuarter - muHalf;
      const RealType Y_2 = -0.5 * ( lambdaQuarter * detTerm - lambdaQuarter + muHalf * Z_1 - mu -
                                    ( lambdaQuarter + muHalf ) * logDetTerm );
      const RealType Y_3 = Y_1 + Y_2;


      // Da
      RealType summand1 = ( -1 / ( a + b - c ) - 1 / ( a - b + c ) + 1 / ( -a + b + c ) - 1 / ( a + b + c )) * Y_3;
      RealType summand2 = muHalf * ( a * X_2 + 2 * a * b_tSqr - X_4 * Z_1 ) / X_1;
      Dest[e0] += a_u * ( summand1 + summand2 );

      // Db
      summand1 = ( -1 / ( a + b - c ) + 1 / ( a - b + c ) - 1 / ( -a + b + c ) - 1 / ( a + b + c )) * Y_3;
      summand2 = muHalf * ( b * X_2 + 2 * b * a_tSqr - X_5 * Z_1 ) / X_1;
      Dest[e1] += a_u * ( summand1 + summand2 );

      // Dc
      summand1 = ( 1. / ( a + b - c ) - 1. / ( a - b + c ) - 1. / ( -a + b + c ) - 1. / ( a + b + c )) * Y_3;
      summand2 = muHalf * ( c * X_3 * Z_1 - c * X_2 ) / X_1;
      Dest[e2] += a_u * ( summand1 + summand2 );
    }

    Dest.array() *= _memWeight;

  }

  /**
   * \brief This constructs the gradient w.r.t. the deformed geometry, i.e. D_2 W[.,.]
   * \param[in] UndeformedGeom lengths and angles of the undeformed geometry
   * \param[in] DeformedGeom lengths and angles of the deformed geometry
   * \param[out] Dest vector which will contain the gradient afterwards
   */
  void applyDefGradient( const VectorType &UndeformedGeom, const VectorType &DeformedGeom, VectorType &Dest ) const {
    if ( UndeformedGeom.size() != 2 * _numEdges || DeformedGeom.size() != 2 * _numEdges )
      throw std::length_error( "NRICMembraneDeformation::applyDefGradient: "
                               "Number of variables doesn't match topology!" );

    Dest.resize( 2 * _numEdges );
    Dest.setZero();

    VectorType undefTriangleAreas( _numFaces );
    _A.apply( UndeformedGeom, undefTriangleAreas );

    VectorType defTriangleAreas( _numFaces );
    _A.apply( DeformedGeom, defTriangleAreas );

    for ( int faceIdx = 0; faceIdx < _numFaces; faceIdx++ ) {
      int e0 = _topology.getEdgeOfTriangle( faceIdx, 0 );
      int e1 = _topology.getEdgeOfTriangle( faceIdx, 1 );
      int e2 = _topology.getEdgeOfTriangle( faceIdx, 2 );

      // Undeformed Lengths
      const RealType &a = UndeformedGeom[e0];
      const RealType &b = UndeformedGeom[e1];
      const RealType &c = UndeformedGeom[e2];

      // Deformed Lengths
      const RealType &a_t = DeformedGeom[e0];
      const RealType &b_t = DeformedGeom[e1];
      const RealType &c_t = DeformedGeom[e2];

      // Squared Lengths
      const RealType aSqr = a * a;
      const RealType bSqr = b * b;
      const RealType cSqr = c * c;
      const RealType a_tSqr = a_t * a_t;
      const RealType b_tSqr = b_t * b_t;
      const RealType c_tSqr = c_t * c_t;

      // Areas
      const RealType &a_u = undefTriangleAreas[faceIdx];
      const RealType &a_d = defTriangleAreas[faceIdx];

      // Relative change of area = determinant term
      const RealType detTerm = a_d * a_d / ( a_u * a_u );
      const RealType logDetTerm = std::log( detTerm );

      // Helper polynomials / rational functions
      const RealType X_3 = ( -aSqr - bSqr + cSqr );
      const RealType X_2 = ( -a_tSqr - b_tSqr + c_tSqr );
      const RealType X_1 = aSqr * bSqr - X_3 * X_3 / 4;
      const RealType X_4 = 2 * a * bSqr + a * X_3;
      const RealType X_5 = 2 * aSqr * b + b * X_3;

      const RealType Z_1 = ( a_tSqr * bSqr - 0.5 * X_3 * X_2 + aSqr * b_tSqr ) / X_1;

      const RealType Y_1 = lambdaQuarter * detTerm - lambdaQuarter - muHalf;
      const RealType Y_2 = -0.5 * ( lambdaQuarter * detTerm - lambdaQuarter + muHalf * Z_1 - mu -
                                    ( lambdaQuarter + muHalf ) * logDetTerm );
      const RealType Y_3 = Y_1 + Y_2;

      RealType summand1, summand2;

      // Da_t
      summand1 = ( 1 / ( a_t + b_t - c_t ) + 1 / ( a_t - b_t + c_t ) - 1 / ( -a_t + b_t + c_t ) +
                   1 / ( a_t + b_t + c_t )) * Y_1;
      summand2 = muHalf * (( 2 * a_t * bSqr + a_t * X_3 ) / X_1 );
      Dest[e0] += a_u * ( summand1 + summand2 );

      // Db_t
      summand1 = ( 1 / ( a_t + b_t - c_t ) - 1 / ( a_t - b_t + c_t ) + 1 / ( -a_t + b_t + c_t ) +
                   1 / ( a_t + b_t + c_t )) * Y_1;
      summand2 = muHalf * (( 2 * aSqr * b_t + b_t * X_3 ) / X_1 );
      Dest[e1] += a_u * ( summand1 + summand2 );

      // Dc_t
      summand1 =
              ( -1. / ( a_t + b_t - c_t ) + 1 / ( a_t - b_t + c_t ) + 1 / ( -a_t + b_t + c_t ) +
                1 / ( a_t + b_t + c_t )) * Y_1;
      summand2 = muHalf * (( -c_t * X_3 ) / X_1 );
      Dest[e2] += a_u * ( summand1 + summand2 );
    }

    Dest.array() *= _memWeight;
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
    VectorType undefTriangleAreas( _numFaces );
    _A.apply( UndeformedGeom, undefTriangleAreas );

    VectorType defTriangleAreas( _numFaces );
    _A.apply( DeformedGeom, defTriangleAreas );

    for ( int faceIdx = 0; faceIdx < _numFaces; faceIdx++ ) {
      int e0 = _topology.getEdgeOfTriangle( faceIdx, 0 );
      int e1 = _topology.getEdgeOfTriangle( faceIdx, 1 );
      int e2 = _topology.getEdgeOfTriangle( faceIdx, 2 );

      // Undeformed Lengths
      const RealType &a = UndeformedGeom[e0];
      const RealType &b = UndeformedGeom[e1];
      const RealType &c = UndeformedGeom[e2];

      // Deformed Lengths
      const RealType &a_t = DeformedGeom[e0];
      const RealType &b_t = DeformedGeom[e1];
      const RealType &c_t = DeformedGeom[e2];

      // Squared Lengths
      const RealType aSqr = a * a;
      const RealType bSqr = b * b;
      const RealType cSqr = c * c;
      const RealType a_tSqr = a_t * a_t;
      const RealType b_tSqr = b_t * b_t;
      const RealType c_tSqr = c_t * c_t;

      // Areas
      const RealType &a_u = undefTriangleAreas[faceIdx];
      const RealType &a_d = defTriangleAreas[faceIdx];

      // Relative change of area = determinant term
      const RealType detTerm = a_d * a_d / ( a_u * a_u );
      const RealType logDetTerm = std::log( detTerm );

      // A terms
      const RealType A_1 = 1 / (( a + b - c ) * ( a_t + b_t + c_t ));
      const RealType A_2 = 1 / (( a - b + c ) * ( a_t + b_t + c_t ));
      const RealType A_3 = 1 / (( -a + b + c ) * ( a_t + b_t + c_t ));
      const RealType A_4 = 1 / (( a + b + c ) * ( a_t + b_t + c_t ));
      const RealType A_5 = 1 / (( a + b - c ) * ( -a_t + b_t + c_t ));
      const RealType A_6 = 1 / (( a - b + c ) * ( -a_t + b_t + c_t ));
      const RealType A_7 = 1 / (( -a + b + c ) * ( -a_t + b_t + c_t ));
      const RealType A_8 = 1 / (( a + b + c ) * ( -a_t + b_t + c_t ));
      const RealType A_9 = 1 / (( a + b - c ) * ( a_t - b_t + c_t ));
      const RealType A_10 = 1 / (( a - b + c ) * ( a_t - b_t + c_t ));
      const RealType A_11 = 1 / (( -a + b + c ) * ( a_t - b_t + c_t ));
      const RealType A_12 = 1 / (( a + b + c ) * ( a_t - b_t + c_t ));
      const RealType A_13 = 1 / (( a + b - c ) * ( a_t + b_t - c_t ));
      const RealType A_14 = 1 / (( a - b + c ) * ( a_t + b_t - c_t ));
      const RealType A_15 = 1 / (( -a + b + c ) * ( a_t + b_t - c_t ));
      const RealType A_16 = 1 / (( a + b + c ) * ( a_t + b_t - c_t ));

      // B terms
      const RealType B_1 = 1 / (( a + b - c ) * ( a + b + c ));
      const RealType B_2 = 1 / (( a - b + c ) * ( a + b + c ));
      const RealType B_3 = 1 / (( -a + b + c ) * ( a + b + c ));
      const RealType B_4 = 1 / (( a + b + c ) * ( a + b + c ));
      const RealType B_5 = 1 / (( a + b - c ) * ( -a + b + c ));
      const RealType B_6 = 1 / (( a - b + c ) * ( -a + b + c ));
      const RealType B_7 = 1 / (( -a + b + c ) * ( -a + b + c ));
      const RealType B_9 = 1 / (( a + b - c ) * ( a - b + c ));
      const RealType B_10 = 1 / (( a - b + c ) * ( a - b + c ));
      const RealType B_13 = 1 / (( a + b - c ) * ( a + b - c ));

      // C terms
      const RealType C_1 = 1 / (( a_t + b_t - c_t ) * ( a_t + b_t + c_t ));
      const RealType C_2 = 1 / (( a_t - b_t + c_t ) * ( a_t + b_t + c_t ));
      const RealType C_3 = 1 / (( -a_t + b_t + c_t ) * ( a_t + b_t + c_t ));
      const RealType C_4 = 1 / (( a_t + b_t + c_t ) * ( a_t + b_t + c_t ));
      const RealType C_5 = 1 / (( a_t + b_t - c_t ) * ( -a_t + b_t + c_t ));
      const RealType C_6 = 1 / (( a_t - b_t + c_t ) * ( -a_t + b_t + c_t ));
      const RealType C_7 = 1 / (( -a_t + b_t + c_t ) * ( -a_t + b_t + c_t ));
      const RealType C_9 = 1 / (( a_t + b_t - c_t ) * ( a_t - b_t + c_t ));
      const RealType C_10 = 1 / (( a_t - b_t + c_t ) * ( a_t - b_t + c_t ));
      const RealType C_13 = 1 / (( a_t + b_t - c_t ) * ( a_t + b_t - c_t ));

      // X polynomials
      const RealType X_2 =
              0.5 *
              (( -aSqr + bSqr + cSqr ) * a_tSqr + ( aSqr + bSqr - cSqr ) * c_tSqr + ( aSqr - bSqr + cSqr ) * b_tSqr );

      // Y polynomials
      const RealType Y_1 = lambdaQuarter * detTerm - lambdaQuarter + ( 2 * X_2 * mu ) / ( 16 * a_u * a_u ) - mu -
                           ( lambdaQuarter + muHalf ) * logDetTerm;

      // Z polynomials
      const RealType Z_1 = (( a + b - c ) * ( a - b + c ) * ( -a + b + c ) +
                            ( a + b - c ) * ( a + b + c ) * ( -a + b + c ) -
                            ( a - b + c ) * ( a + b + c ) * ( -a + b + c ) +
                            ( a + b - c ) * ( a - b + c ) * ( a + b + c ));
      const RealType Z_2 =
              ( a + b - c ) * ( a - b + c ) * ( -a + b + c ) - ( a + b - c ) * ( a + b + c ) * ( -a + b + c ) +
              ( a - b + c ) * ( a + b + c ) * ( -a + b + c ) + ( a + b - c ) * ( a - b + c ) * ( a + b + c );
      const RealType Z_3 =
              ( a + b - c ) * ( a - b + c ) * ( -a + b + c ) + ( a + b - c ) * ( a + b + c ) * ( -a + b + c ) +
              ( a - b + c ) * ( a + b + c ) * ( -a + b + c ) - ( a + b - c ) * ( a - b + c ) * ( a + b + c );

      // S polynomials
      const RealType S_1 = ( +1 / ( a + b - c ) - 1 / ( a - b + c ) - 1 / ( -a + b + c ) - 1 / ( a + b + c ));
      const RealType S_2 = ( -1 / ( a + b - c ) + 1 / ( a - b + c ) - 1 / ( -a + b + c ) - 1 / ( a + b + c ));
      const RealType S_3 = ( -1 / ( a + b - c ) - 1 / ( a - b + c ) + 1 / ( -a + b + c ) - 1 / ( a + b + c ));

      // T polynomials
      RealType T = ( lambdaQuarter * detTerm - lambdaQuarter - muHalf + 2 * X_2 * mu / ( 16 * a_u * a_u ));

      // Da_t Da_t (3,3)
      RealType Cterms = 2 * ( -C_3 + C_2 + C_1 - C_6 - C_5 + C_9 );
      RealType term1 = detTerm * lambdaQuarter * Cterms;
      Cterms = C_13 + C_10 + C_7 + C_4;
      RealType term2 = ( lambdaQuarter + muHalf ) * Cterms;
      RealType term3 = 2 * ( -aSqr + bSqr + cSqr ) * mu / ( 16 * a_u * a_u );

      triplets.emplace_back( e0 + rowOffset, e0 + colOffset, _memWeight * factor * a_u * ( term1 + term2 + term3 ));

      // Db_t Da_t (4,3)
      Cterms = 2 * ( C_1 + C_6 );
      term1 = ( detTerm * lambdaQuarter ) * Cterms;
      Cterms = C_13 - C_10 - C_7 + C_4;
      term2 = ( lambdaQuarter + muHalf ) * Cterms;

      triplets.emplace_back( e1 + rowOffset, e0 + colOffset, _memWeight * factor * a_u * ( term1 + term2 ));
      triplets.emplace_back( e0 + rowOffset, e1 + colOffset, _memWeight * factor * a_u * ( term1 + term2 ));

      // Dc_t Da_t (5,3)
      Cterms = 2 * ( C_2 + C_5 );
      term1 = ( detTerm * lambdaQuarter ) * Cterms;
      Cterms = -C_13 + C_10 - C_7 + C_4;
      term2 = ( lambdaQuarter + muHalf ) * Cterms;

      triplets.emplace_back( e2 + rowOffset, e0 + colOffset, _memWeight * factor * a_u * ( term1 + term2 ));
      triplets.emplace_back( e0 + rowOffset, e2 + colOffset, _memWeight * factor * a_u * ( term1 + term2 ));

      // Db_t Db_t (4,4)
      Cterms = 2 * ( C_3 - C_2 + C_1 - C_6 + C_5 - C_9 );
      term1 = ( detTerm * lambdaQuarter ) * Cterms;
      Cterms = C_13 + C_10 + C_7 + C_4;
      term2 = ( lambdaQuarter + muHalf ) * Cterms;
      term3 = 2 * ( aSqr - bSqr + cSqr ) * mu / ( 16 * a_u * a_u );

      triplets.emplace_back( e1 + rowOffset, e1 + colOffset, _memWeight * factor * a_u * ( term1 + term2 + term3 ));

      // Dc_t Db_t (5,4)
      Cterms = 2 * ( C_3 + C_9 );
      term1 = ( detTerm * lambdaQuarter ) * Cterms;
      Cterms = -C_13 - C_10 + C_7 + C_4;
      term2 = ( lambdaQuarter + muHalf ) * Cterms;

      triplets.emplace_back( e2 + rowOffset, e1 + colOffset, _memWeight * factor * a_u * ( term1 + term2 ));
      triplets.emplace_back( e1 + rowOffset, e2 + colOffset, _memWeight * factor * a_u * ( term1 + term2 ));

      // Dc_t Dc_t (5,5)
      Cterms = 2 * ( C_3 + C_2 - C_1 + C_6 - C_5 - C_9 );
      term1 = ( detTerm * lambdaQuarter ) * Cterms;
      Cterms = C_13 + C_10 + C_7 + C_4;
      term2 = ( lambdaQuarter + muHalf ) * Cterms;
      term3 = 2 * ( aSqr + bSqr - cSqr ) * mu / ( 16 * a_u * a_u );

      triplets.emplace_back( e2 + rowOffset, e2 + colOffset, _memWeight * factor * a_u * ( term1 + term2 + term3 ));
    }
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

    VectorType undefTriangleAreas( _numFaces );
    _A.apply( UndeformedGeom, undefTriangleAreas );

    VectorType defTriangleAreas( _numFaces );
    _A.apply( DeformedGeom, defTriangleAreas );

    for ( int faceIdx = 0; faceIdx < _numFaces; faceIdx++ ) {
      int e0 = _topology.getEdgeOfTriangle( faceIdx, 0 );
      int e1 = _topology.getEdgeOfTriangle( faceIdx, 1 );
      int e2 = _topology.getEdgeOfTriangle( faceIdx, 2 );

      // Undeformed Lengths
      const RealType &a = UndeformedGeom[e0];
      const RealType &b = UndeformedGeom[e1];
      const RealType &c = UndeformedGeom[e2];

      // Deformed Lengths
      const RealType &a_t = DeformedGeom[e0];
      const RealType &b_t = DeformedGeom[e1];
      const RealType &c_t = DeformedGeom[e2];

      // Squared Lengths
      const RealType aSqr = a * a;
      const RealType bSqr = b * b;
      const RealType cSqr = c * c;
      const RealType a_tSqr = a_t * a_t;
      const RealType b_tSqr = b_t * b_t;
      const RealType c_tSqr = c_t * c_t;

      // Areas
      const RealType &a_u = undefTriangleAreas[faceIdx];
      const RealType &a_d = defTriangleAreas[faceIdx];

      // Relative change of area = determinant term
      const RealType detTerm = a_d * a_d / ( a_u * a_u );
      const RealType logDetTerm = std::log( detTerm );

      // A terms
      const RealType A_1 = 1 / (( a + b - c ) * ( a_t + b_t + c_t ));
      const RealType A_2 = 1 / (( a - b + c ) * ( a_t + b_t + c_t ));
      const RealType A_3 = 1 / (( -a + b + c ) * ( a_t + b_t + c_t ));
      const RealType A_4 = 1 / (( a + b + c ) * ( a_t + b_t + c_t ));
      const RealType A_5 = 1 / (( a + b - c ) * ( -a_t + b_t + c_t ));
      const RealType A_6 = 1 / (( a - b + c ) * ( -a_t + b_t + c_t ));
      const RealType A_7 = 1 / (( -a + b + c ) * ( -a_t + b_t + c_t ));
      const RealType A_8 = 1 / (( a + b + c ) * ( -a_t + b_t + c_t ));
      const RealType A_9 = 1 / (( a + b - c ) * ( a_t - b_t + c_t ));
      const RealType A_10 = 1 / (( a - b + c ) * ( a_t - b_t + c_t ));
      const RealType A_11 = 1 / (( -a + b + c ) * ( a_t - b_t + c_t ));
      const RealType A_12 = 1 / (( a + b + c ) * ( a_t - b_t + c_t ));
      const RealType A_13 = 1 / (( a + b - c ) * ( a_t + b_t - c_t ));
      const RealType A_14 = 1 / (( a - b + c ) * ( a_t + b_t - c_t ));
      const RealType A_15 = 1 / (( -a + b + c ) * ( a_t + b_t - c_t ));
      const RealType A_16 = 1 / (( a + b + c ) * ( a_t + b_t - c_t ));

      // B terms
      const RealType B_1 = 1 / (( a + b - c ) * ( a + b + c ));
      const RealType B_2 = 1 / (( a - b + c ) * ( a + b + c ));
      const RealType B_3 = 1 / (( -a + b + c ) * ( a + b + c ));
      const RealType B_4 = 1 / (( a + b + c ) * ( a + b + c ));
      const RealType B_5 = 1 / (( a + b - c ) * ( -a + b + c ));
      const RealType B_6 = 1 / (( a - b + c ) * ( -a + b + c ));
      const RealType B_7 = 1 / (( -a + b + c ) * ( -a + b + c ));
      const RealType B_9 = 1 / (( a + b - c ) * ( a - b + c ));
      const RealType B_10 = 1 / (( a - b + c ) * ( a - b + c ));
      const RealType B_13 = 1 / (( a + b - c ) * ( a + b - c ));

      // C terms
      const RealType C_1 = 1 / (( a_t + b_t - c_t ) * ( a_t + b_t + c_t ));
      const RealType C_2 = 1 / (( a_t - b_t + c_t ) * ( a_t + b_t + c_t ));
      const RealType C_3 = 1 / (( -a_t + b_t + c_t ) * ( a_t + b_t + c_t ));
      const RealType C_4 = 1 / (( a_t + b_t + c_t ) * ( a_t + b_t + c_t ));
      const RealType C_5 = 1 / (( a_t + b_t - c_t ) * ( -a_t + b_t + c_t ));
      const RealType C_6 = 1 / (( a_t - b_t + c_t ) * ( -a_t + b_t + c_t ));
      const RealType C_7 = 1 / (( -a_t + b_t + c_t ) * ( -a_t + b_t + c_t ));
      const RealType C_9 = 1 / (( a_t + b_t - c_t ) * ( a_t - b_t + c_t ));
      const RealType C_10 = 1 / (( a_t - b_t + c_t ) * ( a_t - b_t + c_t ));
      const RealType C_13 = 1 / (( a_t + b_t - c_t ) * ( a_t + b_t - c_t ));

      // X polynomials
      const RealType X_2 =
              0.5 *
              (( -aSqr + bSqr + cSqr ) * a_tSqr + ( aSqr + bSqr - cSqr ) * c_tSqr + ( aSqr - bSqr + cSqr ) * b_tSqr );

      // Y polynomials
      const RealType Y_1 = lambdaQuarter * detTerm - lambdaQuarter + ( 2 * X_2 * mu ) / ( 16 * a_u * a_u ) - mu -
                           ( lambdaQuarter + muHalf ) * logDetTerm;

      // Z polynomials
      const RealType Z_1 = (( a + b - c ) * ( a - b + c ) * ( -a + b + c ) +
                            ( a + b - c ) * ( a + b + c ) * ( -a + b + c ) -
                            ( a - b + c ) * ( a + b + c ) * ( -a + b + c ) +
                            ( a + b - c ) * ( a - b + c ) * ( a + b + c ));
      const RealType Z_2 =
              ( a + b - c ) * ( a - b + c ) * ( -a + b + c ) - ( a + b - c ) * ( a + b + c ) * ( -a + b + c ) +
              ( a - b + c ) * ( a + b + c ) * ( -a + b + c ) + ( a + b - c ) * ( a - b + c ) * ( a + b + c );
      const RealType Z_3 =
              ( a + b - c ) * ( a - b + c ) * ( -a + b + c ) + ( a + b - c ) * ( a + b + c ) * ( -a + b + c ) +
              ( a - b + c ) * ( a + b + c ) * ( -a + b + c ) - ( a + b - c ) * ( a - b + c ) * ( a + b + c );

      // S polynomials
      const RealType S_1 = ( +1 / ( a + b - c ) - 1 / ( a - b + c ) - 1 / ( -a + b + c ) - 1 / ( a + b + c ));
      const RealType S_2 = ( -1 / ( a + b - c ) + 1 / ( a - b + c ) - 1 / ( -a + b + c ) - 1 / ( a + b + c ));
      const RealType S_3 = ( -1 / ( a + b - c ) - 1 / ( a - b + c ) + 1 / ( -a + b + c ) - 1 / ( a + b + c ));

      // T polynomials
      RealType T = ( lambdaQuarter * detTerm - lambdaQuarter - muHalf + 2 * X_2 * mu / ( 16 * a_u * a_u ));

      // Da Da
      RealType term1 = 2 * ( lambdaQuarter + muHalf + Y_1 / 4 ) * ( B_1 + B_2 - B_3 - B_5 - B_6 + B_9 );
      RealType term2 = ( lambdaQuarter + muHalf + T - Y_1 / 4 ) * ( B_4 + B_7 + B_10 + B_13 );
      RealType term3 = ( 2 * ( -a_tSqr + b_tSqr + c_tSqr ) * mu / ( 16 * a_u * a_u )) * ( 1 + a * S_3 );

      RealType term4 = 0;
      RealType term5 = 0;

      triplets.emplace_back( e0 + rowOffset, e0 + colOffset,
                             _memWeight * factor * a_u * ( term1 + term2 + term3 ));


      // Db Da
      term1 = 2 * ( lambdaQuarter + muHalf + Y_1 / 4 ) * ( B_6 + B_1 );
      term2 = ( lambdaQuarter + muHalf + T - Y_1 / 4 ) * ( B_4 - B_7 - B_10 + B_13 );
      term3 = ( a * ( -a_tSqr + b_tSqr + c_tSqr ) * mu / ( 16 * a_u * a_u )) * S_2;
      term4 = ( b * ( a_tSqr - b_tSqr + c_tSqr ) * mu / ( 16 * a_u * a_u )) * S_3;


      triplets.emplace_back( e1 + rowOffset, e0 + colOffset,
                             _memWeight * factor * a_u * ( term1 + term2 + term3 + term4 ));
      triplets.emplace_back( e0 + rowOffset, e1 + colOffset,
                             _memWeight * factor * a_u * ( term1 + term2 + term3 + term4 ));

      // Dc Da
      term1 = 2 * ( lambdaQuarter + muHalf + Y_1 / 4 ) * ( B_5 + B_2 );
      term2 = ( lambdaQuarter + muHalf + T - Y_1 / 4 ) * ( B_4 - B_7 + B_10 - B_13 );
      term3 = ( a * ( -a_tSqr + b_tSqr + c_tSqr ) * mu / ( 16 * a_u * a_u )) * S_1;
      term4 = ( c * ( a_tSqr + b_tSqr - c_tSqr ) * mu / ( 16 * a_u * a_u )) * S_3;


      triplets.emplace_back( e2 + rowOffset, e0 + colOffset,
                             _memWeight * factor * a_u * ( term1 + term2 + term3 + term4 ));
      triplets.emplace_back( e0 + rowOffset, e2 + colOffset,
                             _memWeight * factor * a_u * ( term1 + term2 + term3 + term4 ));


      // Db Db
      term1 = 2 * ( lambdaQuarter + muHalf + Y_1 / 4 ) * ( B_1 - B_2 + B_3 + B_5 - B_6 - B_9 );
      term2 = ( lambdaQuarter + muHalf + T - Y_1 / 4 ) * ( B_4 + B_7 + B_10 + B_13 );
      term3 = ( 2 * ( a_tSqr - b_tSqr + c_tSqr ) * mu / ( 16 * a_u * a_u )) * ( 1 + b * S_2 );

      triplets.emplace_back( e1 + rowOffset, e1 + colOffset,
                             _memWeight * factor * a_u * ( term1 + term2 + term3 ));

      // Dc Db
      term1 = 2 * ( lambdaQuarter + muHalf + Y_1 / 4 ) * ( B_9 + B_3 );
      term2 = ( lambdaQuarter + muHalf + T - Y_1 / 4 ) * ( B_4 + B_7 - B_10 - B_13 );
      term3 = ( b * ( a_tSqr - b_tSqr + c_tSqr ) * mu / ( 16 * a_u * a_u )) * S_1;
      term4 = ( c * ( a_tSqr + b_tSqr - c_tSqr ) * mu / ( 16 * a_u * a_u )) * S_2;


      triplets.emplace_back( e2 + rowOffset, e1 + colOffset,
                             _memWeight * factor * a_u * ( term1 + term2 + term3 + term4 ));
      triplets.emplace_back( e1 + rowOffset, e2 + colOffset,
                             _memWeight * factor * a_u * ( term1 + term2 + term3 + term4 ));

      // Dc Dc
      term1 = 2 * ( lambdaQuarter + muHalf + Y_1 / 4 ) * ( -B_1 + B_2 + B_3 - B_5 + B_6 - B_9 );
      term2 = ( lambdaQuarter + muHalf + T - Y_1 / 4 ) * ( B_4 + B_7 + B_10 + B_13 );
      term3 = ( 2 * ( a_tSqr + b_tSqr - c_tSqr ) * mu / ( 16 * a_u * a_u )) * ( 1 + c * S_1 );

      triplets.emplace_back( e2 + rowOffset, e2 + colOffset,
                             _memWeight * factor * a_u * ( term1 + term2 + term3 ));

    }
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
    VectorType undefTriangleAreas( _numFaces );
    _A.apply( UndeformedGeom, undefTriangleAreas );

    VectorType defTriangleAreas( _numFaces );
    _A.apply( DeformedGeom, defTriangleAreas );

    for ( int faceIdx = 0; faceIdx < _numFaces; faceIdx++ ) {
      int e0 = _topology.getEdgeOfTriangle( faceIdx, 0 );
      int e1 = _topology.getEdgeOfTriangle( faceIdx, 1 );
      int e2 = _topology.getEdgeOfTriangle( faceIdx, 2 );

      // Undeformed Lengths
      const RealType &a = UndeformedGeom[e0];
      const RealType &b = UndeformedGeom[e1];
      const RealType &c = UndeformedGeom[e2];

      // Deformed Lengths
      const RealType &a_t = DeformedGeom[e0];
      const RealType &b_t = DeformedGeom[e1];
      const RealType &c_t = DeformedGeom[e2];

      // Squared Lengths
      const RealType aSqr = a * a;
      const RealType bSqr = b * b;
      const RealType cSqr = c * c;
      const RealType a_tSqr = a_t * a_t;
      const RealType b_tSqr = b_t * b_t;
      const RealType c_tSqr = c_t * c_t;

      // Areas
      const RealType &a_u = undefTriangleAreas[faceIdx];
      const RealType &a_d = defTriangleAreas[faceIdx];

      // Relative change of area = determinant term
      const RealType detTerm = a_d * a_d / ( a_u * a_u );
      const RealType logDetTerm = std::log( detTerm );

      // A terms
      const RealType A_1 = 1 / (( a + b - c ) * ( a_t + b_t + c_t ));
      const RealType A_2 = 1 / (( a - b + c ) * ( a_t + b_t + c_t ));
      const RealType A_3 = 1 / (( -a + b + c ) * ( a_t + b_t + c_t ));
      const RealType A_4 = 1 / (( a + b + c ) * ( a_t + b_t + c_t ));
      const RealType A_5 = 1 / (( a + b - c ) * ( -a_t + b_t + c_t ));
      const RealType A_6 = 1 / (( a - b + c ) * ( -a_t + b_t + c_t ));
      const RealType A_7 = 1 / (( -a + b + c ) * ( -a_t + b_t + c_t ));
      const RealType A_8 = 1 / (( a + b + c ) * ( -a_t + b_t + c_t ));
      const RealType A_9 = 1 / (( a + b - c ) * ( a_t - b_t + c_t ));
      const RealType A_10 = 1 / (( a - b + c ) * ( a_t - b_t + c_t ));
      const RealType A_11 = 1 / (( -a + b + c ) * ( a_t - b_t + c_t ));
      const RealType A_12 = 1 / (( a + b + c ) * ( a_t - b_t + c_t ));
      const RealType A_13 = 1 / (( a + b - c ) * ( a_t + b_t - c_t ));
      const RealType A_14 = 1 / (( a - b + c ) * ( a_t + b_t - c_t ));
      const RealType A_15 = 1 / (( -a + b + c ) * ( a_t + b_t - c_t ));
      const RealType A_16 = 1 / (( a + b + c ) * ( a_t + b_t - c_t ));

      // B terms
      const RealType B_1 = 1 / (( a + b - c ) * ( a + b + c ));
      const RealType B_2 = 1 / (( a - b + c ) * ( a + b + c ));
      const RealType B_3 = 1 / (( -a + b + c ) * ( a + b + c ));
      const RealType B_4 = 1 / (( a + b + c ) * ( a + b + c ));
      const RealType B_5 = 1 / (( a + b - c ) * ( -a + b + c ));
      const RealType B_6 = 1 / (( a - b + c ) * ( -a + b + c ));
      const RealType B_7 = 1 / (( -a + b + c ) * ( -a + b + c ));
      const RealType B_9 = 1 / (( a + b - c ) * ( a - b + c ));
      const RealType B_10 = 1 / (( a - b + c ) * ( a - b + c ));
      const RealType B_13 = 1 / (( a + b - c ) * ( a + b - c ));

      // C terms
      const RealType C_1 = 1 / (( a_t + b_t - c_t ) * ( a_t + b_t + c_t ));
      const RealType C_2 = 1 / (( a_t - b_t + c_t ) * ( a_t + b_t + c_t ));
      const RealType C_3 = 1 / (( -a_t + b_t + c_t ) * ( a_t + b_t + c_t ));
      const RealType C_4 = 1 / (( a_t + b_t + c_t ) * ( a_t + b_t + c_t ));
      const RealType C_5 = 1 / (( a_t + b_t - c_t ) * ( -a_t + b_t + c_t ));
      const RealType C_6 = 1 / (( a_t - b_t + c_t ) * ( -a_t + b_t + c_t ));
      const RealType C_7 = 1 / (( -a_t + b_t + c_t ) * ( -a_t + b_t + c_t ));
      const RealType C_9 = 1 / (( a_t + b_t - c_t ) * ( a_t - b_t + c_t ));
      const RealType C_10 = 1 / (( a_t - b_t + c_t ) * ( a_t - b_t + c_t ));
      const RealType C_13 = 1 / (( a_t + b_t - c_t ) * ( a_t + b_t - c_t ));

      // X polynomials
      const RealType X_2 =
              0.5 *
              (( -aSqr + bSqr + cSqr ) * a_tSqr + ( aSqr + bSqr - cSqr ) * c_tSqr + ( aSqr - bSqr + cSqr ) * b_tSqr );

      // Y polynomials
      const RealType Y_1 = lambdaQuarter * detTerm - lambdaQuarter + ( 2 * X_2 * mu ) / ( 16 * a_u * a_u ) - mu -
                           ( lambdaQuarter + muHalf ) * logDetTerm;

      // Z polynomials
      const RealType Z_1 = (( a + b - c ) * ( a - b + c ) * ( -a + b + c ) +
                            ( a + b - c ) * ( a + b + c ) * ( -a + b + c ) -
                            ( a - b + c ) * ( a + b + c ) * ( -a + b + c ) +
                            ( a + b - c ) * ( a - b + c ) * ( a + b + c ));
      const RealType Z_2 =
              ( a + b - c ) * ( a - b + c ) * ( -a + b + c ) - ( a + b - c ) * ( a + b + c ) * ( -a + b + c ) +
              ( a - b + c ) * ( a + b + c ) * ( -a + b + c ) + ( a + b - c ) * ( a - b + c ) * ( a + b + c );
      const RealType Z_3 =
              ( a + b - c ) * ( a - b + c ) * ( -a + b + c ) + ( a + b - c ) * ( a + b + c ) * ( -a + b + c ) +
              ( a - b + c ) * ( a + b + c ) * ( -a + b + c ) - ( a + b - c ) * ( a - b + c ) * ( a + b + c );

      // S polynomials
      const RealType S_1 = ( +1 / ( a + b - c ) - 1 / ( a - b + c ) - 1 / ( -a + b + c ) - 1 / ( a + b + c ));
      const RealType S_2 = ( -1 / ( a + b - c ) + 1 / ( a - b + c ) - 1 / ( -a + b + c ) - 1 / ( a + b + c ));
      const RealType S_3 = ( -1 / ( a + b - c ) - 1 / ( a - b + c ) + 1 / ( -a + b + c ) - 1 / ( a + b + c ));

      // T polynomials
      RealType T = ( lambdaQuarter * detTerm - lambdaQuarter - muHalf + 2 * X_2 * mu / ( 16 * a_u * a_u ));

      RealType term1, term2, term3;

      if ( FirstDerivWRTDef ) {
        // Da_t Da (3,0)
        RealType Aterms =
                -A_1 - A_2 + A_3 - A_4 + A_5 + A_6 - A_7 + A_8 - A_9 - A_10 + A_11 - A_12 - A_13 - A_14 + A_15 - A_16;

        term1 = 0.5 * a_u * ( detTerm * lambdaQuarter + ( lambdaQuarter + muHalf )) * Aterms;
        term2 = a_t * ( -aSqr + bSqr + cSqr ) * mu / ( 16 * a_u ) * S_3;
        term3 = a * a_t * mu / ( 4 * a_u );

        triplets.emplace_back( e0 + rowOffset, e0 + colOffset, _memWeight * factor * ( term1 + term2 - term3 ));

        // Db_t Da
        Aterms = -A_1 - A_2 + A_3 - A_4 - A_5 - A_6 + A_7 - A_8 + A_9 + A_10 - A_11 + A_12 - A_13 - A_14 + A_15 - A_16;

        term1 = 0.5 * a_u * ( detTerm * lambdaQuarter + ( lambdaQuarter + muHalf )) * Aterms;
        term2 = b_t * ( aSqr - bSqr + cSqr ) * mu / ( 16 * a_u ) * S_3;
        term3 = a * b_t * mu / ( 4 * a_u );

        triplets.emplace_back( e1 + rowOffset, e0 + colOffset, _memWeight * factor * ( term1 + term2 + term3 ));

        // Dc_t Da
        Aterms = -A_1 - A_2 + A_3 - A_4 - A_5 - A_6 + A_7 - A_8 - A_9 - A_10 + A_11 - A_12 + A_13 + A_14 - A_15 + A_16;

        term1 = 0.5 * a_u * ( detTerm * lambdaQuarter + ( lambdaQuarter + muHalf )) * Aterms;
        term2 = c_t * ( aSqr + bSqr - cSqr ) * mu / ( 16 * a_u ) * S_3;
        term3 = a * c_t * mu / ( 4 * a_u );

        triplets.emplace_back( e2 + rowOffset, e0 + colOffset, _memWeight * factor * ( term1 + term2 + term3 ));

        // Da_t Db (3,1)
        Aterms = -A_1 + A_2 - A_3 - A_4 + A_5 - A_6 + A_7 + A_8 - A_9 + A_10 - A_11 - A_12 - A_13 + A_14 - A_15 - A_16;

        term1 = 0.5 * a_u * ( detTerm * lambdaQuarter + ( lambdaQuarter + muHalf )) * Aterms;
        term2 = a_t * ( -aSqr + bSqr + cSqr ) * mu / ( 16 * a_u ) * S_2;
        term3 = b * a_t * mu / ( 4 * a_u );

        triplets.emplace_back( e0 + rowOffset, e1 + colOffset, _memWeight * factor * ( term1 + term2 + term3 ));

        // Db_t Db (4,1)
        Aterms = -A_1 + A_2 - A_3 - A_4 - A_5 + A_6 - A_7 - A_8 + A_9 - A_10 + A_11 + A_12 - A_13 + A_14 - A_15 - A_16;

        term1 = 0.5 * a_u * ( detTerm * lambdaQuarter + ( lambdaQuarter + muHalf )) * Aterms;
        term2 = b_t * ( aSqr - bSqr + cSqr ) * mu / ( 16 * a_u ) * S_2;
        term3 = -b * b_t * mu / ( 4 * a_u );

        triplets.emplace_back( e1 + rowOffset, e1 + colOffset, _memWeight * factor * ( term1 + term2 + term3 ));

        // Dc_t Db (5,1)
        Aterms = -A_1 + A_2 - A_3 - A_4 - A_5 + A_6 - A_7 - A_8 - A_9 + A_10 - A_11 - A_12 + A_13 - A_14 + A_15 + A_16;

        term1 = 0.5 * a_u * ( detTerm * lambdaQuarter + ( lambdaQuarter + muHalf )) * Aterms;
        term2 = c_t * ( aSqr + bSqr - cSqr ) * mu / ( 16 * a_u ) * S_2;
        term3 = b * c_t * mu / ( 4 * a_u );

        triplets.emplace_back( e2 + rowOffset, e1 + colOffset, _memWeight * factor * ( term1 + term2 + term3 ));

        // Da_t Dc (3,2)
        Aterms = A_1 - A_2 - A_3 - A_4 - A_5 + A_6 + A_7 + A_8 + A_9 - A_10 - A_11 - A_12 + A_13 - A_14 - A_15 - A_16;

        term1 = 0.5 * a_u * ( detTerm * lambdaQuarter + ( lambdaQuarter + muHalf )) * Aterms;
        term2 = a_t * ( -aSqr + bSqr + cSqr ) * mu / ( 16 * a_u ) * S_1;
        term3 = c * a_t * mu / ( 4 * a_u );

        triplets.emplace_back( e0 + rowOffset, e2 + colOffset, _memWeight * factor * ( term1 + term2 + term3 ));

        // Db_t Dc (4,2)
        Aterms = A_1 - A_2 - A_3 - A_4 + A_5 - A_6 - A_7 - A_8 - A_9 + A_10 + A_11 + A_12 + A_13 - A_14 - A_15 - A_16;

        term1 = 0.5 * a_u * ( detTerm * lambdaQuarter + ( lambdaQuarter + muHalf )) * Aterms;
        term2 = b_t * ( aSqr - bSqr + cSqr ) * mu / ( 16 * a_u ) * S_1;
        term3 = c * b_t * mu / ( 4 * a_u );

        triplets.emplace_back( e1 + rowOffset, e2 + colOffset, _memWeight * factor * ( term1 + term2 + term3 ));

        // Dc_t Dc (5,2)
        Aterms = A_1 - A_2 - A_3 - A_4 + A_5 - A_6 - A_7 - A_8 + A_9 - A_10 - A_11 - A_12 - A_13 + A_14 + A_15 + A_16;

        term1 = 0.5 * a_u * ( detTerm * lambdaQuarter + ( lambdaQuarter + muHalf )) * Aterms;
        term2 = c_t * ( aSqr + bSqr - cSqr ) * mu / ( 16 * a_u ) * S_1;
        term3 = -c * c_t * mu / ( 4 * a_u );

        triplets.emplace_back( e2 + rowOffset, e2 + colOffset, _memWeight * factor * ( term1 + term2 + term3 ));

      }
      else {
        // Da_t Da (3,0)
        RealType Aterms =
                -A_1 - A_2 + A_3 - A_4 + A_5 + A_6 - A_7 + A_8 - A_9 - A_10 + A_11 - A_12 - A_13 - A_14 + A_15 - A_16;

        term1 = 0.5 * a_u * ( detTerm * lambdaQuarter + ( lambdaQuarter + muHalf )) * Aterms;
        term2 = a_t * ( -aSqr + bSqr + cSqr ) * mu / ( 16 * a_u ) * S_3;
        term3 = a * a_t * mu / ( 4 * a_u );

        triplets.emplace_back( e0 + rowOffset, e0 + colOffset, _memWeight * factor * ( term1 + term2 - term3 ));

        // Db_t Da
        Aterms = -A_1 - A_2 + A_3 - A_4 - A_5 - A_6 + A_7 - A_8 + A_9 + A_10 - A_11 + A_12 - A_13 - A_14 + A_15 - A_16;

        term1 = 0.5 * a_u * ( detTerm * lambdaQuarter + ( lambdaQuarter + muHalf )) * Aterms;
        term2 = b_t * ( aSqr - bSqr + cSqr ) * mu / ( 16 * a_u ) * S_3;
        term3 = a * b_t * mu / ( 4 * a_u );

        triplets.emplace_back( e0 + rowOffset, e1 + colOffset, _memWeight * factor * ( term1 + term2 + term3 ));

        // Dc_t Da
        Aterms = -A_1 - A_2 + A_3 - A_4 - A_5 - A_6 + A_7 - A_8 - A_9 - A_10 + A_11 - A_12 + A_13 + A_14 - A_15 + A_16;

        term1 = 0.5 * a_u * ( detTerm * lambdaQuarter + ( lambdaQuarter + muHalf )) * Aterms;
        term2 = c_t * ( aSqr + bSqr - cSqr ) * mu / ( 16 * a_u ) * S_3;
        term3 = a * c_t * mu / ( 4 * a_u );

        triplets.emplace_back( e0 + rowOffset, e2 + colOffset, _memWeight * factor * ( term1 + term2 + term3 ));

        // Da_t Db (3,1)
        Aterms = -A_1 + A_2 - A_3 - A_4 + A_5 - A_6 + A_7 + A_8 - A_9 + A_10 - A_11 - A_12 - A_13 + A_14 - A_15 - A_16;

        term1 = 0.5 * a_u * ( detTerm * lambdaQuarter + ( lambdaQuarter + muHalf )) * Aterms;
        term2 = a_t * ( -aSqr + bSqr + cSqr ) * mu / ( 16 * a_u ) * S_2;
        term3 = b * a_t * mu / ( 4 * a_u );

        triplets.emplace_back( e1 + rowOffset, e0 + colOffset, _memWeight * factor * ( term1 + term2 + term3 ));

        // Db_t Db (4,1)
        Aterms = -A_1 + A_2 - A_3 - A_4 - A_5 + A_6 - A_7 - A_8 + A_9 - A_10 + A_11 + A_12 - A_13 + A_14 - A_15 - A_16;

        term1 = 0.5 * a_u * ( detTerm * lambdaQuarter + ( lambdaQuarter + muHalf )) * Aterms;
        term2 = b_t * ( aSqr - bSqr + cSqr ) * mu / ( 16 * a_u ) * S_2;
        term3 = -b * b_t * mu / ( 4 * a_u );

        triplets.emplace_back( e1 + rowOffset, e1 + colOffset, _memWeight * factor * ( term1 + term2 + term3 ));

        // Dc_t Db (5,1)
        Aterms = -A_1 + A_2 - A_3 - A_4 - A_5 + A_6 - A_7 - A_8 - A_9 + A_10 - A_11 - A_12 + A_13 - A_14 + A_15 + A_16;

        term1 = 0.5 * a_u * ( detTerm * lambdaQuarter + ( lambdaQuarter + muHalf )) * Aterms;
        term2 = c_t * ( aSqr + bSqr - cSqr ) * mu / ( 16 * a_u ) * S_2;
        term3 = b * c_t * mu / ( 4 * a_u );

        triplets.emplace_back( e1 + rowOffset, e2 + colOffset, _memWeight * factor * ( term1 + term2 + term3 ));

        // Da_t Dc (3,2)
        Aterms = A_1 - A_2 - A_3 - A_4 - A_5 + A_6 + A_7 + A_8 + A_9 - A_10 - A_11 - A_12 + A_13 - A_14 - A_15 - A_16;

        term1 = 0.5 * a_u * ( detTerm * lambdaQuarter + ( lambdaQuarter + muHalf )) * Aterms;
        term2 = a_t * ( -aSqr + bSqr + cSqr ) * mu / ( 16 * a_u ) * S_1;
        term3 = c * a_t * mu / ( 4 * a_u );

        triplets.emplace_back( e2 + rowOffset, e0 + colOffset, _memWeight * factor * ( term1 + term2 + term3 ));

        // Db_t Dc (4,2)
        Aterms = A_1 - A_2 - A_3 - A_4 + A_5 - A_6 - A_7 - A_8 - A_9 + A_10 + A_11 + A_12 + A_13 - A_14 - A_15 - A_16;

        term1 = 0.5 * a_u * ( detTerm * lambdaQuarter + ( lambdaQuarter + muHalf )) * Aterms;
        term2 = b_t * ( aSqr - bSqr + cSqr ) * mu / ( 16 * a_u ) * S_1;
        term3 = c * b_t * mu / ( 4 * a_u );

        triplets.emplace_back( e2 + rowOffset, e1 + colOffset, _memWeight * factor * ( term1 + term2 + term3 ));

        // Dc_t Dc (5,2)
        Aterms = A_1 - A_2 - A_3 - A_4 + A_5 - A_6 - A_7 - A_8 + A_9 - A_10 - A_11 - A_12 - A_13 + A_14 + A_15 + A_16;

        term1 = 0.5 * a_u * ( detTerm * lambdaQuarter + ( lambdaQuarter + muHalf )) * Aterms;
        term2 = c_t * ( aSqr + bSqr - cSqr ) * mu / ( 16 * a_u ) * S_1;
        term3 = -c * c_t * mu / ( 4 * a_u );

        triplets.emplace_back( e2 + rowOffset, e2 + colOffset, _memWeight * factor * ( term1 + term2 + term3 ));
      }
    }
  }

  /**
   * \return the number of nonzero entries of the sparse hessian (twice the number of edges)
   */
  int numOfNonZeroHessianEntries() const {
    return 36 * _numFaces; // for each face we have 6 DOFs, 3 undeformed lengths, 3 deformed lengths
  }

  /**
   * \return the membrane weight of the energy
   */
  RealType getMembraneWeight() const {
    return _memWeight;
  }
};


#endif //NRIC_NONLINEARDEFORMATIONS_H
