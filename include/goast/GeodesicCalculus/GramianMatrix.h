// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Approximation of Gramian matrix.
 * \author Sassen
 * 
 */
#ifndef GRAMIANMATRIX_H
#define GRAMIANMATRIX_H

#include <goast/Core.h>
#include <goast/Optimization.h>

#include "DiscreteGeodesicCalculus.h"

/**
 * \brief Approximation of Gramian matrix.
 * \author Sassen
 *
 *  For \f$ n \f$ input shapes \f$ s_1, \ldots, s_n \f$ and geodesic mean \f$ s \f$,
 *  the Gramian matrix \f$ G \in \R^{n,n} \f$ is defined as \f$ g_{ij} = g( \log_s(s_i), \log_s(s_j) ) \f$
 *  with an appropiate scalar product \f$ g(.,.) \f$.
 *
 *  Here we compute an approximate Gramian matrix as proposed in Heeren et al. 2018, \cite HeZhRu18.
 *  This approximation avoids tangent space operations but works with deformation energies between shapes instead.
 */
template<typename ConfiguratorType>
class GramianMatrixAssembler {
protected:
    typedef typename ConfiguratorType::RealType RealType;

    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;
    typedef std::vector<TripletType> TripletListType;

    const MeshTopologySaver &_topology;
    const DeformationBase<ConfiguratorType> &_W;
    const VectorType &_shapes;
    const VectorType &_mean;
    const unsigned int _numOfShapes;
    std::vector<Eigen::Ref<const VectorType> > _argRefs;
    const int _K;

    bool _quiet;

    const OptimizationParameters<ConfiguratorType>& _optPars;

    std::vector<int> const* _mask;

public:
    GramianMatrixAssembler(const MeshTopologySaver &Topology,
                           const DeformationBase<ConfiguratorType> &W,
                           const VectorType &shapes,
                           const VectorType &mean,
                           const unsigned int numOfShapes,
                           const int K,
                           const OptimizationParameters<ConfiguratorType> &optPars,
                           bool quiet = true)
            : _topology(Topology), _W(W), _shapes(shapes), _mean(mean),
              _numOfShapes(numOfShapes), _K(K), _optPars(optPars), _quiet(quiet) {
        // Create references for this different shapes bec. of convenience
        _argRefs.reserve(numOfShapes);

        const int numLocalDofs = _shapes.size() / _numOfShapes;
        for (int k = 0; k < numOfShapes; k++)
            _argRefs.push_back(_shapes.segment(k * numLocalDofs, numLocalDofs));
    }

    void setBoundaryMask(const std::vector<int> &localMask) {
        _mask = &localMask;
//        fillPathMask(1, 3 * _topology.getNumVertices(), localMask, *_mask);
    }

    // argument can either
    //   * be empty (then the geodesic segments are computed) 
    //   * contain full geodesic segments
    //   * contain only inner shapes
    void assemble( const VectorType& Argument, MatrixType &Dest) {
        if( (Dest.rows() != _numOfShapes) || (Dest.cols() != _numOfShapes) )
            Dest.resize(_numOfShapes, _numOfShapes);

        // fill triplet lists
        TripletListType tripletList;
        tripletList.reserve(_numOfShapes * _numOfShapes);      
            
        pushTriplets( Argument, tripletList);
        Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
    };

    // argument can either
    //   * be empty (then the geodesic segments are computed) 
    //   * contain full geodesic segments
    //   * contain only inner shapes
    void pushTriplets( const VectorType& Argument, TripletListType &tripletList) {

        const int numLocalDofs = _shapes.size() / _numOfShapes;
        const int numGeodesicsDofs = numLocalDofs * (_K - 1);
        
        std::vector<Eigen::Ref<const VectorType> > innerShapesRefs;
        // if K = 1, there are no geodesic segments but only the input shapes!
        if( _K == 1 ){
            if(!_quiet) std::cerr << "K=1: consider input shapes as inner shapes!" << std::endl;
            for (int k = 0; k < _numOfShapes; k++)
                innerShapesRefs.push_back(_shapes.segment(k * numLocalDofs, numLocalDofs));
        }
        // otherwise check whether geodesic segments are given
        else{        
          // nothing is given: compute geodesic segments
          if( Argument.size() == 0 ){
            if(!_quiet) std::cerr << "K>1 && empty argument: compute geodesic segments!" << std::endl;
            VectorType fullSegments( _numOfShapes * numGeodesicsDofs );
            computeGeodesics( fullSegments );
            for (int k = 0; k < _numOfShapes; k++)
              innerShapesRefs.push_back( fullSegments.segment(k * numGeodesicsDofs, numLocalDofs));
  
          }
          // something is given: check whether geodesic segments or only inner shapes
          else{
            // check if argument has right size
            if(!_quiet) std::cerr << "K>1: get inner shapes from argument!" << std::endl;
            if( (Argument.size() != _numOfShapes * numGeodesicsDofs) && (Argument.size() != _numOfShapes * numLocalDofs) )
                throw BasicException("GramianMatrixAssembler::pushTriplets(): argument has unvalid size!");
            
            int offset = (Argument.size() == _numOfShapes * numGeodesicsDofs) ? _K - 1 : 1;
            for (int k = 0; k < _numOfShapes; k++)
                innerShapesRefs.push_back( Argument.segment( k * offset * numLocalDofs, numLocalDofs));
          }
        }


        for (int i = 0; i < _numOfShapes; i++) {
          RealType Energy_i;
            _W.applyEnergy(_mean, innerShapesRefs[i], Energy_i);
            Energy_i *= _K * _K;
            tripletList.push_back(TripletType(i, i, Energy_i / _numOfShapes));
            for (int j = 0; j < i; j++) {
                tripletList.push_back(TripletType(i, j, Energy_i / (2. * _numOfShapes)));
                tripletList.push_back(TripletType(j, i, Energy_i / (2. * _numOfShapes)));

                RealType Energy;
                _W.applyEnergy(_mean, innerShapesRefs[j], Energy);
                Energy *= _K * _K;
                tripletList.push_back(TripletType(i, j, Energy / (2. * _numOfShapes)));
                tripletList.push_back(TripletType(j, i, Energy / (2. * _numOfShapes)));

                _W.applyEnergy(innerShapesRefs[i], innerShapesRefs[j], Energy);
                Energy *= _K * _K;
                tripletList.push_back(TripletType(i, j, -0.5 * Energy / (2. * _numOfShapes)));
                tripletList.push_back(TripletType(j, i, -0.5 * Energy / (2. * _numOfShapes)));

                _W.applyEnergy(innerShapesRefs[j], innerShapesRefs[i], Energy);
                Energy *= _K * _K;
                tripletList.push_back(TripletType(i, j, -0.5 * Energy / (2. * _numOfShapes)));
                tripletList.push_back(TripletType(j, i, -0.5 * Energy / (2. * _numOfShapes)));
            }

        }
    };

protected:

    void computeGeodesics(VectorType &geodesics) {

        const int numLocalDofs = _shapes.size() / _numOfShapes;
        const int numGeodesicsDofs = numLocalDofs * (_K - 1);

        //! \todo Use initialization of GeodesicInterpolation
        if (geodesics.size() != _numOfShapes * numGeodesicsDofs) {
            geodesics.resize(_numOfShapes * numGeodesicsDofs);
            for (int k = 0; k < _numOfShapes; k++) {
                for (int j = 0; j < _K - 1; j++)
                    geodesics.segment(k * numGeodesicsDofs + j * numLocalDofs, numLocalDofs) = _argRefs[k];
            }
        }

//        std::vector<Eigen::Ref<VectorType> > geodesicRef(_numOfShapes);
        std::vector<Eigen::Ref<VectorType> > geodesicRef;
        for (int k = 0; k < _numOfShapes; k++)
            geodesicRef.push_back(geodesics.segment(k * numGeodesicsDofs, numGeodesicsDofs));

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int k = 0; k < _numOfShapes; k++) {
            if(!_quiet) std::cerr << "GramianMatrixAssembler::computeGeodesics: Computing geodesic " << k << std::endl;
            VectorType geodesic(geodesicRef[k]);

            GeodesicInterpolation<ConfiguratorType> interpolationOp(_topology, _mean, _argRefs[k], _W, _K + 1,
                                                                         _optPars, _quiet);
            interpolationOp.setBoundaryMask(*_mask);
            interpolationOp.execute(geodesic, 2);
            geodesicRef[k] = geodesic;
        }
    }
};

#endif //GRAMIANMATRIX_H
