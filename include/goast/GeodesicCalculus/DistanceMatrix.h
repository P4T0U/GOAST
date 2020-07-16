// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Assembling of pairwise geodesic distance matrix.
 * \author Sassen
 * 
 */
#ifndef DISTANCEMATRIX_H
#define DISTANCEMATRIX_H

#include <goast/Core.h>
#include "DiscreteGeodesicCalculus.h"

/**
 * \brief For a selection of shapes \f$ s_1, \ldots, s_n \f$, compute pairwise geodesic distances.
 * \author Sassen
 *
 * In detail, assemble \f$ n \times n \f$-matrix \f$ d \f$ with \f$ d_{ij} = E[(s_0^{ij}, \ldots, s_K^{ij})] \f$,
 * where \f$ (s_0^{ij}, \ldots, s_K^{ij}) \f$ is a discrete geodesic connecting \f$ s_i \f$ and \f$ s_j \f$.
 */
template<typename ConfiguratorType>
class GeodesicDistanceMatrixAssembler {
protected:
    typedef typename ConfiguratorType::RealType RealType;

    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;
    typedef std::vector<TripletType> TripletListType;

    const MeshTopologySaver &_topology;
    const DeformationBase<ConfiguratorType> &_W;
    const VectorType &_shapes;

    const unsigned int _numOfShapes;
    std::vector<Eigen::Ref<const VectorType> > _argRefs;
    const int _K;
    const int _numLocalDofs;

    bool _quiet;

    const OptimizationParameters<ConfiguratorType> &_optPars;

    std::vector<int> const *_mask;

public:
    GeodesicDistanceMatrixAssembler(const MeshTopologySaver &Topology,
                                    const DeformationBase<ConfiguratorType> &W,
                                    const VectorType &shapes,
                                    const unsigned int numOfShapes,
                                    const int K,
                                    const OptimizationParameters<ConfiguratorType> &optPars,
                                    bool quiet = true)
            : _topology(Topology), _W(W), _shapes(shapes), _numOfShapes(numOfShapes), _K(K),
              _optPars(optPars), _quiet(quiet), _mask(NULL), _numLocalDofs(_shapes.size() / _numOfShapes) {
        // Create references for this different shapes bec. of convenience
        _argRefs.reserve(numOfShapes);

            for (int k = 0; k < numOfShapes; k++)
                _argRefs.push_back(_shapes.segment(k * _numLocalDofs, _numLocalDofs));
    }

    ~GeodesicDistanceMatrixAssembler() {
        if (!_quiet) std::cerr << "GeodesicDistanceMatrixAssembler deleted." << std::endl;
    }

    void setBoundaryMask(const std::vector<int> &localMask) {
        _mask = &localMask;
//        fillPathMask(1, 3 * _topology.getNumVertices(), localMask, *_mask);
    }

    void assemble(MatrixType &Dest) {
        if ((Dest.rows() != _numOfShapes) || (Dest.cols() != _numOfShapes))
            Dest.resize(_numOfShapes, _numOfShapes);

        // fill triplet lists
        TripletListType tripletList;
        tripletList.reserve(_numOfShapes * _numOfShapes);

        pushTriplets(tripletList);
        Dest.setFromTriplets(tripletList.cbegin(), tripletList.cend());
    };

    void pushTriplets(TripletListType &tripletList) {

        const int numGeodesicsDofs = _numLocalDofs * (_K - 1);


        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 2)
        #endif
        for (int k = 0; k < _numOfShapes * (_numOfShapes + 1) / 2; k++) {
            int i = k / (_numOfShapes + 1), j = k % (_numOfShapes + 1);
            if (j > i) i = _numOfShapes - i - 1, j = _numOfShapes - j;
            if (i == j) continue;

            if (!_quiet) std::cerr << "Computing geodesic " << i << " -> " << j << std::endl;
            VectorType geodesic, Energies, StartGeom, EndGeom;
          RealType Energy;
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                StartGeom = _argRefs[i];
                EndGeom = _argRefs[j];
            };

            RealType length;

            if (_K > 1) {
                GeodesicInterpolation<ConfiguratorType> interpolationOp(_topology, StartGeom, EndGeom, _W,
                                                                             _K + 1, _optPars, _quiet);
                interpolationOp.setBoundaryMask(*_mask);
                interpolationOp.execute(geodesic, 2);

                DiscretePathEnergy<ConfiguratorType> E(_W, _K, StartGeom, EndGeom);
//                E.evaluateSingleEnergies(geodesic, Energies);
//                length = std::pow(Energies.array().abs().sqrt().sum(), 2);
                E.apply(geodesic, Energy);
                length = Energy;
            }
            else {
                _W.applyEnergy(StartGeom, EndGeom, Energy);
                length = Energy;
            }

            if (!_quiet) std::cerr << "Length(" << i << ", " << j << ") : " << length << std::endl;
            if (!_quiet) std::cerr << "Energies(" << i << ", " << j << ") : " << Energies << std::endl;

#ifdef _OPENMP
#pragma omp critical
#endif
            {
                tripletList.push_back(TripletType(i, j, length));
                tripletList.push_back(TripletType(j, i, length));
            }
        }

    };

protected:


};

#endif //DISTANCEMATRIX_H
