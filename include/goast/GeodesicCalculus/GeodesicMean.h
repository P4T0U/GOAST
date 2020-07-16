// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Geodesic mean computation.
 * \author Heeren, Sassen
 *
 */
// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 [INSERTCOPYRIGHTHOLDER]
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Geodesic mean computation.
 * \author Heeren, Sassen
 *
 */
#ifndef GEODESICMEAN_H
#define GEODESICMEAN_H

#include <vector>

#include <goast/Core.h>
#include <goast/GeodesicCalculus.h>
#include <goast/Optimization.h>


/**
 * \brief Discrete geodesic mean functional
 * \author Heeren, Sassen
 *
 * Given \f$ n \f$ example shapes \f$ s_1,\dots,s_n \f$, weights  \f$ \alpha_1,\ldots,\alpha_n \f$
 * this class implements the geodesic average energy \f[ F[s] = \sum_{i=1}^n \alpha_i \sum_{k=1}^K W[s^i_{k-1}, s^i_k] \, ,\f]
 * where \f$ s^i_0 = s_i \f$ and \f$ s^i_K = s \f$ for \f$ i = 1, \ldots, n \f$.
 */
template<typename ConfiguratorType>
class GeodesicMeanFunctional
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;

    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;

    const DeformationBase<ConfiguratorType> &_W;
    const VectorType &_shapes, _alpha;
    const unsigned int _numOfShapes;
    const unsigned int _numOfFreeShapes;
    const unsigned int _numLocalDofs;
    std::vector<Eigen::Ref<const VectorType> > _argRefs;
    int _K;

public:
    GeodesicMeanFunctional(const DeformationBase<ConfiguratorType> &W,
                           int K,
                           const VectorType &shapes,
                           const VectorType alpha,
                           const unsigned int numOfShapes)
            : _W(W), _shapes(shapes), _alpha(alpha), _numOfShapes(numOfShapes), _K(K),
              _numOfFreeShapes(numOfShapes * (K-1) + 1), _numLocalDofs( shapes.size() / _numOfShapes )
    {
        if (shapes.size() % numOfShapes != 0)
            throw BasicException("GeodesicMeanFunctional: wrong number of dof!");
        if (alpha.size() % numOfShapes != 0)
            throw BasicException("GeodesicMeanFunctional: wrong number of alphas!");

        // Create references for this different shapes bec. of convenience
        _argRefs.reserve(numOfShapes);

        for (int k = 0; k < numOfShapes; k++)
            _argRefs.push_back(_shapes.segment(k * _numLocalDofs, _numLocalDofs));
    }

    //! \param Arg The mean shape s and the different free shapes of the geodesics s_1^1,...,s_{K-1}^1,s_1^1,...,s_{K-1}^n
    void apply(const VectorType &Arg, RealType &Dest) const override {
        if(Arg.size() % _numOfFreeShapes != 0)
            throw BasicException("GeodesicMeanFunctional::apply: wrong number of dofs!");

        if(Arg.size() / _numOfFreeShapes != _argRefs[0].size())
            throw BasicException("GeodesicMeanFunctional::apply: wrong size of dofs!");

        Dest = 0.;
        RealType factor = 1. * _K / _numOfShapes;

        // bring into more convenient form
        std::vector< Eigen::Ref<const VectorType> > argRefsFree;
        argRefsFree.reserve(_numOfFreeShapes);
        const int numLocalDofs = Arg.size() / _numOfFreeShapes;

        for( int k = 0; k < _numOfFreeShapes; k++ )
            argRefsFree.push_back( Arg.segment(k*numLocalDofs, numLocalDofs) );

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int k = 0; k < _numOfShapes; k++) {
            RealType Energy = 0.;

            _W.applyAddEnergy(argRefsFree[0], argRefsFree[1 + k * (_K - 1)], Energy);

            for (int j = 1; j < _K - 1; j++)
                _W.applyAddEnergy(argRefsFree[1 + k * (_K - 1) + j - 1], argRefsFree[1 + k * (_K - 1) + j], Energy);

            _W.applyAddEnergy(argRefsFree[1 + k * (_K - 1) + (_K - 2)], _argRefs[k], Energy);
#pragma omp critical
            {
                Dest += _alpha[k] * factor * Energy;
            }
        }
    }

    //! \param Arg The mean shape s and the different free shapes of the geodesics s_1^1,...,s_{K-1}^1,s_1^1,...,s_{K-1}^n
    void apply(const std::vector<VectorType> &Arg, RealType &Dest) const {
        if( Arg.size() != 1 + _numOfShapes )
            throw BasicException("GeodesicMeanFunctional::apply: wrong size of argument!");
        if( Arg[0].size() != _numLocalDofs)
            throw BasicException("GeodesicMeanFunctional::apply: wrong size of mean!");
        for( int n = 0; n < _numOfShapes; n++ )
            if(Arg[1+n].size() != (_K-1) * _numLocalDofs )
                throw BasicException("GeodesicMeanFunctional::apply: wrong size of segments!");

        Dest = 0.;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int n = 0; n < _numOfShapes; n++) {
            RealType Energy = 0.;
            const Eigen::Ref<const VectorType> endShape = _shapes.segment( n * _numLocalDofs, _numLocalDofs );
            DiscretePathEnergy<ConfiguratorType>( _W, _K, Arg[0], endShape ).apply ( Arg[n+1], Energy );
#pragma omp critical
            {
                Dest += _alpha[n] * Energy / _numOfShapes;
            }
        }
    }
};

//! \brief Gradient of GeodesicMeanFunctional
//! \author Heeren, Sassen
template<typename ConfiguratorType>
class GeodesicMeanFunctionalGradient
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

protected:
    typedef typename ConfiguratorType::RealType RealType;

    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;

    const DeformationBase<ConfiguratorType> &_W;
    const VectorType &_shapes, _alpha;
    const unsigned int _numOfShapes;
    const unsigned int _numOfFreeShapes;
    const unsigned int _numLocalDofs;
    std::vector<Eigen::Ref<const VectorType> > _argRefs;
    int _K;

public:
    GeodesicMeanFunctionalGradient(const DeformationBase<ConfiguratorType> &W,
                                   int K,
                                   const VectorType &shapes,
                                   const VectorType &alpha,
                                   const unsigned int numOfShapes)
            : _W(W), _shapes(shapes), _alpha(alpha), _numOfShapes(numOfShapes), _K(K),
              _numOfFreeShapes(numOfShapes * (K-1) + 1), _numLocalDofs(shapes.size()/numOfShapes) {
        if (shapes.size() % numOfShapes != 0)
            throw BasicException("GeodesicMeanFunctionalGradient: wrong number of dof!");
        if (alpha.size() % numOfShapes != 0)
            throw BasicException("GeodesicMeanFunctionalGradient: wrong number of alphas!");

        // Create references for this different shapes bec. of convenience
        _argRefs.reserve(numOfShapes);

        const int numLocalDofs = _shapes.size() / _numOfShapes;
        for (int k = 0; k < numOfShapes; k++)
            _argRefs.push_back(_shapes.segment(k * numLocalDofs, numLocalDofs));
    }

    void apply(const VectorType &Arg, VectorType &Dest) const override {
        if(Arg.size() % _numOfFreeShapes != 0)
            throw BasicException("GeodesicMeanFunctionalGradient::apply: wrong number of dofs!");

        if(Arg.size() / _numOfFreeShapes != _argRefs[0].size())
            throw BasicException("GeodesicMeanFunctionalGradient::apply: wrong size of dofs!");

        RealType factor = 1. * _K / _numOfShapes;

        Dest.resize(Arg.size());
        Dest.setZero();

        // bring into more convenient form
        std::vector< Eigen::Ref<const VectorType> > argRefsFree;
        argRefsFree.reserve(_numOfFreeShapes);
        for( int k = 0; k < _numOfFreeShapes; k++ )
            argRefsFree.push_back( Arg.segment(k*_numLocalDofs, _numLocalDofs) );

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int n = 0; n < _numOfShapes; n++) {
            VectorType Gradient(Arg.size());
            Gradient.setZero();

            _W.applyAddUndefGradient(argRefsFree[0], argRefsFree[1 + n * (_K - 1)], Gradient.segment(0, _numLocalDofs));
            _W.applyAddDefGradient(  argRefsFree[0], argRefsFree[1 + n * (_K - 1)], Gradient.segment((1 + n * (_K - 1)) * _numLocalDofs, _numLocalDofs));

            for (int k = 1; k < _K - 1; k++) {
                _W.applyAddUndefGradient(argRefsFree[1 + n * (_K - 1) + k - 1], argRefsFree[1 + n * (_K - 1) + k],
                                         Gradient.segment((1 + n * (_K - 1) + k - 1) * _numLocalDofs, _numLocalDofs));
                _W.applyAddDefGradient(argRefsFree[1 + n * (_K - 1) + k - 1], argRefsFree[1 + n * (_K - 1) + k],
                                       Gradient.segment((1 + n * (_K - 1) + k) * _numLocalDofs, _numLocalDofs));
            }

            // _W.applyAddEnergy(argRefsFree[1 + k * (_K - 1) + (_K - 2)], _argRefs[k], Energy);
            _W.applyAddUndefGradient(argRefsFree[1 + n * (_K - 1) + (_K - 2)], _argRefs[n],
                                     Gradient.segment((1 + n * (_K - 1) + (_K - 2)) * _numLocalDofs, _numLocalDofs));


#ifdef _OPENMP
#pragma omp critical
#endif
            Dest += _alpha[n] * factor * Gradient;
        }
    }

    void apply(const std::vector<VectorType> &Arg, VectorType &Dest) const {
        if( Arg.size() != 1 + _numOfShapes )
            throw BasicException("GeodesicMeanFunctionalGradient::apply: wrong size of argument!");
        if( Arg[0].size() != _numLocalDofs)
            throw BasicException("GeodesicMeanFunctionalGradient::apply: wrong size of mean!");
        for( int n = 0; n < _numOfShapes; n++ )
            if(Arg[1+n].size() != (_K-1) * _numLocalDofs )
                throw BasicException("GeodesicMeanFunctionalGradient::apply: wrong size of segments!");

        RealType factor = 1. * _K / _numOfShapes;

        Dest.resize( _numOfFreeShapes * _numLocalDofs );
        Dest.setZero();

//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
        for (int n = 0; n < _numOfShapes; n++) {
            VectorType Gradient( _numOfFreeShapes * _numLocalDofs );
            Gradient.setZero();

            _W.applyAddUndefGradient(Arg[0], Arg[1 + n].segment(0, _numLocalDofs), Gradient.segment(0, _numLocalDofs));
            _W.applyAddDefGradient(  Arg[0], Arg[1 + n].segment(0, _numLocalDofs), Gradient.segment((1 + n * (_K - 1)) * _numLocalDofs, _numLocalDofs));

            for (int k = 1; k < _K-1; k++) {
                _W.applyAddUndefGradient(Arg[1 + n].segment( (k-1)*_numLocalDofs, _numLocalDofs), Arg[1 + n].segment( k*_numLocalDofs, _numLocalDofs), Gradient.segment((1 + n * (_K - 1) + k - 1) * _numLocalDofs, _numLocalDofs));
                _W.applyAddDefGradient(Arg[1 + n].segment( (k-1)*_numLocalDofs, _numLocalDofs), Arg[1 + n].segment( k*_numLocalDofs, _numLocalDofs), Gradient.segment((1 + n * (_K - 1) + k) * _numLocalDofs, _numLocalDofs));
            }

            _W.applyAddUndefGradient(Arg[1 + n].segment( (_K-2)*_numLocalDofs, _numLocalDofs), _shapes.segment( n*_numLocalDofs, _numLocalDofs), Gradient.segment( (1 + n * (_K - 1) + (_K - 2) )* _numLocalDofs, _numLocalDofs));
//#ifdef _OPENMP
//#pragma omp critical
//#endif
            Dest += _alpha[n] * factor * Gradient;
            std::cerr << n << "th shape has partial derivative norm of " << Gradient.norm() << std::endl;
        }
    }
};

/**
 * \brief Compute discrete geodesic mean by minimizing GeodesicMeanFunctional
 * \author Heeren, Sassen
 *
 * Compute discrete geodesic mean by minimizing the geodesic mean functional.
 * For \f$ n \f$ input shapes and chosen \f$ K \f$, we have \f$ n \times (K-1) + 1\f$ free shapes.
 * Two different optimization strategies are distinguished:
 *
 * a) optimize all ( in execute() ): full optimization of all free shapes simultaneously.
 *
 * b) alternating optimization ( in executeAlternating() ): first fixing mean \f$ s and optimize all geodesic segments parallely,
 *    then update mean s based on fixed segments as elastic mean of discrete logar \f$ithms.
 *
 */
template<typename ConfiguratorType>
class GeodesicMean {
protected:
    typedef typename ConfiguratorType::RealType RealType;

    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;

    const MeshTopologySaver &_topology;
    const DeformationBase<ConfiguratorType> &_W;
    const VectorType &_shapes, _alpha;
    const unsigned int _numOfShapes;
    std::vector<Eigen::Ref<const VectorType> > _argRefs;
    const int _K;
    bool _quiet, _innerQuiet;
    bool _overwrite;
    int _saveMode; // 0 = no saving, 1 = mean, 2 = mean and log shapes, 3 = all shapes

    const OptimizationParameters<ConfiguratorType>& _optPars;

    std::vector<int> *_mask;
    std::vector<int> *_singleMask;
    std::vector<int> *_geodMask;

public:
    GeodesicMean(const MeshTopologySaver &Topology,
                 const DeformationBase<ConfiguratorType> &W,
                 const int K,
                 const VectorType &shapes,
                 const VectorType alpha,
                 const unsigned int numOfShapes,
                 const OptimizationParameters<ConfiguratorType>& optPars,
                 bool quiet = true,
                 int saveMode = 1 )
            : _topology(Topology), _W(W), _shapes(shapes), _alpha(alpha), _numOfShapes(numOfShapes),
              _optPars(optPars), _mask(NULL), _singleMask(NULL), _geodMask(NULL), _K(K), _quiet(quiet), _innerQuiet(true), _overwrite(true), _saveMode(saveMode) {

        // Create references for this different shapes bec. of convenience
        _argRefs.reserve(numOfShapes);

        const int numLocalDofs = _shapes.size() / _numOfShapes;
        for (int k = 0; k < numOfShapes; k++)
            _argRefs.push_back(_shapes.segment(k * numLocalDofs, numLocalDofs));
    }

    ~GeodesicMean() {
        if (_mask) {
            if(!_quiet) std::cerr << "Delete boundary mask." << std::endl;
            delete _mask;
        }
        if( _geodMask ){
            if(!_quiet) std::cerr << "Delete geodesic boundary mask." << std::endl;
            delete _geodMask;
        }
        if( _singleMask ){
            if(!_quiet) std::cerr << "Delete single boundary mask." << std::endl;
            delete _singleMask;
        }
    }

    void setBoundaryMask(const std::vector<int> &localMask) {
        if( localMask.size() > 0 ){
            _mask = new std::vector<int>;
            _geodMask = new std::vector<int>;
            _singleMask = new std::vector<int>;
            fillPathMask(_numOfShapes * (_K - 1) + 1, 3 * _topology.getNumVertices(), localMask, *_mask);
            fillPathMask(_numOfShapes * (_K - 1), 3 * _topology.getNumVertices(), localMask, *_geodMask);
            fillPathMask(1, 3 * _topology.getNumVertices(), localMask, *_singleMask);
        }
    }

    //! \param Arg The mean shape s and the different free shapes of the geodesics s_1^1,...,s_{K-1}^1, s_1^1,...,s_{K-1}^n going from the mean to the input shapes
    void executeAlternating(VectorType &geodMeanAndSegments, const unsigned int numIter, std::string savenameStem = "" ) const {

        const int numLocalDofs    = 3 * _topology.getNumVertices();
        const int numOfTotalFreeShapes = 1 + _numOfShapes * (_K-1);
        const int numOfTotalDofs  = numLocalDofs * numOfTotalFreeShapes;
        const int numGeodesicDofs = (_K-1) * numLocalDofs;

        if( geodMeanAndSegments.size() != numOfTotalDofs )
            throw BasicException("GeodesicMean::executeAlternating(): argument has wrong size!");

        bool saving( savenameStem.size() > 0 && _saveMode > 0 );

        GeodesicMeanFunctional<ConfiguratorType> E(_W, _K, _shapes, _alpha, _numOfShapes);
        GeodesicMeanFunctionalGradient<ConfiguratorType> DE(_W, _K, _shapes, _alpha, _numOfShapes);

        // define auxiliary variables
        TriMesh outputMesh ( _topology.getGrid() );
        RealType energy;
        VectorType grad;

        // initial output
        if( !_quiet ){
            E.apply(geodMeanAndSegments, energy);
            std::cerr << "Initial geodesic mean energy  = " << energy << std::endl;
            DE.apply(geodMeanAndSegments, grad);
            std::cerr << "Initial geodesic mean gradient norm = " << grad.norm() << std::endl << std::endl;
        }

        // get geodesic mean
        VectorType geodMean(geodMeanAndSegments.segment(0, numLocalDofs));

        // references to the different geodesics
        std::vector<VectorType> geodRefs;
        geodRefs.reserve(_numOfShapes);
        if(!_quiet) {
            std::cerr << "DOFs per segment: " << _K-1 << " x " << numLocalDofs << " = " << numGeodesicDofs << std::endl;
            std::cerr << "DOFs total = " << numOfTotalDofs << std::endl;
        }

        // shapes closest to the mean for each geodesic, i.e. (s_1^k)_k
        VectorType innerShapes(_numOfShapes * numLocalDofs);
        for( int k = 0; k < _numOfShapes; k++) {
            geodRefs.push_back(geodMeanAndSegments.segment(numLocalDofs + k * numGeodesicDofs, numGeodesicDofs));
            innerShapes.segment(k * numLocalDofs, numLocalDofs) = geodRefs[k].segment(0, numLocalDofs);
        }

        // initialization scheme for initial geodesic (default is linear nodal interpolation)
        int initializationScheme =  _optPars.getInitializationScheme() > 0 ? _optPars.getInitializationScheme() : 2;

        // start alternating optimization
        if(!_quiet) std::cerr << "\n\nSTART ALTERNATING OPTIMIZATION." << std::endl;
        if(!_quiet) std::cerr << "===================================="<< std::endl;
        for (int j = 0; j < numIter; j++) {

            if(!_quiet) std::cerr << "Start " << j+1 << "th outer iteration of " << numIter << std::endl;

            // UPDATE MEAN
            if(!_quiet) std::cerr << "Updating mean..." << std::endl;

            ElasticMean<ConfiguratorType> elasticMeanOp(_topology, _W, innerShapes, _alpha, _numOfShapes, _optPars, _innerQuiet);
            if( _singleMask )
                elasticMeanOp.setBoundaryMask(*_singleMask);
            elasticMeanOp.execute(geodMean);

            // UPDATE GEODESIC SEGMENTS (IN PARALLEL)
            if(!_quiet) std::cerr << "Updating geodesics..." << std::endl;
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for( int k = 0; k < _numOfShapes; k++) {
                if(!_quiet) std::cerr << "Iteration " << j+1 << ", updating geodesic segment " << k << std::endl;

                // Somehow not copying leads to segfaults..
                VectorType startGeom(_argRefs[k]);

                GeodesicInterpolation<ConfiguratorType> interpolationOp(_topology, geodMean, startGeom, _W, _K + 1, _optPars, _innerQuiet);
                if( _singleMask )
                    interpolationOp.setBoundaryMask(*_singleMask);
                if (j == 0)
                    interpolationOp.execute(geodRefs[k], initializationScheme );
                else
                    interpolationOp.execute(geodRefs[k]);

                //interpolationOp.checkIfGeodesicIsRelaxed( geodRefs[k] );
            }

            // update inner shapes
            for( int k = 0; k < _numOfShapes; k++)
                innerShapes.segment(k * numLocalDofs, numLocalDofs) = geodRefs[k].segment(0, numLocalDofs);

            // save results according to save mode
            if( saving ){
                if(!_quiet) std::cerr << "Saving..." << std::endl;
                saveResults( geodMean, geodRefs,  savenameStem, j );
            }

            // COMPUTE GEODESIC MEAN ENERGY
            if(!_quiet){
                // write back
                geodMeanAndSegments.segment(0, numLocalDofs) = geodMean;
                for( int k = 0; k < _numOfShapes; k++)
                    geodMeanAndSegments.segment(numLocalDofs + k * numGeodesicDofs, numGeodesicDofs) = geodRefs[k];

                E.apply(geodMeanAndSegments, energy);
                std::cerr << "Geodesic mean energy = " << energy << std::endl;
                std::cerr << "====================================\n\n";
            }

//            if(!_quiet) {
//              std::cerr << "Check segments..." << std::endl;
//              checkIfAllSegmentsAreRelaxed( geodMeanAndSegments );
//            }
        }

        // write back
        geodMeanAndSegments.segment(0, numLocalDofs) = geodMean;
        for( int k = 0; k < _numOfShapes; k++)
            geodMeanAndSegments.segment(numLocalDofs + k * numGeodesicDofs, numGeodesicDofs) = geodRefs[k];

        // report final energy
        if( !_quiet ){
            E.apply(geodMeanAndSegments, energy);
            std::cerr << "Final geodesic mean energy = " << energy << std::endl;
            DE.apply(geodMeanAndSegments, grad);
            std::cerr << "Final geodesic mean gradient norm = " << grad.norm() << std::endl << std::endl;
        }
    }

    //! \param Arg First entry in vector is the mean shape which is d-dimensional (here d is three times number of nodes),
    //! (n+1)th entry is discrete geodesics (s_1^n,...,s_{K-1}^n) going from the mean to the nth input shape (n = 1, ..., N),
    //! i.e. the (n+1)th entry has dimension (K-1)*d.
    void executeAlternating(std::vector<VectorType> &geodMeanAndSegments, const unsigned int numIter, std::string savenameStem = "" ) const {

        const int numLocalDofs    = 3 * _topology.getNumVertices();
        const int numOfTotalFreeShapes = 1 + _numOfShapes * (_K-1);
        const int numOfTotalDofs  = numLocalDofs * numOfTotalFreeShapes;
        const int numGeodesicDofs = (_K-1) * numLocalDofs;

        if( geodMeanAndSegments.size() != 1 + _numOfShapes )
            throw BasicException("GeodesicMean::executeAlternating: wrong size of argument!");
        if( geodMeanAndSegments[0].size() != numLocalDofs)
            throw BasicException("GeodesicMean::executeAlternating: wrong size of mean!");
        for( int n = 0; n < _numOfShapes; n++ )
            if(geodMeanAndSegments[1+n].size() != numGeodesicDofs )
                throw BasicException("GeodesicMean::executeAlternating: wrong size of segments!");

        bool saving( savenameStem.size() > 0 && _saveMode > 0 );

        // define auxiliary variables
        TriMesh outputMesh ( _topology.getGrid() );
        RealType energy;
        VectorType  grad;
        GeodesicMeanFunctional<ConfiguratorType> E(_W, _K, _shapes, _alpha, _numOfShapes);
        GeodesicMeanFunctionalGradient<ConfiguratorType> DE(_W, _K, _shapes, _alpha, _numOfShapes);

        // initial output
        if( !_quiet ){
            E.apply(geodMeanAndSegments, energy);
            std::cerr << "Initial geodesic mean energy  = " << energy << std::endl;
            DE.apply(geodMeanAndSegments, grad);
            std::cerr << "Initial geodesic mean gradient norm = " << grad.norm() << std::endl << std::endl;
        }

        // discrete logarithm shapes, i.e. shapes closest to the mean for each geodesic, i.e. (s_1^k)_k
        VectorType logShapes(_numOfShapes * numLocalDofs);
        for( int k = 0; k < _numOfShapes; k++)
            logShapes.segment(k * numLocalDofs, numLocalDofs) = geodMeanAndSegments[1+k].segment(0, numLocalDofs);

        // initialization scheme for initial geodesic (default is linear nodal interpolation)
        //int initializationScheme =  _optPars.getInitializationScheme() > 0 ? _optPars.getInitializationScheme() : 2;

        // start alternating optimization
        if(!_quiet) std::cerr << "\n\nSTART ALTERNATING OPTIMIZATION." << std::endl;
        if(!_quiet) std::cerr << "===================================="<< std::endl;
        for (int j = 0; j < numIter; j++) {

            if(!_quiet) std::cerr << "Start " << j+1 << "th outer iteration of " << numIter << std::endl;

            // UPDATE MEAN
            if(!_quiet) std::cerr << "Updating mean..." << std::endl;

            ElasticMean<ConfiguratorType> elasticMeanOp(_topology, _W, logShapes, _alpha, _numOfShapes, _optPars, _innerQuiet);
            if( _singleMask )
                elasticMeanOp.setBoundaryMask(*_singleMask);
            elasticMeanOp.execute(geodMeanAndSegments[0]);

            // UPDATE GEODESIC SEGMENTS (IN PARALLEL)
            if(!_quiet) std::cerr << "Updating geodesics..." << std::endl;
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for( int k = 0; k < _numOfShapes; k++) {
                if(!_innerQuiet) std::cerr << "Iteration " << j+1 << ", updating geodesic segment " << k << std::endl;

                // Somehow not copying leads to segfaults..
                VectorType startGeom(_argRefs[k]);

                GeodesicInterpolation<ConfiguratorType> interpolationOp(_topology, geodMeanAndSegments[0], startGeom, _W, _K + 1, _optPars, true);
                if( _singleMask )
                    interpolationOp.setBoundaryMask(*_singleMask);
                interpolationOp.execute(geodMeanAndSegments[1+k]);

                //interpolationOp.checkIfGeodesicIsRelaxed( geodRefs[k] );
            }

            // update inner shapes
            for( int k = 0; k < _numOfShapes; k++)
                logShapes.segment(k * numLocalDofs, numLocalDofs) = geodMeanAndSegments[1+k].segment(0, numLocalDofs);

            // save results according to save mode
            if( saving ){
                if(!_quiet) std::cerr << "Saving..." << std::endl;
                saveResults( geodMeanAndSegments,  savenameStem, j );
            }

            // COMPUTE GEODESIC MEAN ENERGY
            if(!_quiet){
                E.apply(geodMeanAndSegments, energy);
                std::cerr << "Geodesic mean energy = " << energy << std::endl;
                DE.apply(geodMeanAndSegments, grad);
                std::cerr << "Geodesic mean gradient norm = " << grad.norm() << std::endl << std::endl;
                std::cerr << "====================================\n\n";
            }
        }

        // report final energy
        if( !_quiet ){
            E.apply(geodMeanAndSegments, energy);
            std::cerr << "Final geodesic mean energy = " << energy << std::endl;
            DE.apply(geodMeanAndSegments, grad);
            std::cerr << "Final geodesic mean gradient norm = " << grad.norm() << std::endl << std::endl;
        }

        if(!_quiet) {
            std::cerr << "Check segments..." << std::endl;
            checkIfAllSegmentsAreRelaxed( geodMeanAndSegments );
        }
    }

    //! \param Arg The mean shape s and the different free shapes of the geodesics s_1^1,...,s_{K-1}^1, s_1^1,...,s_{K-1}^n going from the mean to the input shapes
    void execute( VectorType &geodMeanAndSegments ) const {

        // define functional and derivatives
        GeodesicMeanFunctional<ConfiguratorType> E(_W, _K, _shapes, _alpha, _numOfShapes);
        GeodesicMeanFunctionalGradient<ConfiguratorType> DE(_W, _K, _shapes, _alpha, _numOfShapes);
        //GeodesicMeanFunctionalHessian<ConfiguratorType> D2E(_W, _K, _shapes, _alpha, _numOfShapes);

        VectorType initialization(geodMeanAndSegments);
        RealType energy;
        VectorType grad;

        // initial output
        if( !_quiet ){
            E.apply(initialization, energy);
            std::cerr << "Initial functional  = " << energy << std::endl;
            DE.apply(initialization, grad);
            std::cerr << "Initial gradient norm = " << grad.norm() << std::endl << std::endl;
        }

        // optimization with gradient descent
        if (_optPars.getGradientIterations() > 0) {
            if(!_quiet) std::cerr << "Start gradient descent with " << _optPars.getGradientIterations() << " steps... " << std::endl;
            GradientDescent<ConfiguratorType> GD(E, DE, _optPars);
            if (_mask) {
                if(!_quiet) std::cerr << "Set boundary mask... " << std::endl;
                GD.setBoundaryMask(*_mask);
            }
            GD.solve(initialization, geodMeanAndSegments);
            initialization = geodMeanAndSegments;
        }

        // optimization with BFGS
        if (_optPars.getBFGSIterations() > 0) {
            if(!_quiet) std::cerr << "Start Quasi-Newton with " << _optPars.getBFGSIterations() << " steps... " << std::endl;
            QuasiNewtonBFGS<ConfiguratorType> QNM(E, DE, _optPars);
            if (_mask) {
                if(!_quiet) std::cerr << "Set boundary mask... " << std::endl;
                QNM.setBoundaryMask(*_mask);
            }
            QNM.solve(initialization, geodMeanAndSegments);
            initialization = geodMeanAndSegments;
        }

        // final output
        if( !_quiet) {
            E.apply(geodMeanAndSegments, energy);
            std::cerr << "Final energy = " << energy << std::endl;
            DE.apply(geodMeanAndSegments, grad);
            std::cerr << "Final gradient norm = " << grad.norm() << std::endl << std::endl;
        }
    }

    // save geodesic mean and segments
    void saveResults( const VectorType &geodMean, const std::vector<VectorType>& geodSegments, std::string savenameStem, int iter ) const {
        const int numLocalDofs = 3 * _topology.getNumVertices();
        const int numOfTotalFreeShapes = 1 + _numOfShapes * (_K - 1);
        const int numOfTotalDofs = numLocalDofs * numOfTotalFreeShapes;

        if (geodMean.size() != numLocalDofs)
            throw BasicException("GeodesicMean::saveResults(): mean has wrong size!");
        if (geodSegments.size() != _numOfShapes)
            throw BasicException("GeodesicMean::saveResults(): segment vector has wrong size!");

        if( _saveMode == 0 )
            return;

        TriMesh outputMesh(_topology.getGrid());

        std::ostringstream iterString;
        if (!_overwrite)
            iterString << "_iter" << iter;

        // save geodesic mean
        std::ostringstream meanName;
        meanName << savenameStem << "_geodMean" << iterString.str() << ".ply";
        setGeometry(outputMesh, geodMean);
        OpenMesh::IO::write_mesh(outputMesh, meanName.str());

        if( _saveMode == 1 )
            return;

        // only save log shapes
        if( _saveMode == 2 ){
            for (int i = 0; i < _numOfShapes; i++) {
                std::ostringstream saveName;
                saveName << savenameStem << iterString.str() << "_segment" << i << "_shape1.ply";
                setGeometry(outputMesh, geodSegments[i].segment(0, numLocalDofs) );
                OpenMesh::IO::write_mesh(outputMesh, saveName.str());
            }
        }

        // save all segments
        if (_saveMode == 3 ) {
            for (int i = 0; i < _numOfShapes; i++) {
                // mean shape
                std::ostringstream startName;
                startName << savenameStem << iterString.str() << "_segment" << i << "_shape0.ply";
                setGeometry(outputMesh, geodMean);
                OpenMesh::IO::write_mesh(outputMesh, startName.str());
                // segments
                for (int j = 0; j < _K - 1; j++) {
                    std::ostringstream saveName;
                    saveName << savenameStem << iterString.str() << "_segment" << i << "_shape" << j + 1 << ".ply";
                    setGeometry(outputMesh,geodSegments[i].segment(j*numLocalDofs, numLocalDofs));
                    OpenMesh::IO::write_mesh(outputMesh, saveName.str());
                }
                // input shape
                std::ostringstream inputName;
                inputName << savenameStem << iterString.str() << "_segment" << i << "_shape" << _K << ".ply";
                setGeometry(outputMesh, _shapes.segment(i * numLocalDofs, numLocalDofs));
                OpenMesh::IO::write_mesh(outputMesh, inputName.str());
            }
        }
    }

    // save geodesic mean and segments
    void saveResults( const std::vector<VectorType> &geodMeanAndSegments, std::string savenameStem, int iter ) const {
        const int numLocalDofs = 3 * _topology.getNumVertices();
        const int numGeodDofs = numLocalDofs * (_K - 1);

        if (geodMeanAndSegments.size() != 1 + _numOfShapes )
            throw BasicException("GeodesicMean::saveResults(): argument has wrong size!");
        if( geodMeanAndSegments[0].size() != numLocalDofs )
            throw BasicException("GeodesicMean::saveResults(): mean has wrong size!");

        if( _saveMode == 0 )
            return;

        TriMesh meanMesh(_topology.getGrid());

        std::ostringstream iterString;
        if (!_overwrite)
            iterString << "_iter" << iter;

        // save geodesic mean
        std::ostringstream meanName;
        meanName << savenameStem << "_geodMean" << iterString.str() << ".ply";
        setGeometry(meanMesh, geodMeanAndSegments[0]);
        OpenMesh::IO::write_mesh(meanMesh, meanName.str());

        if( _saveMode == 1 )
            return;

        // only save log shapes
        if( _saveMode == 2 ){
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < _numOfShapes; i++) {
                if( geodMeanAndSegments[1+i].size() != numGeodDofs )
                    throw BasicException("GeodesicMean::saveResults(): segment has wrong size!");
                TriMesh currMesh(_topology.getGrid());
                std::ostringstream saveName;
                saveName << savenameStem << iterString.str() << "_segment" << i << "_shape1.ply";
                setGeometry(currMesh, geodMeanAndSegments[1+i].segment(0, numLocalDofs) );
                OpenMesh::IO::write_mesh(currMesh, saveName.str());
            }
        }

        // save all segments
        if (_saveMode == 3 )
            std::cerr << "GeodesicMean::saveResults() WARNING: saveMode 3 has not been implemented for this argument representation! Either choose different mode or different argument representation!\n";
    }

    //
    void checkIfAllSegmentsAreRelaxed( const VectorType &geodMeanAndSegments ) const {

        const int numLocalDofs    = 3 * _topology.getNumVertices();
        const int numOfTotalFreeShapes = 1 + _numOfShapes * (_K-1);
        const int numOfTotalDofs  = numLocalDofs * numOfTotalFreeShapes;
        const int numGeodesicDofs = (_K-1) * numLocalDofs;

        if( geodMeanAndSegments.size() != numOfTotalDofs )
            throw BasicException("GeodesicMean::checkIfAllSegmentsAreRelaxed(): argument has wrong size!");

        // get geodesic mean
        VectorType geodMean(geodMeanAndSegments.segment(0, numLocalDofs));
        VectorType singleEnergies;

        // get geodesic segments
        for( int k = 0; k < _numOfShapes; k++){
            VectorType segment( geodMeanAndSegments.segment(numLocalDofs + k * numGeodesicDofs, numGeodesicDofs) );
            VectorType inputShape( _shapes.segment(k * numLocalDofs, numLocalDofs) );
            DiscretePathEnergy<ConfiguratorType>          E( _W, _K, geodMean, inputShape );
            DiscretePathEnergyGradient<ConfiguratorType> DE( _W, _K, geodMean, inputShape );

            RealType energy;
            VectorType grad;
            E.apply( segment, energy );
            std::cerr << "Path energy of " << k << "th segment is " << energy << std::endl;
            E.evaluateSingleEnergies( segment, singleEnergies );
            std::cerr << "Intermediate deform. energies are " << std::endl;
            for( int i = 0; i < singleEnergies.size(); i++ )
                std::cerr << singleEnergies[i] << ", ";
            std::cerr << std::endl;
            DE.apply( segment, grad );
            std::cerr << "Path energy gradient norm of " << k << "th segment is " << grad.norm() << std::endl << std::endl;
        }

    }

    void checkIfAllSegmentsAreRelaxed( const std::vector<VectorType> &geodMeanAndSegments ) const {

        const int numLocalDofs    = 3 * _topology.getNumVertices();
        const int numGeodesicDofs = (_K-1) * numLocalDofs;

        if( geodMeanAndSegments.size() != 1 + _numOfShapes )
            throw BasicException("GeodesicMean::checkIfAllSegmentsAreRelaxed(): argument has wrong size!");
        if( geodMeanAndSegments[0].size() != numLocalDofs )
            throw BasicException("GeodesicMean::checkIfAllSegmentsAreRelaxed(): mean has wrong size!");

        RealType energy;
        VectorType singleEnergies, gradient;

        // get geodesic segments
        for( int k = 0; k < _numOfShapes; k++){
            if( geodMeanAndSegments[1+k].size() != numGeodesicDofs )
                throw BasicException("GeodesicMean::checkIfAllSegmentsAreRelaxed(): segment has wrong size!");
            VectorType inputShape( _shapes.segment(k * numLocalDofs, numLocalDofs) );
            DiscretePathEnergy<ConfiguratorType>          E( _W, _K, geodMeanAndSegments[0], inputShape );
            DiscretePathEnergyGradient<ConfiguratorType> DE( _W, _K, geodMeanAndSegments[0], inputShape );
            E.apply( geodMeanAndSegments[1+k], energy );
            std::cerr << "==================================================" << std::endl;
            std::cerr << "Path energy of " << k << "th segment is " << energy << std::endl;
            E.evaluateSingleEnergies( geodMeanAndSegments[1+k], singleEnergies );
            std::cerr << "Intermediate deform. energies are ";
            for( int i = 0; i < singleEnergies.size(); i++ )
                std::cerr << singleEnergies[i] << ", ";
            std::cerr << std::endl;
            DE.apply( geodMeanAndSegments[1+k], gradient );
            std::cerr << "Path energy gradient norm of " << k << "th segment is " << gradient.norm() << std::endl << std::endl;
        }
    }

};

#endif //GEODESICMEAN_H
