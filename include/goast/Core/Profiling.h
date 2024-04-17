// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2024 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include <map>
#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>

template<typename Derived>
class TimedClass {
public:
  inline static std::map<std::string, double> m_Runtimes;
  inline static std::map<std::string, int> m_numRuns;

  struct ScopeTimer {
    const std::string m_Name;
    std::chrono::high_resolution_clock::time_point m_StartTime;
  public:
    explicit ScopeTimer( std::string Name ) : m_Name( std::move( Name )),
                                              m_StartTime( std::chrono::high_resolution_clock::now()) {}

    ~ScopeTimer() {
      auto EndTime = std::chrono::high_resolution_clock::now();
      m_Runtimes[m_Name] += std::chrono::duration<double, std::milli>( EndTime - m_StartTime ).count();
      m_numRuns[m_Name]++;
    }
  };


public:

  static void resetTimers() {
    m_Runtimes.clear();
    m_numRuns.clear();
  }

  virtual void printTimings() const {
    std::cout << " | Timings for " << typeid( *this ).name() << std::endl;
    for ( auto const &[Name, _]: m_Runtimes )
      std::cout << " | "
                << std::left << std::setw( 45 ) << Name << " | "
                << std::right << std::setw( 6 ) << m_numRuns[Name] << " | "
                << std::setw( 8 ) << std::fixed << std::setprecision( 2 ) << m_Runtimes[Name] / 1000. << "s | "
                << std::setw( 8 ) << std::fixed << std::setprecision( 2 ) << m_Runtimes[Name] / m_numRuns[Name] << "ms"
                << " | "
                << std::endl;
  }
};