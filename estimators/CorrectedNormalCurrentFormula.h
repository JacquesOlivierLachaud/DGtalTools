/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#pragma once

/**
 * @file CorrectedNormalCurrentFormula.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2017/06/19
 *
 * Header file for module CorrectedNormalCurrentFormula.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(CorrectedNormalCurrentFormula_RECURSES)
#error Recursive header files inclusion detected in CorrectedNormalCurrentFormula.h
#else // defined(CorrectedNormalCurrentFormula_RECURSES)
/** Prevents recursive inclusion of headers. */
#define CorrectedNormalCurrentFormula_RECURSES

#if !defined CorrectedNormalCurrentFormula_h
/** Prevents repeated inclusion of headers. */
#define CorrectedNormalCurrentFormula_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <map>
#include "DGtal/base/Common.h"
#include "DGtal/topology/CDigitalSurfaceContainer.h"
#include "DGtal/topology/IndexedDigitalSurface.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "DGtal/math/linalg/SimpleMatrix.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include "SphericalTriangle.h"

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // class CorrectedNormalCurrentFormula
  /**
     Description of class 'CorrectedNormalCurrentFormula' <p> \brief
     Aim: A helper class that provides static methods to compute
     corrected normal current formulas of curvatures.

     @tparam TDigitalSurfaceContainer any type of digital surface container.
  */
  template < typename TRealPoint, typename TRealVector >
  class CorrectedNormalCurrentFormula
  {
    typedef TRealPoint  RealPoint;
    typedef TRealVector RealVector;
    
  };

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
//#include "CorrectedNormalCurrentFormula.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined CorrectedNormalCurrentFormula_h

#undef CorrectedNormalCurrentFormula_RECURSES
#endif // else defined(CorrectedNormalCurrentFormula_RECURSES)
  
