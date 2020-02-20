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

     Formula for interpolated measures:
     
     MU0=-1/6*((uAz + uBz + uCz)*Bx - (uAz + uBz + uCz)*Cx)*Ay + 1/6*((uAz + uBz + uCz)*Ax - (uAz + uBz + uCz)*Cx)*By - 1/6*((uAz + uBz + uCz)*Ax - (uAz + uBz + uCz)*Bx)*Cy + 1/6*((uAy + uBy + uCy)*Bx - (uAy + uBy + uCy)*Cx - (uAx + uBx + uCx)*By + (uAx + uBx + uCx)*Cy)*Az - 1/6*((uAy + uBy + uCy)*Ax - (uAy + uBy + uCy)*Cx - (uAx + uBx + uCx)*Ay + (uAx + uBx + uCx)*Cy)*Bz + 1/6*((uAy + uBy + uCy)*Ax - (uAy + uBy + uCy)*Bx - (uAx + uBx + uCx)*Ay + (uAx + uBx + uCx)*By)*Cz
     Let UM=uA+uB+uC.
     MU0=-1/6*(uMz*Bx - uMz*Cx)*Ay + 1/6*(uMz*Ax - uMz*Cx)*By - 1/6*(uMz*Ax - uMz*Bx)*Cy + 1/6*(uMy*Bx - uMy*Cx - uMx*By + uMx*Cy)*Az - 1/6*(uMy*Ax - uMy*Cx - uMx*Ay + uMx*Cy)*Bz + 1/6*(uMy*Ax - uMy*Bx - uMx*Ay + uMx*By)*Cz
     We see by simple computations that MU0 can be written as (uM = UM/3)
     MU0=1/2*det( uM, B-A, C-A )

     MU1=1/6*((uBy - uCy)*uAz - (uAy + 2*uCy)*uBz + (uAy + 2*uBy)*uCz)*Ax + 1/6*((uBy + 2*uCy)*uAz - (uAy - uCy)*uBz - (2*uAy + uBy)*uCz)*Bx - 1/6*((2*uBy + uCy)*uAz - (2*uAy + uCy)*uBz - (uAy - uBy)*uCz)*Cx - 1/6*((uBx - uCx)*uAz - (uAx + 2*uCx)*uBz + (uAx + 2*uBx)*uCz)*Ay - 1/6*((uBx + 2*uCx)*uAz - (uAx - uCx)*uBz - (2*uAx + uBx)*uCz)*By + 1/6*((2*uBx + uCx)*uAz - (2*uAx + uCx)*uBz - (uAx - uBx)*uCz)*Cy + 1/6*((uBx - uCx)*uAy - (uAx + 2*uCx)*uBy + (uAx + 2*uBx)*uCy)*Az + 1/6*((uBx + 2*uCx)*uAy - (uAx - uCx)*uBy - (2*uAx + uBx)*uCy)*Bz - 1/6*((2*uBx + uCx)*uAy - (2*uAx + uCx)*uBy - (uAx - uBx)*uCy)*Cz

     This formula can also be written in a clearer form
     6*MU1 = | u_A+u_B+u_C u_C-u_B A | + | u_A+u_B+u_C u_A-u_C B | + | u_A+u_B+u_C u_B-u_A C |
     It follows that 
     MU1=1/2( | uM u_C-u_B A | + | uM u_A-u_C B | + | uM u_B-u_A C |

     Last, Gaussian curvature measure is
     MU2=-1/2*uCx*uBy*uAz + 1/2*uBx*uCy*uAz + 1/2*uCx*uAy*uBz - 1/2*uAx*uCy*uBz - 1/2*uBx*uAy*uCz + 1/2*uAx*uBy*uCz

     which is simply
     MU2=1/2*det( uA, uB, uC )

  */
  template < typename TRealPoint, typename TRealVector >
  class CorrectedNormalCurrentFormula
  {
    typedef TRealPoint                     RealPoint;
    typedef TRealVector                    RealVector;
    typedef typename RealVector::Component Scalar;
    typedef std::vector< RealPoint >       RealPoints;
    typedef std::vector< RealVector >      RealVectors;
    typedef std::size_t                    Index;

    /// Computes mu0 measure (area) of triangle abc given a constant
    /// corrected normal vector \a u.
    /// @param a any point
    /// @param b any point
    /// @param c any point
    /// @param u the constant corrected normal vector to triangle abc
    /// @return the mu0-measure of triangle abc, i.e. its area.
    static
    Scalar mu0ConstantU
    ( const RealPoint& a, const RealPoint& b, const RealPoint& c,
      const RealVector& u )
    {
      return 0.5 * ( b - a ).crossProduct( c - a ).dotProduct( u );
    }

    /// Computes mu0 measure (area) of triangle abc given an interpolated
    /// corrected normal vector \a ua, \a \ub, \a uc.
    /// @param a any point
    /// @param b any point
    /// @param c any point
    /// @param ua the corrected normal vector at point a
    /// @param ub the corrected normal vector at point b
    /// @param uc the corrected normal vector at point c
    /// @param unit_uM when 'true' forces the average normal vector
    /// `(ua+ub+uc)/3` to be unitary, otherwise leave it as is.
    /// @return the mu0-measure of triangle abc, i.e. its area.
    static
    Scalar mu0InterpolatedU
    ( const RealPoint& a, const RealPoint& b, const RealPoint& c,
      const RealVector& ua, const RealVector& ub, const RealVector& uc,
      bool unit_uM = false )
    {
      // MU0=1/2*det( uM, B-A, C-A )
      //    =  1/2 < ( (u_A + u_B + u_C)/3.0 ) | (AB x AC ) >
      RealVector uM = ( ua+ub+uc ) / 3.0;
      if ( uM ) uM /= uM.norm();
      return O.5 * ( b - a ).crossProduct( c - a ).dotProduct( uM );
    }
    
    /// Computes mu0 measure (area) of polygonal face \a pts given a
    /// constant corrected normal vector \a u.
    /// @param pts the (ccw ordered) points forming the vertices of a polygonal face.
    /// @param u the constant corrected normal vector to this polygonal face.
    /// @return the mu0-measure of the given polygonal face, i.e. its area.
    static
    Scalar mu0ConstantU( const RealPoints& pts, const RealVector& u )
    {
      if ( pts.size() <  3 ) return 0.0;
      if ( pts.size() == 3 )
	return mu0ConstantU( pts[ 0 ], pts[ 1 ], pts[ 2 ], u );
      const RealPoint b = barycenter( pts );
      Scalar          a = 0.0;
      for ( Index i = 0; i < pts.size(); i++ )
	a += mu0ConstantU( b, pts[ i ], pts[ (i+1)%pts.size() ], u );
      return a;
    }

    /// Computes area of polygonal face \a pts given an interpolated
    /// corrected normal vector \a ua, \a \ub, \a uc.
    /// @param pts the (ccw ordered) points forming the vertices of a polygonal face.
    /// @param u the (ccw ordered) normal vectors at the corresponding vertices in \a pts.
    /// @param unit_uM when 'true' forces the average normal vectors
    /// to be unitary, otherwise leave it as is.
    /// @return the mu0-measure of the given polygonal face, i.e. its area.
    static
    Scalar mu0InterpolatedU( const RealPoints& pts, const RealVectors& u,
			     bool unit_uM = false )
    {
      ASSERT( pts.size() == u.size() );
      if ( pts.size() <  3 ) return 0.0;
      if ( pts.size() == 3 )
	return mu0InterpolatedU( pts[ 0 ], pts[ 1 ], pts[ 2 ],
				 u[ 0 ], u[ 1 ], u[ 2 ], unit_uM );
      const RealPoint   b = barycenter( pts );
      const RealVector ub = barycenter( pts );
      Scalar          a = 0.0;
      for ( Index i = 0; i < pts.size(); i++ )
	a += mu0InterpolatedU( b,  pts[ i ], pts[ (i+1)%pts.size() ],
			       ub,   u[ i ],   u[ (i+1)%pts.size() ], unit_uM );
      return a;
    }
    
    /// Given a vector of points, returns its barycenter.
    /// @param pts any vector of points
    /// @return the barycenter of these points.
    static 
    RealPoint barycenter( const RealPoints& pts )
    {
      RealPoint b;
      for ( auto p : pts ) b += p;
      b /= pts.size();
      return b;
    }

    /// Given a vector of unit vectors, returns their average unit vector.
    /// @param pts any vector of vectors.
    /// @return the average unit vector.
    static 
    RealVector averageUnitVector( const RealVectors& vecs )
    {
      RealVector avg;
      for ( auto v : vecs ) avg += v;
      return avg.getNormalized();
    }
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
  
