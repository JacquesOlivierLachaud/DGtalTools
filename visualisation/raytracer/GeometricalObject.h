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
 * @file GeometricalObject.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2017/02/28
 *
 * This file is part of the DGtal library.
 */

#if defined(GeometricalObject_RECURSES)
#error Recursive header files inclusion detected in GeometricalObject.h
#else // defined(GeometricalObject_RECURSES)
/** Prevents recursive inclusion of headers. */
#define GeometricalObject_RECURSES

#if !defined GeometricalObject_h
/** Prevents repeated inclusion of headers. */
#define GeometricalObject_h

#include "raytracer/Ray.h"
#include "raytracer/RayIntersection.h"

namespace DGtal {

  namespace rt {

    /// This is an interface specifying methods that any ray-traced
    /// geometrical object should have. 
    struct GeometricalObject {

      /// Default constructor. Nothing to do.
      GeometricalObject() {}

      /// Virtual destructor since object contains virtual methods.
      virtual ~GeometricalObject() {}

      /// @return the normal vector at point \a p on the object (\a p
      /// should be on or close to the sphere).
      virtual Vector3 getNormal( Point3 p ) = 0;

      /// @param[in,out] ray_inter as input the incoming ray, as
      /// output information abour intersection.
      ///
      /// @return true if there was an intersection, false otherwise
      /// (more information is stored in ray_inter)
      virtual bool intersectRay( const Ray& ray, RayIntersection& ray_inter ) = 0;
                    
    };

  } // namespace rt
} // namespace DGtal

#endif // !defined GeometricalObject_h

#undef GeometricalObject_RECURSES
#endif // else defined(GeometricalObject_RECURSES)

