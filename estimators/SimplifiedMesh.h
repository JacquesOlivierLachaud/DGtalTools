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
 * @file SimplifiedMesh.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2020/02/18
 *
 * Header file for module SimplifiedMesh.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(SimplifiedMesh_RECURSES)
#error Recursive header files inclusion detected in SimplifiedMesh.h
#else // defined(SimplifiedMesh_RECURSES)
/** Prevents recursive inclusion of headers. */
#define SimplifiedMesh_RECURSES

#if !defined SimplifiedMesh_h
/** Prevents repeated inclusion of headers. */
#define SimplifiedMesh_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <sstream>
#include <string>
// always include EigenSupport.h before any other Eigen headers
// #include "DGtal/math/linalg/EigenSupport.h"
#include "DGtal/base/Common.h"
#include "DGtal/base/Clone.h"

namespace DGtal
{
  /////////////////////////////////////////////////////////////////////////////
  // template class SimplifiedMesh
  /**
     Description of template class 'SimplifiedMesh' <p> \brief Aim:
     Represents a mesh as faces and a list of vertices. Vertices may
     be shared among faces but no correct topology is required.

     @tparam TSpace an arbitrary model of CSpace.
   */
  
  template < typename TRealPoint, typename TRealVector >
  struct SimplifiedMesh
  {
    typedef TRealPoint                              RealPoint;
    typedef TRealVector                             RealVector;
    typedef SimplifiedMesh< RealPoint, RealVector > Self;
    static const Dimension dimension = RealPoint::dimension;
    BOOST_STATIC_ASSERT( ( dimension == 3 ) );

    /// The type for counting elements.
    typedef std::size_t                             Size;
    /// The type used for numbering vertices
    typedef std::size_t                             Index;
    /// The type that defines the vertices of a face.
    typedef std::vector< Index >                    FaceVertices;
    /// The type that defines the faces of a vertex.
    typedef std::vector< Index >                    VertexFaces;

    /// @name Standard services
    /// @{
    
    /// Default destructor.
    ~SimplifiedMesh() = default;
    /// Default constructor.
    SimplifiedMesh() = default;
    /// Default copy constructor.
    /// @param other the object to clone
    SimplifiedMesh( const Self& other ) = default;
    /// Default move constructor.
    /// @param other the object to move
    SimplifiedMesh( Self&& other ) = default;
    /// Default assignment constructor.
    /// @param other the object to clone
    /// @return a reference to 'this'.
    Self& operator=( const Self& other ) = default;

    /// Builds a mesh from vertex positions and polygonal faces.
    template <typename RealPointIterator, typename FaceIterator>
    SimplifiedMesh( RealPointIterator itPos, RealPointIterator itPosEnd,
                    FaceIterator itFace, FaceIterator itFaceEnd );
    
    /// Clears everything. The object is empty.
    void clear();

    /// @}

    
    // ------------------------- Protected Datas ------------------------------
  protected:

    std::vector< FaceVertices > myFaces;
    std::vector< VertexFaces >  myVertices;
    std::vector< RealPoint >    myPositions;
    std::vector< RealVector >   myVertexNormals;
    std::vector< RealVector >   myFaceNormals;

    
    // ------------------------- Private Datas --------------------------------
  private:


    // ------------------------- Internals ------------------------------------
  private:
    
    
  }; // end of class SimplifiedMesh
  
} // namespace DGtal

///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "SimplifiedMesh.ih"
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined SimplifiedMesh_h

#undef SimplifiedMesh_RECURSES
#endif // else defined(SimplifiedMesh_RECURSES)

