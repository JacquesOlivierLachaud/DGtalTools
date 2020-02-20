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
#include "DGtal/helpers/StdDefs.h"

namespace DGtal
{
  /////////////////////////////////////////////////////////////////////////////
  // template class SimplifiedMesh
  /**
     Description of template class 'SimplifiedMesh' <p> \brief Aim:
     Represents a mesh as faces and a list of vertices. Vertices may
     be shared among faces but no correct topology is required.

     @tparam TRealPoint an arbitrary model of RealPoint.
     @tparam TRealVector an arbitrary model of RealVector.
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

    /// Initializes a mesh from vertex positions and polygonal faces
    /// (clears everything before).
    template <typename RealPointIterator, typename FaceIterator>
    bool init( RealPointIterator itPos, RealPointIterator itPosEnd,
	       FaceIterator itFace, FaceIterator itFaceEnd );

    /// Clears everything. The object is empty.
    void clear();

    /// Given a range of real vectors, sets the normals of every
    /// vertex to the given vectors.
    template <typename RealVectorIterator>
    bool setVertexNormals( RealVectorIterator itN, RealVectorIterator itNEnd );

    /// Given a range of real vectors, sets the normals of every
    /// face to the given vectors.
    template <typename RealVectorIterator>
    bool setFaceNormals( RealVectorIterator itN, RealVectorIterator itNEnd );

    /// Uses the positions of vertices to compute a normal vector to
    /// each face of the mesh. It computes the barycenter,
    /// triangulates implicitly the face to build the normal vector.
    void computeFaceNormalsFromPositions();

    /// Uses the normals associated with vertices to compute a normal
    /// vector to each face of the mesh. It simply averages the
    /// normals at every incident vertex.
    void computeFaceNormalsFromVertexNormals();

    /// Uses the normals associated with faces to compute a normal
    /// vector to each vertex of the mesh. It simply averages the
    /// normals of every incident face.
    void computeVertexNormalsFromFaceNormals();
    
    /// @}

    /// @name Accessors
    /// @{

    /// @return a const reference to the vector giving for each face
    /// its incident vertices.
    const std::vector< FaceVertices >& incidentVertices() const
    { return myFaces; }

    /// @return a const reference to the vector giving for each vertex
    /// its incident faces.
    const std::vector< VertexFaces >& incidentFaces() const
    { return myVertices; }

    /// @return a const reference to the vector of positions (of vertices).
    const std::vector< RealVector >& positions() const
    { return myPositions; }

    /// @return a const reference to the vector of normals to vertices.
    const std::vector< RealVector >& vertexNormals() const
    { return myVertexNormals; }

    /// @return a const reference to the vector of normals to faces.
    const std::vector< RealVector >& faceNormals() const
    { return myFaceNormals; }

    /// @}

    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;
    
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

  /**
   * Overloads 'operator<<' for displaying objects of class 'SimplifiedMesh'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'SimplifiedMesh' to write.
   * @return the output stream after the writing.
   */
  template < typename TRealPoint, typename TRealVector >
  std::ostream&
  operator<< ( std::ostream & out,
	       const SimplifiedMesh<TRealPoint, TRealVector> & object );

  
  /////////////////////////////////////////////////////////////////////////////
  // template class SimplifiedMeshReader
  /**
     Description of template class 'SimplifiedMeshReader' <p> \brief Aim:
     An helper class for reading mesh files and creating a SimplifiedMesh.

     @tparam TRealPoint an arbitrary model of RealPoint.
     @tparam TRealVector an arbitrary model of RealVector.
   */
  template < typename TRealPoint, typename TRealVector >
  struct SimplifiedMeshReader
  {
    typedef TRealPoint                              RealPoint;
    typedef TRealVector                             RealVector;
    typedef SimplifiedMeshReader< RealPoint, RealVector > Self;
    static const Dimension dimension = RealPoint::dimension;
    BOOST_STATIC_ASSERT( ( dimension == 3 ) );

    typedef DGtal::SimplifiedMesh< RealPoint, RealVector > SimplifiedMesh;
    typedef typename SimplifiedMesh::Size           Size;
    typedef typename SimplifiedMesh::Index          Index;
    typedef typename SimplifiedMesh::FaceVertices   FaceVertices;
    typedef typename SimplifiedMesh::VertexFaces    VertexFaces;

    /// Checks that every index in \a indices are different from the others.
    /// @param indices a vector of integer indices
    /// @return 'true' iff the integer indices are all pairwise different.
    static bool verifyIndices( const std::vector< Index > indices );

    /// Splits a string \a str into several strings according to a
    /// delimiter \a delim.
    /// @param[in] str any string.
    /// @param[in] delim any delimiter character.
    /// @return the vector of split strings.
    static
    std::vector< std::string > split( const std::string& str, char delim = ' ');

    /// Reads an input file as an OBJ file format and outputs the
    /// corresponding simplified mesh.
    /// @param[inout] the input stream where the OBJ file is read.
    /// @param[out] the output simplified mesh.
    /// @return 'true' if both reading the input stream was ok and the created mesh is ok.
    static
    bool readOBJ( std::istream & input, SimplifiedMesh & smesh );
  };


  /////////////////////////////////////////////////////////////////////////////
  // template class SimplifiedMeshWriter
  /**
     Description of template class 'SimplifiedMeshWriter' <p> \brief Aim:
     An helper class for writing mesh file formats and creating a SimplifiedMesh.

     @tparam TRealPoint an arbitrary model of RealPoint.
     @tparam TRealVector an arbitrary model of RealVector.
   */
  template < typename TRealPoint, typename TRealVector >
  struct SimplifiedMeshWriter
  {
    typedef TRealPoint                              RealPoint;
    typedef TRealVector                             RealVector;
    typedef SimplifiedMeshWriter< RealPoint, RealVector > Self;
    static const Dimension dimension = RealPoint::dimension;
    BOOST_STATIC_ASSERT( ( dimension == 3 ) );

    typedef DGtal::SimplifiedMesh< RealPoint, RealVector > SimplifiedMesh;
    typedef typename SimplifiedMesh::Size           Size;
    typedef typename SimplifiedMesh::Index          Index;
    typedef typename SimplifiedMesh::FaceVertices   FaceVertices;
    typedef typename SimplifiedMesh::VertexFaces    VertexFaces;

    /// Writes a simplified mesh in an output file (in OBJ file format).
    /// @param[inout] the output stream where the OBJ file is written.
    /// @param[in] the simplified mesh.
    /// @return 'true' if writing in the output stream was ok.
    static
    bool writeOBJ( std::ostream & output, const SimplifiedMesh & smesh );
    
  };

} // namespace DGtal

///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "SimplifiedMesh.ih"
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined SimplifiedMesh_h

#undef SimplifiedMesh_RECURSES
#endif // else defined(SimplifiedMesh_RECURSES)

