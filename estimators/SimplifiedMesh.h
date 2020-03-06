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
#include "DGtal/io/Color.h"

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
    typedef typename RealVector::Component          Scalar;
    typedef std::vector<Scalar>                     Scalars;
    /// The type for counting elements.
    typedef std::size_t                             Size;
    /// The type used for numbering vertices and faces
    typedef std::size_t                             Index;
    typedef Index                                   Face;
    typedef Index                                   Edge;
    typedef Index                                   Vertex;
    typedef std::pair< Face, Scalar >               WeightedFace;
    typedef std::pair< Vertex, Scalar >             WeightedVertex;
    /// The type that defines a range of vertices
    typedef std::vector< Vertex >                   Vertices;
    /// The type that defines a range of faces
    typedef std::vector< Face >                     Faces;
    typedef std::vector< WeightedFace >             WeightedFaces;
    typedef std::pair< Vertex, Vertex >             VertexPair;
    //---------------------------------------------------------------------------
  public:
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
    
    /// Uses the normals associated with faces to compute a normal
    /// vector to each vertex of the mesh. It uses the weights
    /// proposed by [Max, 1999] for combining face information into
    /// vertex information.
    void computeVertexNormalsFromFaceNormalsWithMaxWeights();

    template <typename AnyRing>
    std::vector<AnyRing> computeFaceValuesFromVertexValues
    ( const std::vector<AnyRing>& vvalues ) const;
    
    template <typename AnyRing>
    std::vector<AnyRing> computeVertexValuesFromFaceValues
    ( const std::vector<AnyRing>& fvalues ) const;
    
    /// @}

    //---------------------------------------------------------------------------
  public:
    /// @name Accessors
    /// @{

    /// @return the number of faces of the mesh.
    Size nbFaces() const
    { return myIncidentVertices.size(); }

    /// @return the number of (unordered) edges of the mesh.
    Size nbEdges() const
    { return myEdgeVertices.size(); }

    /// @return the number of vertices of the mesh.
    Size nbVertices() const
    { return myIncidentFaces.size(); }

    /// @param i any vertex of the mesh
    /// @param j any vertex of the mesh
    /// @return the edge index of edge (i,j) or `nbEdges()` if this
    /// edge does not exist.
    /// @note O(log E) time complexity.
    Edge makeEdge( Vertex i, Vertex j ) const;
    
    /// @return a const reference to the vector giving for each face
    /// its incident vertices.
    const std::vector< Vertices >& incidentVertices() const
    { return myIncidentVertices; }

    /// @return a const reference to the vector giving for each vertex
    /// its incident faces.
    const std::vector< Faces >& incidentFaces() const
    { return myIncidentFaces; }

    /// @return a const reference to the vector of positions (of vertices).
    const std::vector< RealVector >& positions() const
    { return myPositions; }

    /// @return a const reference to the vector of normals to vertices.
    const std::vector< RealVector >& vertexNormals() const
    { return myVertexNormals; }

    /// @return a const reference to the vector of normals to faces.
    const std::vector< RealVector >& faceNormals() const
    { return myFaceNormals; }

    /// @return a const reference to the vector of neighbor faces for each face.
    const std::vector< Faces >& neighborFaces() const
    { return myNeighborFaces; }

    /// @return a const reference to the vector of neighbor vertices for each vertex.
    const std::vector< Vertices >& neighborVertices() const
    { return myNeighborVertices; }

    /// @return a const reference to the vector giving for each edge
    /// its two vertices (as a pair (i,j), i<j).
    /// @note edges are sorted in increasing order.
    const std::vector< VertexPair >& edgeVertices() const
    { return myEdgeVertices; }
    
    /// @return a const reference to the vector giving for each edge
    /// its incident faces (one, two, or more if non manifold)
    const std::vector< Faces >& edgeFaces() const
    { return myEdgeFaces; }

    
    
    /// @}

    //---------------------------------------------------------------------------
  public:
    /// @name Geometric services
    /// @{

    /// @return the average of the length of edges.
    Scalar averageEdgeLength() const;

    /// Perturbate the positions with a uniform random noise of 'p *
    /// averageEdgeLength' along arbitrary directions.
    /// @param p any positive real value.
    void perturbateWithUniformRandomNoise( Scalar p );
    void perturbateWithAdaptiveUniformRandomNoise( Scalar p );

    /// @param f any valid face index.
    /// @return the centroid (or barycenter) of face \a f.
    RealPoint faceCentroid( Index f ) const;

    /// @param f any valid face index.
    /// @return the area of face \a f.
    Scalar faceArea( Index f ) const;

    /// @param v any valid vertex index.
    /// @return the Max's weights for each incident face to \a v, in the same order as `myIncidentFaces[ v ]`.
    /// @note Used in computeVertexNormalsFromFaceNormalsWithMaxWeights
    Scalars getMaxWeights( Index v ) const;
    
    WeightedFaces
    computeFacesInclusionsInBall( Scalar r, Index f ) const;
    
    /// Computes an approximation of the inclusion ratio of a given
    /// face \a f with a ball of radius \a r and center \a p.
    ///
    /// @param p the center of the ball.
    /// @param r the radius of the ball.
    /// @param f any index of face.
    ///
    /// @return the inclusion ratio as a scalar between 0 (no
    /// intersection) and 1 (inclusion).
    Scalar inclusionRatio( RealPoint p, Scalar r, Index f ) const;

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
    /// For each face, its range of incident vertices
    std::vector< Vertices >     myIncidentVertices;
    /// For each vertex, its range of incident faces
    std::vector< Faces >        myIncidentFaces;
    /// For each vertex, its position
    std::vector< RealPoint >    myPositions;
    /// For each vertex, its normal vector
    std::vector< RealVector >   myVertexNormals;
    /// For each face, its normal vector
    std::vector< RealVector >   myFaceNormals;
    /// For each face, its range of neighbor faces (no particular order)
    std::vector< Faces >        myNeighborFaces;
    /// For each vertex, its range of neighbor vertices (no particular order)
    std::vector< Vertices >     myNeighborVertices;
    /// For each edge, its two vertices
    std::vector< VertexPair >   myEdgeVertices;
    /// For each edge, its faces (one, two, or more if non manifold)
    std::vector< Faces >        myEdgeFaces;
    
    // ------------------------- Private Datas --------------------------------
  private:


    // ------------------------- Internals ------------------------------------
  protected:

    /// Computes neighboring information.
    void computeNeighbors();
    /// Computes edge information.
    void computeEdges();

    /// @return a random number between 0.0 and 1.0
    static Scalar rand01()
    { return (Scalar) rand() / (Scalar) RAND_MAX; }
    
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
    typedef typename SimplifiedMesh::Vertices       Vertices;
    typedef typename SimplifiedMesh::Faces          Faces;

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
    typedef typename SimplifiedMesh::Vertices       Vertices;
    typedef typename SimplifiedMesh::Faces          Faces;
    typedef std::vector< Color >                    Colors;

    /// Writes a simplified mesh in an output file (in OBJ file format).
    /// @param[inout] the output stream where the OBJ file is written.
    /// @param[in] the simplified mesh.
    /// @return 'true' if writing in the output stream was ok.
    static
    bool writeOBJ( std::ostream & output, const SimplifiedMesh & smesh );

    static
    bool writeOBJ( std::string            objfile,
                   const SimplifiedMesh & smesh, 
                   const Colors&          diffuse_colors = Colors(),
                   const Color&           ambient_color  = Color( 32, 32, 32 ),
                   const Color&           diffuse_color  = Color( 200, 200, 255 ),
                   const Color&           specular_color = Color::White );

    template <typename EdgePredicate>
    static
    bool writeEdgeLinesOBJ( std::string            objfile,
			    const SimplifiedMesh & smesh,
			    const EdgePredicate &  edge_predicate,
			    const double           relative_thickness = 0.05,
			    const Color&           ambient_color = Color::Black,
			    const Color&           diffuse_color = Color::Black,
			    const Color&           specular_color= Color::Black );
    
  };


  /////////////////////////////////////////////////////////////////////////////
  // template class SimplifiedMeshHelper
  /**
     Description of template class 'SimplifiedMeshHelper' <p> \brief Aim:
     An helper class for building classical meshes.

     @tparam TRealPoint an arbitrary model of RealPoint.
     @tparam TRealVector an arbitrary model of RealVector.
   */
  template < typename TRealPoint, typename TRealVector >
  struct SimplifiedMeshHelper
  {
    typedef TRealPoint                              RealPoint;
    typedef TRealVector                             RealVector;
    typedef SimplifiedMeshHelper< RealPoint, RealVector > Self;
    static const Dimension dimension = RealPoint::dimension;
    BOOST_STATIC_ASSERT( ( dimension == 3 ) );

    typedef DGtal::SimplifiedMesh< RealPoint, RealVector > SimplifiedMesh;
    typedef typename RealVector::Component          Scalar;
    typedef std::vector<Scalar>                     Scalars;
    typedef typename SimplifiedMesh::Size           Size;
    typedef typename SimplifiedMesh::Index          Index;
    typedef typename SimplifiedMesh::Vertices       Vertices;
    typedef typename SimplifiedMesh::Faces          Faces;

    enum class Normals { NO_NORMALS, VERTEX_NORMALS, FACE_NORMALS };
    static
    SimplifiedMesh makeSphere( Scalar radius, RealPoint center, Size m, Size n, Normals normals );
    
    static
    Scalars sphereMeanCurvatures( Scalar radius, Size m, Size n );
    
    static
    Scalars sphereGaussianCurvatures( Scalar radius, Size m, Size n );
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

