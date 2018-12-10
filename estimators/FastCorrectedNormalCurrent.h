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
 * @file FastCorrectedNormalCurrent.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2017/06/19
 *
 * Header file for module FastCorrectedNormalCurrent.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(FastCorrectedNormalCurrent_RECURSES)
#error Recursive header files inclusion detected in FastCorrectedNormalCurrent.h
#else // defined(FastCorrectedNormalCurrent_RECURSES)
/** Prevents recursive inclusion of headers. */
#define FastCorrectedNormalCurrent_RECURSES

#if !defined FastCorrectedNormalCurrent_h
/** Prevents repeated inclusion of headers. */
#define FastCorrectedNormalCurrent_h

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
  // class FastCorrectedNormalCurrent
  /**
     Description of class 'FastCorrectedNormalCurrent' <p> \brief Aim:
     Represent a current over a digital surface, whose tangent plane
     is described by an additional normal vector field. It is useful
     to define geometric measures onto digital surface, which provide
     area, mean and gaussian curvature information.

     @note it is tagged "fast" since it used the new
     IndexedDigitalSurface to represent the digital surface.
     
     @tparam TDigitalSurfaceContainer any type of digital surface container.
  */
  template <typename TDigitalSurfaceContainer>
  class FastCorrectedNormalCurrent
  {
  public:
    typedef TDigitalSurfaceContainer                              DigitalSurfaceContainer;
    typedef FastCorrectedNormalCurrent< DigitalSurfaceContainer > Self;
    BOOST_CONCEPT_ASSERT(( concepts::CDigitalSurfaceContainer<DigitalSurfaceContainer> ));

    typedef IndexedDigitalSurface< DigitalSurfaceContainer > Surface;
    typedef typename Surface::KSpace                  KSpace;
    typedef typename Surface::Cell                    Cell;
    typedef typename Surface::SCell                   SCell;
    typedef typename Surface::Surfel                  Surfel;
    typedef typename Surface::Vertex                  Vertex; ///< an index
    typedef typename Surface::Arc                     Arc;    ///< an index
    typedef typename Surface::Face                    Face;   ///< an index
    typedef typename Surface::ConstIterator           ConstIterator;
    typedef typename Surface::VertexRange             VertexRange;
    typedef typename Surface::ArcRange                ArcRange;
    typedef typename Surface::FaceRange               FaceRange;
    typedef typename KSpace::Space                    Space;
    typedef typename Space::Point                     Point;
    typedef typename Space::Vector                    Vector;
    typedef typename Space::RealPoint                 RealPoint;
    typedef typename Space::RealVector                RealVector;
    typedef typename RealVector::Component            Scalar;
    typedef std::vector< RealVector >                 NormalVectorField;
    typedef std::vector< RealPoint >                  CentroidMap;
    typedef std::vector< Scalar >                     MeasureMap;
    typedef SimpleMatrix< Scalar, 3, 3 >              RealTensor;
    typedef std::vector< RealTensor >                 TensorMeasureMap;
    
    // Checks that dimension is 3.
    BOOST_STATIC_ASSERT(( KSpace::dimension == 3 ));
    
    // ----------------------- Standard services ------------------------------
  public:
  
    /**
     * Destructor.
     */
    ~FastCorrectedNormalCurrent() {}

    /// Default constructor. The object is invalid.
    FastCorrectedNormalCurrent()
      : theSurface( nullptr ) {}

    /**
       Constructor from surface. The surface is shared by the
       current. Computes also its natural tangent bundle
       (i.e. using the plane defined by each surfel).
 
       @param surface the digital surface which defines the current.

       @param h the digitization gridstep of the surface.

       @param crisp when 'false', measures on a ball are computed with
       an estimation of the relative intersection with each cell (more
       precise slower), otherwise the intersection is approximated
       either by 0 or 1 (less accurate, 30% faster).
    */
    FastCorrectedNormalCurrent( Alias<Surface> surface,
				Scalar h = 1.0, bool crisp = false )
      : theSurface( surface )
    {
      setParams( h, crisp );
      trace.info() << "computeTrivialNormals" << std::endl;
      computeTrivialNormals();
      trace.info() << "setCorrectedNormals" << std::endl;
      setCorrectedNormals( myTrivialNormals );
      trace.info() << "computeCellCentroids" << std::endl;
      computeCellCentroids();
      trace.info() << "...done" << std::endl;
    }
  
    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    FastCorrectedNormalCurrent ( const FastCorrectedNormalCurrent & other ) = default;
  
    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    FastCorrectedNormalCurrent & operator= ( const FastCorrectedNormalCurrent & other ) = default;

    /// Sets the parameters.
    ///
    /// @param h the digitization gridstep of the surface.
    ///
    /// @param crisp when 'false', measures on a ball are computed with
    /// an estimation of the relative intersection with each cell (more
    /// precise slower), otherwise the intersection is approximated
    /// either by 0 or 1 (less accurate, 30% faster).
    void setParams( Scalar h, bool crisp = false )
    {
      myH     = h;
      myCrisp = crisp;
      myEmbedder.init( myH );
    }
    
    /// Attaches the given surface to this current. Computes also its
    /// natural tangent bundle.
    /// @param surface the digital surface which defines the current.
    void attach( Alias<Surface> surface )
    {
      theSurface = surface;
      myTrivialNormals.clear();
      computeTrivialNormals();
      setCorrectedNormals( myTrivialNormals );
      computeCellCentroids();
    }
    
    /// @return the Khalimsky space associated with this current.
    const KSpace& space() const
    {
      ASSERT( theSurface != 0 );
      return theSurface->space();
    }
    
    /// @return the surface associated with this current.
    Surface& surface()
    {
      ASSERT( theSurface != 0 );
      return *theSurface;
    }

    /// @return an iterator on the first vertex of the surface.
    ConstIterator begin() const
    {
      ASSERT( theSurface != 0 );
      return theSurface->begin();
    }

    /**
       @return a ConstIterator after the last vertex of the surface.
    */
    ConstIterator end() const
    {
      ASSERT( theSurface != 0 );
      return theSurface->end();
    }

    // ----------------------- Normal services --------------------------------------
  public:

    /// Computes the trivial normals of the attached surface.
    void computeTrivialNormals()
    {
      const KSpace& K = space();
      myTrivialNormals.resize( theSurface->nbVertices() );
      for ( auto s : *this )
	{
	  SCell aSurfel = theSurface->surfel( s );
	  Dimension  k  = K.sOrthDir( aSurfel );
	  bool  direct  = K.sDirect( aSurfel, k );
	  RealVector t  = RealVector::zero;
	  t[ k ]        = direct ? -1.0 : 1.0;
	  myTrivialNormals[ s ] = t ;
	}
    }

    /// Computes the centroids of the cells of the attached surface.
    void computeCellCentroids()
    {
      const KSpace&   K = space();
      myVertexCentroids.resize( theSurface->nbVertices() );
      myArcCentroids   .resize( theSurface->nbArcs()     );
      myFaceCentroids  .resize( theSurface->nbFaces()    );
      for ( unsigned int i = 0; i < myVertexCentroids.size(); ++i )
	{
	  const SCell    aSurfel = theSurface->surfel( i );
	  myVertexCentroids[ i ] = computeCentroid( space().unsigns( aSurfel ) );
	}
      for ( unsigned int i = 0; i < myArcCentroids.size(); ++i )
	{
	  const SCell  aLinel = theSurface->linel( i );
	  myArcCentroids[ i ] = computeCentroid( space().unsigns( aLinel ) );
	}
      for ( unsigned int i = 0; i < myFaceCentroids.size(); ++i )
	{
	  const SCell aPointel = theSurface->pointel( i );
	  myFaceCentroids[ i ] = computeCentroid( space().unsigns( aPointel ) );
	}
    }
    
    /// Sets the corrected normal vector field of the current.
    /// @param nvf a map Surfel -> RealVector
    void setCorrectedNormals( const NormalVectorField& nvf )
    {
      myCorrectedNormals = nvf;
    }

    /// Sets the corrected normal vector field of the current.
    template <typename SurfelIterator, typename RealVectorIterator>
    void setCorrectedNormals( SurfelIterator itS, SurfelIterator itSEnd,
			      RealVectorIterator itRV )
    {
      myCorrectedNormals.resize( theSurface->nbVertices() );
      for ( ; itS != itSEnd; ++itS )
	{
	  const Vertex v = theSurface->getVertex( *itS );
	  if ( v != theSurface->INVALID_FACE )
	    myCorrectedNormals[ v ] = *itRV++;
	  else
	    trace.warning() << "[FastCorrectedNormalCurrent::setCorrectedNormals]"
			    << " Surfel " << *itS << " is not in the surface." << std::endl;
	}
    }

    
    /// Computes all measures mu0 per vertex.
    void computeAllMu0()
    {
      myMu0.resize( theSurface->nbVertices() );
      auto vtcs = theSurface->allVertices();
      for ( auto v : vtcs ) myMu0[ v ] = computeMu0( v );
    }
    /// Computes all measures mu1 per arc.
    void computeAllMu1()
    {
      myMu1.resize( theSurface->nbArcs() );
      auto arcs = theSurface->allArcs();
      for ( auto a : arcs ) myMu1[ a ] = computeMu1( a );
    }
    /// Computes all measures mu2 per face.
    void computeAllMu2()
    {
      myMu2.resize( theSurface->nbFaces() );
      auto faces = theSurface->allFaces();
      for ( auto f : faces ) myMu2[ f ] = computeMu2( f );
    }
    /// Computes all measures muOmega per face.
    void computeAllMuOmega()
    {
      myMuOmega.resize( theSurface->nbArcs() );
      auto arcs = theSurface->allArcs();
      for ( auto a : arcs ) myMuOmega[ a ] = computeMuOmega( a );
    }

    /// Computes all anisotropic measures per arc.
    void computeAllAnisotropicMu()
    {
      myAnisotropicMu.resize( theSurface->nbArcs() );
      auto arcs = theSurface->allArcs();
      for ( auto a : arcs ) myAnisotropicMu[ a ] = computeAnisotropicMu( a );
    }


    // ----------------------- Indexed Digital Surface services --------------------------
  public:

    /// @param[in] v any vertex index.
    /// @return the corresponding surfel.
    const SCell& surfel( Vertex v ) const
    {
      return theSurface->surfel( v );
    }

    /// @param[in] a any arc (index).
    /// @return the corresponding separator linel.
    const SCell& linel( Arc a ) const
    {
      return theSurface->linel( a );
    }

    /// @param[in] f any face index.
    /// @return the corresponding pivot pointel.
    const SCell& pointel( Face f ) const
    {
      return theSurface->pointel( f );
    }
    
    /// @param[in] aSurfel any surfel of the surface
    ///
    /// @return the vertex (ie an index) corresponding to this surfel,
    /// or INVALID_FACE if it does not exist.
    Vertex getVertex( const SCell& aSurfel ) const
    {
      return theSurface->getVertex( aSurfel );
    }

    /// @param[in] aLinel any linel that is a separator on the surface (orientation is important).
    ///
    /// @return the arc (ie an index) corresponding to this separator linel,
    /// or INVALID_FACE if it does not exist.
    Arc getArc( const SCell& aLinel ) const
    {
      return theSurface->getArc( aLinel );
    }

    /// @param[in] aPointel any pointel that is a pivot on the surface (orientation is positive).
    ///
    /// @return the face (ie an index) corresponding to this pivot pointel,
    /// or INVALID_FACE if it does not exist.
    Face getFace( const SCell& aPointel ) const
    {
      return theSurface->getFace( aPointel );
    }

    
    // ----------------------- Measure services -------------------------------
  public:

    /// @param[in] c any signed cell.
    /// @return the centroid of the signed cell c in real space;
    RealPoint vertexCentroid( Vertex v )
    {
      return myVertexCentroids[ v ];
    }

    /// @param[in] c any cell.
    /// @return the centroid of the cell c in real space;
    RealPoint computeCentroid( Cell c ) const
    {
      Point    kp = space().uKCoords( c );
      return   0.5 * myEmbedder( kp );
    }
    
    
    /// Computes the d-dimensional Hausdorff measure of a d-dimensional cell.
    Scalar Hmeasure( Dimension d ) const
    {
      Scalar H = 1.0;
      for ( Dimension k = d; k != 0; --k ) H *= myH;
      return H;
    }

    Scalar sRelativeIntersection( RealPoint p, Scalar r, SCell c )
    {
      return relativeIntersection( p, r, space().unsigns( c ) );
    }

    Scalar sCrispIntersection( RealPoint p, Scalar r, SCell c )
    {
      return crispIntersection( p, r, space().unsigns( c ) );
    }

    /// Computes an approximation of the relative intersection of the
    /// ball of radius \a r and center \a p with the given cell. Cells
    /// are embedded naturally in the grid of step h.
    ///
    /// @return the relative intersection as a scalar between 0 (no
    /// intersection) and 1 (inclusion).
    Scalar relativeIntersection( RealPoint p, Scalar r, Cell c )
    {
      const KSpace & K = space();
      auto       faces = K.uFaces( c );
      faces.push_back( c );
      bool       first = true;
      Scalar     d_max = 0.0;
      Scalar     d_min = 0.0;
      for ( Cell f : faces )
	{
	  RealPoint x = computeCentroid( f );
	  Scalar    d = ( x - p ).norm();
	  if ( first ) { d_min = d; first = false; }
	  d_max = std::max( d_max, d );
	  d_min = std::min( d_min, d );
	}
      if      ( d_max <= r     ) return 1.0;
      else if ( r     <= d_min ) return 0.0;
      return ( r - d_min ) / ( d_max - d_min );
    }

    /// Computes an approximation of the crisp intersection of the
    /// ball of radius \a r and center \a p with the given cell. Cells
    /// are embedded naturally in the grid of step h.
    ///
    /// @return the crisp intersection as the scalar 0 (no
    /// intersection) and 1 (intersection).
    Scalar crispIntersection( RealPoint p, Scalar r, Cell c )
    {
      const KSpace & K = space();
      RealPoint x = computeCentroid( c );
      Scalar    d = ( x - p ).norm();
      return d <= r ? 1.0 : 0.0;
    }

    /// \f$ \mu_0 \f$ Lipschitz-Killing measure. It corresponds to a
    /// corrected area measure, and is non null only on 2-cells (or Vertex).
    ///
    /// @param[in] c any 2-dimensional cell (or a Vertex in 3D digital
    /// surfaces).
    ///
    /// @return the corrected area measure \f$ \mu_0 := \cos \alpha
    /// d\mathcal{H}^2 \f$.
    Scalar mu0( Vertex c )
    {
      return myMu0[ c ];
    }
    
    /// \f$ \mu_0 \f$ Lipschitz-Killing measure. It corresponds to a
    /// corrected area measure, and is non null only on 2-cells (or Vertex).
    ///
    /// @param[in] c any 2-dimensional cell (or a Vertex in 3D digital
    /// surfaces).
    ///
    /// @return the corrected area measure \f$ \mu_0 := \cos \alpha
    /// d\mathcal{H}^2 \f$.
    Scalar computeMu0( Vertex c )
    {
      return Hmeasure( 2 ) * myTrivialNormals[ c ].dot( myCorrectedNormals[ c ] );
    }

    /// \f$ \mu_1 \f$ Lipschitz-Killing measure. It corresponds to a
    /// measure of mean curvature, since normal vectors may turn
    /// around an edge, and is non null only on 1-cells (or Arc). The
    /// measure is oriented, but gives the same result of an arc and
    /// its opposite.
    ///
    /// @param[in] a any arc (ie a 1-cell in-between two 2-cells on 3D digital
    /// surfaces).
    ///
    /// @return the corrected mean curvature measure \f$ \mu_1 := \Psi
    /// \vec{e} \cdot \vec{e}_1 d\mathcal{H}^1 \f$.
    Scalar mu1( Arc arc )
    {
      return myMu1[ arc ];
    }

    /// Contains all the information necessary to compute
    /// lipschitz-killing measures mu_k along a given arc.
    struct LocalArcInformation
    {
      Arc        arc; ///< the arc under interest
      RealVector ui;  ///< the corrected normal vecteur to the left of the arc
      RealVector uj;  ///< the corrected normal vecteur to the right of the arc
      RealVector e;   ///< the unit vector along the linel
      RealVector e1;  ///< u_i x u_j / || u_i x u_j ||
      Scalar     psi; ///< the dihedral angle from u_i to u_j
    };

    /// Computes all the information necessary to compute
    /// lipschitz-killing measures mu_k along a given arc.
    /// @param[out] lai the useful curvature information for the given \a arc
    /// @param[in]  arc any arc.
    void computeLocalArcInformation( LocalArcInformation& lai, Arc arc )
    {
      const KSpace & K = space();
      Vertex   si_plus = theSurface->tail( arc );
      Vertex  si_minus = theSurface->head( arc );
      Surfel    s_plus = theSurface->surfel( si_plus );
      Surfel   s_minus = theSurface->surfel( si_minus );
      Cell       linel = K.unsigns( theSurface->linel( arc ) ); 
      Dimension      l = *( K.uDirs( linel ) );
      Cell         pta = K.uIncident( linel, l, true );
      Cell         ptb = K.uIncident( linel, l, false );
      auto       faces = theSurface->facesAroundArc( arc );
      if ( faces.size() != 1 )
	faces = theSurface->facesAroundArc( theSurface->opposite( arc ) );	
      Cell       pivot = K.unsigns( theSurface->pointel( faces[ 0 ] ) );
      if ( pivot != pta ) std::swap( pta, ptb );
      RealPoint      a = computeCentroid( pta );
      lai.e            = ( computeCentroid( ptb ) - a ).getNormalized();
      lai.ui           = myCorrectedNormals[ si_plus ];
      lai.uj           = myCorrectedNormals[ si_minus ];
      RealVector  u_e1 = lai.ui.crossProduct( lai.uj );
      Scalar       ne1 = std::min( u_e1.norm(), 1.0 );
      lai.e1           = (ne1 != 0.0 ) ? ( u_e1 / ne1 ) : lai.e;
      lai.psi          = asin( ne1 );
    }

    /// \f$ \mu_1 \f$ Lipschitz-Killing measure. It corresponds to a
    /// measure of mean curvature, since normal vectors may turn
    /// around an edge, and is non null only on 1-cells (or Arc). The
    /// measure is oriented, but gives the same result of an arc and
    /// its opposite.
    ///
    /// @param[in] a any arc (ie a 1-cell in-between two 2-cells on 3D digital
    /// surfaces).
    ///
    /// @return the corrected mean curvature measure \f$ \mu_1 := \Psi
    /// \vec{e} \cdot \vec{e}_1 d\mathcal{H}^1 \f$.
    Scalar computeMu1( Arc arc )
    {
      LocalArcInformation info;
      computeLocalArcInformation( info, arc );
      return  Hmeasure( 1 ) * info.psi * info.e.dot( info.e1 );
    }

    /// \f$ \mu_2 \f$ Lipschitz-Killing measure. It corresponds to a
    /// measure of Gaussian curvature, since normal vectors form a cone
    /// around a vertex, and is non null only on 0-cells (or Face). The
    /// measure is oriented, but gives the same result on a face or its opposite.
    ///
    /// @param[in] f any face (ie a 0-cell in-between several 2-cells on 3D digital
    /// surfaces).
    ///
    /// @return the corrected Gaussian curvature measure \f$ \mu_0
    /// \f$, which is the area of the spherical polygon made by the
    /// normals onto the Gauss sphere.
    Scalar mu2( Face face )
    {
      return myMu2[ face ];
    }

    /// \f$ \mu_2 \f$ Lipschitz-Killing measure. It corresponds to a
    /// measure of Gaussian curvature, since normal vectors form a cone
    /// around a vertex, and is non null only on 0-cells (or Face). The
    /// measure is oriented, but gives the same result on a face or its opposite.
    ///
    /// @param[in] f any face (ie a 0-cell in-between several 2-cells on 3D digital
    /// surfaces).
    ///
    /// @return the corrected Gaussian curvature measure \f$ \mu_0
    /// \f$, which is the area of the spherical polygon made by the
    /// normals onto the Gauss sphere.
    Scalar computeMu2( Face face )
    {
      VertexRange     vtcs = theSurface->verticesAroundFace( face );
      const unsigned int n = vtcs.size();
      if ( n < 3 )
	{
	  trace.warning() << "[FastCorrectedNormalCurrent::computeMu2]"
			  << " Only " << n << " vertices around face "
			  << face << ", ie pivot " << theSurface->pointel( face ) << std::endl;
	  return 0.0;
	}
      // std::cout << n << std::endl;
      std::vector< RealVector > normals( n );
      for ( unsigned int i = 0; i < n; ++i )
	normals[ i ] = myCorrectedNormals[ vtcs[ i ] ];
      Scalar S = 0.0;
      for ( unsigned int i = 0; i <= n-3; ++i )
	{
	  SphericalTriangle<Space> ST( normals[ 0 ], normals[ i+2 ], normals[ i+1 ] );
	  S += ST.algebraicArea();
	}
      return S;
    }

    /// \f$ \mu_Omega \f$ Lipschitz-Killing measure, i.e. the measure
    /// associated to the symplectic form. It corresponds to a measure
    /// of the inconcistency between the given normal field and the
    /// geometry normal field. It may be non null only on 1-cells (or
    /// Arc). The measure is oriented, but gives the same result of an
    /// arc and its opposite.
    ///
    /// @param[in] a any arc (ie a 1-cell in-between two 2-cells on 3D digital
    /// surfaces).
    ///
    /// @return the corrected symplectic form measure \f$ \mu_Omega := \vec{e}
    /// \cdot (\vec{u}_i - \vec{u}_j) d\mathcal{H}^1 \f$.
    Scalar muOmega( Arc arc )
    {
      return myMuOmega[ arc ];
    }
    
    /// \f$ \mu_Omega \f$ Lipschitz-Killing measure, i.e. the measure
    /// associated to the symplectic form. It corresponds to a measure
    /// of the inconcistency between the given normal field and the
    /// geometry normal field. It may be non null only on 1-cells (or
    /// Arc). The measure is oriented, but gives the same result of an
    /// arc and its opposite.
    ///
    /// @param[in] a any arc (ie a 1-cell in-between two 2-cells on 3D digital
    /// surfaces).
    ///
    /// @return the corrected  symplectic form measure \f$ \mu_Omega := \vec{e}
    /// \cdot (\vec{u}_i - \vec{u}_j) d\mathcal{H}^1 \f$.
    Scalar computeMuOmega( Arc arc )
    {
      const KSpace & K = space();
      Vertex   si_plus = theSurface->tail( arc );
      Vertex  si_minus = theSurface->head( arc );
      Surfel    s_plus = theSurface->surfel( si_plus );
      Surfel   s_minus = theSurface->surfel( si_minus );
      Cell       linel = K.unsigns( theSurface->linel( arc ) ); 
      Dimension      l = *( K.uDirs( linel ) );
      Cell         pta = K.uIncident( linel, l, true );
      Cell         ptb = K.uIncident( linel, l, false );
      auto       faces = theSurface->facesAroundArc( arc );
      if ( faces.size() != 1 )
	faces = theSurface->facesAroundArc( theSurface->opposite( arc ) );	
      Cell       pivot = K.unsigns( theSurface->pointel( faces[ 0 ] ) );
      if ( pivot != pta ) std::swap( pta, ptb );
      RealPoint      a = computeCentroid( pta );
      RealVector     e = ( computeCentroid( ptb ) - a ).getNormalized();
      RealPoint     s0 = myVertexCentroids[ si_plus ];
      RealPoint     s1 = myVertexCentroids[ si_minus ];
      // s_plus must be to the left of e.
      // if ( e.crossProduct( s0 - a ).dot( myTrivialNormals[ s_plus ] ) < 0.0 )
      // 	   std::swap( s_plus, s_minus );
      // Computes u_+ and u_-, then their cross product.
      RealVector    u_p = myCorrectedNormals[ si_plus ];
      RealVector    u_m = myCorrectedNormals[ si_minus ];
      return  Hmeasure( 1 ) * e.dot( u_p - u_m );
    }

    /// Computes the anisotropic measure mu_XY according to the
    /// information \a lai and to the directions \a X and \a Y.
    ///
    /// @param lai the local curvature measure of an arc.
    /// @param X any direction
    /// @param Y any direction
    Scalar getMuXY( const LocalArcInformation& lai, RealVector X, RealVector Y )
    {
      RealVector       X_x_e = X.crossProduct( lai.e );
      RealVector     e1_x_ui = lai.e1.crossProduct( lai.ui );
      Scalar          Y_d_ui = Y.dot( lai.ui );
      Scalar      X_x_e_d_ui = X_x_e.dot( lai.ui );
      Scalar     Y_d_e1_x_ui = Y.dot( e1_x_ui );
      Scalar X_x_e_d_e1_x_ui = X_x_e.dot( e1_x_ui );
      Scalar         sin_psi = sin( lai.psi );
      Scalar        sin_2psi = sin( 2.0*lai.psi );
      // v1
      // Scalar m = ( lai.psi + sin_2psi ) * X_x_e_d_e1_x_ui*Y_d_ui
      // 	- ( sin_2psi - lai.psi ) * X_x_e_d_ui*Y_d_e1_x_ui
      // 	+ sin_psi * sin_psi * ( X_x_e_d_ui*Y_d_ui - X_x_e_d_e1_x_ui*Y_d_e1_x_ui );
      // return 0.5 * m;
      // v2
      // Scalar m = (1.0 + sin_2psi ) * X_x_e_d_e1_x_ui*Y_d_ui
      // 	- ( 1.0 - sin_2psi ) * X_x_e_d_ui*Y_d_e1_x_ui
      // 	+ sin_psi * sin_psi * ( X_x_e_d_ui*Y_d_ui - X_x_e_d_e1_x_ui*Y_d_e1_x_ui );
      // return 0.5 * lai.psi * m;
      Scalar m = (-lai.psi - 0.5*sin_2psi)*X_x_e_d_ui*Y_d_e1_x_ui
	+(sin_psi*sin_psi) * ( X_x_e_d_ui*Y_d_ui - X_x_e_d_e1_x_ui*Y_d_e1_x_ui )
	+(lai.psi - 0.5*sin_2psi) *  X_x_e_d_e1_x_ui*Y_d_ui;
      return 0.5 * m;
    }

    /// Anisotropic curvature measure. It is a 3x3 tensor that more or
    /// less describes the curvature tensor.
    ///
    /// @param[in] arc any arc (ie a 1-cell in-between two 2-cells on 3D digital
    /// surfaces).
    ///
    /// @return the anisotropic curvature measure along the \a arc.
    RealTensor computeAnisotropicMu( Arc arc )
    {
      LocalArcInformation info;
      computeLocalArcInformation( info, arc );
      RealTensor T;
      Scalar h = Hmeasure( 1 ); 
      for ( Dimension i = 0; i < 3; ++i ) {
	RealVector X = RealVector::base( i, 1.0 );
	for ( Dimension j = 0; j < 3; ++j ) {
	  RealVector Y = RealVector::base( j, 1.0 );
	  T.setComponent( i, j, h * getMuXY( info, X, Y ) );
	}
      }
      return T;
    }
    
    /// Anisotropic curvature measure. It is a 3x3 tensor that more or
    /// less describes the curvature tensor.
    ///
    /// @param[in] arc any arc (ie a 1-cell in-between two 2-cells on 3D digital
    /// surfaces).
    ///
    /// @return the anisotropic curvature measure along the \a arc.
    RealTensor anisotropicMu( Arc arc )
    {
      return myAnisotropicMu[ arc ];
    }
      
    
    Scalar mu0( RealPoint p, Scalar r, Vertex c )
    {
      SCell aSurfel = theSurface->surfel( c );
      Scalar     ri = myCrisp
	? sCrispIntersection( p, r, aSurfel )
	: sRelativeIntersection( p, r, aSurfel );
      return ri != 0.0 ? ri * mu0( c ) : 0.0;
    }

    Scalar mu1( RealPoint p, Scalar r, Arc a )
    {
      SCell aLinel = theSurface->linel( a ); // oriented 1-cell
      Scalar    ri = myCrisp
	? sCrispIntersection( p, r, aLinel )
	: sRelativeIntersection( p, r, aLinel );
      return ri != 0.0 ? ri * mu1( a ) : 0.0;
    }

    Scalar mu2( RealPoint p, Scalar r, const Face& f )
    {
      SCell aPointel = theSurface->pointel( f ); 
      Scalar      ri = myCrisp
	? sCrispIntersection( p, r, aPointel )
	: sRelativeIntersection( p, r, aPointel );
      return ri != 0.0 ? ri * mu2( f ) : 0.0;
    }

    Scalar muOmega( RealPoint p, Scalar r, Arc a )
    {
      SCell aLinel = theSurface->linel( a ); // oriented 1-cell
      Scalar    ri = myCrisp
	? sCrispIntersection( p, r, aLinel )
	: sRelativeIntersection( p, r, aLinel );
      return ri != 0.0 ? ri * muOmega( a ) : 0.0;
    }

    RealTensor anisotropicMu( RealPoint p, Scalar r, Arc a )
    {
      SCell aLinel = theSurface->linel( a ); // oriented 1-cell
      Scalar    ri = myCrisp
	? sCrispIntersection( p, r, aLinel )
	: sRelativeIntersection( p, r, aLinel );
      return ( ri != 0.0 ) ? ( ri * anisotropicMu( a ) ) : RealTensor();
    }

    
    struct SquaredDistance2Point {
      typedef Scalar Value;
      FastCorrectedNormalCurrent& current;
      RealPoint center;
      SquaredDistance2Point( FastCorrectedNormalCurrent& aCurrent,
			     const RealPoint&              aPoint )
	: current( aCurrent ), center( aPoint ) {}
      Value operator() ( Vertex v ) const
      {
	RealPoint  x = current.vertexCentroid( v );
	RealVector w = x - center;
	return w.dot( w );
      }
    };
    
    VertexRange getVerticesInBall( Vertex vc, Scalar r )
    {
      typedef DistanceBreadthFirstVisitor
	< Surface, SquaredDistance2Point >     DistanceVisitor;
      typedef typename DistanceVisitor::Node   MyNode;
      typedef typename DistanceVisitor::Scalar MySize;
      
      VertexRange output;
      RealPoint   center = myVertexCentroids[ vc ];
      Scalar       limit = r*r;
      SquaredDistance2Point d2pfct( *this, center );
      DistanceVisitor       visitor( *theSurface, d2pfct, vc );
      while ( ! visitor.finished() )
	{
	  MyNode node = visitor.current();
	  Vertex    v = node.first;
	  Scalar    d = node.second;
	  output.push_back( v );
	  if ( d <= limit ) visitor.expand();
	  else              visitor.ignore();
	}
      return output;
    }

    ArcRange getArcsInBall( Vertex vc, Scalar r )
    {
      VertexRange scells = getVerticesInBall( vc, r );
      std::set<Arc> arcs;
      for ( auto s : scells )
	{
	  auto l_arcs = theSurface->outArcs( s );
	  for ( auto a : l_arcs )
	    {
	      auto b = theSurface->opposite( a );
	      if ( arcs.find( a ) == arcs.end() && arcs.find( b ) == arcs.end() )
		arcs.insert( a );
	    }
	}
      return ArcRange( arcs.begin(), arcs.end() );
    }

    FaceRange getFacesInBall( Vertex vc, Scalar r )
    {
      VertexRange scells = getVerticesInBall( vc, r );
      std::set<Face> faces;
      for ( auto s : scells )
	{
	  auto l_faces = theSurface->facesAroundVertex( s );
	  faces.insert( l_faces.begin(), l_faces.end() );
	}
      return FaceRange( faces.begin(), faces.end() );
    }

    Scalar mu0Ball( Vertex vc, Scalar r )
    {
      VertexRange vtcs = getVerticesInBall( vc, r );
      Scalar        m0 = 0.0;
      RealPoint      x = myVertexCentroids[ vc ];
      for ( auto v : vtcs ) {
	m0      += mu0( x, r, v ); 
      }
      return m0;
    }

    Scalar mu1Ball( Vertex vc, Scalar r )
    {
      ArcRange    arcs = getArcsInBall( vc, r );
      Scalar        m1 = 0.0;
      RealPoint      x = myVertexCentroids[ vc ];
      for ( auto a : arcs ) {
	m1      += mu1( x, r, a ); 
      }
      return m1;
    }

    Scalar mu2Ball( Vertex vc, Scalar r )
    {
      // std::cout << "mu2ball c=" << c << " r=" << r << std::endl;
      FaceRange  faces = getFacesInBall( vc, r );
      // std::cout << "mu2ball #faces=" << faces.size() << std::endl;
      Scalar        m2 = 0.0;
      RealPoint      x = myVertexCentroids[ vc ];
      for ( auto f : faces ) {
	m2      += mu2( x, r, f ); 
      }
      return m2;
    }

    Scalar muOmegaBall( Vertex vc, Scalar r )
    {
      ArcRange    arcs = getArcsInBall( vc, r );
      Scalar        m1 = 0.0;
      RealPoint      x = myVertexCentroids[ vc ];
      for ( auto a : arcs ) {
	m1      += muOmega( x, r, a ); 
      }
      return m1;
    }

    RealTensor anisotropicMuBall( Vertex vc, Scalar r )
    {
      ArcRange    arcs = getArcsInBall( vc, r );
      RealTensor    mT = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, };
      RealPoint      x = myVertexCentroids[ vc ];
      for ( auto a : arcs ) {
	mT      += anisotropicMu( x, r, a ); 
      }
      return mT;
    }
    
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

    // ------------------------- Private Datas --------------------------------
  private:

    /// A smart securable pointer onto the digital surface.
    CountedPtrOrPtr<Surface> theSurface;
    /// The digitization grid step.
    Scalar                   myH;
    /// Specifies how intersection are computed.
    bool                     myCrisp;
    /// The standard embedding with gridstep h.
    RegularPointEmbedder<Space> myEmbedder;
    /// The natural normal vector field.
    NormalVectorField        myTrivialNormals;
    /// The corrected normal vector field.
    NormalVectorField        myCorrectedNormals;
    /// The map vertex -> centroid (to limit computations).
    CentroidMap              myVertexCentroids;
    /// The map arc -> centroid (to limit computations).
    CentroidMap              myArcCentroids;
    /// The map face -> centroid (to limit computations).
    CentroidMap              myFaceCentroids;
    /// The map Vertex -> mu0
    MeasureMap               myMu0;
    /// The map Arc -> mu1
    MeasureMap               myMu1;
    /// The map Face -> mu2
    MeasureMap               myMu2;
    /// The map Arc -> muOmega
    MeasureMap               myMuOmega;
    /// The map Arc -> AnisotropicMu
    TensorMeasureMap         myAnisotropicMu;
    
    // ------------------------- Hidden services ------------------------------
  protected:


    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class FastCorrectedNormalCurrent


  /**
   * Overloads 'operator<<' for displaying objects of class 'FastCorrectedNormalCurrent'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'FastCorrectedNormalCurrent' to write.
   * @return the output stream after the writing.
   */
  template <typename Surface>
  std::ostream&
  operator<< ( std::ostream & out, const FastCorrectedNormalCurrent<Surface> & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "FastCorrectedNormalCurrent.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined FastCorrectedNormalCurrent_h

#undef FastCorrectedNormalCurrent_RECURSES
#endif // else defined(FastCorrectedNormalCurrent_RECURSES)
