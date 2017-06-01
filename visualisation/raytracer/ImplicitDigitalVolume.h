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
 * @file ImplicitDigitalVolume.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2017/02/28
 *
 * This file is part of the DGtal library.
 */

#if defined(ImplicitDigitalVolume_RECURSES)
#error Recursive header files inclusion detected in ImplicitDigitalVolume.h
#else // defined(ImplicitDigitalVolume_RECURSES)
/** Prevents recursive inclusion of headers. */
#define ImplicitDigitalVolume_RECURSES

#if !defined ImplicitDigitalVolume_h
/** Prevents repeated inclusion of headers. */
#define ImplicitDigitalVolume_h

#include <map>
#include "raytracer/GeometricalObject.h"
#include "raytracer/Parallelogram.h"
#include "raytracer/StandardDSL3d.h"
#include "DGtal/base/Clone.h"
#include "DGtal/topology/SetOfSurfels.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"

namespace DGtal {
  namespace rt {

    /// A triangle is a model of GraphicalObject.
    template <typename TImage>
    struct ImplicitDigitalVolume : public virtual GeometricalObject {

      typedef TImage                           Image;
      typedef KSpace3                          KSpace;
      typedef typename TImage::Value           Value;
      typedef functors::SimpleThresholdForegroundPredicate<Image> ThresholdedImage;
      typedef KSpace3::SCell                   SCell;
      typedef KSpace3::Cell                    Cell;
      typedef KSpace3::SCellSet                SCellSet;
      typedef SetOfSurfels<KSpace3, SCellSet>  SurfaceStorage;
      typedef DigitalSurface< SurfaceStorage > Surface;
      typedef std::map<SCell,Vector3>          VectorField;
      typedef VectorField                      NormalMap;

      typedef KSpace::Space                         Space;
      typedef Space::Point                          Point;
      typedef Space::RealVector                     RealVector;
      typedef RealVector::Component                 Scalar;
      typedef KSpace::Surfel                        Surfel;
      typedef HyperRectDomain<Space>                Domain;
      typedef ImplicitPolynomial3Shape<Space>       PolynomialShape;              
      typedef typename PolynomialShape::Polynomial3 Polynomial;
      typedef std::map<Point,Polynomial>            PolynomialStore;
      typedef ImageContainerBySTLVector<Domain, Vector3> VectorImage;
      BOOST_STATIC_ASSERT(( KSpace::dimension == 3 ));

      
      inline
      ImplicitDigitalVolume( Clone<TImage> anImage, Value threshold )
        : myImage( anImage ), myThresholdedImage( myImage, threshold ),
          myCfgImage( myImage.domain() ), myGradientImage( myImage.domain() ),
          myT( threshold ), ptrSurface( 0 )
      {
        // // We need to extend the domain since we are not considering voxels, but 8-cubes.
        // Domain extDomain( myImage.domain().lowerBound() - Point3i::diagonal( 1 ),
        //                   myImage.domain().upperBound() + Point3i::diagonal( 1 ) );
        K.init( myImage.domain().lowerBound(),
                myImage.domain().upperBound(), true );
        Point3i lo   = K.uCoords( K.lowerCell() ); 
        Point3i hi   = K.uCoords( K.upperCell() );
        Point3 A000 = Point3( lo ); // - Point3::diagonal( 0.5 );
        Point3 A111 = Point3( hi ); // - Point3::diagonal( 0.5 );
        Point3 A100( A111[ 0 ], A000[ 1 ], A000[ 2 ] );
        Point3 A010( A000[ 0 ], A111[ 1 ], A000[ 2 ] );
        Point3 A001( A000[ 0 ], A000[ 1 ], A111[ 2 ] );
        Point3 A110( A111[ 0 ], A111[ 1 ], A000[ 2 ] );
        Point3 A011( A000[ 0 ], A111[ 1 ], A111[ 2 ] );
        Point3 A101( A111[ 0 ], A000[ 1 ], A111[ 2 ] );
        sides.push_back( Parallelogram( A000, A010, A100 ) ); // bottom (xy0)
        sides.push_back( Parallelogram( A000, A001, A010 ) ); // left   (0yz)
        sides.push_back( Parallelogram( A000, A100, A001 ) ); // front  (x0z)
        sides.push_back( Parallelogram( A111, A011, A101 ) ); // top    (xy1)
        sides.push_back( Parallelogram( A111, A101, A110 ) ); // right  (1yz)
        sides.push_back( Parallelogram( A111, A110, A011 ) ); // back   (x1z)
	trace.beginBlock( "Compute boundary" );
        SCellSet boundary;
        Surfaces<KSpace3>::sMakeBoundary( boundary, K, myThresholdedImage,
                                          K.lowerBound(), K.upperBound() );
        ptrSurface = new Surface( new SurfaceStorage( K, true, boundary ) );
	trace.endBlock();
        // interpolateVectorField( trivialNormals, true );
        // Init polynomials
        Polynomial X[ 2 ][ 3 ];
        X[ 1 ][ 0 ] = mmonomial<Scalar>( 1, 0, 0 );
        X[ 1 ][ 1 ] = mmonomial<Scalar>( 0, 1, 0 );
        X[ 1 ][ 2 ] = mmonomial<Scalar>( 0, 0, 1 );
        for ( Dimension i = 0; i < 3; ++i )
          X[ 0 ][ i ] = mmonomial<Scalar>( 0, 0, 0 ) - X[ 1 ][ i ];
        for ( Dimension i = 0; i < 8; ++i )
          {
            P[ i ] = X[ (i & 1) ? 1 : 0 ][ 0 ]
              *      X[ (i & 2) ? 1 : 0 ][ 1 ]
              *      X[ (i & 4) ? 1 : 0 ][ 2 ];
            //std::cout << "P[ " << i << " ] = " << P[ i ] << std::endl;
          }
	trace.beginBlock( "Compute gradients" );
        auto it = myCfgImage.begin();
        for ( auto p : myImage.domain() )
	  {
	    *it++ = getConfiguration( p );
	    Vector3   g = computeGradientAtVoxel( p );
	    myGradientImage.setValue( p,  g );
	    //std::cout << "At p=" << p << " g=" << g << std::endl;
	  }
	for ( auto surfel : surface() )
	  {
	    Dimension k = space().sOrthDir( surfel );
	    Point3i  p0 = space().sCoords( space().sIncident( surfel, k, false ) );
	    Point3i  p1 = p0 + Point3i::base( k );
	    Vector3   g = computeGradientAtSurfel( surfel );
	    myGradientImage.setValue( p0,  g );
	    myGradientImage.setValue( p1,  g );
	    //std::cout << "At surfel=" << surfel << " g=" << g << std::endl;
	  }
	trace.endBlock();
      }
    
      /// Virtual destructor since object contains virtual methods.
      ~ImplicitDigitalVolume() {}

      
      Image&            image()         { return myImage; }
      ThresholdedImage& thresholdedImage()   { return myThresholdedImage; }
      KSpace3&          space()         { return K; }
      const KSpace3&    space() const   { return K; }
      Surface&          surface()       { return *ptrSurface; }
      const Surface&    surface() const { return *ptrSurface; }

      /// Secured accessor to voxel value.
      Value at( Point3i voxel ) const
      { return myImage.domain().isInside( voxel ) ? myImage( voxel ) : 0; }

      /// @return the gradient at the given surfel.
      Vector3 computeGradientAtSurfel( const SCell& surfel ) const
      {
	Vector3 grad;
	Dimension k = space().sOrthDir( surfel );
	Point3i  p0 = space().sCoords( space().sIncident( surfel, k, false ) );
	Point3i  p1 = p0 + Point3i::base( k );
	Dimension i = (k+1)%3;
	Dimension j = (k+2)%3;
	grad[ k ]   = (Scalar) at( p1 ) - (Scalar) at( p0 );
	grad[ i ]   =
	  0.25*( (Scalar) at( p0 + Point3i::base( i ) ) - (Scalar) at( p0 - Point3i::base( i ) )
		 + (Scalar) at( p1 + Point3i::base( i ) ) - (Scalar) at( p1 - Point3i::base( i ) ) );
	grad[ j ]   =
	  0.25*( (Scalar) at( p0 + Point3i::base( j ) ) - (Scalar) at( p0 - Point3i::base( j ) )
		 + (Scalar) at( p1 + Point3i::base( j ) ) - (Scalar) at( p1 - Point3i::base( j ) ) );
	return grad != Point3::zero ? grad.getNormalized() : grad;
      }

      /// @return the gradient at the given surfel.
      Vector3 computeGradientAtVoxel( const Point3i& p ) const
      {
	Vector3 grad;
	for ( Dimension i = 0; i < 3; i++ )
	  grad[ i ]   =
	    0.5*( (Scalar) at( p + Point3i::base( i ) ) - (Scalar) at( p - Point3i::base( i ) ) );
	return grad != Point3::zero ? grad.getNormalized() : grad;
      }
      
      Vector3 interpolateGradient( const Point3i& base, const Point3& p ) const
      {
	const Scalar  a = p[ 0 ];
	const Scalar  b = p[ 1 ];
	const Scalar  c = p[ 2 ];
	const Scalar ma = 1.0 - a;
	const Scalar mb = 1.0 - b;
	const Scalar mc = 1.0 - c;
	Vector3 grad = (ma*mb*mc) * myGradientImage( base )
	  +            ( a*mb*mc) * myGradientImage( base + Point3i(1,0,0) )
	  +            (ma* b*mc) * myGradientImage( base + Point3i(0,1,0) )
	  +            ( a* b*mc) * myGradientImage( base + Point3i(1,1,0) )
	  +            (ma*mb* c) * myGradientImage( base + Point3i(0,0,1) )
	  +            ( a*mb* c) * myGradientImage( base + Point3i(1,0,1) )
	  +            (ma* b* c) * myGradientImage( base + Point3i(0,1,1) )
	  +            ( a* b* c) * myGradientImage( base + Point3i(1,1,1) );
	return grad.getNormalized();
      }

      /// @return the normal vector at point \a p on the object (\a p
      /// should be on or close to the sphere).
      virtual Vector3 getNormal( Point3 p )
      {
        if ( p == last_p )
          return last_n;
        trace.warning() << "[ImplicitDigitalVolume::getNormal] Movement after ray ?"
                        << " p=" << p << " last_p=" << last_p
                        << " last_n=" << last_n << std::endl;
        return last_n;
      }

      /// @param[in,out] ray_inter as input the incoming ray, as
      /// output information abour intersection.
      ///
      /// @return true if there was an intersection, false otherwise
      /// (more information is stored in ray_inter)
      virtual bool intersectRay( const Ray& ray, RayIntersection& ray_inter )
      {
        // Checks first that it intersects the bounding box.
        Point3 p  [ 6 ];
        bool   hit[ 6 ];
        int    nb_hit  = 0;
        Real   dist    = 100000000.0;
        Real   alpha   = 0.0;
        unsigned int j = 6;
        for ( unsigned int i = 0; i < 6; i++ )
          {
            hit[ i ] = sides[ i ].intersectRay( ray, ray_inter );
            Real d   = ray_inter.distance;
            p[ i ]   = ray_inter.intersection;
            dist     = std::min( d, dist );
            hit[ i ] = d < 0.0;
            if ( hit[ i ] )
              {
                nb_hit++;
                Real beta = ( p[ i ] - ray.origin ).dot( ray.direction );
                if ( ( beta > RT_EPSILON )
                     && ( ( j == 6 ) || ( beta < alpha ) ) )
                  {
                    j = i;
                    alpha = beta;
                  }
              }
          }
        ray_inter.distance = 1.0f;
        // ray_inter.distance = dist;
        if ( nb_hit == 0 ) return false;
        if ( j == 6 ) {
          trace.warning() << "No closest point !" << std::endl;
          return false;
        }
        // To check that the bounding box is correct.
        // p_intersect = p[ j ];
        // last_n = sides[ j ].getNormal( p_intersect );
        // return -1.0;

        // (0) Checks if origin is in the volume
        Point3  s_origin = ray.origin + 0.0001* ray.direction;
        Point3i origin ( (Integer) floor( s_origin[ 0 ] ),
                         (Integer) floor( s_origin[ 1 ] ),
                         (Integer) floor( s_origin[ 2 ] ) );
        bool origin_in = myImage.domain().isInside( origin );
        Point3  q      = origin_in ? s_origin : ( p[ j ] + 0.001* ray.direction );
        // trace.info() << "Inside: q=" << q << std::endl;
        // Now casts the ray in the volume
        // (1) build a digital ray
        Point3 shift_origin = q - Point3(0.5,0.5,0.5);
        Ray ray2( shift_origin,
                  ray.direction, ray.depth ); //p[ j ]
        StandardDSL3d D( ray2,
                         RT_PRECISION*( K.upperBound() - K.lowerBound() ).norm1() );
        // (2) sort points along ray
        Point3i first  = origin_in
          ? origin
          : Point3i( (Integer) floor( q[ 0 ] ),
                     (Integer) floor( q[ 1 ] ),
                     (Integer) floor( q[ 2 ] ) );
        StandardDSL3d::ConstIterator it  = D.begin( first );
        if ( ! D.isInDSL( *it  ) )
          {
            trace.warning() << " *it=" << *it << " not inside"
                            << " Dxy=" << D.xy
                            << " Dxz=" << D.xz
                            << " Dyz=" << D.yz << std::endl;
            return false;
          }
        if ( ! myImage.domain().isInside( *it  ) ) return false;
        if ( origin_in )   it++;
        // std::cout << "- DSL:";
        while ( true )
          {
            // This is the current 8-cube
            Point3i p  = *it;
            // std::cout << " (" << p[0] << "," << p[1] << "," << p[2] << ")";
            if ( ! myImage.domain().isInside( *it  ) )
              break;
            // Get configuration
            int cfg = myCfgImage( p ); // getConfiguration( p );
            if ( cfg > 0 && cfg < 255 ) // intersection possible
              {
                // std::cout << " *";
                if ( intersect8Cube( ray, ray_inter, p, cfg ) )
                  {
                    Scalar proj  =
                      ( ray_inter.intersection - ray.origin )
                      .dot( ray.direction );
                    const Point3& xi = ray_inter.intersection;
                    // std::cout << "(pi=" << proj << ")"
                    //           << "/(" << xi[0]<< "," << xi[1] << "," << xi[2] << ")"
                    //           << std::endl;
                    if ( proj > 0.01 )
                      return true;
                  }
              }
            ++it;
          }
        // std::cout << "[NO]" << std::endl;
        return false;
      }
      int getConfiguration( Point3i cube ) const
      {
        int cfg   = 0;
        int shift = 1;
        for ( Dimension l = 0; l < 8; l++, shift <<= 1 )
          {
            Point3i voxel( cube[ 0 ] + ( (l & 1) ? 1 : 0 ),
                           cube[ 1 ] + ( (l & 2) ? 1 : 0 ),
                           cube[ 2 ] + ( (l & 4) ? 1 : 0 ) );
            if ( myImage.domain().isInside( voxel ) && myThresholdedImage( voxel ) )
              cfg = cfg | shift;
          }
        return cfg;
      }

      // Computes and memorizes polynomial for the given base.
      void getPolynomial( Polynomial& bP, const Point3i& base )
      {
        if ( base == last_base )
          {
            bP = last_poly;
            return;
          }
        auto it = poly_mem.find( base );
        if ( it != poly_mem.end() )
          {
            bP = last_poly = it->second;
          }
        else
          {
            Polynomial& L = poly_mem[ base ];
            for ( Dimension l = 0; l < 8; l++ )
              {
                Point3i voxel( base[ 0 ] + ( (l & 1) ? 1 : 0 ),
                               base[ 1 ] + ( (l & 2) ? 1 : 0 ),
                               base[ 2 ] + ( (l & 4) ? 1 : 0 ) );
                Value v = myImage.domain().isInside( voxel ) ? myImage( voxel ) : 0;
                // std::cout << " v=" << (Scalar) v << endl;
                L      += ( ( (Scalar) v - (Scalar) myT ) / 256.0 ) * P[ l ];
              }
            bP = last_poly = L;
            // std::cout << "base=" << base << " L=" << L << endl;
          }
        last_base = base;
      }
      
      bool intersect8Cube( const Ray& ray, RayIntersection& ray_inter, Point3i base, int cfg )
      {
        // Make the local polynomial isosurface
        Polynomial L;
        getPolynomial( L, base );
        PolynomialShape PS( L );
        // std::cout << " PS=" << PS << std::endl;
        // Checks for intersection along ray
        Point3         p = Point3( 0.5, 0.5, 0.5 );
        // Project onto ray.
        ray_inter.distance = 1.0;
        const Point3&  o = ray.origin; // p + ray.origin;
        const Vector3& u = ray.direction;
        Point3  B = Point3( base[ 0 ], base[ 1 ], base[ 2 ] );
        Point3 rO = o - B;
        Point3 q  = rO + (p - rO).dot( u ) * u;
	Scalar d2 = (q-p).dot(q-p);
	for ( Dimension l = 0; l < 8; l++ )
	  if ( cfg & (1<<l) )
	    {
	      Point3 pp( (l & 1) ? 1 : 0, (l & 2) ? 1 : 0, (l & 4) ? 1 : 0 );
	      Point3 qq  = rO + (pp - rO).dot( u ) * u;
	      Scalar dd2 = (qq-pp).dot(qq-pp);
	      if ( dd2 < d2 )
		{
		  d2 = dd2;
		  q = qq;
		}
	    }
        bool found = false;
        found = intersectPolynomialShape( PS, u, q );
        if ( ! found ) return false;
        // if ( ! ( fabs( PS( q ) ) < 0.005 ) ) return false;
        if ( ( q[ 0 ] < -0.001 ) || ( q[ 0 ] > 1.001 )
             || ( q[ 1 ] < -0.001 ) || ( q[ 1 ] > 1.001 )
             || ( q[ 2 ] < -0.001 ) || ( q[ 2 ] > 1.001 ) ) return false;

	// GOURAUD type shading
        /* last_n     = -PS.gradient( q ); */
        /* if ( last_n == Vector3::zero ) return false; */
        /* last_n    /= last_n.norm(); */
	// PHONG shading
        last_n     = - interpolateGradient( base, q );
        // std::cout << "- base=" << base << " P=" << L << std::endl;
        // std::cout << "  => found q=" << q << " n=" << last_n << " pi=" << (q + B) << std::endl;
        q         += B + Point3::diagonal(0.5); // + p; // B - p;
        last_p     = q;
        ray_inter.normal       = last_n;
        ray_inter.intersection = q; // + 0.001*last_n;
        ray_inter.reflexion    = q + 0.001*last_n; // ray_inter.intersection;
        ray_inter.refraction   = q - 0.001*last_n;
        ray_inter.distance     = -1.0;
        return true;
      }


      bool intersectPolynomialShapeE( const PolynomialShape& PS, const Vector3& u,
                                      Point3& q )
      {
        Scalar min = PS( q );
        Scalar max = min;
        Vector3 t  = u * 0.25;
        Point3 bmin= q;
        Point3 bmax= q;
        Point3 q0  = q;
        Point3 q1  = q;
        for ( int i = 0; i < 4; i++ )
          {
            q0 -= t;
            Scalar v0 = PS( q0 );
            if ( v0 > max )      { bmax = q0; max = v0; }
            else if ( v0 < min ) { bmin = q0; min = v0; }
            q1 += t;
            Scalar v1 = PS( q1 );
            if ( v1 > max )      { bmax = q1; max = v1; }
            else if ( v1 < min ) { bmin = q1; min = v1; }
            if ( min <= 0.0 && max > 0.0 )
              {
                bool ok = intersectPolynomialShapeDicho( PS, q0, q1 );
                q = q0;
                return ok;
              }
          }
        return false;
      }

      bool intersectPolynomialShapeDicho( const PolynomialShape& PS,
                                          Point3& q0, Point3& q1 )
      {
        const Scalar eps = 0.0025;
        const int   iter = 12;
        for ( int n = 0; n < iter; n++ )
          {
            Point3 qm = 0.5*(q0+q1);
            Scalar  v = PS( qm );
            if ( fabs( v ) < eps ) { q0 = qm; return true; }
            if ( v > 0.0 ) q1 = qm;
            else           q0 = qm;
          }
        Scalar  v = PS( q0 );
        return ( fabs( v ) < eps );
      }

      bool intersectPolynomialShapeR( const PolynomialShape& PS, const Vector3& u,
                                      Point3& q )
      {
        const int   iter = 8;
        Point3 q0 = q-u; 
        Point3 q1 = q+u;
        Point3 qm = q;
        bool sq0 = PS( q0 ) >= 0.0;
        bool sq1 = PS( q1 ) >= 0.0;
        bool sqm = PS( qm ) >= 0.0;
        if ( ( sq0 == sqm ) && ( sq1 == sqm ) ) return false;
        if ( sq0 == sqm ) q0 = qm;
        else              q1 = qm;
        for ( int n = 0; n < iter; n++ )
          {
            qm       = 0.5*(q0+q1);
            bool sqm = PS( qm ) >= 0.0;
            if ( sqm == sq0 ) q0 = qm;
            else              q1 = qm;
          }
        q = qm;
        return true;
      }

      bool intersectPolynomialShape( const PolynomialShape& PS, const Vector3& u,
                                     Point3& q )
      {
        const int   iter = 12;
        //const Scalar att = 0.5;
        const Scalar eps = 0.005;
        Scalar      val = PS( q );
        for ( int n = 0; n < iter; n++ )
          {
            if ( fabs( val ) <= eps ) return true;
            Vector3  g = PS.gradient( q );
	    Scalar  ng = g.norm();
	    if ( ng < 0.00001 ) break;
            Scalar   t = - val * ng / g.dot( u );
	    q         += t * u;
            if ( q.normInfinity() > 2.0 ) return false;
            val        = PS( q );
            // std::cout << "PS(q)=" << val << std::endl;
          }
        return ( fabs( val ) <= eps );
      }
      
      // ----------------------------------------------------------------------
      
      /// A clone of the image.
      Image            myImage;
      /// The thresholded image.
      ThresholdedImage myThresholdedImage;
      /// The configuration image.
      Image            myCfgImage;
      /// The gradient image (for Phong shading).
      VectorImage      myGradientImage;
      /// The threshold
      Value            myT;
      /// the Khalimsky space
      KSpace3       K;
      /// Polynomial associated with each vertex of the 8-cube;
      Polynomial    P[ 8 ];
      /// Sides of the bounding box of the digital volume.
      std::vector<Parallelogram> sides;
      /// The digital surface corresponding to the boundary of bimage.
      Surface*      ptrSurface;
      /// Last point of intersection
      Point3        last_p;
      /// Last normal at intersection
      Vector3       last_n;
      /// Last surfel at intersection
      SCell         last_surfel;
      /// Last parallelogram at intersection
      Parallelogram   last_square;
      Point3i         last_base;
      Polynomial      last_poly;
      PolynomialStore poly_mem;
    };

  } // namespace rt
} // namespace DGtal

#endif // !defined ImplicitDigitalVolume_h

#undef ImplicitDigitalVolume_RECURSES
#endif // else defined(ImplicitDigitalVolume_RECURSES)
