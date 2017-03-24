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
 * @file StandardDSL3d.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2017/02/28
 *
 * This file is part of the DGtal library.
 */

#if defined(StandardDSL3d_RECURSES)
#error Recursive header files inclusion detected in StandardDSL3d.h
#else // defined(StandardDSL3d_RECURSES)
/** Prevents recursive inclusion of headers. */
#define StandardDSL3d_RECURSES

#if !defined StandardDSL3d_h
/** Prevents repeated inclusion of headers. */
#define StandardDSL3d_h

#include "DGtal/geometry/curves/ArithmeticalDSL.h"
#include "raytracer/Ray.h"

namespace DGtal {
  namespace rt {

    struct StandardDSL3d {
      typedef Integer                                      Coordinate;
      typedef DGtal::int64_t                               LongInteger;
      typedef ArithmeticalDSL<Coordinate, LongInteger, 4>  StandardDSL;
      typedef StandardDSL::Point                           Point2i;
      typedef DGtal::rt::Point3                            Point3;
      typedef DGtal::rt::Point3i                           Point3i;

      /**
       * \brief Aim: This class aims at representing an iterator
       * that provides a way to scan the points of a 3d DSL.
       * It is both a model of readable iterator and of
       * bidirectional iterator.
       *
       * Instances of this class are returned by begin and end methods.
       */
      class ConstIterator
        : public boost::iterator_facade< ConstIterator, //derived type
                                         Point3i const,   //value type
                                         boost::bidirectional_traversal_tag,
                                         Point3i const,   //reference type
                                         Coordinate >     //difference type
      {

        // ------------------------- Private data -----------------------
      private:
        
        StandardDSL::ConstIterator myXY;
        StandardDSL::ConstIterator myXZ;
        StandardDSL::ConstIterator myYZ;
        
        /// The current point
        Point3i  myP;

        // ------------------------- Standard services -----------------------
      public:
        
        /**
         * Default constructor (not valid).
         */
        ConstIterator() {}

        /**
         * Constructor.
         * @param D an arithmetical 3D DSL
         * @param p a point of the DSL containing @a aDSL
         */
        ConstIterator( const StandardDSL3d* D, const Point3i& p )
          : myP( p )
        {
          myXY = D->xy.begin( Point2i( p[ 0 ], p[ 1 ] ) );
          myXZ = D->xz.begin( Point2i( p[ 0 ], p[ 2 ] ) );
          myYZ = D->yz.begin( Point2i( p[ 1 ], p[ 2 ] ) );
        }
        
        /**
         * Copy constructor.
         * @param aOther the iterator to clone.
         */
        ConstIterator( const ConstIterator & other )
          : myXY( other.myXY ), myXZ( other.myXZ ), myYZ( other.myYZ ),
            myP( other.myP )
        {}
        
        /**
         * Assignment.
         * @param other the iterator to copy.
         * @return a reference on 'this'.
         */
        ConstIterator& operator= ( const ConstIterator & other )
        {
          if ( this != &other )
            {
              myXY = other.myXY;
              myXZ = other.myXZ;
              myYZ = other.myYZ;
              myP  = other.myP;
            }
          return *this;
        }
        
        /**
         * Destructor. Does nothing.
         */
        ~ConstIterator() {}

        // ------------------------- iteration services -------------------------
      private:
        friend class boost::iterator_core_access;
        
        /**
         * Dereference operator
         * @return the current point
         */
        Point3i const dereference() const
        { return myP; }
        
        /**
         * Moves @a myCurrentPoint to the next point of the DSL
         */
        void increment()
        {
          myXY++;
          bool xy_x = (*myXY)[ 0 ] != myP[ 0 ];
          myXZ++;
          bool xz_x = (*myXZ)[ 0 ] != myP[ 0 ];
          myYZ++;
          bool yz_y = (*myYZ)[ 0 ] != myP[ 1 ];
          if ( xy_x && xz_x )
            { // displacement along x
              myYZ--;
              myP[ 0 ] = (*myXY)[ 0 ];
            }
          else if ( ( ! xy_x ) && yz_y )
            { // displacement along y
              myXZ--;
              myP[ 1 ] = (*myYZ)[ 0 ];
            }
          else
            { // displacement along z
              myXY--;
              myP[ 2 ] = (*myXZ)[ 1 ];
            }
        }
        
        /**
         * Moves @a myCurrentPoint to the previous point of the DSL
         */
        void decrement()
        {
          myXY--;
          bool xy_x = (*myXY)[ 0 ] != myP[ 0 ];
          myXZ--;
          bool xz_x = (*myXZ)[ 0 ] != myP[ 0 ];
          myYZ--;
          bool yz_y = (*myYZ)[ 0 ] != myP[ 1 ];
          if ( xy_x && xz_x )
            { // displacement along x
              myYZ++;
              myP[ 0 ] = (*myXY)[ 0 ];
            }
          else if ( ( ! xy_x ) && yz_y )
            { // displacement along y
              myXZ++;
              myP[ 1 ] = (*myYZ)[ 0 ];
            }
          else
            { // displacement along z
              myXY++;
              myP[ 2 ] = (*myXZ)[ 1 ];
            }
        }
        
        /**
         * Equality operator.
         *
         * @param aOther the iterator to compare with
         * (must be defined on the same DSL).
         *
         * @return 'true' if their current points coincide.
         */
        bool equal(const ConstIterator& aOther) const
        {
          return myP == aOther.myP;
        }
                
      };
      
      /**
       * Type of const reverse iterator, adapted from ConstIterator.
       */
      typedef DGtal::ReverseIterator<ConstIterator> ConstReverseIterator;
      
      /// Constructor from ray and digitization precision.
      StandardDSL3d( const Ray& aRay, Coordinate precision )
        : xy(0,0,0), xz(0,0,0), yz(0,0,0),
          xyValid( false ), xzValid( false ), yzValid( false )
      {
        Coordinate ux = (Coordinate) round( aRay.direction[ 0 ] * (Real) precision );
        Coordinate uy = (Coordinate) round( aRay.direction[ 1 ] * (Real) precision );
        Coordinate uz = (Coordinate) round( aRay.direction[ 2 ] * (Real) precision );
        Coordinate gxy = IntegerComputer<Coordinate>::staticGcd( ux, uy );
        Coordinate gxz = IntegerComputer<Coordinate>::staticGcd( ux, uz );
        Coordinate gyz = IntegerComputer<Coordinate>::staticGcd( uy, uz );
        // trace.info() << "or=" << aRay.origin << std::endl;
        // trace.info() << "di=" << aRay.direction << std::endl;
        // trace.info() << "ux=" << ux << " uy=" << uy << " uz=" << uz << std::endl;
        // trace.info() << "gxy=" << gxy << " gxz=" << gxz << " gyz=" << gyz << std::endl;
        // trace.info() << "po=" << getPoint( aRay.origin ) << std::endl;
        Point3 q = aRay.origin;
        if ( gxy != 0 )
          {
            Coordinate a  = uy / gxy;
            Coordinate b  = ux / gxy;
            Coordinate mu = (Coordinate) round( a * q[ 0 ] - b * q[ 1 ]
                                                - 0.5*( std::abs(a)+std::abs(b)-1 ) );
            xy = StandardDSL( a, b, mu );
            xyValid = true;
          }
        if ( gxz != 0 )
          {
            Coordinate a  = uz / gxz;
            Coordinate b  = ux / gxz;
            Coordinate mu = (Coordinate) round( a * q[ 0 ] - b * q[ 2 ]
                                                - 0.5*( std::abs(a)+std::abs(b)-1 ) );
            xz = StandardDSL( a, b, mu );
            xzValid = true;
          }
        if ( gyz != 0 )
          {
            Coordinate a  = uz / gyz;
            Coordinate b  = uy / gyz;
            Coordinate mu = (Coordinate) round( a * q[ 1 ] - b * q[ 2 ]
                                                - 0.5*( std::abs(a)+std::abs(b)-1 ) );
            yz = StandardDSL( a, b, mu );
            yzValid = true;
          }
        // trace.info() << "xy=" << xy << std::endl;
        // trace.info() << "xz=" << xz << std::endl;
        // trace.info() << "yz=" << yz << std::endl;
        if ( ! isInDSL( getPoint( aRay.origin ) ) )
          trace.warning() << "Ray origin is not in DSL 3d." << std::endl;
      }

      Point3i getPoint( Point3 p ) const
      {
        return Point3i( (Coordinate) round( p[ 0 ] ),
                        (Coordinate) round( p[ 1 ] ),
                        (Coordinate) round( p[ 2 ] ) );
      }

      ConstIterator begin( Point3i p ) const
      { return ConstIterator( this, p ); }

      ConstIterator end( Point3i p ) const
      { return ConstIterator( this, p ); }

      bool isInDSL(const Point3i& p) const
      {
        return xy.isInDSL( Point2i( p[ 0 ], p[ 1 ] ) )
          &&   xz.isInDSL( Point2i( p[ 0 ], p[ 2 ] ) )
          &&   yz.isInDSL( Point2i( p[ 1 ], p[ 2 ] ) );
      }

      static Point2i projectXY( const Point3i& p )
      {
        return Point2i( p[ 0 ], p[ 1 ] );
      }
      static Point2i projectXZ( const Point3i& p )
      {
        return Point2i( p[ 0 ], p[ 2 ] );
      }
      static Point2i projectYZ( const Point3i& p )
      {
        return Point2i( p[ 1 ], p[ 2 ] );
      }
      bool before( const Point3i& p1, const Point3i& p2 ) const
      {
        return ( p1 != p2 ) && beforeOrEqual( p1, p2 );
      }
      bool beforeOrEqual( const Point3i& p1, const Point3i& p2 ) const
      {
        return ( !xyValid || ( xy.beforeOrEqual( projectXY( p1 ), projectXY( p2 ) ) ) )
          &&   ( !xzValid || ( xz.beforeOrEqual( projectXZ( p1 ), projectXZ( p2 ) ) ) )
          &&   ( !yzValid || ( yz.beforeOrEqual( projectYZ( p1 ), projectYZ( p2 ) ) ) );
      }
      
      /// The three projections of the 3D standard digital straight line
      StandardDSL xy, yz, xz;
      bool xyValid, yzValid, xzValid;
    };
  } //  namespace rt {
} // namespace DGtal {

///////////////////////////////////////////////////////////////////////////////
#endif // !defined StandardDSL3d_h

#undef StandardDSL3d_RECURSES
#endif // else defined(StandardDSL3d_RECURSES)