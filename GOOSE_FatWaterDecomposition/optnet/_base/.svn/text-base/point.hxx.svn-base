/*
 ==========================================================================
 |   
 |   $Id: point.hxx 2137 2007-07-02 03:26:31Z kangli $
 |
 |   Written by Kang Li <kangl@cmu.edu>
 |   Department of Electrical and Computer Engineering
 |   Carnegie Mellon University
 |   
 ==========================================================================
 |   This file is a part of the OptimalNet library.
 ==========================================================================
 | Copyright (c) 2003-2007 Kang Li <kangl@cmu.edu>. All rights reserved.
 | 
 | Version: MPL 1.1/GPL 2.0/LGPL 2.1
 | 
 | The contents of  this file are  subject to the  Mozilla Public License
 | Version  1.1 (the  "License"); you  may not  use this  file except  in
 | compliance with the License. You may obtain a copy of the License at
 | 
 | http://www.mozilla.org/MPL/
 | 
 | Software distributed under  the License is  distributed on an  "AS IS"
 | basis, WITHOUT WARRANTY  OF ANY KIND,  either express or  implied. See
 | the License for the specific language governing rights and limitations
 | under the License.
 | 
 | The Original Code is OptimalNet (optnet) Library code.
 | 
 | The  Initial  Developer of  the  Original Code  is  Kang Li.  Portions
 | created  by  the Initial  Developer  are Copyright  (C)  2003-2007 the
 | Initial Developer. All Rights Reserved.
 | 
 | Contributor(s): None
 | 
 | Alternatively, the contents of this  file may be used under  the terms
 | of either of the  GNU General Public License  Version 2 or later  (the
 | "GPL"), or the GNU Lesser General Public License Version 2.1 or  later
 | (the "LGPL"), in which case the provisions of the GPL or the LGPL  are
 | applicable instead of those  above. If you wish  to allow use of  your
 | version of this  file only under  the terms of  either the GPL  or the
 | LGPL, and not to allow others  to use your version of this  file under
 | the  terms  of  the  MPL,  indicate  your  decision  by  deleting  the
 | provisions above and replace them with the notice and other provisions
 | required by the GPL or the  LGPL. If you do not delete  the provisions
 | above, a recipient may use your  version of this file under the  terms
 | of any one of the MPL, the GPL or the LGPL.
 ==========================================================================
 */

#ifndef ___POINT_HXX___
#   define ___POINT_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <cmath>
#   include <cassert>
#   include <cstddef>
#   ifdef min       //
#       undef min   //
#   endif           // These would interfere with
#   ifdef max       //  std::numerical_limits<_T>::max
#       undef max   //  std::numerical_limits<_T>::min
#   endif           //
#   include <limits>

/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class point
///  @brief Multidimensional point type.
///////////////////////////////////////////////////////////////////////////
template <typename _T, int _Dim>
struct point
{
    typedef _T          value_type;
    typedef _T&         reference;
    typedef const _T&   const_reference;
    typedef _T*         pointer;
    typedef const _T*   const_pointer;
    typedef size_t      size_type;

    point()
    {
        clear();
    }

#   ifdef __OPTNET_CRAPPY_MSC__
    template<typename _Point>
    point(const _Point& rhs)
    {
        assert(_Dim == rhs.dim());
        for (int i = 0; i < _Dim; ++i)
            v[i] = static_cast<value_type>(rhs.v[i]);
    }
#   else // normal compilers
    point(const point& rhs)
    {
        for (int i = 0; i < _Dim; ++i)
            v[i] = rhs.v[i];
    }

    template<typename _T2>
    point(const point<_T2, _Dim>& rhs)
    {
        for (int i = 0; i < _Dim; ++i)
            v[i] = static_cast<value_type>(rhs.v[i]);
    }
#   endif

    //
    // operator =
    //
    ///////////////////////////////////////////////////////////////////////
#   ifdef __OPTNET_CRAPPY_MSC__
    template<typename _Point>
    inline point& operator=(const _Point& rhs)
    {
        assert(_Dim == rhs.dim());
        for (int i = 0; i < _Dim; ++i)
            v[i] = static_cast<value_type>(rhs.v[i]);
        return *this;
    }
#   else // normal compilers
    inline point& operator=(const point& rhs)
    {
        for (int i = 0; i < _Dim; ++i)
            v[i] = rhs.v[i];
        return *this;
    }

    template<typename _T2>
    inline point& operator=(const point<_T2, _Dim>& rhs)
    {
        for (int i = 0; i < _Dim; ++i)
            v[i] = static_cast<value_type>(rhs.v[i]);
        return *this;
    }
#   endif

    ///////////////////////////////////////////////////////////////////////
    ///  Point component access.
    ///
    ///  @param i The index of the component.
    ///
    ///  @return A constant reference to the i-th component.
    ///////////////////////////////////////////////////////////////////////
    inline const_reference operator[](const size_type& i) const
    {
        return v[i];
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Point component access.
    ///
    ///  @param i The index of the component.
    ///
    ///  @return A reference to the i-th component.
    ///////////////////////////////////////////////////////////////////////
    inline reference operator[](const size_type& i)
    {
        return v[i];
    }

    //
    // operator +=
    //
    ///////////////////////////////////////////////////////////////////////
#   ifndef __OPTNET_CRAPPY_MSC__
    inline point& operator+=(const point& rhs)
    {
        for (int i = 0; i < _Dim; ++i)
            v[i] += rhs.v[i];
        return *this;
    }
#   endif

    template<typename _T2>
    inline point& operator+=(const point<_T2, _Dim>& rhs)
    {
        for (int i = 0; i < _Dim; ++i)
            v[i] = static_cast<value_type>(v[i] + rhs.v[i]);
        return *this;
    }

#   define POINT_OPERATOR_AE_PROTO(type)                    \
        inline point& operator+=(const type& rhs)           \
        {                                                   \
            for (int i = 0; i < _Dim; ++i)                  \
                v[i] = static_cast<value_type>(v[i] + rhs); \
            return *this;                                   \
        }
    POINT_OPERATOR_AE_PROTO(float           )
    POINT_OPERATOR_AE_PROTO(double          )
    POINT_OPERATOR_AE_PROTO(long double     )
    POINT_OPERATOR_AE_PROTO(char            )
    POINT_OPERATOR_AE_PROTO(int             )
    POINT_OPERATOR_AE_PROTO(long            )
    POINT_OPERATOR_AE_PROTO(short           )
    POINT_OPERATOR_AE_PROTO(unsigned char   )
    POINT_OPERATOR_AE_PROTO(unsigned int    )
    POINT_OPERATOR_AE_PROTO(unsigned long   )
    POINT_OPERATOR_AE_PROTO(unsigned short  )

    //
    // operator -=
    //
    ///////////////////////////////////////////////////////////////////////
#   ifndef __OPTNET_CRAPPY_MSC__
    inline point& operator-=(const point& rhs)
    {
        for (int i = 0; i < _Dim; ++i)
            v[i] -= rhs.v[i];
        return *this;
    }
#   endif

    template<typename _T2>
    inline point& operator-=(const point<_T2, _Dim>& rhs)
    {
        for (int i = 0; i < _Dim; ++i)
            v[i] = static_cast<value_type>(v[i] - rhs.v[i]);
        return *this;
    }
    
#   define POINT_OPERATOR_SE_PROTO(type)                    \
        inline point& operator-=(const type& rhs)           \
        {                                                   \
            for (int i = 0; i < _Dim; ++i)                  \
                v[i] = static_cast<value_type>(v[i] - rhs); \
            return *this;                                   \
        }
    POINT_OPERATOR_SE_PROTO(float           )
    POINT_OPERATOR_SE_PROTO(double          )
    POINT_OPERATOR_SE_PROTO(long double     )
    POINT_OPERATOR_SE_PROTO(char            )
    POINT_OPERATOR_SE_PROTO(int             )
    POINT_OPERATOR_SE_PROTO(long            )
    POINT_OPERATOR_SE_PROTO(short           )
    POINT_OPERATOR_SE_PROTO(unsigned char   )
    POINT_OPERATOR_SE_PROTO(unsigned int    )
    POINT_OPERATOR_SE_PROTO(unsigned long   )
    POINT_OPERATOR_SE_PROTO(unsigned short  )

    //
    // operator *=
    //
    ///////////////////////////////////////////////////////////////////////
#   ifndef __OPTNET_CRAPPY_MSC__
    inline point& operator*=(const point& rhs)
    {
        for (int i = 0; i < _Dim; ++i)
            v[i] *= rhs.v[i];
        return *this;
    }
#   endif
    
    template<typename _T2>
    inline point& operator*=(const point<_T2, _Dim>& rhs)
    {
        for (int i = 0; i < _Dim; ++i)
            v[i] = static_cast<value_type>(v[i] * rhs.v[i]);
        return *this;
    }
    
#   define POINT_OPERATOR_ME_PROTO(type)                    \
        inline point& operator*=(const type& rhs)           \
        {                                                   \
            for (int i = 0; i < _Dim; ++i)                  \
                v[i] = static_cast<value_type>(v[i] * rhs); \
            return *this;                                   \
        }
    POINT_OPERATOR_ME_PROTO(float           )
    POINT_OPERATOR_ME_PROTO(double          )
    POINT_OPERATOR_ME_PROTO(long double     )
    POINT_OPERATOR_ME_PROTO(char            )
    POINT_OPERATOR_ME_PROTO(int             )
    POINT_OPERATOR_ME_PROTO(long            )
    POINT_OPERATOR_ME_PROTO(short           )
    POINT_OPERATOR_ME_PROTO(unsigned char   )
    POINT_OPERATOR_ME_PROTO(unsigned int    )
    POINT_OPERATOR_ME_PROTO(unsigned long   )
    POINT_OPERATOR_ME_PROTO(unsigned short  )
    
    //
    // operator /=
    //
    ///////////////////////////////////////////////////////////////////////
#   ifndef __OPTNET_CRAPPY_MSC__
    inline point& operator/=(const point& rhs)
    {
        for (int i = 0; i < _Dim; ++i)
            v[i] /= rhs.v[i];
        return *this;
    }
#   endif

    template<typename _T2>
    inline point& operator/=(const point<_T2, _Dim>& rhs)
    {
        for (int i = 0; i < _Dim; ++i)
            v[i] = static_cast<value_type>(v[i] / rhs.v[i]);
        return *this;
    }
    
#   define POINT_OPERATOR_DE_PROTO(type)                    \
        inline point& operator/=(const type& rhs)           \
        {                                                   \
            for (int i = 0; i < _Dim; ++i)                  \
                v[i] = static_cast<value_type>(v[i] / rhs); \
            return *this;                                   \
        }
    POINT_OPERATOR_DE_PROTO(float           )
    POINT_OPERATOR_DE_PROTO(double          )
    POINT_OPERATOR_DE_PROTO(long double     )
    POINT_OPERATOR_DE_PROTO(char            )
    POINT_OPERATOR_DE_PROTO(int             )
    POINT_OPERATOR_DE_PROTO(long            )
    POINT_OPERATOR_DE_PROTO(short           )
    POINT_OPERATOR_DE_PROTO(unsigned char   )
    POINT_OPERATOR_DE_PROTO(unsigned int    )
    POINT_OPERATOR_DE_PROTO(unsigned long   )
    POINT_OPERATOR_DE_PROTO(unsigned short  )

    //
    // operator <
    //
    ///////////////////////////////////////////////////////////////////////
    ///  Compares two points.
    ///
    ///  @return Returns true if the point is less than the point on the
    ///          right-hand side, false otherwise.
    ///////////////////////////////////////////////////////////////////////
    inline bool operator< (const point& rhs) const
    {
        return memcmp(v, rhs.v, sizeof(value_type) * _Dim) < 0;
    }

    //
    // operator ==
    //
#   ifndef __OPTNET_CRAPPY_MSC__
    ///////////////////////////////////////////////////////////////////////
    ///  Compares two points.
    ///
    ///  @return Returns true if the two points are equal.
    ///////////////////////////////////////////////////////////////////////
    inline bool operator==(const point& rhs) const
    {
        for (int i = 0; i < _Dim; ++i)
            if (v[i] != rhs.v[i]) return false;
        return true;
    }
#   endif

    ///////////////////////////////////////////////////////////////////////
    ///  Compares two points.
    ///
    ///  @return Returns true if the two points are equal.
    ///////////////////////////////////////////////////////////////////////
    template <typename _T2>
    inline bool operator==(const point<_T2, _Dim>& rhs) const
    {
        for (int i = 0; i < _Dim; ++i)
            if (v[i] != static_cast<value_type>(rhs.v[i])) return false;
        return true;
    }

    //
    // operator !=
    //
#   ifndef __OPTNET_CRAPPY_MSC__
    ///////////////////////////////////////////////////////////////////////
    ///  Compares two points.
    ///
    ///  @return Returns true if the two points are not equal.
    ///////////////////////////////////////////////////////////////////////
    inline bool operator!=(const point& rhs) const
    {
        return !(*this == rhs);
    }
#   endif

    ///////////////////////////////////////////////////////////////////////
    ///  Compares two points.
    ///
    ///  @return Returns true if the two points are not equal.
    ///////////////////////////////////////////////////////////////////////
    template <typename _T2>
    inline bool operator!=(const point<_T2, _Dim>& rhs) const
    {
        return !(*this == rhs);
    }

    //
    // operator -
    //
    ///////////////////////////////////////////////////////////////////////
    ///  Negate the coordinates of the point.
    ///
    ///  @return The negated point.
    ///////////////////////////////////////////////////////////////////////
    inline point operator-() const
    {
        point pt;
        for (int i = 0; i < _Dim; ++i)
            pt.v[i] = -v[i];
        return pt;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Computes the dot-product of two points (vectors).
    ///
    ///  @return The computed dot-product.
    ///////////////////////////////////////////////////////////////////////
#   ifndef __OPTNET_CRAPPY_MSC__
    inline value_type dot(const point& rhs)
    {
        value_type ans = value_type();
        for (int i = 0; i < _Dim; ++i)
            ans += v[i] * rhs.v[i];
        return ans;
    }
#   endif

    ///////////////////////////////////////////////////////////////////////
    ///  Computes the dot-product of two points (vectors).
    ///
    ///  @return The computed dot-product.
    ///////////////////////////////////////////////////////////////////////
    template<typename _T2>
    inline value_type dot(const point<_T2, _Dim>& rhs)
    {
        value_type ans = value_type();
        for (int i = 0; i < _Dim; ++i)
            ans += v[i] * static_cast<value_type>(rhs.v[i]);
        return ans;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the squared Euclidean distance between the point and
    ///  the origin.
    ///
    ///  @return The computed squared length.
    ///////////////////////////////////////////////////////////////////////
    inline value_type length2() const
    {
        value_type ans = value_type();
        for (int i = 0; i < _Dim; ++i)
            ans += v[i] * v[i];
        return ans;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the Euclidean distance between the point and the origin.
    ///
    ///  @return The computed length.
    ///////////////////////////////////////////////////////////////////////
    inline value_type length() const
    {
        value_type ans = length2();
        return static_cast<value_type>(sqrt((double)(ans)));
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Normalize the point such that its distance to the origin is one.
    ///////////////////////////////////////////////////////////////////////
    inline void normalize()
    {
        int i;
        value_type len = length();
        
        if (len > std::numeric_limits<value_type>::epsilon()) {
            double invlen = 1.0 / (double)len;
            for (i = 0; i < _Dim; ++i)
                v[i] = static_cast<value_type>(v[i] * invlen);
        }
        else {
            for (i = 0; i < _Dim; ++i)
                v[i] = static_cast<value_type>(0);
        }
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Sets the point to zero.
    ///////////////////////////////////////////////////////////////////////
    inline void clear()
    {
        for (value_type* p = v; p != v + _Dim; ++p)
            *p = value_type(); // Clear value.
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the dimension of the point.
    ///
    ///  @return The dimension of the point. Same as the _Dim template
    ///          parameter.
    ///////////////////////////////////////////////////////////////////////
    inline size_type dim() const
    {
        return size_type(_Dim);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Point components.
    ///
    ///  @remarks Note that point[x] = point.v[x].
    ///////////////////////////////////////////////////////////////////////
    value_type v[_Dim];
};



//
// operator +
//
///////////////////////////////////////////////////////////////////////////
template <typename _T, int _Dim>
inline point<_T, _Dim>
    operator+ (const point<_T, _Dim>& lhs, const point<_T, _Dim>& rhs)
{
    point<_T, _Dim> obj(lhs); obj += rhs; return obj;
}

#   ifndef __OPTNET_CRAPPY_MSC__
template <typename _T, typename _T2, int _Dim>
inline point<_T, _Dim>
    operator+ (const point<_T, _Dim>& lhs, const point<_T2, _Dim>& rhs)
{
    point<_T, _Dim> obj(lhs); obj += rhs; return obj;
}
#   endif

///////////////////////////////////////////////////////////////////////////
#   define POINT_OPERATOR_AR_PROTO(type)                            \
        template <typename _T, int _Dim>                            \
        inline point<_T, _Dim>                                      \
            operator+ (const point<_T, _Dim>& lhs, const type& rhs) \
        {                                                           \
            point<_T, _Dim> obj(lhs); obj += rhs; return obj;       \
        }
    POINT_OPERATOR_AR_PROTO(float           )
    POINT_OPERATOR_AR_PROTO(double          )
    POINT_OPERATOR_AR_PROTO(long double     )
    POINT_OPERATOR_AR_PROTO(char            )
    POINT_OPERATOR_AR_PROTO(int             )
    POINT_OPERATOR_AR_PROTO(long            )
    POINT_OPERATOR_AR_PROTO(short           )
    POINT_OPERATOR_AR_PROTO(unsigned char   )
    POINT_OPERATOR_AR_PROTO(unsigned int    )
    POINT_OPERATOR_AR_PROTO(unsigned long   )
    POINT_OPERATOR_AR_PROTO(unsigned short  )

///////////////////////////////////////////////////////////////////////////
#   define POINT_OPERATOR_AL_PROTO(type)                            \
        template <typename _T, int _Dim>                            \
        inline point<_T, _Dim>                                      \
            operator+ (const type& lhs, const point<_T, _Dim>& rhs) \
        {                                                           \
            point<_T, _Dim> obj(rhs); obj += lhs; return obj;       \
        }
    POINT_OPERATOR_AL_PROTO(float           )
    POINT_OPERATOR_AL_PROTO(double          )
    POINT_OPERATOR_AL_PROTO(long double     )
    POINT_OPERATOR_AL_PROTO(char            )
    POINT_OPERATOR_AL_PROTO(int             )
    POINT_OPERATOR_AL_PROTO(long            )
    POINT_OPERATOR_AL_PROTO(short           )
    POINT_OPERATOR_AL_PROTO(unsigned char   )
    POINT_OPERATOR_AL_PROTO(unsigned int    )
    POINT_OPERATOR_AL_PROTO(unsigned long   )
    POINT_OPERATOR_AL_PROTO(unsigned short  )

//
// operator -
//
///////////////////////////////////////////////////////////////////////////

template <typename _T, int _Dim>
inline point<_T, _Dim>
    operator- (const point<_T, _Dim>& lhs, const point<_T, _Dim>& rhs)
{
    point<_T, _Dim> obj(lhs); obj -= rhs; return obj;
}

#   ifndef __OPTNET_CRAPPY_MSC__
template <typename _T, typename _T2, int _Dim>
inline point<_T, _Dim>
    operator- (const point<_T, _Dim>& lhs, const point<_T2, _Dim>& rhs)
{
    point<_T, _Dim> obj(lhs); obj -= rhs; return obj;
}
#   endif //__OPTNET_CRAPPY_MSC__

///////////////////////////////////////////////////////////////////////////
#   define POINT_OPERATOR_SR_PROTO(type)                            \
        template <typename _T, int _Dim>                            \
        inline point<_T, _Dim>                                      \
            operator- (const point<_T, _Dim>& lhs, const type& rhs) \
        {                                                           \
            point<_T, _Dim> obj(lhs); obj -= rhs; return obj;       \
        }
    POINT_OPERATOR_SR_PROTO(float           )
    POINT_OPERATOR_SR_PROTO(double          )
    POINT_OPERATOR_SR_PROTO(long double     )
    POINT_OPERATOR_SR_PROTO(char            )
    POINT_OPERATOR_SR_PROTO(int             )
    POINT_OPERATOR_SR_PROTO(long            )
    POINT_OPERATOR_SR_PROTO(short           )
    POINT_OPERATOR_SR_PROTO(unsigned char   )
    POINT_OPERATOR_SR_PROTO(unsigned int    )
    POINT_OPERATOR_SR_PROTO(unsigned long   )
    POINT_OPERATOR_SR_PROTO(unsigned short  )

//
// operator *
//
///////////////////////////////////////////////////////////////////////////

template <typename _T, int _Dim>
inline point<_T, _Dim>
    operator* (const point<_T, _Dim>& lhs, const point<_T, _Dim>& rhs)
{
    point<_T, _Dim> obj(lhs); obj *= rhs; return obj;
}

#   ifndef __OPTNET_CRAPPY_MSC__
template <typename _T, typename _T2, int _Dim>
inline point<_T, _Dim>
    operator* (const point<_T, _Dim>& lhs, const point<_T2, _Dim>& rhs)
{
    point<_T, _Dim> obj(lhs); obj *= rhs; return obj;
}
#   endif //__OPTNET_CRAPPY_MSC__

///////////////////////////////////////////////////////////////////////////
#   define POINT_OPERATOR_MR_PROTO(type)                            \
        template <typename _T, int _Dim>                            \
        inline point<_T, _Dim>                                      \
            operator* (const point<_T, _Dim>& lhs, const type& rhs) \
        {                                                           \
            point<_T, _Dim> obj(lhs); obj *= rhs; return obj;       \
        }
    POINT_OPERATOR_MR_PROTO(float           )
    POINT_OPERATOR_MR_PROTO(double          )
    POINT_OPERATOR_MR_PROTO(long double     )
    POINT_OPERATOR_MR_PROTO(char            )
    POINT_OPERATOR_MR_PROTO(int             )
    POINT_OPERATOR_MR_PROTO(long            )
    POINT_OPERATOR_MR_PROTO(short           )
    POINT_OPERATOR_MR_PROTO(unsigned char   )
    POINT_OPERATOR_MR_PROTO(unsigned int    )
    POINT_OPERATOR_MR_PROTO(unsigned long   )
    POINT_OPERATOR_MR_PROTO(unsigned short  )

///////////////////////////////////////////////////////////////////////////
#   define POINT_OPERATOR_ML_PROTO(type)                            \
        template <typename _T, int _Dim>                            \
        inline point<_T, _Dim>                                      \
            operator* (const type& lhs, const point<_T, _Dim>& rhs) \
        {                                                           \
            point<_T, _Dim> obj(rhs); obj *= lhs; return obj;       \
        }
    POINT_OPERATOR_ML_PROTO(float           )
    POINT_OPERATOR_ML_PROTO(double          )
    POINT_OPERATOR_ML_PROTO(long double     )
    POINT_OPERATOR_ML_PROTO(char            )
    POINT_OPERATOR_ML_PROTO(int             )
    POINT_OPERATOR_ML_PROTO(long            )
    POINT_OPERATOR_ML_PROTO(short           )
    POINT_OPERATOR_ML_PROTO(unsigned char   )
    POINT_OPERATOR_ML_PROTO(unsigned int    )
    POINT_OPERATOR_ML_PROTO(unsigned long   )
    POINT_OPERATOR_ML_PROTO(unsigned short  )

//
// operator /
//
///////////////////////////////////////////////////////////////////////////
template <typename _T, int _Dim>
inline point<_T, _Dim>
    operator/ (const point<_T, _Dim>& lhs, const point<_T, _Dim>& rhs)
{
    point<_T, _Dim> obj(lhs); obj /= rhs; return obj;
}

#   ifndef __OPTNET_CRAPPY_MSC__
template <typename _T, typename _T2, int _Dim>
inline point<_T, _Dim>
    operator/ (const point<_T, _Dim>& lhs, const point<_T2, _Dim>& rhs)
{
    point<_T, _Dim> obj(lhs); obj /= rhs; return obj;
}
#   endif //__OPTNET_CRAPPY_MSC__

///////////////////////////////////////////////////////////////////////////
#   define POINT_OPERATOR_DR_PROTO(type)                            \
        template <typename _T, int _Dim>                            \
        inline point<_T, _Dim>                                      \
            operator/ (const point<_T, _Dim>& lhs, const type& rhs) \
        {                                                           \
            point<_T, _Dim> obj(lhs); obj /= rhs; return obj;       \
        }
    POINT_OPERATOR_DR_PROTO(float           )
    POINT_OPERATOR_DR_PROTO(double          )
    POINT_OPERATOR_DR_PROTO(long double     )
    POINT_OPERATOR_DR_PROTO(char            )
    POINT_OPERATOR_DR_PROTO(int             )
    POINT_OPERATOR_DR_PROTO(long            )
    POINT_OPERATOR_DR_PROTO(short           )
    POINT_OPERATOR_DR_PROTO(unsigned char   )
    POINT_OPERATOR_DR_PROTO(unsigned int    )
    POINT_OPERATOR_DR_PROTO(unsigned long   )
    POINT_OPERATOR_DR_PROTO(unsigned short  )

} // namespace

#endif // ___POINT_HXX___
