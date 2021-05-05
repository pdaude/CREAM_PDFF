/*
 ==========================================================================
 |   
 |   $Id: array2_base.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___ARRAY2_BASE_HXX___
#   define ___ARRAY2_BASE_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   include <assert.h>
#   include <memory.h>
#   include <optnet/_base/memory.hxx>
#   include <optnet/_base/array_traits.hxx>
#   include <optnet/_base/point.hxx>

#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200)
#       pragma warning(disable: 4018)
#       pragma warning(disable: 4146)
#   endif

#   include <algorithm>
#   include <vector>

///  @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class array2_base array2_base.hxx "_base/array2_base.hxx"
///  @brief Base array template class.
///////////////////////////////////////////////////////////////////////////
template <typename _Ty>
class array2_base
{
    typedef array2_base<_Ty>            _Self;
    typedef array_traits<_Ty, _Self>    _Traits;

public:
    
    /// @brief See array_traits::value_type
    typedef typename _Traits::value_type        value_type;

    /// @brief See array_traits::reference
    typedef typename _Traits::reference         reference;

    /// @brief See array_traits::const_reference
    typedef typename _Traits::const_reference   const_reference;

    /// @brief See array_traits::iterator
    typedef typename _Traits::iterator          iterator;
    
    /// @brief See array_traits::const_iterator
    typedef typename _Traits::const_iterator    const_iterator;

    /// @brief See array_traits::pointer
    typedef typename _Traits::pointer           pointer;

    /// @brief See array_traits::const_pointer
    typedef typename _Traits::const_pointer     const_pointer;

    /// @brief See array_traits::difference_type
    typedef typename _Traits::difference_type   difference_type;

    /// @brief See array_traits::size_type
    typedef typename _Traits::size_type         size_type;

    // constructor/destructors
    array2_base() :
        m_p(0), m_sz(0)
    {
        m_asz[0] = 0;
        m_asz[1] = 0;
    }
    
    // queries

    ///////////////////////////////////////////////////////////////////////
    ///  Returns a pointer to the first element in the array.
    ///
    ///  @return A pointer to the first element in the array.
    ///////////////////////////////////////////////////////////////////////
    inline pointer          data()          { return m_p; }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns a const_pointer to the first element in the array.
    ///
    ///  @return A const_pointer to the first element in the array.
    ///////////////////////////////////////////////////////////////////////
    inline const_pointer    data() const    { return m_p; }

    // iterators

    ///////////////////////////////////////////////////////////////////////
    ///  Returns an iterator to the first element in the array.
    ///
    ///  @return An iterator to the first element in the array.
    ///////////////////////////////////////////////////////////////////////
    inline iterator         begin()         { return iterator(m_p); }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns a const_iterator to the first element in the array.
    ///
    ///  @return A const_iterator to the first element in the array.
    ///////////////////////////////////////////////////////////////////////
    inline const_iterator   begin()  const  { return const_iterator(m_p); }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns an iterator pointing just beyond the end of the array.
    ///
    ///  @return An iterator pointing just beyond the end of the array.
    ///////////////////////////////////////////////////////////////////////
    inline iterator         end()           { return iterator(m_p + m_sz); }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns a const_iterator pointing just beyond the end of the array.
    ///
    ///  @return A const_iterator pointing just beyond the end of the array.
    ///////////////////////////////////////////////////////////////////////
    inline const_iterator   end()    const  { return const_iterator(m_p + m_sz); }

    ///////////////////////////////////////////////////////////////////////
    ///  Tests if the array container is empty. (New in version 0.9.16)
    ///
    ///  @return True if the array is empty; false if it is not empty.
    ///////////////////////////////////////////////////////////////////////
    inline bool             empty()  const  { return (0 == m_sz); }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the total number of elements in the array.
    ///////////////////////////////////////////////////////////////////////
    inline size_type        size()   const  { return m_sz; }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns a pointer to the vector that contains the sizes of all
    ///  the 5 dimensions of the array.
    ///////////////////////////////////////////////////////////////////////
    inline const size_type* sizes()  const  { return m_asz; }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the size of the first  dimension of the array.
    ///////////////////////////////////////////////////////////////////////
    inline size_type        size_0() const  { return m_asz[0]; }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the size of the second dimension of the array.
    ///////////////////////////////////////////////////////////////////////
    inline size_type        size_1() const  { return m_asz[1]; }

    ///////////////////////////////////////////////////////////////////////
    ///  Similar to size_0, but returns size as an int.
    ///////////////////////////////////////////////////////////////////////
    inline int              xdim() const    { return (int)m_asz[0]; }

    ///////////////////////////////////////////////////////////////////////
    ///  Similar to size_1, but returns size as an int.
    ///////////////////////////////////////////////////////////////////////
    inline int              ydim() const    { return (int)m_asz[1]; }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the offset value of the specified element from the
    ///  begining of the array.
    ///
    ///  @param  i0 size_type  The first  index of the array element. 
    ///  @param  i1 size_type  The second index of the array element. 
    ///
    ///  @return The offset of the specified element.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline difference_type  offset(size_type i0,
                                   size_type i1
                                   ) const
    {
        return (difference_type)(i1 * m_asz[0] + i0);
    }
    
    // modifier

    ///////////////////////////////////////////////////////////////////////
    ///  Fill the array with a certain value.
    ///
    ///  @param value The value to be filled in the array.
    ///
    ///  @see clear
    ///////////////////////////////////////////////////////////////////////
    void fill(const value_type& value)
    {
        std::fill(begin(), end(), value);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Clear the content of the array.
    ///
    ///  @remarks Equivalent to fill(value_type()).
    ///
    ///  @see fill
    ///////////////////////////////////////////////////////////////////////
    void clear()
    {
        std::fill(begin(), end(), value_type());
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Gets the minimum and maximum elements of the array.
    ///
    ///  @param  pmin_val  A pointer to the variable to store the minimum
    ///                    element.
    ///  @param  pmax_val  A pointer to the variable to store the maximum
    ///                    element.
    ///
    ///  @return Returns false if the array is empty, true otherwise.
    ///
    ///////////////////////////////////////////////////////////////////////
    bool min_max_element(value_type* pmin_val,
                         value_type* pmax_val
                         ) const
    {
        if (m_sz > 0) {
            *pmin_val = *pmax_val = m_p[0];
            for (size_type i = 1; i < m_sz; ++i) {
                if (m_p[i] > *pmax_val)
                    *pmax_val = m_p[i];
                else
                if (m_p[i] < *pmin_val)
                    *pmin_val = m_p[i];
            }
            return true;
        }
        return false;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Scales the values of the array elements to a certain range.
    ///
    ///  @param  lo  The minimum value in the output range.
    ///  @param  hi  The maximum value in the output range.
    ///
    ///  @return Returns false if the array is empty, true otherwise.
    ///
    ///////////////////////////////////////////////////////////////////////
    bool scale_to_range(const value_type& lo, const value_type& hi)
    {
        value_type min_val, max_val, d1, d2;
        double     dd;

        if (!this->min_max_element(&min_val, &max_val))
            return false;

        d1 = max_val - min_val;
        d2 = hi - lo;

        if (d1 != 0 && d2 != 0) {
            dd = 1.0 / (double)d1;
            for (size_type i = 1; i < m_sz; ++i) {
                m_p[i] = (value_type)((m_p[i] - min_val) * d2 * dd) + lo;
            }
        }
        else {
            this->fill(lo);
        }

        return true;
    }

    // element access

    ///////////////////////////////////////////////////////////////////////
    ///  Returns a reference to the array element at a specified position.
    ///
    ///  @param  index  The index of the requested element.
    ///
    ///  @return A reference to the requested element.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline reference operator()(const point<size_type, 2>& index)
    {
        return m_p[index.v[1] * m_asz[0] + index.v[0]];
    }
    ///////////////////////////////////////////////////////////////////////
    ///  Returns a reference to the array element at a specified position.
    ///
    ///  @param  index  The index of the requested element.
    ///
    ///  @return A const_reference to the requested element.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline const_reference operator()(const point<size_type, 2>& index) const
    {
        return m_p[index.v[1] * m_asz[0] + index.v[0]];
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns a reference to the array element at a specified position.
    ///
    ///  @param  i0 size_type  The first  index of the array element. 
    ///  @param  i1 size_type  The second index of the array element. 
    ///
    ///  @return A reference to the requested element.
    ///
    ///  @remarks If the position specified is greater than the size of
    ///           the container, the result is undefined.
    ///////////////////////////////////////////////////////////////////////
    inline reference        operator()(size_type i0,
                                       size_type i1
                                       )
    {
        return m_p[i1 * m_asz[0] + i0];
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns a reference to the array element at a specified position.
    ///
    ///  @param  i0 size_type  The first  index of the array element. 
    ///  @param  i1 size_type  The second index of the array element. 
    ///  @param  i2 size_type  The third  index of the array element. 
    ///
    ///  @return A const_reference to the requested element.
    ///
    ///  @remarks If the position specified is greater than the size of
    ///           the container, the result is undefined.
    ///////////////////////////////////////////////////////////////////////
    inline const_reference  operator()(size_type i0,
                                       size_type i1
                                       ) const
    {
        return m_p[i1 * m_asz[0] + i0];
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns a reference to the array element at a specified index.
    ///
    ///  @param  i size_type  The index of the array element. 
    ///
    ///  @return A reference to the requested element.
    ///
    ///  @remarks If the position specified is greater than the size of
    ///           the container, the result is undefined.
    ///////////////////////////////////////////////////////////////////////
    inline reference        operator[](size_type i)
    {
        return m_p[i];
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns a reference to the array element at a specified index.
    ///
    ///  @param  i size_type  The index of the array element. 
    ///
    ///  @return A const_reference to the requested element.
    ///
    ///  @remarks If the position specified is greater than the size of
    ///           the container, the result is undefined.
    ///////////////////////////////////////////////////////////////////////
    inline const_reference  operator[](size_type i) const
    {
        return m_p[i];
    }

    // comparators

    ///////////////////////////////////////////////////////////////////////
    ///  Determine if two array2_base objects are of the same size.
    ///
    ///  @param  rhs  array2_base& The right-hand-side array2_base object.
    ///
    ///  @return Returns true if the two arrays are of the same size.
    ///
    ///////////////////////////////////////////////////////////////////////
#   ifndef __OPTNET_CRAPPY_MSC__
    inline bool is_same_size(const array2_base& rhs) const
    {
        return (m_asz[0] == rhs.m_asz[0]) &&
               (m_asz[1] == rhs.m_asz[1]);
    }
#   endif

    ///////////////////////////////////////////////////////////////////////
    ///  Determine if two array2_base objects are of the same size.
    ///
    ///  @param  rhs  array2_base& The right-hand-side array2_base object.
    ///
    ///  @return Returns true if the two arrays are of the same size.
    ///
    ///////////////////////////////////////////////////////////////////////
    template <typename _Ta>
    inline bool is_same_size(const array2_base<_Ta>& rhs) const
    {
        return (m_asz[0] == rhs.size_0()) &&
               (m_asz[1] == rhs.size_1());
    }
    
    ///////////////////////////////////////////////////////////////////////
    ///  Determines if the data stored in two array2_base objects are equal.
    ///
    ///  @param  rhs  array2_base& The right-hand-side array2_base object.
    ///
    ///  @return Returns true if the two array2_base objects are not equal.
    ///////////////////////////////////////////////////////////////////////
    bool operator!=(const array2_base& rhs) const
    {
        return
            (m_asz[0] != rhs.m_asz[0]) ||
            (m_asz[1] != rhs.m_asz[1]) ||
            (::memcmp(m_p,
                      rhs.m_p,
                      m_sz * sizeof(value_type)) != 0);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Determines if the data stored in two array2_base objects are equal.
    ///
    ///  @param  rhs  array2_base& The right-hand-side array2_base object.
    ///
    ///  @return Returns true if the two array2_base objects are equal.
    ///////////////////////////////////////////////////////////////////////
    bool operator==(const array2_base& rhs) const
    {
        return !(*this != rhs);
    }


protected:

    // protected member variables
    value_type*     m_p;
    size_type       m_asz[2];
    size_type       m_sz;

private:
    
    // not implemented
    array2_base(const array2_base&);
    array2_base& operator=(const array2_base&);

};

} // namespace

#endif
