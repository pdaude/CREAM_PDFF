/*
 ==========================================================================
 |   
 |   $Id: iterator.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___ITERATOR_HXX___
#   define ___ITERATOR_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <iterator>
#   include <optnet/config.h>

/// @namespace optnet
namespace optnet {

namespace detail {

//
// class normal_iterator
//
#   ifndef __OPTNET_CRAPPY_MSC__

/* 
 * 
 * Copyright (c) 1994 Hewlett-Packard Company
 * 
 * Permission to use,  copy, modify, distribute  and sell this  software and
 * its documentation for any purpose is hereby granted without fee, provided
 * that the above copyright notice appear  in all copies and that both  that
 * copyright  notice  and  this  permission  notice  appear  in   supporting
 * documentation.  Hewlett-Packard  Company makes  no representations  about
 * the suitability of this software for any purpose.  It is provided "as is"
 * without express or implied warranty.
 * 
 * 
 * Copyright (c) 1996-1998 Silicon Graphics Computer Systems, Inc.
 * 
 * Permission to use,  copy, modify, distribute  and sell this  software and
 * its documentation for any purpose is hereby granted without fee, provided
 * that the above copyright notice appear  in all copies and that both  that
 * copyright  notice  and  this  permission  notice  appear  in   supporting
 * documentation.   Silicon  Graphics  makes  no  representations  about the
 * suitability of  this software  for any  purpose.  It  is provided "as is"
 * without express or implied warranty.
 */

    ///////////////////////////////////////////////////////////////////////
    ///  @class normal_iterator
    ///  @brief Standard iterator type.
    ///////////////////////////////////////////////////////////////////////
    template<typename _Iterator, typename _Container>
    class normal_iterator
        : public std::iterator<
            typename std::iterator_traits<_Iterator>::iterator_category,
            typename std::iterator_traits<_Iterator>::value_type,
            typename std::iterator_traits<_Iterator>::difference_type,
            typename std::iterator_traits<_Iterator>::pointer,
            typename std::iterator_traits<_Iterator>::reference>
    {
    protected:
        _Iterator m_current;

    public:
        typedef typename std::iterator_traits<_Iterator>::difference_type    
            difference_type;
        typedef typename std::iterator_traits<_Iterator>::reference reference;
        typedef typename std::iterator_traits<_Iterator>::pointer   pointer;

        normal_iterator() : m_current(_Iterator()) { }

        explicit 
            normal_iterator(const _Iterator& __i) : m_current(__i) { }

        // Allow iterator to const_iterator conversion
        template<typename _Iter>
        inline normal_iterator(const normal_iterator<_Iter, _Container>& __i)
        : m_current(__i.base()) { }

        // Forward iterator requirements
        reference operator* () const { return *m_current; }
        pointer   operator->() const { return  m_current; }

        normal_iterator& operator++() { ++m_current; return *this; }
        normal_iterator  operator++(int)
        { return normal_iterator(m_current++); }

        // Bidirectional iterator requirements
        normal_iterator& operator--() { --m_current; return *this; }
        normal_iterator  operator--(int)
        { return normal_iterator(m_current--); }

        // Random access iterator requirements
        reference        operator[](const difference_type& __n) const
        { return m_current[__n]; }

        normal_iterator& operator+=(const difference_type& __n)
        { m_current += __n; return *this; }

        normal_iterator  operator+ (const difference_type& __n) const
        { return normal_iterator(m_current + __n); }

        normal_iterator& operator-=(const difference_type& __n)
        { m_current -= __n; return *this; }

        normal_iterator  operator- (const difference_type& __n) const
        { return normal_iterator(m_current - __n); }

        const _Iterator& base() const { return m_current; }
    };

    // Forward iterator requirements
    template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline bool
    operator==(const normal_iterator<_IteratorL, _Container>& __lhs,
               const normal_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() == __rhs.base(); }

    template<typename _Iterator, typename _Container>
    inline bool
    operator==(const normal_iterator<_Iterator,  _Container>& __lhs,
               const normal_iterator<_Iterator,  _Container>& __rhs)
    { return __lhs.base() == __rhs.base(); }

    template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline bool
    operator!=(const normal_iterator<_IteratorL, _Container>& __lhs,
               const normal_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() != __rhs.base(); }

    template<typename _Iterator, typename _Container>
    inline bool
    operator!=(const normal_iterator<_Iterator,  _Container>& __lhs,
               const normal_iterator<_Iterator,  _Container>& __rhs)
    { return __lhs.base() != __rhs.base(); }

    // Random access iterator requirements
    template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline bool 
    operator< (const normal_iterator<_IteratorL, _Container>& __lhs,
               const normal_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() < __rhs.base(); }

    template<typename _Iterator, typename _Container>
    inline bool
    operator< (const normal_iterator<_Iterator,  _Container>& __lhs,
               const normal_iterator<_Iterator,  _Container>& __rhs)
    { return __lhs.base() < __rhs.base(); }

    template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline bool
    operator> (const normal_iterator<_IteratorL, _Container>& __lhs,
               const normal_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() > __rhs.base(); }

    template<typename _Iterator, typename _Container>
    inline bool
    operator> (const normal_iterator<_Iterator,  _Container>& __lhs,
               const normal_iterator<_Iterator,  _Container>& __rhs)
    { return __lhs.base() > __rhs.base(); }

    template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline bool
    operator<=(const normal_iterator<_IteratorL, _Container>& __lhs,
               const normal_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() <= __rhs.base(); }

    template<typename _Iterator, typename _Container>
    inline bool
    operator<=(const normal_iterator<_Iterator,  _Container>& __lhs,
               const normal_iterator<_Iterator,  _Container>& __rhs)
    { return __lhs.base() <= __rhs.base(); }

    template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline bool
    operator>=(const normal_iterator<_IteratorL, _Container>& __lhs,
               const normal_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() >= __rhs.base(); }

    template<typename _Iterator, typename _Container>
    inline bool
    operator>=(const normal_iterator<_Iterator,  _Container>& __lhs,
               const normal_iterator<_Iterator,  _Container>& __rhs)
    { return __lhs.base() >= __rhs.base(); }

    template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline typename normal_iterator<_IteratorL, _Container>::difference_type
    operator- (const normal_iterator<_IteratorL, _Container>& __lhs,
               const normal_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() - __rhs.base(); }

    template<typename _Iterator, typename _Container>
    inline normal_iterator<_Iterator, _Container>
    operator+ (
        typename normal_iterator<_Iterator, _Container>::difference_type __n,
        const    normal_iterator<_Iterator, _Container>& __i)
    { return normal_iterator<_Iterator, _Container>(__i.base() + __n); }

#   endif // __OPTNET_CRAPPY_MSC__

} // namespace

} // namespace


#endif 
