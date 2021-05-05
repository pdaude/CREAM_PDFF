/*
 ==========================================================================
 |   
 |   $Id: array_traits.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___ARRAY_TRAITS_HXX___
#   define ___ARRAY_TRAITS_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <optnet/_base/type.hxx>
#   include <optnet/_base/iterator.hxx>

///  @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class array_traits array_traits.hxx "_base/array_traits.hxx"
///  @brief Array character traits class template.
///////////////////////////////////////////////////////////////////////////
template <typename _T, typename _Container>
class array_traits
{
public:

    ///////////////////////////////////////////////////////////////////////
    /// @brief A type that represents the data type used in an array.
    ///////////////////////////////////////////////////////////////////////
    typedef _T                  value_type;

    ///////////////////////////////////////////////////////////////////////
    /// @brief A type that provides a pointer to an element in an array.
    ///////////////////////////////////////////////////////////////////////
    typedef value_type*         pointer;

    ///////////////////////////////////////////////////////////////////////
    /// @brief A type that provides a pointer to a const element in an
    ///        array.
    ///////////////////////////////////////////////////////////////////////
    typedef const value_type*   const_pointer;

    ///////////////////////////////////////////////////////////////////////
    /// @brief A type that provides a random-access iterator that can read
    ///        or modify any element in a vector.
    ///////////////////////////////////////////////////////////////////////
#   ifndef __OPTNET_CRAPPY_MSC__
    typedef optnet::detail::normal_iterator<pointer, _Container>
                                iterator;
#   else
    typedef value_type*         iterator;
#   endif

    ///////////////////////////////////////////////////////////////////////
    /// @brief A type that provides a random-access iterator that can read
    ///        a const element in an array.
    ///////////////////////////////////////////////////////////////////////
#   ifndef __OPTNET_CRAPPY_MSC__
    typedef optnet::detail::normal_iterator<const_pointer, _Container>
                                const_iterator;
#   else
    typedef const value_type*   const_iterator;
#   endif

    ///////////////////////////////////////////////////////////////////////
    /// @brief A type that provides a reference to an element stored in an
    ///        array.
    ///////////////////////////////////////////////////////////////////////
    typedef value_type&         reference;

    ///////////////////////////////////////////////////////////////////////
    /// @brief A type that provides a reference to a const element stored
    ///        in an array for reading and performing const operations.
    ///////////////////////////////////////////////////////////////////////
    typedef const value_type&   const_reference;

    ///////////////////////////////////////////////////////////////////////
    /// @brief A type that provides the difference between the addresses of
    ///        two elements in an array.
    ///////////////////////////////////////////////////////////////////////
    typedef ptrdiff_t           difference_type;

    ///////////////////////////////////////////////////////////////////////
    /// @brief A type that counts the number of elements in an array.
    ///////////////////////////////////////////////////////////////////////
    typedef size_t              size_type;

};

} // namespace


#endif 
