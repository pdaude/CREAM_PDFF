/*
 ==========================================================================
 |   
 |   $Id: array.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___ARRAY_HXX___
#   define ___ARRAY_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   include <optnet/_base/array_base.hxx>

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma warning(disable: 4244)
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class array array.hxx "_base/array.hxx"
///  @brief Array template class.
///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg = net_f_xy>
class array: public array_base<_Ty, _Tg>
{
    typedef array_base<_Ty, _Tg>    _Base;

public:

    typedef typename _Base::value_type      value_type;
    typedef typename _Base::reference       reference;
    typedef typename _Base::const_reference const_reference;
    typedef typename _Base::iterator        iterator;
    typedef typename _Base::const_iterator  const_iterator;
    typedef typename _Base::pointer         pointer;
    typedef typename _Base::const_pointer   const_pointer;
    typedef typename _Base::difference_type difference_type;
    typedef typename _Base::size_type       size_type;

    // constructor/destructors
    ///////////////////////////////////////////////////////////////////////
    /// Default constructor.
    ///////////////////////////////////////////////////////////////////////
    array();

    ///////////////////////////////////////////////////////////////////////
    ///  Construct and allocate memory for the array.
    ///
    ///  @param asz Array of size_type specifying the array dimensions.
    ///
    ///  @see create
    ///////////////////////////////////////////////////////////////////////
    array(const size_type* asz);

    ///////////////////////////////////////////////////////////////////////
    ///  Construct and allocate memory for the array.
    ///
    ///  @param s0 The size of the first  dimension.
    ///  @param s1 The size of the second dimension.
    ///  @param s2 The size of the third  dimension.
    ///  @param s3 The size of the fourth dimension (default: 1).
    ///  @param s4 The size of the fifth  dimension (default: 1).
    ///
    ///  @see create
    ///////////////////////////////////////////////////////////////////////
    array(size_type s0, 
          size_type s1, 
          size_type s2, 
          size_type s3 = 1, 
          size_type s4 = 1
          );

    ///////////////////////////////////////////////////////////////////////
    /// Copy constructor.
    ///////////////////////////////////////////////////////////////////////
    array(const array& robj);

    ///////////////////////////////////////////////////////////////////////
    /// Copy constructor.
    ///////////////////////////////////////////////////////////////////////
    array(const array_base<_Ty, _Tg>& robj);

    ///////////////////////////////////////////////////////////////////////
    /// Copy constructor.
    ///////////////////////////////////////////////////////////////////////
    array(const array2_base<_Ty>& robj);

    ~array();

    template <typename _Ta>
    ///////////////////////////////////////////////////////////////////////
    ///  Copy to another array.
    ///
    ///  @param copy The destination array.
    ///////////////////////////////////////////////////////////////////////
    void copy_to(array<_Ta>& copy)
    {
        copy.create(_Base::m_asz[0], _Base::m_asz[1],
                    _Base::m_asz[2], _Base::m_asz[3],
                    _Base::m_asz[4]);

        const_iterator it_from = this->begin();
        typename array<_Ta>::iterator it_to = copy.begin();

        for (; it_from != this->end(); ++it_from, ++it_to) {
            *it_to = static_cast<_Ta>(*it_from);
        }
    }

    // modifiers
    array& operator=(const array& robj);
    array& operator=(const array_base<_Ty, _Tg>& robj);
    array& operator=(const array2_base<_Ty>& robj);

    // I have to put this function inline because of the buggy
    // Visual C++ 6.0.
    template <typename _Ta> array& operator=(_Ta& robj)
    {
        if (&robj == this) return *this;
        
        size_type sz = robj.size();
        
        if (_Base::m_sz != sz) {
            // Free allocated resources first.
            SAFE_DELETE_ARRAY(_Base::m_p);
            if (0 != sz) _Base::m_p = new value_type[sz];
            _Base::m_sz = sz;
        }

        const size_type* asz = robj.sizes();
        
        if (asz[0] != _Base::m_asz[0] ||
            asz[1] != _Base::m_asz[1] ||
            asz[2] != _Base::m_asz[2] ||
            asz[3] != _Base::m_asz[3] ||
            asz[4] != _Base::m_asz[4]) {
        
            _Base::m_asz[0] = asz[0];
            _Base::m_asz[1] = asz[1];
            _Base::m_asz[2] = asz[2];
            _Base::m_asz[3] = asz[3];
            _Base::m_asz[4] = asz[4];
              
            _Base::free_lut();
            if (0 != sz) _Base::init_lut();
        }
        
        // Perform deep copy. Use std::copy since the elements may not
        // be simple type.
        std::copy(robj.begin(), robj.end(), _Base::begin());
        
        return *this;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Allocate memory for the array.
    ///
    ///  @param asz Array of size_type specifying the array dimensions.
    ///
    ///  @see create_and_fill
    ///////////////////////////////////////////////////////////////////////
    bool            create (const size_type* asz);

    ///////////////////////////////////////////////////////////////////////
    ///  Allocate memory for the array.
    ///
    ///  @param s0 The size of the first  dimension.
    ///  @param s1 The size of the second dimension.
    ///  @param s2 The size of the third  dimension.
    ///  @param s3 The size of the fourth dimension (default: 1).
    ///  @param s4 The size of the fifth  dimension (default: 1).
    ///
    ///  @see create_and_fill
    ///////////////////////////////////////////////////////////////////////
    bool            create (size_type s0, 
                            size_type s1, 
                            size_type s2, 
                            size_type s3 = 1, 
                            size_type s4 = 1
                            );

    ///////////////////////////////////////////////////////////////////////
    ///  Allocate memory for the array and fill it with some value.
    ///
    ///  @param value The value to fill in the array.
    ///  @param s0    The size of the first  dimension.
    ///  @param s1    The size of the second dimension.
    ///  @param s2    The size of the third  dimension.
    ///  @param s3    The size of the fourth dimension.
    ///  @param s4    The size of the fifth  dimension.
    ///
    ///  @see create
    ///////////////////////////////////////////////////////////////////////
    bool            create_and_fill(const value_type& value,
                                    size_type         s0, 
                                    size_type         s1, 
                                    size_type         s2, 
                                    size_type         s3 = 1, 
                                    size_type         s4 = 1
                                    );
    ///////////////////////////////////////////////////////////////////////
    ///  Free the memory occupied by the array.
    ///
    ///  @see create
    ///////////////////////////////////////////////////////////////////////
    void            release();

    ///////////////////////////////////////////////////////////////////////
    ///  Reshape the array, keeping the total size unchanged.
    ///////////////////////////////////////////////////////////////////////
    bool            reshape(size_type s0, 
                            size_type s1, 
                            size_type s2,
                            size_type s3 = 1, 
                            size_type s4 = 1
                            );
};

} // namespace

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_base/array.cxx>
#   endif

#endif 
