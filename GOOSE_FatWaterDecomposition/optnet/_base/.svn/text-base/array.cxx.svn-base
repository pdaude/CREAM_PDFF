/*
 ==========================================================================
 |   
 |   $Id: array.cxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___ARRAY_CXX___
#   define ___ARRAY_CXX___

#   include <optnet/_base/array.hxx>

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4244)
#       if (_MSC_VER <= 1200)
#           pragma warning(disable: 4018)
#           pragma warning(disable: 4146)
#       endif
#   endif

#   include <algorithm>


namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>::array() :
    array_base<_Ty, _Tg>()
{
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>::array(const size_type* asz)
{
    assert(asz[0] > 0);
    assert(asz[1] > 0);
    assert(asz[2] > 0);
    assert(asz[3] > 0);
    assert(asz[4] > 0);

    size_type sz = asz[0] * asz[1] * asz[2] * asz[3] * asz[4];

    if (0 != sz) {
        _Base::m_p = new value_type[sz];

        // FIXME: Just to guard against non-standard new operator,
        //        which does not throw.
        assert(0 != _Base::m_p); 
        
        std::fill(_Base::m_p, _Base::m_p + sz, value_type());
    }

    _Base::m_sz     = sz;
    _Base::m_asz[0] = asz[0];
    _Base::m_asz[1] = asz[1];
    _Base::m_asz[2] = asz[2];
    _Base::m_asz[3] = asz[3];
    _Base::m_asz[4] = asz[4];

    if (0 != sz) _Base::init_lut();
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>::array(size_type s0, 
                       size_type s1, 
                       size_type s2, 
                       size_type s3, 
                       size_type s4
                       )
{
    size_type sz = s0 * s1 * s2 * s3 * s4;

    if (0 != sz) {
        _Base::m_p = new value_type[sz];
        assert(0 != _Base::m_p); // Note: see above.
        std::fill(_Base::m_p, _Base::m_p + sz, value_type());
    }

    _Base::m_sz     = sz;
    _Base::m_asz[0] = s0;
    _Base::m_asz[1] = s1;
    _Base::m_asz[2] = s2;
    _Base::m_asz[3] = s3;
    _Base::m_asz[4] = s4;

    if (0 != sz) _Base::init_lut();
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>::array(const array& robj)
{
    if (0 != robj.m_sz) {
        _Base::m_p = new value_type[robj.m_sz];
        
        assert(0 != _Base::m_p); // Note: see above.
        
        // Perform deep copy. Use std::copy since the elements may not
        // be simple type.
        std::copy(robj.begin(), robj.end(), _Base::begin());
    }
    
    _Base::m_sz = robj.m_sz;
    _Base::m_asz[0] = robj.m_asz[0];
    _Base::m_asz[1] = robj.m_asz[1];
    _Base::m_asz[2] = robj.m_asz[2];
    _Base::m_asz[3] = robj.m_asz[3];
    _Base::m_asz[4] = robj.m_asz[4];
      
    if (0 != robj.m_sz) _Base::init_lut();
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>::array(const array_base<_Ty, _Tg>& robj)
{
    size_type sz = robj.size();
    const size_type* asz = robj.sizes();
    
    if (0 != sz) {
        _Base::m_p = new value_type[sz];
        
        assert(0 != _Base::m_p); // Note: see above.
        
        // Perform deep copy. Use std::copy since the elements may not
        // be simple type.
        std::copy(robj.begin(), robj.end(), _Base::begin());
    }
    
    _Base::m_sz = sz;
    _Base::m_asz[0] = asz[0];
    _Base::m_asz[1] = asz[1];
    _Base::m_asz[2] = asz[2];
    _Base::m_asz[3] = asz[3];
    _Base::m_asz[4] = asz[4];
    
    if (0 != sz) _Base::init_lut();
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>::array(const array2_base<_Ty>& robj)
{
    size_type sz = robj.size();
    
    if (0 != sz) {
        _Base::m_p = new value_type[sz];
        
        assert(0 != _Base::m_p); // Note: see above.
        
        // Perform deep copy. Use std::copy since the elements may not
        // be simple type.
        std::copy(robj.begin(), robj.end(), _Base::begin());
    }

    const size_type* asz = robj.sizes();
        
    _Base::m_sz = sz;
    _Base::m_asz[0] = asz[0];
    _Base::m_asz[1] = asz[1];
    _Base::m_asz[2] = 1;
    _Base::m_asz[3] = 1;
    _Base::m_asz[4] = 1;
    
    if (0 != sz) _Base::init_lut();
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>::~array()
{
    release();
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>&
array<_Ty, _Tg>::operator=(const array& robj)
{
    if (&robj == this) return *this;
    
    if (_Base::m_sz != robj.m_sz) {
        // Free allocated resources first.
        SAFE_DELETE_ARRAY(_Base::m_p);
        if (0 != robj.m_sz) _Base::m_p = new value_type[robj.m_sz];
        _Base::m_sz = robj.m_sz;
    }

    if (robj.m_asz[0] != _Base::m_asz[0] ||
        robj.m_asz[1] != _Base::m_asz[1] ||
        robj.m_asz[2] != _Base::m_asz[2] ||
        robj.m_asz[3] != _Base::m_asz[3] ||
        robj.m_asz[4] != _Base::m_asz[4]) {
    
        _Base::m_asz[0] = robj.m_asz[0];
        _Base::m_asz[1] = robj.m_asz[1];
        _Base::m_asz[2] = robj.m_asz[2];
        _Base::m_asz[3] = robj.m_asz[3];
        _Base::m_asz[4] = robj.m_asz[4];
          
        _Base::free_lut();
        if (0 != robj.m_sz) _Base::init_lut();
    }
    
    // Perform deep copy. Use std::copy since the elements may not
    // be simple type.
    std::copy(robj.begin(), robj.end(), _Base::begin());
          
    return *this;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>&
array<_Ty, _Tg>::operator=(const array_base<_Ty, _Tg>& robj)
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

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>&
array<_Ty, _Tg>::operator=(const array2_base<_Ty>& robj)
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

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
bool
array<_Ty, _Tg>::create (const size_type* asz)
{
    size_type sz = asz[0] * asz[1] * asz[2] * asz[3] * asz[4];

    // Relocate memory only if size is changed.
    if (_Base::m_sz != sz) {

        // Free allocated resources first.
        SAFE_DELETE_ARRAY(_Base::m_p);

        _Base::m_p = new value_type[sz];
        if (0 == _Base::m_p) // Should not happen for standard new.
            return false;

        _Base::m_sz = sz;
    }

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
    
    std::fill(_Base::m_p, _Base::m_p + sz, value_type());

    return true;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
bool
array<_Ty, _Tg>::create (size_type s0, 
                         size_type s1, 
                         size_type s2, 
                         size_type s3,
                         size_type s4
                         )
{
    size_type sz = s0 * s1 * s2 * s3 * s4;

    // Relocate memory only if size is changed.
    if (_Base::m_sz != sz) {

        // Free allocated resources first.
        SAFE_DELETE_ARRAY(_Base::m_p);

        _Base::m_p = new value_type[sz];
        if (0 == _Base::m_p) // Should not happen for standard new.
            return false;

        _Base::m_sz = sz;
    }

    if (s0 != _Base::m_asz[0] ||
        s1 != _Base::m_asz[1] ||
        s2 != _Base::m_asz[2] ||
        s3 != _Base::m_asz[3] ||
        s4 != _Base::m_asz[4]) {
    
        _Base::m_asz[0] = s0;
        _Base::m_asz[1] = s1;
        _Base::m_asz[2] = s2;
        _Base::m_asz[3] = s3;
        _Base::m_asz[4] = s4;
          
        _Base::free_lut();
        if (0 != sz) _Base::init_lut();
    }
    
    std::fill(_Base::m_p, _Base::m_p + sz, value_type());

    return true;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
bool
array<_Ty, _Tg>::create_and_fill(const value_type& value,
                                 size_type         s0, 
                                 size_type         s1, 
                                 size_type         s2, 
                                 size_type         s3,
                                 size_type         s4
                                 )
{
    size_type sz = s0 * s1 * s2 * s3 * s4;

    // Relocate memory only if size is changed.
    if (_Base::m_sz != sz) {

        // Free allocated resources first.
        SAFE_DELETE_ARRAY(_Base::m_p);

        _Base::m_p = new value_type[sz];
        if (0 == _Base::m_p) // Should not happen for standard new.
            return false;

        _Base::m_sz = sz;
    }

    if (s0 != _Base::m_asz[0] ||
        s1 != _Base::m_asz[1] ||
        s2 != _Base::m_asz[2] ||
        s3 != _Base::m_asz[3] ||
        s4 != _Base::m_asz[4]) {
    
        _Base::m_asz[0] = s0;
        _Base::m_asz[1] = s1;
        _Base::m_asz[2] = s2;
        _Base::m_asz[3] = s3;
        _Base::m_asz[4] = s4;
          
        _Base::free_lut();
        if (0 != sz) _Base::init_lut();
    }
    
    std::fill(_Base::m_p, _Base::m_p + sz, value);
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
void
array<_Ty, _Tg>::release()
{
    // Free array data buffer.
    SAFE_DELETE_ARRAY(_Base::m_p);

    // Free look-up-table.
    _Base::free_lut();
    
    // Reset sizes.
    _Base::m_sz     = 0;
    _Base::m_asz[0] = 0;
    _Base::m_asz[1] = 0;
    _Base::m_asz[2] = 0;
    _Base::m_asz[3] = 0;
    _Base::m_asz[4] = 0;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
bool
array<_Ty, _Tg>::reshape(size_type s0, 
                         size_type s1, 
                         size_type s2, 
                         size_type s3,
                         size_type s4
                         )
{
    if (0 == _Base::m_p) return false;

    size_type sz = s0 * s1 * s2 * s3 * s4;

    // The array sizes must be equal.
    if (_Base::m_sz == sz) {

        _Base::free_lut();

        _Base::m_asz[0] = s0;
        _Base::m_asz[1] = s1;
        _Base::m_asz[2] = s2;
        _Base::m_asz[3] = s3;
        _Base::m_asz[4] = s4;

        return _Base::init_lut();
    }

    return false;
}

} // namespace

#endif
