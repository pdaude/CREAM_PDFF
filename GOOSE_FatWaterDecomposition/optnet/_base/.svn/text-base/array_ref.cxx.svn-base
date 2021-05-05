/*
 ==========================================================================
 |   
 |   $Id: array_ref.cxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___ARRAY_REF_CXX___
#   define ___ARRAY_REF_CXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4244)
#   endif

#   include <optnet/_base/array_ref.hxx>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array_ref<_Ty, _Tg>::array_ref()
{
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array_ref<_Ty, _Tg>::array_ref(_Ty* p, 
                               size_type s0, 
                               size_type s1, 
                               size_type s2, 
                               size_type s3, 
                               size_type s4
                               )
{
    assert(p != 0);
    assert(s0 > 0);
    assert(s1 > 0);
    assert(s2 > 0);
    assert(s3 > 0);
    assert(s4 > 0);

    _Base::m_p      = p;
    _Base::m_sz     = s0 * s1 * s2 * s3 * s4;
    _Base::m_asz[0] = s0;
    _Base::m_asz[1] = s1;
    _Base::m_asz[2] = s2;
    _Base::m_asz[3] = s3;
    _Base::m_asz[4] = s4;

    _Base::init_lut();
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array_ref<_Ty, _Tg>::array_ref(const array_ref<_Ty, _Tg>& robj)
{
    assert(&robj != this);

    _Base::m_p  = robj.m_p;
    _Base::m_sz = robj.m_sz;
    ::memcpy(_Base::m_asz, 
             robj.m_asz, 
             OPTNET_ARRAY_DIMS * sizeof(size_type));

    // Initialize the look-up-table.
    _Base::init_lut();
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array_ref<_Ty, _Tg>::~array_ref()
{
    _Base::free_lut();
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array_ref<_Ty, _Tg>&
array_ref<_Ty, _Tg>::operator=(const array_ref<_Ty, _Tg>& robj)
{
    if (&robj != this) {
        _Base::m_p  = robj.m_p;
        _Base::m_sz = robj.m_sz;
        ::memcpy(_Base::m_asz, 
                 robj.m_asz, 
                 OPTNET_ARRAY_DIMS * sizeof(size_type));

        // Reinitialize the look-up-table.
        _Base::free_lut();
        if (0 != _Base::m_sz) _Base::init_lut();
    }
    return *this;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
bool
array_ref<_Ty, _Tg>::assign (const _Ty* p,
                             size_type s0, 
                             size_type s1, 
                             size_type s2, 
                             size_type s3, 
                             size_type s4
                             )
{
    if (0 == p) return false;

    _Base::m_p = p;

    if (_Base::m_asz[0] != s0 ||
        _Base::m_asz[1] != s1 || 
        _Base::m_asz[2] != s2 ||
        _Base::m_asz[3] != s3 || 
        _Base::m_asz[4] != s4) {

        _Base::free_lut();

        _Base::m_sz     = s0 * s1 * s2 * s3 * s4;
        _Base::m_asz[0] = s0;
        _Base::m_asz[1] = s1;
        _Base::m_asz[2] = s2;
        _Base::m_asz[3] = s3;
        _Base::m_asz[4] = s4;
        
        return _Base::init_lut();
    }

    return true;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
bool
array_ref<_Ty, _Tg>::reshape(size_type s0, 
                             size_type s1, 
                             size_type s2, 
                             size_type s3, 
                             size_type s4
                             )
{
    assert(s0 > 0);
    assert(s1 > 0);
    assert(s2 > 0);
    assert(s3 > 0);
    assert(s4 > 0);

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
