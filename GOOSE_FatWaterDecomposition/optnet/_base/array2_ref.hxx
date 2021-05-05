/*
 ==========================================================================
 |   
 |   $Id: array2_ref.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___ARRAY2_REF_HXX___
#   define ___ARRAY2_REF_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   include <optnet/_base/array2_base.hxx>

///  @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class array2_ref array2_ref.hxx "_base/array2_ref.hxx"
///  @brief Array template class.
///////////////////////////////////////////////////////////////////////////
template <typename _Ty>
class array2_ref: public array2_base<_Ty>
{
    typedef array2_base<_Ty>    _Base;

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
    array2_ref() :
        array2_base<_Ty>()
    {
    }
    
    explicit array2_ref(_Ty* p, size_type s0, size_type s1)
    {
        assert(p != 0);
        assert(s0 > 0);
        assert(s1 > 0);

        _Base::m_p      = p;
        _Base::m_sz     = s0 * s1;
        _Base::m_asz[0] = s0;
        _Base::m_asz[1] = s1;
    }

    array2_ref(const array2_base<_Ty>& robj)
    {
        assert(&robj != this);
        
        _Base::m_p = const_cast<_Ty*>(robj.data());

        _Base::m_asz[0] = robj.sizes()[0];
        _Base::m_asz[1] = robj.sizes()[1];
        _Base::m_sz     = robj.size();
    }
    
    array2_ref(const array2_ref& robj)
    {
        assert(&robj != this);

        _Base::m_p      = const_cast<_Ty*>(robj._Base::m_p);
        _Base::m_asz[0] = robj._Base::m_asz[0];
        _Base::m_asz[1] = robj._Base::m_asz[1];
        _Base::m_sz     = robj._Base::m_sz;
    }
    
    ~array2_ref()
    {
    }

    // modifiers
    array2_ref& operator=(const array2_ref& robj)
    {
        if (&robj != this) {
            _Base::m_p      = const_cast<_Ty*>(robj._Base::m_p);
            _Base::m_asz[0] = robj._Base::m_asz[0];
            _Base::m_asz[1] = robj._Base::m_asz[1];
            _Base::m_sz     = robj._Base::m_sz;
        }
        return *this;
    }

    array2_ref& operator=(const array2_base<_Ty>& robj)
    {
        _Base::m_p  = const_cast<_Ty*>(robj.data());

        _Base::m_asz[0] = robj.sizes()[0];
        _Base::m_asz[1] = robj.sizes()[1];
        _Base::m_sz     = robj.size();

        return *this;
    }

    bool assign(const _Ty* p, size_type s0, size_type s1)
    {
        assert(s0 > 0);
        assert(s1 > 0);

        if (0 == p) return false;

        _Base::m_p      = p;
        _Base::m_sz     = s0 * s1;
        _Base::m_asz[0] = s0;
        _Base::m_asz[1] = s1;
        
        return true;
    }    
};

} // namespace

#endif
