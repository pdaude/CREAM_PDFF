/*
 ==========================================================================
 |   
 |   $Id: array_ref.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___ARRAY_REF_HXX___
#   define ___ARRAY_REF_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4244)
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <optnet/_base/array_base.hxx>

///  @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class array_ref array_ref.hxx "_base/array_ref.hxx"
///  @brief Array template class.
///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg = net_f_xy>
class array_ref: public array_base<_Ty, _Tg>
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
    array_ref();

    ///////////////////////////////////////////////////////////////////////
    ///  Construct and assign reference to the array_ref.
    ///
    ///  @param p     Pointer to the block of memory to be referenced.
    ///  @param s0    The size of the first  dimension.
    ///  @param s1    The size of the second dimension.
    ///  @param s2    The size of the third  dimension.
    ///  @param s3    The size of the fourth dimension.
    ///  @param s4    The size of the fifth  dimension.
    ///
    ///  @see assign.
    ///////////////////////////////////////////////////////////////////////
    explicit array_ref(_Ty* p, 
              size_type s0, 
              size_type s1, 
              size_type s2, 
              size_type s3 = 1, 
              size_type s4 = 1
              );

    // I have to put this function inline because of the buggy
    // Visual C++ 6.0.
    template <typename _Ta>
    array_ref(const array_base<_Ty, _Ta>& robj)
    {
        _Base::m_sz = robj.size();
        ::memcpy(_Base::m_asz, 
                 robj.sizes(), 
                 OPTNET_ARRAY_DIMS * sizeof(size_type));

        _Base::m_p = const_cast<_Ty*>(robj.data());

        // Initialize the look-up-table.
        _Base::init_lut();
    }
    
    array_ref(const array_ref& robj);
    ~array_ref();

    // modifiers
    array_ref& operator=(const array_ref& robj);

    // I have to put this function inline because of the buggy
    // Visual C++ 6.0.
    template <typename _Ta>
    array_ref&
    operator=(const array_base<_Ty, _Ta>& robj)
    {
        if (&robj != this) {
            _Base::m_sz = robj.size();
            _Base::m_p  = const_cast<_Ty*>(robj.data());

            ::memcpy(_Base::m_asz, 
                     robj.sizes(), 
                     OPTNET_ARRAY_DIMS * sizeof(size_type));
                         
            // Reinitialize the look-up-table.
            _Base::free_lut();
            if (0 != _Base::m_sz) _Base::init_lut();
        }

        return *this;
    }


    ///////////////////////////////////////////////////////////////////////
    ///  Assign reference to the array_ref.
    ///
    ///  @param p     Pointer to the block of memory to be referenced.
    ///  @param s0    The size of the first  dimension.
    ///  @param s1    The size of the second dimension.
    ///  @param s2    The size of the third  dimension.
    ///  @param s3    The size of the fourth dimension.
    ///  @param s4    The size of the fifth  dimension.
    ///
    ///  @return Returns true if the memory is sucessfully referenced.
    ///////////////////////////////////////////////////////////////////////
    bool            assign (const _Ty* p,
                            size_type s0, 
                            size_type s1, 
                            size_type s2 = 1, 
                            size_type s3 = 1, 
                            size_type s4 = 1
                            );

    ///////////////////////////////////////////////////////////////////////
    ///  Reshape the array, keeping the total size unchanged.
    ///////////////////////////////////////////////////////////////////////
    bool            reshape(size_type s0, 
                            size_type s1, 
                            size_type s2 = 1, 
                            size_type s3 = 1, 
                            size_type s4 = 1
                            );

};


} // namespace


#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_base/array_ref.cxx>
#   endif

#endif
