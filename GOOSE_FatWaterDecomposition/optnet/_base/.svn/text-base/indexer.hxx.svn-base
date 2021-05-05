/*
 ==========================================================================
 |   
 |   $Id: indexer.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___INDEXER_HXX___
#   define ___INDEXER_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <optnet/_base/tags.hxx>
#   include <cassert>

//
// Note: Microsoft VC6 does not support class template partial
//       specialization. For compatibility, we are passing all
//       parameters as function arguments. This makes the im-
//       plementation much more clumsy.
//

/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class indexer
///  @brief Array indexer template.
///
///  This template should never be instantiated in real programs.
///  Use its specializations instead.
///////////////////////////////////////////////////////////////////////////
template <typename _Tg>
class indexer
{
public:
    template <typename _Ty> 
    inline _Ty& size_0(const _Ty* s) const 
    { 
        assert(false); 
        return _Ty(); 
    }
    template <typename _Ty> 
    inline _Ty& size_1(const _Ty* s) const
    {
        assert(false);
        return _Ty();
    }
    template <typename _Ty> 
    inline _Ty& size_2(const _Ty* s) const
    {
        assert(false);
        return _Ty();
    }
    template <typename _Ty> 
    inline _Ty& size_3(const _Ty* s) const
    { 
        assert(false);
        return _Ty();
    }
    template <typename _Ty> 
    inline _Ty& size_4(const _Ty* s) const
    {
        assert(false);
        return _Ty();
    }
    
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***   p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     const _Ti* s
                     ) const
    {
        assert(false);
        return _Ty();
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty****  p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     const _Ti* s
                     ) const
    {
        assert(false);
        return _Ty();
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***** p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     _Ti i4,
                     const _Ti* s
                     ) const
    {
        assert(false);
        return _Ty();
    }
};

///////////////////////////////////////////////////////////////////////////
///  @class indexer<net_f_xy>
///  @brief Template specialization of indexer.
///
///  Array indexer for the net in the form Net = f(i0, i1).
///////////////////////////////////////////////////////////////////////////
template <>
class indexer<net_f_xy>
{
public:
    template <typename _Ty> 
        inline _Ty& size_0(_Ty* s) const { return s[0]; }
    template <typename _Ty> 
        inline _Ty& size_1(_Ty* s) const { return s[1]; }
    template <typename _Ty> 
        inline _Ty& size_2(_Ty* s) const { return s[2]; }
    template <typename _Ty> 
        inline _Ty& size_3(_Ty* s) const { return s[3]; }
    template <typename _Ty> 
        inline _Ty& size_4(_Ty* s) const { return s[4]; }

    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***   p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     const _Ti*
                     ) const
    {
       return p[i2][i1][i0];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty****  p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     const _Ti*
                     ) const
    {
       return p[i3][i2][i1][i0];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***** p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     _Ti i4,
                     const _Ti*
                     ) const
    {
       return p[i4][i3][i2][i1][i0];
    }
};

///////////////////////////////////////////////////////////////////////////
///  @class indexer<net_f_yz>
///  @brief Template specialization of indexer.
///
///  Array indexer for the net in the form Net = f(i1, i2).
///////////////////////////////////////////////////////////////////////////
template <>
class indexer<net_f_yz>
{ 
public:
    template <typename _Ty>
        inline _Ty& size_0(_Ty* s) const { return s[1]; }
    template <typename _Ty>
        inline _Ty& size_1(_Ty* s) const { return s[2]; }
    template <typename _Ty>
        inline _Ty& size_2(_Ty* s) const { return s[0]; }
    template <typename _Ty>
        inline _Ty& size_3(_Ty* s) const { return s[3]; }
    template <typename _Ty>
        inline _Ty& size_4(_Ty* s) const { return s[4]; }

    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***   p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     const _Ti*
                     ) const
    {
       return p[i1][i0][i2];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty****  p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     const _Ti*
                     ) const
    {
       return p[i3][i1][i0][i2];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***** p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     _Ti i4,
                     const _Ti*
                     ) const
    {
       return p[i4][i3][i1][i0][i2];
    }
};

///////////////////////////////////////////////////////////////////////////
///  @class indexer<net_f_zx>
///  @brief Template specialization of indexer.
///
///  Array indexer for the net in the form Net = f(z, x).
///////////////////////////////////////////////////////////////////////////
template <>
class indexer<net_f_zx>
{ 
public:
    template <typename _Ty> 
        inline _Ty& size_0(_Ty* s) const { return s[2]; }
    template <typename _Ty> 
        inline _Ty& size_1(_Ty* s) const { return s[0]; }
    template <typename _Ty> 
        inline _Ty& size_2(_Ty* s) const { return s[1]; }
    template <typename _Ty> 
        inline _Ty& size_3(_Ty* s) const { return s[3]; }
    template <typename _Ty> 
        inline _Ty& size_4(_Ty* s) const { return s[4]; }

    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***   p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     const _Ti*
                     ) const
    {
       return p[i0][i2][i1];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty****  p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     const _Ti*
                     ) const
    {
       return p[i3][i0][i2][i1];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***** p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     _Ti i4,
                     const _Ti*
                     ) const
    {
       return p[i4][i3][i0][i2][i1];
    }
};

///////////////////////////////////////////////////////////////////////////
///  @class indexer<net_f_xy_zflipped>
///  @brief Template specialization of indexer.
///
///  Array indexer for the net in the form Net = s2 - 1 - f(i0, i1).
///////////////////////////////////////////////////////////////////////////
template <>
class indexer<net_f_xy_zflipped>
{
public:
    template <typename _Ty> 
        inline _Ty& size_0(_Ty* s) const { return s[0]; }
    template <typename _Ty> 
        inline _Ty& size_1(_Ty* s) const { return s[1]; }
    template <typename _Ty> 
        inline _Ty& size_2(_Ty* s) const { return s[2]; }
    template <typename _Ty> 
        inline _Ty& size_3(_Ty* s) const { return s[3]; }
    template <typename _Ty> 
        inline _Ty& size_4(_Ty* s) const { return s[4]; }

    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***   p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     const _Ti* s
                     ) const
    {
       return p[s[2] - i2 - 1][i1][i0];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty****  p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     const _Ti* s
                     ) const
    {
       return p[i3][s[2] - i2 - 1][i1][i0];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***** p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     _Ti i4,
                     const _Ti* s
                     ) const
    {
       return p[i4][i3][s[2] - i2 - 1][i1][i0];
    }
};

///////////////////////////////////////////////////////////////////////////
///  @class indexer<net_f_yz_xflipped>
///  @brief Template specialization of indexer.
///
///  Array indexer for the net in the form Net = s0 - 1 - f(i1, i2).
///////////////////////////////////////////////////////////////////////////
template <>
class indexer<net_f_yz_xflipped>
{ 
public:
    template <typename _Ty>
        inline _Ty& size_0(_Ty* s) const { return s[1]; }
    template <typename _Ty>
        inline _Ty& size_1(_Ty* s) const { return s[2]; }
    template <typename _Ty>
        inline _Ty& size_2(_Ty* s) const { return s[0]; }
    template <typename _Ty>
        inline _Ty& size_3(_Ty* s) const { return s[3]; }
    template <typename _Ty>
        inline _Ty& size_4(_Ty* s) const { return s[4]; }

    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***   p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     const _Ti* s
                     ) const
    {
       return p[i1][i0][s[0] - i2 - 1];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty****  p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     const _Ti* s
                     ) const
    {
       return p[i3][i1][i0][s[0] - i2 - 1];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***** p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     _Ti i4,
                     const _Ti* s
                     ) const
    {
       return p[i4][i3][i1][i0][s[0] - i2 - 1];
    }
};

///////////////////////////////////////////////////////////////////////////
///  @class indexer<net_f_zx_yflipped>
///  @brief Template specialization of indexer.
///
///  Array indexer for the net in the form Net = s1 - 1 - f(z, x).
///////////////////////////////////////////////////////////////////////////
template <>
class indexer<net_f_zx_yflipped>
{ 
public:
    template <typename _Ty> 
        inline _Ty& size_0(_Ty* s) const { return s[2]; }
    template <typename _Ty> 
        inline _Ty& size_1(_Ty* s) const { return s[0]; }
    template <typename _Ty> 
        inline _Ty& size_2(_Ty* s) const { return s[1]; }
    template <typename _Ty> 
        inline _Ty& size_3(_Ty* s) const { return s[3]; }
    template <typename _Ty> 
        inline _Ty& size_4(_Ty* s) const { return s[4]; }

    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***   p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     const _Ti* s
                     ) const
    {
       return p[i0][s[1] - i2 - 1][i1];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty****  p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     const _Ti* s
                     ) const
    {
       return p[i3][i0][s[1] - i2 - 1][i1];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***** p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     _Ti i4,
                     const _Ti* s
                     ) const
    {
       return p[i4][i3][i0][s[1] - i2 - 1][i1];
    }
};

} // namespace


#endif 
