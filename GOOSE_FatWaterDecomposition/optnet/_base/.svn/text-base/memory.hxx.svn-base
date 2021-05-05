/*
 ==========================================================================
 |   
 |   $Id: memory.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___MEMORY_HXX___
#   define ___MEMORY_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <memory.h>
#   include <cstddef>

/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
/// Allocate an m-by-n (row major) matrix.
///////////////////////////////////////////////////////////////////////////
template <typename _T>
_T** new_matrix(size_t m, size_t n)
{
    _T** pp = new _T*[m];
    *pp     = new _T [m * n];
    for (_T** pp1 = pp; pp1 < pp + m - 1; ++pp1)
        *(pp1 + 1) = *pp1 + n;
    return pp;
}

///////////////////////////////////////////////////////////////////////////
/// Allocate an m-by-n (row major) matrix with initial value.
///////////////////////////////////////////////////////////////////////////
template <typename _T>
_T** new_matrix(size_t m, size_t n, const _T& v)
{
    size_t i;
    _T** pp = new _T*[m];
    *pp     = new _T [m * n];
    for (i = 0; i < m * n; ++i) (*pp)[i] = v;
    for (i = 0; i + 1 < m; ++i) pp[i + 1] = pp[i] + n;
    return pp;
}

///////////////////////////////////////////////////////////////////////////
/// Free a 2-D matrix.
///////////////////////////////////////////////////////////////////////////
template <typename _T>
inline void delete_matrix(_T** pp)
{
    delete[] *pp;
    delete[]  pp;
}

} // namespace

///////////////////////////////////////////////////////////////////////////
//  macros

#   ifndef SAFE_DELETE
/// Delete a variable pointed by p, and set p to NULL.
#       define SAFE_DELETE(p)       { if (p) { delete   (p); (p) = 0; } }
#   endif // SAFE_DELETE

#   ifndef SAFE_DELETE_ARRAY
/// Delete an array   pointed by p, and set p to NULL.
#       define SAFE_DELETE_ARRAY(p) { if (p) { delete[] (p); (p) = 0; } }
#   endif // SAFE_DELETE_ARRAY

#   ifndef SAFE_DELETE_MATRIX
/// Delete a matrix pointed by pp, and set pp to NULL.
#       define SAFE_DELETE_MATRIX(pp) \
                 { if (pp) { optnet::delete_matrix(pp); (pp) = 0; } }
#   endif // SAFE_DELETE_MATRIX

#endif 
