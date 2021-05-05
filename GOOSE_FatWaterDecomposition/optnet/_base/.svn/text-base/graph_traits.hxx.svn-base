/*
 ==========================================================================
 |   
 |   $Id: graph_traits.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___GRAPH_TRAITS_HXX___
#   define ___GRAPH_TRAITS_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/_base/type.hxx>


/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class graph_traits graph_traits.hxx "_base/graph_traits.hxx"
///  @brief Graph character traits class template.
///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _NodeCnt,
          typename _FwdArcCnt = dummy_type,
          typename _RevArcCnt = dummy_type>
class graph_traits
{
public:

    typedef _Cap                    capacity_type;

    typedef typename 
        _NodeCnt::value_type        node_type;
    typedef typename 
        _NodeCnt::reference         node_reference;
    typedef typename 
        _NodeCnt::const_reference   node_const_reference;
    typedef typename 
        _NodeCnt::iterator          node_iterator;
    typedef typename 
        _NodeCnt::const_iterator    node_const_iterator;
    typedef typename 
        _NodeCnt::pointer           node_pointer;
    typedef typename 
        _NodeCnt::const_pointer     node_const_pointer;

    typedef typename 
        _FwdArcCnt::value_type      forward_arc_type;
    typedef typename 
        _FwdArcCnt::reference       forward_arc_reference;
    typedef typename 
        _FwdArcCnt::const_reference forward_arc_const_reference;
    typedef typename 
        _FwdArcCnt::iterator        forward_arc_iterator;
    typedef typename 
        _FwdArcCnt::const_iterator  forward_arc_const_iterator;
    typedef typename 
        _FwdArcCnt::pointer         forward_arc_pointer;
    typedef typename 
        _FwdArcCnt::const_pointer   forward_arc_const_pointer;

    typedef typename
        _RevArcCnt::value_type      reverse_arc_type;
    typedef typename
        _RevArcCnt::reference       reverse_arc_reference;
    typedef typename
        _RevArcCnt::const_reference reverse_arc_const_reference;
    typedef typename 
        _RevArcCnt::iterator        reverse_arc_iterator;
    typedef typename 
        _RevArcCnt::const_iterator  reverse_arc_const_iterator;
    typedef typename 
        _RevArcCnt::pointer         reverse_arc_pointer;
    typedef typename 
        _RevArcCnt::const_pointer   reverse_arc_const_pointer;

    typedef ptrdiff_t               difference_type;
    typedef size_t                  size_type;
};

} // namespace


#endif 
