/*
 ==========================================================================
 |   
 |   $Id: tags.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___TAGS_HXX___
#   define ___TAGS_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class net_f_xy
///  @brief Surface direction tag.
///
///  This indicates that the surface can be expressed as z = f(x, y).
///////////////////////////////////////////////////////////////////////////
class net_f_xy {};

///////////////////////////////////////////////////////////////////////////
///  @class net_f_yz
///  @brief Surface direction tag.
///
///  This indicates that the surface can be expressed as x = f(y, z).
///////////////////////////////////////////////////////////////////////////
class net_f_yz {};

///////////////////////////////////////////////////////////////////////////
///  @class net_f_zx
///  @brief Surface direction tag.
///
///  This indicates that the surface can be expressed as y = f(z, x).
///////////////////////////////////////////////////////////////////////////
class net_f_zx {};

///////////////////////////////////////////////////////////////////////////
///  @class net_f_xy_zflipped
///  @brief Surface direction tag.
///
///  This indicates that the surface can be written as z = Z - 1 - f(x, y).
///////////////////////////////////////////////////////////////////////////
class net_f_xy_zflipped {};

///////////////////////////////////////////////////////////////////////////
///  @class net_f_yz_xflipped
///  @brief Surface direction tag.
///
///  This indicates that the surface can be written as x = X - 1 - f(y, z).
///////////////////////////////////////////////////////////////////////////
class net_f_yz_xflipped {};

///////////////////////////////////////////////////////////////////////////
///  @class net_f_zx_yflipped
///  @brief Surface direction tag.
///
///  This indicates that the surface can be written as y = Y - 1 - f(z, x).
///////////////////////////////////////////////////////////////////////////
class net_f_zx_yflipped {};

} // namespace


#endif 
