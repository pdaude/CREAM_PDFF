/*
 ==========================================================================
 |   
 |   $Id: define.h 2137 2007-07-02 03:26:31Z kangli $
 |
 |   OptimalNet Library Standard C Interface Header
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


#ifndef ___DEFINE_H___
#   define ___DEFINE_H___

/** \file define.h
 *  The header file containing special type definitions for the library.
 */

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   include <stddef.h>

/* Define missing math constants. */
#   ifndef M_E
#       define M_E              2.7182818284590452353602874713527
#   endif
#   ifndef M_LOG10E
#       define M_LOG10E         0.4342944819032518276511289189165
#   endif
#   ifndef M_LN2
#       define M_LN2            0.6931471805599453094172321214582
#   endif
#   ifndef M_LN10
#       define M_LN10           2.3025850929940456840179914546844
#   endif
#   ifndef M_PI
#       define M_PI             3.1415926535897932384626433832795
#   endif
#   ifndef M_TWOPI
#       define M_TWOPI          (M_PI * 2.0)
#   endif
#   ifndef M_PI_2
#       define M_PI_2           1.5707963267948966192313216916398
#   endif
#   ifndef M_PI_4
#       define M_PI_4           0.7853981633974483096156608458199
#   endif
#   ifndef M_3PI_4
#       define M_3PI_4          2.3561944901923449288469825374596
#   endif
#   ifndef M_1_PI
#       define M_1_PI           0.3183098861837906715377675267450
#   endif
#   ifndef M_2_PI
#       define M_2_PI           0.6366197723675813430755350534901
#   endif
#   ifndef M_SQRTPI
#       define M_SQRTPI         1.7724538509055160272981674833411
#   endif
#   ifndef M_SQRTTWOPI
#       define M_SQRTTWOPI      2.5066282746310005024157652848110
#   endif
#   ifndef M_SQRTPI_2
#       define M_SQRTPI_2       1.2533141373155002512078826424055
#   endif
#   ifndef M_2_SQRTPI
#       define M_2_SQRTPI       1.1283791670955125738961589031215
#   endif
#   ifndef M_SQRT1_2
#       define M_SQRT1_2        0.7071067811865475244008443621049
#   endif
#   ifndef M_SQRT2
#       define M_SQRT2          1.4142135623730950488016887242097
#   endif
#   ifndef M_SQRT3
#       define M_SQRT3          1.7320508075688772935274463415059
#   endif

/**
  * \brief The definition of an ROI node, specifying the range of voxels
  *        in the region-of-interest for a column of voxels.
  *
  *        For example, assume the columns are parallel to the z-axis,
  *        and are denoted Col(x, y). The roi_node corresponding to
  *        Col(x, y) is Rn(x, y). Then a voxel I(x, y, z) is in
  *        the ROI if Rn(x, y).lower <= z < Rn(x, y).upper.
  *
  */
typedef struct roi_node
{
    size_t  lower;  /** Lower bound of the column */
    size_t  upper;  /** Upper bound of the column */
} roi_node_t,
*proi_node_t;

#endif
