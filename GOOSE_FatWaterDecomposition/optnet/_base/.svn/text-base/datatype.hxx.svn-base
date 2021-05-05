/*
 ==========================================================================
 |   
 |   $Id: datatype.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___DATATYPE_HXX___
#   define ___DATATYPE_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <cstddef>

/// @namespace optnet
namespace optnet {

/// @enum datatype_id
enum datatype_id
{
    OPTNET_DT_UNKNOWN           = 0x00,
    // integer datatypes
    OPTNET_DT_CHAR              = 0x01,
    OPTNET_DT_UNSIGNED_CHAR     = 0x02,
    OPTNET_DT_SHORT             = 0x03,
    OPTNET_DT_UNSIGNED_SHORT    = 0x04,
    OPTNET_DT_INT               = 0x05,
    OPTNET_DT_UNSIGNED_INT      = 0x06,
    OPTNET_DT_LONG              = 0x07,
    OPTNET_DT_UNSIGNED_LONG     = 0x08,
    OPTNET_DT_INT64             = 0x09, //
    OPTNET_DT_UNSIGNED_INT64    = 0x0A, //
    OPTNET_DT_LONGLONG          = 0x0B, // not formally supported
    OPTNET_DT_UNSIGNED_LONGLONG = 0x0C, //
    // floating point datatypes
    OPTNET_DT_FLOAT             = 0x11,
    OPTNET_DT_DOUBLE            = 0x12,
    OPTNET_DT_LONG_DOUBLE       = 0x13
};

///////////////////////////////////////////////////////////////////////////
///  A template class to obtain the numerical id of datatypes.
///////////////////////////////////////////////////////////////////////////
template <typename _T>
struct datatype // This should never be called.
{ static long id() { assert(0); return OPTNET_DT_UNKNOWN; }};

template <> struct datatype<char>
{ static long id() { return OPTNET_DT_CHAR; }};

template <> struct datatype<unsigned char>
{ static long id() { return OPTNET_DT_UNSIGNED_CHAR; }};

template <> struct datatype<short>
{ static long id() { return OPTNET_DT_SHORT; }};

template <> struct datatype<unsigned short>
{ static long id() { return OPTNET_DT_UNSIGNED_SHORT; }};

template <> struct datatype<int>
{ static long id() { return OPTNET_DT_INT; }};

template <> struct datatype<unsigned int>
{ static long id() { return OPTNET_DT_UNSIGNED_INT; }};

template <> struct datatype<long>
{ static long id() { return OPTNET_DT_LONG; }};

template <> struct datatype<unsigned long>
{ static long id() { return OPTNET_DT_UNSIGNED_LONG; }};

template <> struct datatype<float>
{ static long id() { return OPTNET_DT_FLOAT; }};

template <> struct datatype<double>
{ static long id() { return OPTNET_DT_DOUBLE; }};

template <> struct datatype<long double>
{ static long id() { return OPTNET_DT_LONG_DOUBLE; }};

} // namespace

#endif 
