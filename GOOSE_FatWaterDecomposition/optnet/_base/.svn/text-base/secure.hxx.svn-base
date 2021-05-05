/*
 ==========================================================================
 |   
 |   $Id: secure.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___SECURE_HXX___
#   define ___SECURE_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   include <optnet/config.h>
#   include <assert.h>      // for assert
#   include <string.h>      // for strtok, etc.
#   include <stdlib.h>      // for wcstombs, etc.
#   include <stdarg.h>      // for va_list, etc.
#   include <stdio.h>       // for sprintf, etc.
#   include <errno.h>       // for EINVAL, etc.

/// @namespace optnet
namespace optnet {

//
// Macros
//
#   ifdef __OPTNET_UNICODE__
#       if defined(__OPTNET_SECURE_CRT__) // secure version
#           define secure_sscanf swscanf_s
#       else // insecure version
#           define secure_sscanf swscanf
#       endif
#       define secure_strchr wcschr
#       define secure_strcmp wcscmp
#   else // !Unicode
#       if defined(__OPTNET_SECURE_CRT__) // secure version
#           define secure_sscanf sscanf_s
#       else // insecure version
#           define secure_sscanf sscanf
#       endif
#       define secure_strchr strchr
#       define secure_strcmp strcmp
#   endif

//
// Functions
//

///////////////////////////////////////////////////////////////////////////
///  More secure version of fopen.
///////////////////////////////////////////////////////////////////////////
inline int secure_fopen(FILE**      ppfile,
                        const char* filename,
                        const char* mode
                        )
{
#   if defined(__OPTNET_SECURE_CRT__) // secure version
    return (int)fopen_s(ppfile, filename, mode);
#   else // insecure version
    if (NULL == ppfile) return EINVAL;
    *ppfile = fopen(filename, mode);
    return (NULL == *ppfile) ? 
        EINVAL : 0;
#   endif
}

///////////////////////////////////////////////////////////////////////////
///  More secure version of strtok.
///////////////////////////////////////////////////////////////////////////
inline char* secure_strtok(char*       token,
                           const char* delimit,
                           char**      context
                           )
{
#   if defined(__OPTNET_SECURE_CRT__) // secure version
    return strtok_s(token, delimit, context);
#   else // insecure version
    OPTNET_UNUSED(context); // unused
    return strtok(token, delimit);
#   endif
}

///////////////////////////////////////////////////////////////////////////
///  More secure version of vsprintf.
///////////////////////////////////////////////////////////////////////////
inline int secure_vsprintf(char*       buffer,
                           size_t      size_in_bytes,
                           const char* format,
                           va_list     argptr
                           )
{
#   if defined(__OPTNET_SECURE_CRT__) // secure version
    return vsprintf_s(buffer, size_in_bytes, format, argptr);
#   elif defined(__OPTNET_HAS_VSNPRINTF__) // secure version
    return vsnprintf(buffer, size_in_bytes, format, argptr);
#   else // insecure version
    OPTNET_UNUSED(size_in_bytes);
    return vsprintf(buffer, format, argptr);
#   endif
}

///////////////////////////////////////////////////////////////////////////
///  More secure version of vsprintf.
///////////////////////////////////////////////////////////////////////////
inline int secure_vswprintf(wchar_t*       buffer,
                            size_t         size_in_bytes,
                            const wchar_t* format,
                            va_list        argptr
                            )
{
#   if defined(__OPTNET_SUPPORT_UTF16__)
#       if defined(__OPTNET_SECURE_CRT__) // secure version
            return vswprintf_s(buffer, size_in_bytes, format, argptr);
#       elif defined(__OPTNET_HAS_VSNPRINTF__) // secure version
            return vsnwprintf(buffer, size_in_bytes, format, argptr);
#       else // insecure version
            OPTNET_UNUSED(size_in_bytes);
            return vswprintf(buffer, format, argptr);
#       endif
#   else // __OPTNET_SUPPORT_UTF16__
        assert(false);
        return 0;
#   endif
}

///////////////////////////////////////////////////////////////////////////
///  More secure version of sprintf.
///////////////////////////////////////////////////////////////////////////
inline int secure_sprintf(char*       buffer,
                          size_t      size_in_bytes,
                          const char* format,
                          ...
                          )
{
    va_list argptr;
    va_start(argptr, format);
    int ret = secure_vsprintf(buffer, size_in_bytes, format, argptr);
    va_end(argptr);
    return ret;
}

///////////////////////////////////////////////////////////////////////////
///  More secure version of sprintf.
///////////////////////////////////////////////////////////////////////////
inline int secure_swprintf(wchar_t*       buffer,
                           size_t         size_in_bytes,
                           const wchar_t* format,
                           ...
                           )
{
    va_list argptr;
    va_start(argptr, format);
    int ret = secure_vswprintf(buffer, size_in_bytes, format, argptr);
    va_end(argptr);
    return ret;
}

///////////////////////////////////////////////////////////////////////////
///  More secure version of wcstombs.
///////////////////////////////////////////////////////////////////////////
inline size_t secure_wcstombs(char*          mbstr,
                              const wchar_t* wcstr,
                              size_t         count
                              )
{
#   if defined(__OPTNET_SECURE_CRT__) // secure version
    size_t converted = 0;
    wcstombs_s(&converted, mbstr, count, wcstr, count);
    return converted;
#   else // insecure version
    return wcstombs(mbstr, wcstr, count);
#   endif
}

///////////////////////////////////////////////////////////////////////////
///  More secure version of mbstowcs.
///////////////////////////////////////////////////////////////////////////
inline size_t secure_mbstowcs(wchar_t*    wcstr,
                              const char* mbstr,
                              size_t      count
                              )
{
#   if defined(__OPTNET_SECURE_CRT__) // secure version
    size_t converted = 0;
    mbstowcs_s(&converted, wcstr, count, mbstr, count);
    return converted;
#   else // insecure version
    return mbstowcs(wcstr, mbstr, count);
#   endif
}

#   ifdef __OPTNET_UNICODE__
#       define secure_stprintf secure_swprintf
#       define secure_vstprintf secure_vswprintf
#   else
#       define secure_stprintf secure_sprintf
#       define secure_vstprintf secure_vsprintf
#   endif

} // namespace

#endif // ___SECURE_HXX___
