/*
 ==========================================================================
 |   
 |   $Id: config.h 2137 2007-07-02 03:26:31Z kangli $
 |
 |   OptimalNet Library Platform-Specific Configurations
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

#ifndef ___CONFIG_H___
#   define ___CONFIG_H___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   if (defined(__unix__) || defined(unix)) && !defined(USG)
#       include <sys/param.h>
#   endif

/* operating system recognition */
#   if defined(linux) || defined(__linux) || defined(__linux__)
        /* Linux */
#       define __OPTNET_OS_LINUX__
#   elif !defined(SAG_COM) \
    && (defined(WIN64) || defined(_WIN64) || defined(__WIN64__))
        /* Win32 or Win64 */
#       define __OPTNET_OS_WIN32__
#       define __OPTNET_OS_WIN64__
#       define __OPTNET_OS_WINNT__
#   elif !defined(SAG_COM) \
    && ( \
        defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined \
            (__NT__) \
        )
#       define __OPTNET_OS_WIN32__
#       define __OPTNET_OS_WINNT__
#   elif defined(__hpux)
        /* HP Unix */
#       define __OPTNET_OS_HPUX__
#   elif defined(BSD) || defined(__FreeBSD__)
        /* BSD */
#       define __OPTNET_OS_BSD__
#   elif defined(_MAC) || defined(__APPLE__)
        /* Mac */
#       define __OPTNET_OS_MAC__
#   endif

#   if defined(WIN64) || defined(_WIN64)
        /* Win64 */
#       define __OPTNET_OS_WIN64__
#   endif

/* compiler specific configuration */
#   if defined(__GNUC__)
#       define __OPTNET_CC_GNUC__
#       define __OPTNET_CC_GNUC_VER__
#   elif defined(__INTEL_COMPILER) || defined(__ICL) || defined(__ECL)
#       define __OPTNET_CC_INTEL__
#       define __OPTNET_CC_INTEL_VER__
#   elif defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#       define __OPTNET_CC_MSC__
#       define __OPTNET_CC_MSC_VER__   _MSC_VER
#   elif defined(__LCC__)
#       define __OPTNET_CC_LCC__
#       define __OPTNET_CC_LCC_VER__
#   else
#       define __OPTNET_CC_UNKNOWN__
#   endif

#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200) // VC 6.0
#       define __OPTNET_CRAPPY_MSC__
#   endif

#   if defined(_MSC_VER) && (_MSC_VER >= 1400) // VC 8.0
#       define __OPTNET_SECURE_CRT__
#       if defined(_OPENMP)
#           define __OPTNET_PRAGMA_OMP__
#           define __OPTNET_OMP_NUM_THREADS__ 4
#       endif
#   endif

#   ifndef OPTNET_IMPEXP
#      ifdef OPTNET_EXPORTS
#         define OPTNET_IMPEXP __declspec(dllexport)
#      else
#         define OPTNET_IMPEXP __declspec(dllimport)
#      endif
#   endif

/* API calling convention */
#   if defined(__OPTNET_CC_MSC__) || defined(__OPTNET_CC_INTEL__) \
        || defined(__OPTNET_CC_LCC__)
#       define OPTNET_API /* OPTNET_IMPEXP */
#   else
#       define OPTNET_API
#   endif

/* debugging facilities */
#   if defined(_DEBUG) || defined(DEBUG) && !defined(NDEBUG) \
    || (defined(__LCC__) && __LCCDEBUGLEVEL >= 2)
#       define __OPTNET_DEBUG__ 1
#   endif

/* compiler features */
#   if defined(__OPTNET_CC_MSC__)
#       if (__OPTNET_CC_MSC_VER__ > 1200)
#           define __OPTNET_MEMBER_TEMPLATES__
#       endif
#   else
#       define __OPTNET_MEMBER_TEMPLATES__
#       define __OPTNET_REMOVE_UNUSED_ARG__
#   endif

#   if defined(__OPTNET_REMOVE_UNUSED_ARG__)
#       define OPTNET_UNUSED(arg) /* arg */
#   else  // stupid, broken compiler
#       define OPTNET_UNUSED(arg) arg
#   endif

/* library features */
#   define __OPTNET_SUPPORT_ROI__

#   if defined(__OPTNET_OS_WINNT__) && \
            (defined(__OPTNET_CC_MSC__) || defined(__OPTNET_CC_INTEL__))
#       define __OPTNET_USE_MSXML__      // use MSXML
#       define __OPTNET_SUPPORT_UTF16__  // support unicode
#   endif

#endif
