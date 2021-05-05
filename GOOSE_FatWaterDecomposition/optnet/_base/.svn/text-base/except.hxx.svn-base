/*
 ==========================================================================
 |   
 |   $Id: except.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___EXCEPT_HXX___
#   define ___EXCEPT_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <optnet/config.h>

#   ifndef ___OPTNET_NO_EXCEPTIONS__
#       include <stdexcept>
#       include <string>
#   endif

///  @namespace optnet
namespace optnet {

#   ifdef ___OPTNET_NO_EXCEPTIONS__

    template<class _E>
        inline void throw_exception(const _E& e) { /*dummy*/}

#   else // !___OPTNET_NO_EXCEPTIONS__

    template<class _E>
        inline void throw_exception(const _E& e) { throw e; }

#   endif

#   ifndef ___OPTNET_NO_EXCEPTIONS__
    /// @namespace optnet::io
    namespace io {

        ///////////////////////////////////////////////////////////////////
        ///  @class io_error
        ///  @brief The class serves as the base class for all exceptions
        ///         thrown to report an I/O error.
        ///////////////////////////////////////////////////////////////////
        class io_error: public std::runtime_error
        {
        public:
            io_error(const std::string& message) :
                std::runtime_error(message)
            {
            }
        };
    } // namespace
#   endif

inline void
runtime_assert(const char* expr, const char* file, unsigned int line)
{
    std::string errmsg("Assertion failed: ");
    errmsg.append(expr);
    if (NULL != file) {
        const char* format = "\nFile: %s\nLine: %d";
        char  buffer[4096];
#   if defined(__OPTNET_SECURE_CRT__) // secure version
        sprintf_s(buffer, sizeof(buffer), format, file, line);
#   else
        sprintf(buffer, format, file, line);
#   endif
        errmsg.append(buffer);
    }
#   ifndef ___OPTNET_NO_EXCEPTIONS__
    throw_exception(std::runtime_error(errmsg.c_str()));
#   else
    //FIXME
#   endif
}

} // namespace

/* ===================================================================== */
/*   EXCEPTION HANDLING MACROS                                           */
/* ===================================================================== */

#   ifdef RUNTIME_ASSERT
#       undef RUNTIME_ASSERT
#   endif

#   ifdef  NDEBUG
#       define RUNTIME_ASSERT(expr) \
        (void)((expr) || (runtime_assert(#expr, NULL, 0), 0))
#   else
#       define RUNTIME_ASSERT(expr) \
        (void)((expr) || (runtime_assert(#expr, __FILE__, __LINE__), 0))
#   endif

#   ifndef ___OPTNET_NO_EXCEPTIONS__
#       define __OPTNET_TRY try
#       define __OPTNET_CATCH_AND_RETURN            \
            catch (optnet::io::io_error&) {         \
                return OPTNET_E_IO_ERROR;           \
            }                                       \
            catch (std::invalid_argument&) {        \
                return OPTNET_E_INVALID_ARG;        \
            }                                       \
            catch (std::runtime_error&) {           \
                return OPTNET_E_RUNTIME_ERROR;      \
            }                                       \
            catch (std::exception&) {               \
                return OPTNET_E_FAILED;             \
            }
#   else
#       define __OPTNET_TRY
#       define __OPTNET_CATCH_AND_RETURN
#   endif

#endif 
