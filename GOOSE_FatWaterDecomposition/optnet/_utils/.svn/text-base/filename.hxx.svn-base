/*
 ==========================================================================
 |   
 |   $Id: filename.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___FILENAME_HXX__
#   define ___FILENAME_HXX__

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4018)
#       pragma warning(disable: 4146)
#       pragma warning(disable: 4786)
#   endif

#   include <string>
#   include <optnet/config.h>
#   include <optnet/_base/secure.hxx>
#   include <optnet/_utils/stringutils.hxx>

/// @namespace optnet
namespace optnet { 
    /// @namespace optnet::utils
    namespace utils {

///////////////////////////////////////////////////////////////////////////
///  Extract file name without directory.
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline std::basic_string<_CharType>
extract_file_name(const _CharType* path)
{
    return extract_file_name(std::basic_string<_CharType>(path));
}

///////////////////////////////////////////////////////////////////////////
///  Extract file name without directory.
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline std::basic_string<_CharType>
extract_file_name(const std::basic_string<_CharType>& path)
{
    typename std::basic_string<_CharType>::size_type pos1, pos2;

#   ifdef __OPTNET_OS_WINNT__
    pos1 = path.rfind((_CharType)('\\'));
    pos2 = path.rfind((_CharType)('/' ));
    if (
      (std::basic_string<_CharType>::npos != pos2)
    && ((pos2 > pos1) || (std::basic_string<_CharType>::npos == pos1)))
        pos1 = pos2;
    pos2 = path.rfind((_CharType)(':' ));
    if (
      (std::basic_string<_CharType>::npos != pos2)
    && ((pos2 > pos1) || (std::basic_string<_CharType>::npos == pos1)))
        pos1 = pos2;
#   else // Unixes
    pos1 = path.rfind((_CharType)('/' ));
    pos2 = 0; // unused
#   endif
    if (std::basic_string<_CharType>::npos != pos1)
        return path.substr(pos1 + 1);
    // else
    return path;
}

///////////////////////////////////////////////////////////////////////////
///  Extract file basename, including directory
///    (e.g., ./myfile.txt becomes ./myfile).
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline std::basic_string<_CharType>
extract_file_base(const _CharType* path)
{
    return extract_file_base(std::basic_string<_CharType>(path));
}

///////////////////////////////////////////////////////////////////////////
///  Extract file basename, including directory
///    (e.g., ./myfile.txt becomes ./myfile).
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline std::basic_string<_CharType>
extract_file_base(const std::basic_string<_CharType>& path)
{
    typename std::basic_string<_CharType>::size_type pos, len;

#   ifdef __OPTNET_OS_WINNT__
    len = path.length();
    if (len == 0 || path[len-1] == '\\' ||
        path[len-1] == '/' || path[len-1] == ':')
        return std::basic_string<_CharType>();
#   else // Unixes
    len = path.length();
    if (len == 0 || path[len-1] == '/')
        return std::basic_string<_CharType>();
#   endif

    pos = path.rfind((_CharType)('.'));
    if (std::basic_string<_CharType>::npos != pos) {
        return path.substr(0, pos);
    }
    // else
    return path;
}

///////////////////////////////////////////////////////////////////////////
///  Extract file directory.
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline std::basic_string<_CharType>
extract_file_dir (const _CharType* path)
{
    return extract_file_dir (std::basic_string<_CharType>(path));
}

///////////////////////////////////////////////////////////////////////////
///  Extract file directory.
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline std::basic_string<_CharType>
extract_file_dir (const std::basic_string<_CharType>& path)
{
    typename std::basic_string<_CharType>::size_type pos1, pos2;

#   ifdef __OPTNET_OS_WINNT__
    pos1 = path.rfind((_CharType)('\\'));
    pos2 = path.rfind((_CharType)('/' ));
    if (
      (std::basic_string<_CharType>::npos != pos2)
    && ((pos2 > pos1) || (std::basic_string<_CharType>::npos == pos1)))
        pos1 = pos2;
    pos2 = path.rfind((_CharType)(':' ));
    if (
      (std::basic_string<_CharType>::npos != pos2)
    && ((pos2 > pos1) || (std::basic_string<_CharType>::npos == pos1)))
        pos1 = pos2;
#   else // Unixes
    pos1 = path.rfind((_CharType)('/' ));
    pos2 = 0; // unused
#   endif
    if (std::basic_string<_CharType>::npos != pos1)
        return path.substr(0, pos1 + 1);
    // else
    return std::basic_string<_CharType>();
}

///////////////////////////////////////////////////////////////////////////
///  Extract file extension.
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline std::basic_string<_CharType>
extract_file_ext (const _CharType* path)
{
    return extract_file_ext (std::basic_string<_CharType>(path));
}

///////////////////////////////////////////////////////////////////////////
///  Extract file extension.
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline std::basic_string<_CharType>
extract_file_ext (const std::basic_string<_CharType>& path)
{
    typename std::basic_string<_CharType>::size_type pos;

    pos = path.rfind((_CharType)('.'));
    if (std::basic_string<_CharType>::npos != pos)
        return path.substr(pos + 1);
    // else
    return std::basic_string<_CharType>();
}

///////////////////////////////////////////////////////////////////////////
/// Generate a serially numbered filename, e.g. pix001.png
///////////////////////////////////////////////////////////////////////////
inline std::string
make_serial_filename(const std::string& filename,
                     unsigned int       index,
                     unsigned int       total,
                     const char*        ext
                     )
{
    std::string name;

    if ((total > 0) && (index <= total)) {
        char tmp1[32], tmp2[32];

        if (str_ends_with_no_case(filename, ext))
            name = extract_file_base(filename);
        else
            name = filename;

        if (total > 1) {
            secure_sprintf(tmp1, sizeof(tmp1), "%d", total);
            secure_sprintf(tmp1, sizeof(tmp1), "%%0%dd", strlen(tmp1));
            secure_sprintf(tmp2, sizeof(tmp2), tmp1, index); 
            name += tmp2;
        }

        name += ext;
    }

    return name;
}


    } // namespace
} // namespace

#endif
