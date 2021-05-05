/*
 ==========================================================================
 |   
 |   $Id: stringutils.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___STRINGUTILS_HXX__
#   define ___STRINGUTILS_HXX__

#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200)
#       pragma once
#       pragma warning(disable: 4018)
#       pragma warning(disable: 4146)
#   endif

#   include <optnet/_base/except.hxx>
#   include <optnet/_base/secure.hxx>
#   include <algorithm>
#   include <string>
#   include <locale>


/// @namespace optnet
namespace optnet { 
    /// @namespace optnet::utils
    namespace utils {

///////////////////////////////////////////////////////////////////////////
///  @class char_equal_no_case
///  @brief Case-insensitive char comparing functional.
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
class char_equal_no_case
{
public:
    bool operator()(const _CharType& c1, const _CharType& c2) const
    {
        return std::tolower(c1, std::locale::classic()) 
            == std::tolower(c2, std::locale::classic());
    }
};

///////////////////////////////////////////////////////////////////////////
///  Determines whether a string ends with another string.
///
///  @param[in] str1 A string.
///  @param[in] pc2  A pointer to a c-string whose value is to be searched
///                  for at the end of str1.
///
///  @return Returns true if the string str1 ends with the c-string
///          pointed to by pc2, false otherwise.
///
///  @remarks The function is case-sensitive.
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline bool
str_ends_with(const std::basic_string<_CharType>& str1,
              const _CharType* pc2)
{
    std::basic_string<_CharType> str2(pc2);
    return (str2.length() <= str1.length())
        && (std::equal(str1.end() - str2.length(), 
            str1.end(), str2.begin()));
}

///////////////////////////////////////////////////////////////////////////
///  Determines whether a string ends with another string.
///
///  @param[in] str1 A string.
///  @param[in] str2 A substring to be searched for at the end of str1.
///
///  @return Returns true if the string str1 ends with str2,
///          false otherwise.
///
///  @remarks The function is case-sensitive.
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline bool
str_ends_with(const std::basic_string<_CharType>& str1,
              const std::basic_string<_CharType>& str2)
{
    return (str2.length() <= str1.length())
        && (std::equal(str1.end() - str2.length(), 
            str1.end(), str2.begin()));
}

///////////////////////////////////////////////////////////////////////////
///  Determines whether a string ends with another string, ignoring case.
///
///  @param[in] pc1  A null-terminated string.
///  @param[in] pc2  A pointer to a c-string whose value is to be searched
///                  for at the end of str1.
///
///  @return Returns true if the string str1 ends with the c-string
///          pointed to by pc2, false otherwise.
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline bool
str_ends_with_no_case(const _CharType* pc1,
                      const _CharType* pc2)
{
    std::basic_string<_CharType> str1(pc1);
    std::basic_string<_CharType> str2(pc2);
    return (str2.length() <= str1.length())
        && (std::equal(str1.end() - str2.length(), 
            str1.end(), str2.begin(), char_equal_no_case<_CharType>()));
}

///////////////////////////////////////////////////////////////////////////
///  Determines whether a string ends with another string, ignoring case.
///
///  @param[in] str1 A string.
///  @param[in] pc2  A pointer to a c-string whose value is to be searched
///                  for at the end of str1.
///
///  @return Returns true if the string str1 ends with the c-string
///          pointed to by pc2, false otherwise.
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline bool
str_ends_with_no_case(const std::basic_string<_CharType>& str1,
                      const _CharType* pc2)
{
    std::basic_string<_CharType> str2(pc2);
    return (str2.length() <= str1.length())
        && (std::equal(str1.end() - str2.length(), 
            str1.end(), str2.begin(), char_equal_no_case<_CharType>()));
}

///////////////////////////////////////////////////////////////////////////
///  Determines whether a string ends with another string, ignoring case.
///
///  @param[in] str1 A string.
///  @param[in] str2 A substring to be searched for at the end of str1.
///
///  @return Returns true if the string str1 ends with str2,
///          false otherwise.
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline bool
str_ends_with_no_case(const std::basic_string<_CharType>& str1,
                      const std::basic_string<_CharType>& str2)
{
    return (str2.length() <= str1.length())
        && (std::equal(str1.end() - str2.length(), 
            str1.end(), str2.begin(), char_equal_no_case<_CharType>()));
}

///////////////////////////////////////////////////////////////////////////
///  Convert a string to a long integer.
///
///  @param str  The null-terminated string to convert.
///  @param base Number base to use (default: 10). 
///
///  @exception  std::invalid_argument Conversion unsuccessful.
///
///  @return The converted long integer.
///
///////////////////////////////////////////////////////////////////////////
inline long str_to_long(const char* str, int base = 10)
{
    char* endptr;
    long ans = strtol(str, &endptr, base);
    if (endptr == str || *endptr != '\0') {
        throw_exception(std::invalid_argument(
            "str_to_long: Could not convert the string to long."
            ));
    }
    return ans;
}

///////////////////////////////////////////////////////////////////////////
///  Convert a string to an unsigned long integer.
///
///  @param str  The null-terminated string to convert.
///  @param base Number base to use (default: 10). 
///
///  @exception  std::invalid_argument Conversion unsuccessful.
///
///  @return The converted unsigned long integer.
///
///////////////////////////////////////////////////////////////////////////
inline unsigned long str_to_unsigned_long(const char* str, int base = 10)
{
    char* endptr;
    unsigned long ans = strtoul(str, &endptr, base);
    if (endptr == str || *endptr != '\0') {
        throw_exception(std::invalid_argument(
            "str_to_long: Could not convert the string to unsigned long."
            ));
    }
    return ans;
}

///////////////////////////////////////////////////////////////////////////
///  Convert a string to a double-precision value.
///
///  @param str  The null-terminated string to convert.
///  @param base Number base to use. 
///
///  @exception  std::invalid_argument Conversion unsuccessful.
///
///  @return The converted double-precision value.
///
///////////////////////////////////////////////////////////////////////////
inline double str_to_double(const char* str)
{
    char* endptr;
    double ans = strtod(str, &endptr);
    if (endptr == str || *endptr != '\0') {
        throw_exception(std::invalid_argument(
            "str_to_double: Could not convert the string to double."
            ));
    }
    return ans;
}

#ifndef __CYGWIN__ // Too bad, my Cygwin does not support wcsto* functions.

///////////////////////////////////////////////////////////////////////////
///  Convert a string to a long integer.
///
///  @param str  The null-terminated string to convert.
///  @param base Number base to use (default: 10). 
///
///  @exception  std::invalid_argument Conversion unsuccessful.
///
///  @return The converted long integer.
///
///////////////////////////////////////////////////////////////////////////
inline long str_to_long(const wchar_t* str, int base = 10)
{
    wchar_t* endptr;
    long ans = wcstol(str, &endptr, base);
    if (endptr == str || *endptr != '\0') {
        throw_exception(std::invalid_argument(
            "str_to_long: Could not convert the string to long."
            ));
    }
    return ans;
}

///////////////////////////////////////////////////////////////////////////
///  Convert a string to an unsigned long integer.
///
///  @param str  The null-terminated string to convert.
///  @param base Number base to use (default: 10). 
///
///  @exception  std::invalid_argument Conversion unsuccessful.
///
///  @return The converted unsigned long integer.
///
///////////////////////////////////////////////////////////////////////////
inline unsigned long str_to_unsigned_long(const wchar_t* str, int base = 10)
{
    wchar_t* endptr;
    unsigned long ans = wcstoul(str, &endptr, base);
    if (endptr == str || *endptr != '\0') {
        throw_exception(std::invalid_argument(
            "str_to_long: Could not convert the string to unsigned long."
            ));
    }
    return ans;
}
///////////////////////////////////////////////////////////////////////////
///  Convert a string to a double-precision value.
///
///  @param str  The null-terminated string to convert.
///  @param base Number base to use. 
///
///  @exception  std::invalid_argument Conversion unsuccessful.
///
///  @return The converted double-precision value.
///
///////////////////////////////////////////////////////////////////////////
inline double str_to_double(const wchar_t* str)
{
    wchar_t* endptr;
    double ans = wcstod(str, &endptr);
    if (endptr == str || *endptr != '\0') {
        throw_exception(std::invalid_argument(
            "str_to_double: Could not convert the string to double."
            ));
    }
    return ans;
}
#endif

///////////////////////////////////////////////////////////////////////////
///  Converts any string to an Ansi (multibyte) string.
///
///  @param str  The null-terminated string to convert.
///
///  @return The converted Ansi string.
///
///////////////////////////////////////////////////////////////////////////
inline std::basic_string<char> str_to_ansi(const char* str)
{
    return std::basic_string<char>(str);
}

///////////////////////////////////////////////////////////////////////////
///  Converts any string to an Ansi (multibyte) string.
///
///  @param str  The null-terminated string to convert.
///
///  @return The converted Ansi string.
///
///////////////////////////////////////////////////////////////////////////
inline std::basic_string<char>
    str_to_ansi(const std::basic_string<char>& str)
{
    return str;
}

///////////////////////////////////////////////////////////////////////////
///  Converts any string to an Ansi (multibyte) string.
///
///  @param str  The wide (UTF-16) string to convert.
///
///  @exception  std::invalid_argument Conversion unsuccessful.
///
///  @return The converted Ansi string.
///
///////////////////////////////////////////////////////////////////////////
inline std::basic_string<char> str_to_ansi(const wchar_t* str)
{
    size_t size = secure_wcstombs(NULL, str, 0) + 2;
    char* mbsbuf = new char[size];
    secure_wcstombs(mbsbuf, str, size);
    std::basic_string<char> ansi(mbsbuf);
    delete [] mbsbuf;
    return ansi;
}

///////////////////////////////////////////////////////////////////////////
///  Converts any string to an Ansi (multibyte) string.
///
///  @param str  The wide (UTF-16) string to convert.
///
///  @exception  std::invalid_argument Conversion unsuccessful.
///
///  @return The converted Ansi string.
///
///////////////////////////////////////////////////////////////////////////
inline std::basic_string<char>
    str_to_ansi(const std::basic_string<wchar_t>& str)
{
    size_t size = secure_wcstombs(NULL, str.c_str(), 0) + 2;
    char* mbstr = new char[size];
    secure_wcstombs(mbstr, str.c_str(), size);
    std::basic_string<char> ansi(mbstr);
    delete [] mbstr;
    return ansi;
}

    } // namespace
} // namespace

#endif //___STRINGUTILS_HXX__
