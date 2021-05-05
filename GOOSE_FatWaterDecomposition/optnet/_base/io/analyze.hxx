/*
 ==========================================================================
 |   
 |   $Id: analyze.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___ANALYZE_HXX__
#   define ___ANALYZE_HXX__

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   include <optnet/_base/except.hxx>
#   include <optnet/_base/io/analyze_core.hxx>

/// @namespace optnet
namespace optnet {
    
    /// @namespace optnet::io
    /// @brief     The namespace that contains data I/O facilities.
    namespace io {


///////////////////////////////////////////////////////////////////////////
///  A template class to obtain the numerical id used by Analyze 7.5
///  format for the specified datatype.
///////////////////////////////////////////////////////////////////////////
template <typename _T>
struct analyze_datatype // This should never be called.
{ static long id() { assert(0); return DT_NONE; }};

///////////////////////////////////////////////////////////////////////////
///  Specialization for type unsigned char.
///////////////////////////////////////////////////////////////////////////
template <>
struct analyze_datatype<unsigned char>
{ static long id() { return DT_UNSIGNED_CHAR;   }};

///////////////////////////////////////////////////////////////////////////
///  Specialization for type signed short.
///////////////////////////////////////////////////////////////////////////
template <>
struct analyze_datatype<signed short>
{ static long id() { return DT_SIGNED_SHORT;    }};

///////////////////////////////////////////////////////////////////////////
///  Specialization for type signed short.
///////////////////////////////////////////////////////////////////////////
template <>
struct analyze_datatype<unsigned short>
{ static long id() { return DT_UNSIGNED_SHORT;  }};

///////////////////////////////////////////////////////////////////////////
///  Specialization for type signed int.
///////////////////////////////////////////////////////////////////////////
template <>
struct analyze_datatype<signed int>
{ static long id() { return DT_SIGNED_INT;      }};

///////////////////////////////////////////////////////////////////////////
///  Specialization for type float.
///////////////////////////////////////////////////////////////////////////
template <>
struct analyze_datatype<float>
{ static long id() { return DT_FLOAT;           }};

///////////////////////////////////////////////////////////////////////////
///  Specialization for type double.
///////////////////////////////////////////////////////////////////////////
template <>
struct analyze_datatype<double>
{ static long id() { return DT_DOUBLE;          }};


///////////////////////////////////////////////////////////////////////////
///  Determines the datatype that the specified Analyze(TM) file uses.
///
///  @param[in]  name   The name of the image file.
///
///  @return     The numerical identifier of the datatype of the image.
///
///  @exception  optnet::io_error Failed reading image header.
///
///////////////////////////////////////////////////////////////////////////
inline long
analyze_get_datatype(const char* name)
{
    using namespace optnet::io::detail;

    analyze_volume  volume;
    std::string     hdr_name, img_name;

    assert(0 != name);
    
    analyze_make_filename_pair(name, hdr_name, img_name);

    // Load Analyze image header.
    if (!analyze_load_info(&volume, hdr_name, img_name)) {
        throw_exception(io_error(
            "analyze_get_datatype: Failed reading image header."
            ));
    }

    return (long)volume.datatype;
}

///////////////////////////////////////////////////////////////////////////
///  Loads an Analyze(TM) 7.5 format image file into memory.
///
///  @param[out] a      The a for storing the image data.
///  @param[in]  name   The name of the image file.
///
///  @exception  optnet::io_error Failed reading image header or data, or
///                               the data type of a does not match
///                               that used in the file.
///
///  @remarks This function outputs the image data only, ignoring any
///           extended information such as patient identity, medical
///           history or scanning protocols.
///
///////////////////////////////////////////////////////////////////////////
template <typename _Array>
void
analyze_load(_Array& a, const char* name)
{
    using namespace optnet::io::detail;

    analyze_volume  volume;
    std::string     hdr_name, img_name;

    assert(0 != name);
    
    analyze_make_filename_pair(name, hdr_name, img_name);

    // Load Analyze image header.
    if (!analyze_load_info(&volume, hdr_name, img_name)) {
        throw_exception(io_error(
            "analyze_load: Failed reading image header."
            ));
    }

    // Check data type of the image.
#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200) // VC 6.0

    if (volume.datatype != 
        analyze_datatype</*typename*/ _Array::value_type>::id()) {
        throw_exception(io_error(
            "analyze_load: Data type does not match."
            ));
    }
    
#   else

    if (volume.datatype != 
        analyze_datatype<  typename   _Array::value_type>::id()) {
        throw_exception(io_error(
            "analyze_load: Data type does not match."
            ));
    }
    
#   endif

    // Allocate storage for the image data.
    a.create(volume.image_size_x,
             volume.image_size_y,
             volume.image_size_z,
             volume.image_size_w,
             1
             );
    a.extension().voxel_size[0] = volume.voxel_size_x;
    a.extension().voxel_size[1] = volume.voxel_size_y;
    a.extension().voxel_size[2] = volume.voxel_size_z;
    volume.data = reinterpret_cast<char*>(a.data());

    
    // Load image data into the allocated memory.
    if (!analyze_load_data(&volume, hdr_name, img_name)) {
        throw_exception(io_error(
            "analyze_load: Failed reading image."
            ));
    }

}

///////////////////////////////////////////////////////////////////////////
///  Saves an image into an Analyze(TM) 7.5 format image file.
///
///  @param[in] a      The a for storing the image data.
///  @param[in] name   The name of the image file.
///  @param[in] multi  If true, a multi-phase (3D+time) image will be
///                    saved as multiple separate single-phase files;
///                    if false, the image will be saved in a single
///                    (4-D) file.
///
///  @exception optnet::io_error Failed writing image header or data.
///
///  @remarks This function saves the image data only, ignoring any
///           extended information such as patient identity, voxel
///           dimensions, medical history or scanning properties.
///
///////////////////////////////////////////////////////////////////////////
template <typename _Array>
void
analyze_save(const _Array&  a,
             const char*    name,
             bool           multi = false,
             int            compression = 6
             )
{
    using namespace optnet::io::detail;

    if (a.size() == 0) return; // Nothing to write.

    analyze_volume  volume;
    std::string     hdr_name, img_name;

    memset(&volume, 0, sizeof(analyze_volume));

#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200) // VC 6.0

    volume.datatype 
        = analyze_datatype</*typename*/ _Array::value_type>::id();
        
#   else

    volume.datatype 
        = analyze_datatype<  typename   _Array::value_type>::id();
        
#   endif

    volume.voxel_bytes
        = sizeof(typename _Array::value_type);

    volume.image_size_x = a.sizes()[0];
    volume.image_size_y = a.sizes()[1];
    volume.image_size_z = a.sizes()[2];

    volume.voxel_size_x = a.extension().voxel_size[0];
    volume.voxel_size_y = a.extension().voxel_size[1];
    volume.voxel_size_z = a.extension().voxel_size[2];

    if (multi) {

        volume.image_size_w = 1;
        volume.image_size   = volume.image_size_x * 
                              volume.image_size_y * 
                              volume.image_size_z;

        for (size_t i = 0; i < a.sizes()[3]; ++i) {
            
            analyze_make_filename_pair(name, hdr_name, img_name, (int)(i + 1));
            
            if (!analyze_save_info(&volume, hdr_name, img_name)) {
                throw_exception(io_error(
                    "analyze_save: Failed writing image header."
                    ));
            }

            volume.data = const_cast<char*>
                            (reinterpret_cast<const char*>
                                (a.data() + volume.image_size * i));
            
            if (!analyze_save_data(&volume, hdr_name, img_name, compression)) {
                throw_exception(io_error(
                    "analyze_save: Failed writing image."
                    ));
            }

        } // for
    }
    else {

        volume.image_size_w = a.sizes()[3];
        volume.image_size   = volume.image_size_x * 
                              volume.image_size_y * 
                              volume.image_size_z *
                              volume.image_size_w;

        analyze_make_filename_pair(name, hdr_name, img_name);

        if (!analyze_save_info(&volume, hdr_name, img_name)) {
            throw_exception(io_error(
                "analyze_save: Failed writing image header."
                ));
        }

        volume.data = const_cast<char*>
                        (reinterpret_cast<const char*>(a.data()));
        
        if (!analyze_save_data(&volume, hdr_name, img_name, compression)) {
            throw_exception(io_error(
                "analyze_save: Failed writing image."
                ));
        }

    }
    
}

    } // namespace
} // namespace

#endif
