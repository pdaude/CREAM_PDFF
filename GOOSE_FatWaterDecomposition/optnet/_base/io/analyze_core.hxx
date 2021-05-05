/*
 ==========================================================================
 |   
 |   $Id: analyze_core.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___ANALYZE_CORE_HXX__
#   define ___ANALYZE_CORE_HXX__

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       if defined(_WIN64)  // Win64
#           if defined(_DEBUG)
#               if defined(_DLL)
#                   pragma comment(lib, "zlib64MDd.lib")
#               else
#                   pragma comment(lib, "zlib64MTd.lib")
#               endif
#           else // !_DEBUG
#               if defined(_DLL)
#                   pragma comment(lib, "zlib64MD.lib")
#               else
#                   pragma comment(lib, "zlib64MT.lib")
#               endif
#           endif
#       else // Win32
#           if defined(_DEBUG)
#               if defined(_MT)
#                   if defined(_DLL)
#                       pragma comment(lib, "zlibMDd.lib")
#                   else
#                       pragma comment(lib, "zlibMTd.lib")
#                   endif
#               else
#                   pragma comment(lib, "zlibd.lib")
#               endif
#           else // !_DEBUG
#               if defined(_MT)
#                   if defined(_DLL)
#                       pragma comment(lib, "zlibMD.lib")
#                   else
#                       pragma comment(lib, "zlibMT.lib")
#                   endif
#               else
#                   pragma comment(lib, "zlib.lib")
#               endif
#           endif
#       endif
#   endif

#include <cstdio>
#include <cassert>
#include <optnet/_base/io/detail/zlib.h>
#include <optnet/_base/io/detail/analyze.h>
#include <optnet/_utils/stringutils.hxx>
#include <optnet/_utils/endian.hxx>

//
// namespace optnet::io::detail
//
namespace optnet { namespace io { namespace detail {

///////////////////////////////////////////////////////////////////////////
//  A simple struct to store image data and information.
///////////////////////////////////////////////////////////////////////////
struct analyze_volume
{
    long    datatype;
    bool    swap_endian;
    double  voxel_size_x;
    double  voxel_size_y;
    double  voxel_size_z;
    size_t  voxel_bytes;
    size_t  image_size_x;
    size_t  image_size_y;
    size_t  image_size_z;
    size_t  image_size_w;
    size_t  image_size;
    char*   data;
};

///////////////////////////////////////////////////////////////////////////
inline void
analyze_make_filename_pair(const std::string& name,
                           std::string&       hdr_name,
                           std::string&       img_name,
                           int                n = -1
                           )
{
    using namespace optnet::utils;

    std::string base_name;

    // Check if the supplied file name is ended with '.hdr' or '.img'.
    if (str_ends_with_no_case(name, ".hdr") ||
        str_ends_with_no_case(name, ".img")) {
        base_name = name.substr(0, name.length() - 4);
    }
    else if (str_ends_with_no_case(name, ".img.gz")) {
        base_name = name.substr(0, name.length() - 7);
    }
    else {
        base_name = name;
    }

    if (n >= 0) {
        char number[32];
        secure_sprintf(number, sizeof(number), ".%04d", n);
        base_name += number;
    }

    hdr_name = base_name + ".hdr";
    img_name = base_name + ".img";
}


///////////////////////////////////////////////////////////////////////////
inline void
analyze_header_swap_endian(struct dsr* p)
{
    using namespace optnet::utils;

    swap_endian32(&(p->hk.sizeof_hdr)); 
    swap_endian32(&(p->hk.extents)); 
    swap_endian16(&(p->hk.session_error)); 
    swap_endian16(&(p->dime.dim[0]));
    swap_endian16(&(p->dime.dim[1]));
    swap_endian16(&(p->dime.dim[2]));
    swap_endian16(&(p->dime.dim[3]));
    swap_endian16(&(p->dime.dim[4]));
    swap_endian16(&(p->dime.dim[5]));
    swap_endian16(&(p->dime.dim[6]));
    swap_endian16(&(p->dime.dim[7]));
    swap_endian16(&(p->dime.bitpix));
    swap_endian16(&(p->dime.datatype));
    swap_endian32(&(p->dime.pixdim[0])); 
    swap_endian32(&(p->dime.pixdim[1]));
    swap_endian32(&(p->dime.pixdim[2]));
    swap_endian32(&(p->dime.pixdim[3]));
    swap_endian32(&(p->dime.pixdim[4])); 
    swap_endian32(&(p->dime.pixdim[5])); 
    swap_endian32(&(p->dime.pixdim[6])); 
    swap_endian32(&(p->dime.pixdim[7])); 
    swap_endian32(&(p->dime.vox_offset));
    swap_endian32(&(p->dime.funused1)); 
    swap_endian32(&(p->dime.funused2)); 
    swap_endian32(&(p->dime.cal_max)); 
    swap_endian32(&(p->dime.cal_min)); 
    swap_endian32(&(p->dime.compressed)); 
    swap_endian32(&(p->dime.verified)); 
    swap_endian16(&(p->dime.dim_un0)); 
    swap_endian32(&(p->dime.glmax));
    swap_endian32(&(p->dime.glmin));
}


///////////////////////////////////////////////////////////////////////////
inline bool
analyze_load_info(analyze_volume*    volume, 
                  const std::string& hdr_name,
                  const std::string& img_name
                  )
{
    dsr     header; // defined in detail/analyze.h
    size_t  file_size, bytes_per_voxel;
    FILE*   pfile;
    int     ret;
    
    assert(0 != volume);

    // Clear the analyze_volume struct.
    memset(volume, 0, sizeof(analyze_volume));
    
    // Let's read in the header struct first.
    //
    ret = secure_fopen(&pfile, hdr_name.c_str(), "rb");
    if (0 != ret) return false; // Failed opening file.
    if (fread(&header, sizeof(char), sizeof(dsr), pfile) 
        != sizeof(dsr)) {
        fclose(pfile); // Failed reading file.
        return false;
    }
    fclose(pfile);

    // Change the endian type of the header if necessary.
    if (header.hk.sizeof_hdr == 348 || header.hk.sizeof_hdr == 384) {
        volume->swap_endian = false;
    }
    else {
        analyze_header_swap_endian(&header);
        volume->swap_endian = true;
    }

    // Some .hdr files have the dim[0] corrupted, let's find out otherwise.
    if (header.dime.dim[0] < 1 || header.dime.dim[0] > 15) {
        if (header.dime.dim[4] > 0)
            header.dime.dim[0] = 5;
        else if (header.dime.dim[3] > 0)
            header.dime.dim[0] = 4;
        else if (header.dime.dim[2] > 0)
            header.dime.dim[0] = 3;
        else
            header.dime.dim[0] = 2;
    }

    // Save voxel and image dimensions and the data type used.
    volume->voxel_size_x = header.dime.pixdim[1];
    volume->voxel_size_y = header.dime.pixdim[2];
    volume->voxel_size_z = header.dime.pixdim[3];

    volume->image_size_x = header.dime.dim[1];
    volume->image_size_y = header.dime.dim[2];
    volume->image_size_z = header.dime.dim[3];
    volume->image_size_w = header.dime.dim[4];

    if (0 == volume->image_size_x || 0 == volume->image_size_y || 
        0 == volume->image_size_z) return false;
    if (0 == volume->image_size_w) volume->image_size_w = 1;

    volume->image_size = volume->image_size_x * 
                         volume->image_size_y * 
                         volume->image_size_z * 
                         volume->image_size_w;


    switch (header.dime.datatype) {
    case DT_UNSIGNED_CHAR:  volume->datatype   = DT_UNSIGNED_CHAR; 
                            volume->voxel_bytes = 1;
                            break;
    case DT_UNSIGNED_SHORT: volume->datatype   = DT_UNSIGNED_SHORT;
                            volume->voxel_bytes = 2;
                            break;
    case DT_SIGNED_SHORT:   volume->datatype   = DT_SIGNED_SHORT;
                            volume->voxel_bytes = 2;
                            break;
    case DT_FLOAT:          volume->datatype   = DT_FLOAT;
                            volume->voxel_bytes = 4;
                            break;
    default:                volume->datatype   = DT_NONE;
                            volume->voxel_bytes = 0;
    }

    if (volume->datatype == DT_UNSIGNED_CHAR)
        volume->swap_endian = false;
    if (volume->datatype != DT_NONE)
        return true;

    // If the data type information in the header is invalid,
    // let's try to figure out the correct data type.
    //
    ret = secure_fopen(&pfile, img_name.c_str(), "rb");
    if (0 != ret) return false;
    if (fseek(pfile, 0, SEEK_END) != 0) {
        fclose(pfile);
        return false;
    }
    file_size = (unsigned int)::ftell(pfile);
    fclose(pfile);

    bytes_per_voxel = file_size / volume->image_size;
    switch (bytes_per_voxel) {
    case 1: volume->datatype = DT_UNSIGNED_CHAR; break;
    case 2: volume->datatype = DT_SIGNED_SHORT;  break;
    default:
        return false;
    }

    return true;
}


///////////////////////////////////////////////////////////////////////////
inline bool
analyze_load_data(analyze_volume*      volume,
                  const std::string& /*hdr_name*/,
                  const std::string&   img_name
                  )
{
    using namespace optnet::utils;

    std::string img_gz_name;
    size_t      file_size;
    FILE*       pfile;

    assert(0 != volume && 0 != volume->data);

    file_size = volume->image_size * volume->voxel_bytes;
    if (file_size == 0) return false;

    // Open image file for reading.
    // First try to open an (uncompressed) file with a '.img' extension.
    // If failed, try to open a (compressed) file with a '.img.gz' extension.
    int ret = secure_fopen(&pfile, img_name.c_str(), "rb");

    if (0 == ret) {
        if (fread(volume->data, sizeof(char), file_size, pfile) 
            != file_size) {
            fclose(pfile);
            return false;
        }
        fclose(pfile);
    }
    else {
        img_gz_name = img_name + ".gz";
        gzFile pgzfile = gzopen(img_gz_name.c_str(), "rb");
        if (0 == pgzfile) return false;

        if (gzread(pgzfile, volume->data, (unsigned int)file_size)
            != (int)file_size) {
            gzclose(pgzfile);
            return false;
        }
        gzclose(pgzfile);
    }

    // Adjust endian of the image data if necessary.
    if (volume->swap_endian) {
        switch (volume->voxel_bytes) {
        case 2: {
            for (short* p = (short*)volume->data; 
                        p < (short*)volume->data + volume->image_size;
                        ++p)
                swap_endian16(p);
            break;
        }
        case 4: {
            for (long*  p = (long *)volume->data; 
                        p < (long *)volume->data + volume->image_size;
                        ++p)
                swap_endian32(p);
            break;
        }
        default:
            ;
        }
    }

    return true;
}


///////////////////////////////////////////////////////////////////////////
inline bool
analyze_save_info(const analyze_volume* volume,
                  const std::string&    hdr_name,
                  const std::string&  /*img_name*/
                  )
{
    dsr     header;
    FILE*   pfile;
    
    assert(0 != volume);

    // Clear header struct.
    memset(&header, 0, sizeof(dsr));

    // Fill in known fields in the header struct.
    header.hk.sizeof_hdr  = sizeof(dsr);
    header.hk.extents     = 16384;

    header.dime.dim[0]    = 5;
    header.dime.dim[1]    = (short)volume->image_size_x;
    header.dime.dim[2]    = (short)volume->image_size_y;
    header.dime.dim[3]    = (short)volume->image_size_z;
    header.dime.dim[4]    = (short)volume->image_size_w;

    header.dime.pixdim[1] = (float)volume->voxel_size_x;
    header.dime.pixdim[2] = (float)volume->voxel_size_y;
    header.dime.pixdim[3] = (float)volume->voxel_size_z;

    header.dime.bitpix    = (short)(volume->voxel_bytes * 8);
    header.dime.datatype  = (short)(volume->datatype);

    // Save header data in to a file.
    int ret = secure_fopen(&pfile, hdr_name.c_str(), "wb");
    if (0 != ret) return false;

    if (fwrite(&header, sizeof(char), sizeof(dsr), pfile) 
        != sizeof(dsr)) {
        fclose(pfile); // Failed writing file.
        return false;
    }
    fclose(pfile);

    return true;
}


///////////////////////////////////////////////////////////////////////////
inline bool
analyze_save_data(const analyze_volume* volume,
                  const std::string&  /*hdr_name*/,
                  const std::string&    img_name,
                  int                   compression
                  )
{
    std::string img_gz_name;
    size_t      file_size;
    FILE*       pfile;
    
    assert(0 != volume);

    file_size = volume->image_size * volume->voxel_bytes;
    if (0 == file_size || 0 == volume->data)
        return false; // Nothing to save.

    // Save image data, compress if requested.
    if (0 == compression) {

        int ret = secure_fopen(&pfile, img_name.c_str(), "wb");
        if (0 != ret) return false;

        if (fwrite(volume->data, sizeof(char), file_size, pfile)
            != file_size) {
            // Failed writing the image.
            fclose(pfile);
            return false;
        }
        fclose(pfile);
    }
    else {

        char flags[16];
        secure_sprintf(flags, sizeof(flags), "wb%d", compression);

        // Open a file with a ".img.gz" extension.
        img_gz_name = img_name + ".gz";
        gzFile pgzfile = gzopen(img_gz_name.c_str(), flags);
        if (0 == pgzfile) return false;

        if (gzwrite(pgzfile, (const voidp)volume->data, (unsigned int)file_size) 
            != (int)file_size) {
            // Failed writing the image.
            gzclose(pgzfile);
            return false;
        }
        gzclose(pgzfile);
    }

    return true;
}


} } } // namespace

#endif
