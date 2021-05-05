/**************************************************************************
 *
 *  ANALYZE TM Header File Format
 * 
 *  Copyright (C) 1986-1995
 *  Biomedical Imaging Resource
 *  Mayo Foundation 
 * 
 *  http://www.mayo.edu/bir/Software/Analyze/AnalyzeTechInfo.html#formats
 *  http://www.mrc-cbu.cam.ac.uk/Imaging/analyze_fmt.htm
 *
 **************************************************************************
 */ 

#ifndef __ANALYZE_H__
#define __ANALYZE_H__

/*------------------------------------------------------------------------*/
#define DT_NONE 0
#define DT_UNKNOWN 0        /*  Unknown data type                       */ 
#define DT_BINARY 1         /*  Binary (1 bit per voxel)                */ 
#define DT_UNSIGNED_CHAR 2  /*  Unsigned character (8 bits per voxel)   */ 
#define DT_SIGNED_SHORT 4   /*  Signed short (16 bits per voxel)        */ 
#define DT_UNSIGNED_SHORT 5 /*  Unsigned short (16 bits per voxel)      */ 
#define DT_SIGNED_INT 8     /*  Signed integer (32 bits per voxel)      */ 
#define DT_UNSIGNED_INT 9   /*  Unsigned integer (32 bits per voxel)    */ 
#define DT_FLOAT 16         /*  Floating point (32 bits per voxel)      */ 
#define DT_COMPLEX 32       /*  Complex (64 bits per voxel)             */ 
#define DT_DOUBLE 64        /*  Double precision (64 bits per voxel)    */ 
#define DT_RGB 128 
#define DT_ALL 255
/*------------------------------------------------------------------------*/

enum Orientation
{
    TRANSVERSE_UNFLIPPED = 0, 
    CORONAL_UNFLIPPED,
    SAGITTAL_UNFLIPPED,
    TRANSVERSE_FLIPPED,
    CORONAL_FLIPPED,
    SAGITAL_FLIPPED
};


struct header_key
{                               /* offset + size    */ 
  int sizeof_hdr;               /*      0 +  4      */ 
  char data_type[10];           /*      4 + 10      */ 
  char db_name[18];             /*     14 + 18      */ 
  int extents;                  /*     32 +  4      */ 
  short int session_error;      /*     36 +  2      */ 
  char regular;                 /*     38 +  1      */ 
  char hkey_un0;                /*     39 +  1      */ 
};                              /* total: 40 bytes  */


struct image_dimension 
{                               /* offset + size    */ 
  short int dim[8];             /*      0 + 16      */ 
  /*
    dim[0] Number of dimensions in database; usually 4 
    dim[1] Volume X dimension; number of pixels in an image row 
    dim[2] Volume Y dimension; number of pixel rows in slice 
    dim[3] Volume Z dimension; number of slices in a volume 
    dim[4] Time points, number of volumes in database.
  */
  short int unused8;            /*     16 +  2      */ 
  short int unused9;            /*     18 +  2      */ 
  short int unused10;           /*     20 +  2      */ 
  short int unused11;           /*     22 +  2      */ 
  short int unused12;           /*     24 +  2      */ 
  short int unused13;           /*     26 +  2      */ 
  short int unused14;           /*     28 +  2      */ 
  short int datatype;           /*     30 +  2      */ 
  short int bitpix;             /*     32 +  2      */ 
  short int dim_un0;            /*     34 +  2      */ 
  float pixdim[8];              /*     36 + 32      */ 
  /* 
     pixdim specifies the voxel dimensitons: 
     pixdim[1] - voxel width 
     pixdim[2] - voxel height 
     pixdim[3] - interslice distance 
     ...etc 
  */ 
  float vox_offset;             /*     68 +  4      */ 
  float funused1;               /*     72 +  4      */ 
  float funused2;               /*     76 +  4      */ 
  float funused3;               /*     80 +  4      */ 
  float cal_max;                /*     84 +  4      */ 
  float cal_min;                /*     88 +  4      */ 
  float compressed;             /*     92 +  4      */ 
  float verified;               /*     96 +  4      */ 
  int glmax,glmin;              /*    100 +  8      */ 
};                              /* total: 108 bytes */


struct data_history 
{                               /* offset + size    */ 
  char descrip[80];             /*      0 + 80      */ 
  char aux_file[24];            /*     80 + 24      */ 
  char orient;                  /*    104 +  1      */ 
  /*
    orient                slice orientation for this dataset
    --------------------------------------------------------
    TRANSVERSE_UNFLIPPED  transverse unflipped 
    CORONAL_UNFLIPPED     coronal unflipped 
    SAGITTAL_UNFLIPPED    sagittal unflipped 
    TRANSVERSE_FLIPPED    transverse flipped 
    CORONAL_FLIPPED       coronal flipped 
    SAGITAL_FLIPPED       sagittal flipped 
  */
  char originator[10];          /*    105 + 10      */ 
  char generated[10];           /*    115 + 10      */ 
  char scannum[10];             /*    125 + 10      */ 
  char patient_id[10];          /*    135 + 10      */ 
  char exp_date[10];            /*    145 + 10      */ 
  char exp_time[10];            /*    155 + 10      */ 
  char hist_un0[3];             /*    165 +  3      */ 
  int views;                    /*    168 +  4      */ 
  int vols_added;               /*    172 +  4      */ 
  int start_field;              /*    176 +  4      */ 
  int field_skip;               /*    180 +  4      */ 
  int omax, omin;               /*    184 +  8      */ 
  int smax, smin;               /*    192 +  8      */ 
};


struct dsr 
{ 
  struct header_key hk;         /*     0 + 40       */ 
  struct image_dimension dime;  /*    40 + 108      */ 
  struct data_history hist;     /*   148 + 200      */ 
};                              /* total: 348 bytes */ 


typedef struct 
{ 
  float real; 
  float imag; 
} COMPLEX;


#endif  /*  __ANALYZE_H__ */
