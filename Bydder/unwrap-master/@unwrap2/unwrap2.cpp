// 2D phase unwrapping, modified for inclusion in scipy by Gregor Thalhammer
// Templated and mex wrapped here https://github.com/marcsous/unwrap
// Original code at https://github.com/geggo/phase-unwrap

#include "mex.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include "mxShowCriticalErrorMessage.h"

// forward declaration of unwrap2D function (below)
template <typename T>
        void unwrap2D(T* wrapped_image, T* UnwrappedImage, unsigned char* input_mask,
        int image_width, int image_height,
        int wrap_around_x, int wrap_around_y);

// gateway function to MATLAB
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    if(nrhs < 1 || nrhs > 2)
        mxShowCriticalErrorMessage("wrong number of input arguments",nrhs);
    
    if(nlhs > 1)
        mxShowCriticalErrorMessage("wrong number of output arguments",nlhs);
    
    switch(mxGetClassID(prhs[0]))
    {
        case mxSINGLE_CLASS:
        case mxDOUBLE_CLASS: break;
        default: mxShowCriticalErrorMessage("argument 1 type incompatible");
    }
    
    if(mxIsComplex(prhs[0]))
        mxShowCriticalErrorMessage("argument 1 must be real valued");
    
    int ndims = mxGetNumberOfDimensions(prhs[0]);
    
    if(ndims != 2)
        mxShowCriticalErrorMessage("argument 1 must be a 2D array");
    
    const mwSize *dims = mxGetDimensions(prhs[0]);
    int nx = dims[0], ny = dims[1];
    
    if(nx <= 1 || ny <= 1)
        mxShowCriticalErrorMessage("argument 1 must be a 2D array");
    
    // mask (needs to be unsigned char): 0=keep 1=reject (default to 0)
    unsigned char* input_mask = (unsigned char*)mxCalloc(nx*ny,sizeof(unsigned char));
    
    if (nrhs==2)
    {
        switch(mxGetClassID(prhs[1]))
        {
            case mxLOGICAL_CLASS:
            case mxINT8_CLASS:
            case mxUINT8_CLASS:
            case mxINT16_CLASS:
            case mxUINT16_CLASS:
            case mxINT32_CLASS:
            case mxUINT32_CLASS:
            case mxSINGLE_CLASS:
            case mxINT64_CLASS:
            case mxUINT64_CLASS:
            case mxDOUBLE_CLASS: break;
            default: mxShowCriticalErrorMessage("argument 2 type incompatible");
        }
        
        if(mxGetNumberOfDimensions(prhs[1]) != 2)
            mxShowCriticalErrorMessage("argument 2 must be a 2D array");
        
        dims = mxGetDimensions(prhs[1]);
        if(dims[0] != nx || dims[1] != ny)
            mxShowCriticalErrorMessage("argument 1 and 2 must be same size");
        
        if(mxIsComplex(prhs[1]))
            mxShowCriticalErrorMessage("argument 2 must be real valued");
        
        // handle all valid types: 1 bit, 8 bit, 16 bit, 32 bit, 64 bit
        for (int i=0; i<nx*ny; i++)
        {
            switch(mxGetClassID(prhs[1]))
            {
                case mxLOGICAL_CLASS: input_mask[i] = ((bool*)mxGetData(prhs[1]))[i] == 0;
                break;
                
                case mxINT8_CLASS:
                case mxUINT8_CLASS:   input_mask[i] = ((uint8_t*)mxGetData(prhs[1]))[i] == 0;
                break;
                
                case mxINT16_CLASS:
                case mxUINT16_CLASS:  input_mask[i] = ((uint16_t*)mxGetData(prhs[1]))[i] == 0;
                break;
                
                case mxINT32_CLASS:
                case mxUINT32_CLASS:
                case mxSINGLE_CLASS:  input_mask[i] = ((uint32_t*)mxGetData(prhs[1]))[i] == 0;
                break;
                
                case mxINT64_CLASS:
                case mxUINT64_CLASS:
                case mxDOUBLE_CLASS:  input_mask[i] = ((uint64_t*)mxGetData(prhs[1]))[i] == 0;
                break;
            }
        }
    }
    
    // check for NaNs in the data and set mask to 1 (reject)
    for (int i=0; i<nx*ny; i++)
    {
        float *f; double *d;
        switch(mxGetClassID(prhs[0]))
        {
            case mxSINGLE_CLASS:
                f = (float*)mxGetData(prhs[0]);
                if (!mxIsFinite(f[i])) input_mask[i] = 1;
                break;
                
            case mxDOUBLE_CLASS:
                d = (double*)mxGetData(prhs[0]);
                if (!mxIsFinite(d[i])) input_mask[i] = 1;
                break;
        }
    }
    
    // circular wrap (1=yes 0=no)
    int wrap_around_x = 0;
    int wrap_around_y = 0;
    int wrap_around_z = 0;
    
    // output array - handle both floating point types
    switch(mxGetClassID(prhs[0]))
    {
        case mxSINGLE_CLASS: plhs[0] = mxCreateNumericArray(ndims, dims, mxSINGLE_CLASS, mxREAL);
        unwrap2D((float*)mxGetData(prhs[0]), (float*)mxGetData(plhs[0]),
                input_mask, nx, ny, wrap_around_x, wrap_around_y);
        break;
        
        case mxDOUBLE_CLASS: plhs[0] = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL);
        unwrap2D((double*)mxGetData(prhs[0]), (double*)mxGetData(plhs[0]),
                input_mask, nx, ny, wrap_around_x, wrap_around_y);
        break;
    }
    
}

//This program was written by Munther Gdeisat and Miguel Arevallilo Herraez to program the two-dimensional unwrapper
//entitled "Fast two-dimensional phase-unwrapping algorithm based on sorting by
//reliability following a noncontinuous path"
//by  Miguel Arevallilo Herraez, David R. Burton, Michael J. Lalor, and Munther A. Gdeisat
//published in the Journal Applied Optics, Vol. 41, No. 35, pp. 7437, 2002.
//This program was written by Munther Gdeisat, Liverpool John Moores University, United Kingdom.
//Date 26th August 2007
//The wrapped phase map is assumed to be of floating point data type. The resultant unwrapped phase map is also of floating point type.
//The mask is of byte data type.
//When the mask is 1 this means that the pixel is valid
//When the mask is 0 this means that the pixel is invalid (noisy or corrupted pixel)
//This program takes into consideration the image wrap around problem encountered in MRI imaging.

#define NOMASK 0
#define MASK 1

#define PI static_cast<T>(M_PI)
#define TWOPI static_cast<T>(2*M_PI)

template <typename T>
        struct params_t
{
    T mod;
    int x_connectivity;
    int y_connectivity;
    int no_of_edges;
};

//PIXELM information
template <typename T>
        struct PIXELM
{
    int increment;		//No. of 2*pi to add to the pixel to unwrap it
    int number_of_pixels_in_group;//No. of pixel in the pixel group
    T value;			//value of the pixel
    T reliability;
    unsigned char input_mask;	//0 pixel is masked. NOMASK pixel is not masked
    unsigned char extended_mask;	//0 pixel is masked. NOMASK pixel is not masked
    int group;			//group No.
    int new_group;
    struct PIXELM<T> *head;		//pointer to the first pixel in the group in the linked list
    struct PIXELM<T> *last;		//pointer to the last pixel in the group
    struct PIXELM<T> *next;		//pointer to the next pixel in the group
};

//the EDGE is the line that connects two pixels.
//if we have S pixels, then we have S horizontal edges and S vertical edges
template <typename T>
        struct EDGE
{
    T reliab;			//reliabilty of the edge and it depends on the two pixels
    PIXELM<T> *pointer_1;		//pointer to the first pixel
    PIXELM<T> *pointer_2;		//pointer to the second pixel
    int increment;		//No. of 2*pi to add to one of the pixels to
    //unwrap it with respect to the second
};

//---------------start quicker_sort algorithm --------------------------------
#define EDGEswap(x,y) {EDGE<T> t; t=x; x=y; y=t;}
#define EDGEorder(x,y) if (x.reliab > y.reliab) EDGEswap(x,y)
#define EDGEo2(x,y) EDGEorder(x,y)
#define EDGEo3(x,y,z) EDGEo2(x,y); EDGEo2(x,z); EDGEo2(y,z)

typedef enum {yes, no} yes_no;

template <typename T>
        yes_no find_pivot(EDGE<T> *left, EDGE<T> *right, T *pivot_ptr)
{
    EDGE<T> a, b, c, *p;
    
    a = *left;
    b = *(left + (right - left) /2 );
    c = *right;
    EDGEo3(a,b,c);
    
    if (a.reliab < b.reliab)
    {
        *pivot_ptr = b.reliab;
        return yes;
    }
    
    if (b.reliab < c.reliab)
    {
        *pivot_ptr = c.reliab;
        return yes;
    }
    
    for (p = left + 1; p <= right; ++p)
    {
        if (p->reliab != left->reliab)
        {
            *pivot_ptr = (p->reliab < left->reliab) ? left->reliab : p->reliab;
            return yes;
        }
        return no;
    }
    
    return no;
}

template <typename T>
        EDGE<T> *partition(EDGE<T> *left, EDGE<T> *right, T pivot)
{
    while (left <= right)
    {
        while (left->reliab < pivot)
            ++left;
        while (right->reliab >= pivot)
            --right;
        if (left < right)
        {
            EDGEswap (*left, *right);
            ++left;
            --right;
        }
    }
    return left;
}

template <typename T>
        void quicker_sort(EDGE<T> *left, EDGE<T> *right)
{
    EDGE<T> *p;
    T pivot;
    
    if (find_pivot(left, right, &pivot) == yes)
    {
        p = partition(left, right, pivot);
        quicker_sort(left, p - 1);
        quicker_sort(p, right);
    }
}
//--------------end quicker_sort algorithm -----------------------------------

//--------------------start initialize pixels ----------------------------------
//initialize pixels. See the explination of the pixel class above.
//initially every pixel is assumed to belong to a group consisting of only itself
template <typename T>
        void initialisePIXELs(T *wrapped_image, unsigned char *input_mask, unsigned char *extended_mask, PIXELM<T> *pixel, int image_width, int image_height)
{
    PIXELM<T> *pixel_pointer = pixel;
    T *wrapped_image_pointer = wrapped_image;
    unsigned char *input_mask_pointer = input_mask;
    unsigned char *extended_mask_pointer = extended_mask;
    int i, j;
    
    for (i=0; i < image_height; i++)
    {
        for (j=0; j < image_width; j++)
        {
            pixel_pointer->increment = 0;
            pixel_pointer->number_of_pixels_in_group = 1;
            pixel_pointer->value = *wrapped_image_pointer;
            pixel_pointer->reliability = 9999999999. + rand(); // MB rand?
            pixel_pointer->input_mask = *input_mask_pointer;
            pixel_pointer->extended_mask = *extended_mask_pointer;
            pixel_pointer->head = pixel_pointer;
            pixel_pointer->last = pixel_pointer;
            pixel_pointer->next = NULL;
            pixel_pointer->new_group = 0;
            pixel_pointer->group = -1;
            pixel_pointer++;
            wrapped_image_pointer++;
            input_mask_pointer++;
            extended_mask_pointer++;
        }
    }
}
//-------------------end initialize pixels -----------

//gamma function in the paper
template <typename T>
        T wrap(T pixel_value)
{
    T wrapped_pixel_value;
    if (pixel_value > PI)	wrapped_pixel_value = pixel_value - TWOPI;
    else if (pixel_value < -PI) wrapped_pixel_value = pixel_value + TWOPI;
    else wrapped_pixel_value = pixel_value;
    return wrapped_pixel_value;
}

// pixelL_value is the left pixel,	pixelR_value is the right pixel
template <typename T>
        int find_wrap(T pixelL_value, T pixelR_value)
{
    T difference;
    int wrap_value;
    difference = pixelL_value - pixelR_value;
    
    if (difference > PI)	wrap_value = -1;
    else if (difference < -PI)	wrap_value = 1;
    else wrap_value = 0;
    
    return wrap_value;
}

template <typename T>
        void extend_mask(unsigned char *input_mask, unsigned char *extended_mask,
        int image_width, int image_height,
        params_t<T> *params)
{
    int i,j;
    int image_width_plus_one = image_width + 1;
    int image_width_minus_one = image_width - 1;
    unsigned char *IMP;	//input mask pointer
    unsigned char *EMP = extended_mask;	//extended mask pointer
    
    for (i=0; i < image_height; ++i)
    {
        for (j=0; j < image_width; ++j)
        {
            *EMP = MASK;
            ++EMP;
        }
    }
    
    IMP = input_mask    + image_width + 1;
    EMP = extended_mask + image_width + 1;
    //extend the mask for the image except borders
    for (i=1; i < image_height - 1; ++i)
    {
        for (j=1; j < image_width - 1; ++j)
        {
            if ( (*IMP) == NOMASK && (*(IMP + 1) == NOMASK) && (*(IMP - 1) == NOMASK) &&
                    (*(IMP + image_width) == NOMASK) && (*(IMP - image_width) == NOMASK) &&
                    (*(IMP - image_width_minus_one) == NOMASK) && (*(IMP - image_width_plus_one) == NOMASK) &&
                    (*(IMP + image_width_minus_one) == NOMASK) && (*(IMP + image_width_plus_one) == NOMASK) )
            {
                *EMP = NOMASK;
            }
            ++EMP;
            ++IMP;
        }
        EMP += 2;
        IMP += 2;
    }
    
    if (params->x_connectivity == 1)
    {
        //extend the mask for the right border of the image
        IMP = input_mask    + 2 * image_width - 1;
        EMP = extended_mask + 2 * image_width -1;
        for (i=1; i < image_height - 1; ++ i)
        {
            if ( (*IMP) == NOMASK && (*(IMP - 1) == NOMASK) &&  (*(IMP + 1) == NOMASK) &&
                    (*(IMP + image_width) == NOMASK) && (*(IMP - image_width) == NOMASK) &&
                    (*(IMP - image_width - 1) == NOMASK) && (*(IMP - image_width + 1) == NOMASK) &&
                    (*(IMP + image_width - 1) == NOMASK) && (*(IMP - 2 * image_width + 1) == NOMASK) )
            {
                *EMP = NOMASK;
            }
            EMP += image_width;
            IMP += image_width;
        }
        
        //extend the mask for the left border of the image
        IMP = input_mask    + image_width;
        EMP = extended_mask + image_width;
        for (i=1; i < image_height - 1; ++i)
        {
            if ( (*IMP) == NOMASK && (*(IMP - 1) == NOMASK) && (*(IMP + 1) == NOMASK) &&
                    (*(IMP + image_width) == NOMASK) && (*(IMP - image_width) == NOMASK) &&
                    (*(IMP - image_width + 1) == NOMASK) && (*(IMP + image_width + 1) == NOMASK) &&
                    (*(IMP + image_width - 1) == NOMASK) && (*(IMP + 2 * image_width - 1) == NOMASK) )
            {
                *EMP = NOMASK;
            }
            EMP += image_width;
            IMP += image_width;
        }
    }
    
    if (params->y_connectivity == 1)
    {
        //extend the mask for the top border of the image
        IMP = input_mask    + 1;
        EMP = extended_mask + 1;
        for (i=1; i < image_width - 1; ++i)
        {
            if ( (*IMP) == NOMASK && (*(IMP - 1) == NOMASK) && (*(IMP + 1) == NOMASK) &&
                    (*(IMP + image_width) == NOMASK) && (*(IMP + image_width * (image_height - 1)) == NOMASK) &&
                    (*(IMP + image_width + 1) == NOMASK) && (*(IMP + image_width - 1) == NOMASK) &&
                    (*(IMP + image_width * (image_height - 1) - 1) == NOMASK) && (*(IMP + image_width * (image_height - 1) + 1) == NOMASK) )
            {
                *EMP = NOMASK;
            }
            EMP++;
            IMP++;
        }
        
        //extend the mask for the bottom border of the image
        IMP = input_mask    + image_width * (image_height - 1) + 1;
        EMP = extended_mask + image_width * (image_height - 1) + 1;
        for (i=1; i < image_width - 1; ++i)
        {
            if ( (*IMP) == NOMASK && (*(IMP - 1) == NOMASK) && (*(IMP + 1) == NOMASK) &&
                    (*(IMP - image_width) == NOMASK) && (*(IMP - image_width - 1) == NOMASK) && (*(IMP - image_width + 1) == NOMASK) &&
                    (*(IMP - image_width * (image_height - 1)    ) == NOMASK) &&
                    (*(IMP - image_width * (image_height - 1) - 1) == NOMASK) &&
                    (*(IMP - image_width * (image_height - 1) + 1) == NOMASK) )
            {
                *EMP = NOMASK;
            }
            EMP++;
            IMP++;
        }
    }
}

template <typename T>
        void calculate_reliability(T *wrappedImage, PIXELM<T> *pixel,
        int image_width, int image_height,
        params_t<T> *params)
{
    int image_width_plus_one = image_width + 1;
    int image_width_minus_one = image_width - 1;
    PIXELM<T> *pixel_pointer = pixel + image_width_plus_one;
    T *WIP = wrappedImage + image_width_plus_one; //WIP is the wrapped image pointer
    T H, V, D1, D2;
    int i, j;
    
    for (i = 1; i < image_height -1; ++i)
    {
        for (j = 1; j < image_width - 1; ++j)
        {
            if (pixel_pointer->extended_mask == NOMASK)
            {
                H = wrap(*(WIP - 1) - *WIP) - wrap(*WIP - *(WIP + 1));
                V = wrap(*(WIP - image_width) - *WIP) - wrap(*WIP - *(WIP + image_width));
                D1 = wrap(*(WIP - image_width_plus_one) - *WIP) - wrap(*WIP - *(WIP + image_width_plus_one));
                D2 = wrap(*(WIP - image_width_minus_one) - *WIP) - wrap(*WIP - *(WIP + image_width_minus_one));
                pixel_pointer->reliability = H*H + V*V + D1*D1 + D2*D2;
            }
            pixel_pointer++;
            WIP++;
        }
        pixel_pointer += 2;
        WIP += 2;
    }
    
    if (params->x_connectivity == 1)
    {
        //calculating the reliability for the left border of the image
        PIXELM<T> *pixel_pointer = pixel + image_width;
        T *WIP = wrappedImage + image_width;
        
        for (i = 1; i < image_height - 1; ++i)
        {
            if (pixel_pointer->extended_mask == NOMASK)
            {
                H = wrap(*(WIP + image_width - 1) - *WIP) - wrap(*WIP - *(WIP + 1));
                V = wrap(*(WIP - image_width) - *WIP) - wrap(*WIP - *(WIP + image_width));
                D1 = wrap(*(WIP - 1) - *WIP) - wrap(*WIP - *(WIP + image_width_plus_one));
                D2 = wrap(*(WIP - image_width_minus_one) - *WIP) - wrap(*WIP - *(WIP + 2* image_width - 1));
                pixel_pointer->reliability = H*H + V*V + D1*D1 + D2*D2;
            }
            pixel_pointer += image_width;
            WIP += image_width;
        }
        
        //calculating the reliability for the right border of the image
        pixel_pointer = pixel + 2 * image_width - 1;
        WIP = wrappedImage + 2 * image_width - 1;
        
        for (i = 1; i < image_height - 1; ++i)
        {
            if (pixel_pointer->extended_mask == NOMASK)
            {
                H = wrap(*(WIP - 1) - *WIP) - wrap(*WIP - *(WIP - image_width_minus_one));
                V = wrap(*(WIP - image_width) - *WIP) - wrap(*WIP - *(WIP + image_width));
                D1 = wrap(*(WIP - image_width_plus_one) - *WIP) - wrap(*WIP - *(WIP + 1));
                D2 = wrap(*(WIP - 2 * image_width - 1) - *WIP) - wrap(*WIP - *(WIP + image_width_minus_one));
                pixel_pointer->reliability = H*H + V*V + D1*D1 + D2*D2;
            }
            pixel_pointer += image_width;
            WIP += image_width;
        }
    }
    
    if (params->y_connectivity == 1)
    {
        //calculating the reliability for the top border of the image
        PIXELM<T> *pixel_pointer = pixel + 1;
        T *WIP = wrappedImage + 1;
        
        for (i = 1; i < image_width - 1; ++i)
        {
            if (pixel_pointer->extended_mask == NOMASK)
            {
                H =  wrap(*(WIP - 1) - *WIP) - wrap(*WIP - *(WIP + 1));
                V =  wrap(*(WIP + image_width*(image_height - 1)) - *WIP) - wrap(*WIP - *(WIP + image_width));
                D1 = wrap(*(WIP + image_width*(image_height - 1) - 1) - *WIP) - wrap(*WIP - *(WIP + image_width_plus_one));
                D2 = wrap(*(WIP + image_width*(image_height - 1) + 1) - *WIP) - wrap(*WIP - *(WIP + image_width_minus_one));
                pixel_pointer->reliability = H*H + V*V + D1*D1 + D2*D2;
            }
            pixel_pointer++;
            WIP++;
        }
        
        //calculating the reliability for the bottom border of the image
        pixel_pointer = pixel + (image_height - 1) * image_width + 1;
        WIP = wrappedImage + (image_height - 1) * image_width + 1;
        
        for (i = 1; i < image_width - 1; ++i)
        {
            if (pixel_pointer->extended_mask == NOMASK)
            {
                H =  wrap(*(WIP - 1) - *WIP) - wrap(*WIP - *(WIP + 1));
                V =  wrap(*(WIP - image_width) - *WIP) - wrap(*WIP - *(WIP -(image_height - 1) * (image_width)));
                D1 = wrap(*(WIP - image_width_plus_one) - *WIP) - wrap(*WIP - *(WIP - (image_height - 1) * (image_width) + 1));
                D2 = wrap(*(WIP - image_width_minus_one) - *WIP) - wrap(*WIP - *(WIP - (image_height - 1) * (image_width) - 1));
                pixel_pointer->reliability = H*H + V*V + D1*D1 + D2*D2;
            }
            pixel_pointer++;
            WIP++;
        }
    }
}

//calculate the reliability of the horizontal edges of the image
//it is calculated by adding the reliability of pixel and the relibility of
//its right-hand neighbour
//edge is calculated between a pixel and its next neighbour
template <typename T>
        void horizontalEDGEs(PIXELM<T> *pixel, EDGE<T> *edge,
        int image_width, int image_height,
        params_t<T> *params)
{
    int i, j;
    EDGE<T> *edge_pointer = edge;
    PIXELM<T> *pixel_pointer = pixel;
    int no_of_edges = params->no_of_edges;
    
    for (i = 0; i < image_height; i++)
    {
        for (j = 0; j < image_width - 1; j++)
        {
            if (pixel_pointer->input_mask == NOMASK && (pixel_pointer + 1)->input_mask == NOMASK)
            {
                edge_pointer->pointer_1 = pixel_pointer;
                edge_pointer->pointer_2 = (pixel_pointer+1);
                edge_pointer->reliab = pixel_pointer->reliability + (pixel_pointer + 1)->reliability;
                edge_pointer->increment = find_wrap(pixel_pointer->value, (pixel_pointer + 1)->value);
                edge_pointer++;
                no_of_edges++;
            }
            pixel_pointer++;
        }
        pixel_pointer++;
    }
    
    //construct edges at the right border of the image
    if (params->x_connectivity == 1)
    {
        pixel_pointer = pixel + image_width - 1;
        for (i = 0; i < image_height; i++)
        {
            if (pixel_pointer->input_mask == NOMASK && (pixel_pointer - image_width + 1)->input_mask == NOMASK)
            {
                edge_pointer->pointer_1 = pixel_pointer;
                edge_pointer->pointer_2 = (pixel_pointer - image_width + 1);
                edge_pointer->reliab = pixel_pointer->reliability + (pixel_pointer - image_width + 1)->reliability;
                edge_pointer->increment = find_wrap(pixel_pointer->value, (pixel_pointer  - image_width + 1)->value);
                edge_pointer++;
                no_of_edges++;
            }
            pixel_pointer+=image_width;
        }
    }
    params->no_of_edges = no_of_edges;
}

//calculate the reliability of the vertical edges of the image
//it is calculated by adding the reliability of pixel and the relibility of
//its lower neighbour in the image.
template <typename T>
        void verticalEDGEs(PIXELM<T> *pixel, EDGE<T> *edge,
        int image_width, int image_height,
        params_t<T> *params)
{
    int i, j;
    int no_of_edges = params->no_of_edges;
    PIXELM<T> *pixel_pointer = pixel;
    EDGE<T> *edge_pointer = edge + no_of_edges;
    
    for (i=0; i < image_height - 1; i++)
    {
        for (j=0; j < image_width; j++)
        {
            if (pixel_pointer->input_mask == NOMASK && (pixel_pointer + image_width)->input_mask == NOMASK)
            {
                edge_pointer->pointer_1 = pixel_pointer;
                edge_pointer->pointer_2 = (pixel_pointer + image_width);
                edge_pointer->reliab = pixel_pointer->reliability + (pixel_pointer + image_width)->reliability;
                edge_pointer->increment = find_wrap(pixel_pointer->value, (pixel_pointer + image_width)->value);
                edge_pointer++;
                no_of_edges++;
            }
            pixel_pointer++;
        } 
    } 
    
    //construct edges that connect at the bottom border of the image
    if (params->y_connectivity == 1)
    {
        pixel_pointer = pixel + image_width *(image_height - 1);
        for (i = 0; i < image_width; i++)
        {
            if (pixel_pointer->input_mask == NOMASK && (pixel_pointer - image_width *(image_height - 1))->input_mask == NOMASK)
            {
                edge_pointer->pointer_1 = pixel_pointer;
                edge_pointer->pointer_2 = (pixel_pointer - image_width *(image_height - 1));
                edge_pointer->reliab = pixel_pointer->reliability + (pixel_pointer - image_width *(image_height - 1))->reliability;
                edge_pointer->increment = find_wrap(pixel_pointer->value, (pixel_pointer - image_width *(image_height - 1))->value);
                edge_pointer++;
                no_of_edges++;
            }
            pixel_pointer++;
        }
    }
    params->no_of_edges = no_of_edges;
}

//gather the pixels of the image into groups
template <typename T>
        void  gatherPIXELs(EDGE<T> *edge, params_t<T> *params)
{
    int k;
    PIXELM<T> *PIXEL1;
    PIXELM<T> *PIXEL2;
    PIXELM<T> *group1;
    PIXELM<T> *group2;
    EDGE<T> *pointer_edge = edge;
    int incremento;
    
    for (k = 0; k < params->no_of_edges; k++)
    {
        PIXEL1 = pointer_edge->pointer_1;
        PIXEL2 = pointer_edge->pointer_2;
        
        //PIXELM 1 and PIXELM 2 belong to different groups
        //initially each pixel is a group by it self and one pixel can construct a group
        //no else or else if to this if
        if (PIXEL2->head != PIXEL1->head)
        {
            //PIXELM 2 is alone in its group
            //merge this pixel with PIXELM 1 group and find the number of 2 pi to add
            //to or subtract to unwrap it
            if ((PIXEL2->next == NULL) && (PIXEL2->head == PIXEL2))
            {
                PIXEL1->head->last->next = PIXEL2;
                PIXEL1->head->last = PIXEL2;
                (PIXEL1->head->number_of_pixels_in_group)++;
                PIXEL2->head=PIXEL1->head;
                PIXEL2->increment = PIXEL1->increment-pointer_edge->increment;
            }
            
            //PIXELM 1 is alone in its group
            //merge this pixel with PIXELM 2 group and find the number of 2 pi to add
            //to or subtract to unwrap it
            else if ((PIXEL1->next == NULL) && (PIXEL1->head == PIXEL1))
            {
                PIXEL2->head->last->next = PIXEL1;
                PIXEL2->head->last = PIXEL1;
                (PIXEL2->head->number_of_pixels_in_group)++;
                PIXEL1->head = PIXEL2->head;
                PIXEL1->increment = PIXEL2->increment+pointer_edge->increment;
            }
            
            //PIXELM 1 and PIXELM 2 both have groups
            else
            {
                group1 = PIXEL1->head;
                group2 = PIXEL2->head;
                //if the no. of pixels in PIXELM 1 group is larger than the
                //no. of pixels in PIXELM 2 group.  Merge PIXELM 2 group to
                //PIXELM 1 group and find the number of wraps between PIXELM 2
                //group and PIXELM 1 group to unwrap PIXELM 2 group with respect
                //to PIXELM 1 group.  the no. of wraps will be added to PIXELM 2
                //group in the future
                if (group1->number_of_pixels_in_group > group2->number_of_pixels_in_group)
                {
                    //merge PIXELM 2 with PIXELM 1 group
                    group1->last->next = group2;
                    group1->last = group2->last;
                    group1->number_of_pixels_in_group = group1->number_of_pixels_in_group + group2->number_of_pixels_in_group;
                    incremento = PIXEL1->increment-pointer_edge->increment - PIXEL2->increment;
                    //merge the other pixels in PIXELM 2 group to PIXELM 1 group
                    while (group2 != NULL)
                    {
                        group2->head = group1;
                        group2->increment += incremento;
                        group2 = group2->next;
                    }
                }
                
                //if the no. of pixels in PIXELM 2 group is larger than the
                //no. of pixels in PIXELM 1 group.  Merge PIXELM 1 group to
                //PIXELM 2 group and find the number of wraps between PIXELM 2
                //group and PIXELM 1 group to unwrap PIXELM 1 group with respect
                //to PIXELM 2 group.  the no. of wraps will be added to PIXELM 1
                //group in the future
                else
                {
                    //merge PIXELM 1 with PIXELM 2 group
                    group2->last->next = group1;
                    group2->last = group1->last;
                    group2->number_of_pixels_in_group = group2->number_of_pixels_in_group + group1->number_of_pixels_in_group;
                    incremento = PIXEL2->increment + pointer_edge->increment - PIXEL1->increment;
                    //merge the other pixels in PIXELM 2 group to PIXELM 1 group
                    while (group1 != NULL)
                    {
                        group1->head = group2;
                        group1->increment += incremento;
                        group1 = group1->next;
                    } // while
                    
                } // else
            } //else
        } //if
        pointer_edge++;
    }
}

//unwrap the image
template <typename T>
        void  unwrapImage(PIXELM<T> *pixel, int image_width, int image_height)
{
    int i;
    int image_size = image_width * image_height;
    PIXELM<T> *pixel_pointer=pixel;
    
    for (i = 0; i < image_size; i++)
    {
        pixel_pointer->value += TWOPI * (T)(pixel_pointer->increment);
        pixel_pointer++;
    }
}

//the input to this unwrapper is an array that contains the wrapped
//phase map.  copy the image on the buffer passed to this unwrapper to
//over-write the unwrapped phase map on the buffer of the wrapped
//phase map.
template <typename T>
        void  returnImage(PIXELM<T> *pixel, T *unwrapped_image, int image_width, int image_height)
{
    int i;
    int image_size = image_width * image_height;
    T *unwrapped_image_pointer = unwrapped_image;
    PIXELM<T> *pixel_pointer = pixel;
    
    for (i=0; i < image_size; i++)
    {
        *unwrapped_image_pointer = pixel_pointer->value;
        pixel_pointer++;
        unwrapped_image_pointer++;
    }
}

//the main function of the unwrapper
template <typename T>
        void unwrap2D(T* wrapped_image, T* UnwrappedImage, unsigned char* input_mask,
        int image_width, int image_height,
        int wrap_around_x, int wrap_around_y)
{
    params_t<T> params = {TWOPI, wrap_around_x, wrap_around_y, 0};
    unsigned char *extended_mask;
    PIXELM<T> *pixel;
    EDGE<T> *edge;
    int image_size = image_height * image_width;
    int No_of_Edges_initially = 2 * image_width * image_height;
    
    extended_mask = (unsigned char *) calloc(image_size, sizeof(unsigned char));
    pixel = (PIXELM<T> *) calloc(image_size, sizeof(PIXELM<T>));
    edge = (EDGE<T> *) calloc(No_of_Edges_initially, sizeof(EDGE<T>));
    
    //srand(time(NULL)); // MB new seed every time we initialize?
    
    extend_mask(input_mask, extended_mask, image_width, image_height, &params);
    initialisePIXELs(wrapped_image, input_mask, extended_mask, pixel, image_width, image_height);
    calculate_reliability(wrapped_image, pixel, image_width, image_height, &params);
    horizontalEDGEs(pixel, edge, image_width, image_height, &params);
    verticalEDGEs(pixel, edge, image_width, image_height, &params);
    
    if (params.no_of_edges > 0) //avoid crash MB
    {
        //sort the EDGEs depending on their reiability. The PIXELs with higher
        //relibility (small value) first
        quicker_sort(edge, edge + params.no_of_edges - 1);
        
        //std::cout << edge[image_size-1].reliab << std::endl;
        //std::cout << edge[0].pointer_1 << std::endl;
        //std::cout << edge[0].pointer_2 << std::endl;
        
        //gather PIXELs into groups
        gatherPIXELs(edge, &params);
        
        unwrapImage(pixel, image_width, image_height);
        
        //copy the image from PIXELM structure to the unwrapped phase array
        //passed to this function
        //TODO: replace by (cython?) function to directly write into numpy array ?
        returnImage(pixel, UnwrappedImage, image_width, image_height);
    }
    
    free(edge);
    free(pixel);
    free(extended_mask);
}
