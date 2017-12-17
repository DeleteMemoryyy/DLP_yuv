#ifndef CONVERT_H_
#define CONVERT_H_
#include <emmintrin.h>
#include <immintrin.h>
#include <mmintrin.h>

#include <assert.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <string>

using namespace std;

typedef unsigned char BYTE;
typedef short SHWORD;

#define MODE_NORM 0
#define MODE_MMX 1
#define MODE_SSE2 2
#define MODE_AVX 3

// typedef struct SRGB_C
// {
//     BYTE r;
//     BYTE g;
//     BYTE b;
// } RGB_C;

class YUV420
{
  public:
    BYTE *Y;
    BYTE *U;
    BYTE *V;
    SHWORD *YH;
    SHWORD *UH;
    SHWORD *VH;

    int h;
    int w;

    YUV420(int _h, int _w) : h(_h), w(_w)
    {
        assert(h % 2 == 0 && w % 2 == 0);
        int size = h * w;
        Y = new BYTE[size];
        U = new BYTE[size / 4];
        V = new BYTE[size / 4];
        YH = new SHWORD[size];
        UH = new SHWORD[size / 4];
        VH = new SHWORD[size / 4];
    }

    ~YUV420()
    {
        if (Y != NULL)
            delete[] Y;
        if (U != NULL)
            delete[] U;
        if (V != NULL)
            delete[] V;
        if (YH != NULL)
            delete[] YH;
        if (UH != NULL)
            delete[] UH;
        if (VH != NULL)
            delete[] VH;
    }
    void load_from_block(BYTE *block);
    void save_block(BYTE *block);
    void byte_to_shword();
    void shword_to_byte();
};

class RGB
{
  public:
    BYTE *R;
    BYTE *G;
    BYTE *B;
    SHWORD *RH;
    SHWORD *GH;
    SHWORD *BH;

    int h;
    int w;

    RGB(int _h, int _w) : h(_h), w(_w)
    {
        assert(h % 2 == 0 && w % 2 == 0);
        int size = h * w;
        R = new BYTE[size];
        G = new BYTE[size];
        B = new BYTE[size];
        RH = new SHWORD[size];
        GH = new SHWORD[size];
        BH = new SHWORD[size];
    }

    ~RGB()
    {
        if (R != NULL)
            delete[] R;
        if (G != NULL)
            delete[] G;
        if (B != NULL)
            delete[] B;
        if (RH != NULL)
            delete[] RH;
        if (GH != NULL)
            delete[] GH;
        if (BH != NULL)
            delete[] BH;
    }

    void byte_to_shword();
    void shword_to_byte();
    void desaturation();
};

#endif