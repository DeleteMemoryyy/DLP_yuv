#include "Convert.h"

int mode = MODE_NORM;

inline BYTE CONVERT_ADJUST(double tmp)
{
    return (BYTE)((tmp >= 0 && tmp <= 255) ? tmp : (tmp < 0 ? 0 : 255));
}

inline BYTE CONVERT_CUTOFF(long long tmp)
{
    return (BYTE)((tmp >= 0 && tmp <= 255) ? tmp : (tmp < 0 ? 0 : 255));
}

void RGB::byte_to_shword()
{
    int size = h * w;

    for (int i = 0; i < size; ++i)
        RH[i] = (SHWORD)R[i];
    for (int i = 0; i < size; ++i)
        GH[i] = (SHWORD)G[i];
    for (int i = 0; i < size; ++i)
        BH[i] = (SHWORD)B[i];
}

void RGB::shword_to_byte()
{
    int size = h * w;

    for (int i = 0; i < size; ++i)
        R[i] = CONVERT_CUTOFF(RH[i]);
    for (int i = 0; i < size; ++i)
        G[i] = CONVERT_CUTOFF(GH[i]);
    for (int i = 0; i < size; ++i)
        B[i] = CONVERT_CUTOFF(BH[i]);
}

void RGB::desaturation()
{
    int size = h * w;

    for (int i = 0; i < size; ++i)
        RH[i] = CONVERT_CUTOFF(RH[i]);
    for (int i = 0; i < size; ++i)
        GH[i] = CONVERT_CUTOFF(GH[i]);
    for (int i = 0; i < size; ++i)
        BH[i] = CONVERT_CUTOFF(BH[i]);
}

void YUV420::byte_to_shword()
{
    int size = h * w, uvsize = (h * w) / 4;

    for (int i = 0; i < size; ++i)
        YH[i] = (SHWORD)Y[i];
    for (int i = 0; i < uvsize; ++i)
        UH[i] = (SHWORD)U[i];
    for (int i = 0; i < uvsize; ++i)
        VH[i] = (SHWORD)V[i];
}

void YUV420::shword_to_byte()
{
    int size = h * w, uvsize = (h * w) / 4;

    for (int i = 0; i < size; ++i)
        Y[i] = CONVERT_CUTOFF(YH[i]);
    for (int i = 0; i < uvsize; ++i)
        U[i] = CONVERT_CUTOFF(UH[i]);
    for (int i = 0; i < uvsize; ++i)
        V[i] = CONVERT_CUTOFF(VH[i]);
}

void YUV420::load_from_block(BYTE *block)
{
    int ysize = h * w, uvsize = h * w / 4;
    for (int i = 0; i < ysize; ++i)
        {
            Y[i] = *block;
            block++;
        }
    for (int i = 0; i < uvsize; ++i)
        {
            U[i] = *block;
            block++;
        }
    for (int i = 0; i < uvsize; ++i)
        {
            V[i] = *block;
            block++;
        }
}

void YUV420::save_block(BYTE *block)
{
    int ysize = h * w, uvsize = (h / 2) * (w / 2);
    for (int i = 0; i < ysize; ++i)
        {
            *block = Y[i];
            block++;
        }
    for (int i = 0; i < uvsize; ++i)
        {
            *block = U[i];
            block++;
        }
    for (int i = 0; i < uvsize; ++i)
        {
            *block = V[i];
            block++;
        }
}

void YUV2RGB(YUV420 *yuv, RGB *rgb)
{
    int height = yuv->h;
    int width = yuv->w;

    SHWORD YUV_R[3] = {SHWORD(0.164383 * (1 << 16)), SHWORD(0.017232 * (1 << 16)), 0 * (1 << 16)};
    SHWORD YUV_G[3] = {SHWORD(0.164383 * (1 << 16)), SHWORD(-0.391762 * (1 << 16)),
                       SHWORD(-0.312968 * (1 << 16))};
    SHWORD YUV_B[3] = {SHWORD(0.164383 * (1 << 16)), 0 * (1 << 16), SHWORD(0.096027 * (1 << 16))};

    SHWORD *tmp_u = new SHWORD[height * width];
    SHWORD *tmp_v = new SHWORD[height * width];

    yuv->byte_to_shword();
    for (int i = 0, k = 0; i < height; ++i)
        for (int j = 0; j < width; ++j, ++k)
            {
                int uvi = (i / 2) * (width / 2) + (j / 2);
                tmp_u[k] = yuv->UH[uvi];
                tmp_v[k] = yuv->VH[uvi];
            }

    switch (mode)
        {
            case MODE_NORM:
                {
                    for (int i = 0; i < height; ++i)
                        for (int j = 0; j < width; ++j)
                            {
                                int yi = i * width + j, uvi = (i / 2) * (width / 2) + j / 2;
                                BYTE y = yuv->Y[yi], u = yuv->U[uvi], v = yuv->V[uvi];
                                rgb->R[yi] = CONVERT_ADJUST(y + 1.407 * (v - 128));
                                rgb->G[yi] =
                                    CONVERT_ADJUST(y - 0.344 * (u - 128) - 0.714 * (v - 128));
                                rgb->B[yi] = CONVERT_ADJUST(y + 1.772 * (u - 128));
                            }
                    rgb->byte_to_shword();
                }
                break;
            case MODE_MMX:
                {
                    __m64 tmp, tmp_data;

                    _mm_empty();

                    const __m64 OFFSET_128 = _mm_set_pi16(128, 128, 128, 128);
                    const __m64 OFFSET_16 = _mm_set_pi16(16, 16, 16, 16);
                    const __m64 Y_R = _mm_set_pi16(YUV_R[0], YUV_R[0], YUV_R[0], YUV_R[0]);
                    const __m64 U_R = _mm_set_pi16(YUV_R[1], YUV_R[1], YUV_R[1], YUV_R[1]);
                    // const __m64 V_R = _mm_set_pi16(YUV_R[2], YUV_R[2], YUV_R[2], YUV_R[2]);
                    const __m64 Y_G = _mm_set_pi16(YUV_G[0], YUV_G[0], YUV_G[0], YUV_G[0]);
                    const __m64 U_G = _mm_set_pi16(YUV_G[1], YUV_G[1], YUV_G[1], YUV_G[1]);
                    const __m64 V_G = _mm_set_pi16(YUV_G[2], YUV_G[2], YUV_G[2], YUV_G[2]);
                    const __m64 Y_B = _mm_set_pi16(YUV_B[0], YUV_B[0], YUV_B[0], YUV_B[0]);
                    // const __m64 U_B = _mm_set_pi16(YUV_B[1], YUV_B[1], YUV_B[1], YUV_B[1]);
                    const __m64 V_B = _mm_set_pi16(YUV_B[2], YUV_B[2], YUV_B[2], YUV_B[2]);

                    __m64 *dst = (__m64 *)rgb->RH;
                    __m64 *src_y = (__m64 *)yuv->YH;
                    __m64 *src_u = (__m64 *)tmp_u;
                    __m64 *src_v = (__m64 *)tmp_v;
                    for (int i = 0; i < (height * width) / 4; ++i)
                        {
                            tmp_data = _m_psubw(*src_y, OFFSET_16);  // (Y - 16)
                            tmp = _m_pmulhw(tmp_data, Y_R);          // R = (Y - 16) * 0.164383
                            *dst = _m_paddsw(tmp, tmp_data);         // R += Y - 16

                            tmp_data = _m_psubw(*src_u, OFFSET_128);  // (U - 128)
                            tmp = _m_pmulhw(tmp_data, U_R);           // (U - 128) * 0.017232
                            *dst = _m_paddsw(*dst, tmp);              // R += (U - 128) * 0.017232
                            tmp = _m_psllwi(tmp_data, 1);
                            *dst = _m_paddsw(*dst, tmp);  // R += (U - 128) << 1;

                            dst++;
                            src_y++;
                            src_u++;
                            src_v++;
                        }

                    dst = (__m64 *)rgb->GH;
                    src_y = (__m64 *)yuv->YH;
                    src_u = (__m64 *)tmp_u;
                    src_v = (__m64 *)tmp_v;
                    for (int i = 0; i < (height * width) / 4; ++i)
                        {
                            tmp_data = _m_psubw(*src_y, OFFSET_16);  // (Y - 16)
                            tmp = _m_pmulhw(tmp_data, Y_G);          // G = (Y - 16) * 0.164383
                            *dst = _m_paddsw(tmp, tmp_data);         // G += Y - 16

                            tmp_data = _m_psubw(*src_u, OFFSET_128);  // (U - 128)
                            tmp = _m_pmulhw(tmp_data, U_G);           // (U - 128) * (-0.391762)
                            *dst = _m_paddsw(*dst, tmp);  // G += (U - 128) * (-0.391762)

                            tmp_data = _m_psubw(*src_v, OFFSET_128);  // (V - 128)
                            tmp = _m_pmulhw(tmp_data, V_G);           // (V - 128) * (-0.312968)
                            *dst = _m_paddsw(*dst, tmp);  // G += (V - 128) * (-0.312968)
                            tmp = _m_psrawi(tmp_data, 1);
                            *dst = _m_psubsw(*dst, tmp);  // G -= (V - 128) >> 1;

                            dst++;
                            src_y++;
                            src_u++;
                            src_v++;
                        }

                    dst = (__m64 *)rgb->BH;
                    src_y = (__m64 *)yuv->YH;
                    src_u = (__m64 *)tmp_u;
                    src_v = (__m64 *)tmp_v;
                    for (int i = 0; i < (height * width) / 4; ++i)
                        {
                            tmp_data = _m_psubw(*src_y, OFFSET_16);  // (Y - 16)
                            tmp = _m_pmulhw(tmp_data, Y_B);          // B = (Y - 16) * 0.164383
                            *dst = _m_paddsw(tmp, tmp_data);         // B += Y - 16

                            tmp_data = _m_psubw(*src_v, OFFSET_128);  // (V - 128)
                            tmp = _m_pmulhw(tmp_data, V_B);           // (V - 128) * 0.096027
                            *dst = _m_paddsw(*dst, tmp);              // B += (V - 128) * 0.096027
                            tmp = _m_psrawi(tmp_data, 1);
                            *dst = _m_paddsw(*dst, tmp);  // B += (V - 128) >> 1;
                            tmp = _m_psllwi(tmp_data, 1);
                            *dst = _m_paddsw(*dst, tmp);  // B += (V - 128) << 1;

                            dst++;
                            src_y++;
                            src_u++;
                            src_v++;
                        }
                    _mm_empty();

                    rgb->desaturation();
                }
                break;
            case MODE_SSE2:
                {
                    __m128i tmp, tmp_data;

                    const __m128i OFFSET_128 =
                        _mm_set_epi16(128, 128, 128, 128, 128, 128, 128, 128);
                    const __m128i OFFSET_16 = _mm_set_epi16(16, 16, 16, 16, 16, 16, 16, 16);
                    const __m128i Y_R = _mm_set_epi16(YUV_R[0], YUV_R[0], YUV_R[0], YUV_R[0],
                                                      YUV_R[0], YUV_R[0], YUV_R[0], YUV_R[0]);
                    const __m128i U_R = _mm_set_epi16(YUV_R[1], YUV_R[1], YUV_R[1], YUV_R[1],
                                                      YUV_R[1], YUV_R[1], YUV_R[1], YUV_R[1]);
                    // const __m128i V_R = _mm_set_epi16(YUV_R[2], YUV_R[2], YUV_R[2], YUV_R[2],
                    //                                   YUV_R[2], YUV_R[2], YUV_R[2], YUV_R[2]);
                    const __m128i Y_G = _mm_set_epi16(YUV_G[0], YUV_G[0], YUV_G[0], YUV_G[0],
                                                      YUV_G[0], YUV_G[0], YUV_G[0], YUV_G[0]);
                    const __m128i U_G = _mm_set_epi16(YUV_G[1], YUV_G[1], YUV_G[1], YUV_G[1],
                                                      YUV_G[1], YUV_G[1], YUV_G[1], YUV_G[1]);
                    const __m128i V_G = _mm_set_epi16(YUV_G[2], YUV_G[2], YUV_G[2], YUV_G[2],
                                                      YUV_G[2], YUV_G[2], YUV_G[2], YUV_G[2]);
                    const __m128i Y_B = _mm_set_epi16(YUV_B[0], YUV_B[0], YUV_B[0], YUV_B[0],
                                                      YUV_B[0], YUV_B[0], YUV_B[0], YUV_B[0]);
                    // const __m128i U_B = _mm_set_epi16(YUV_B[1], YUV_B[1], YUV_B[1], YUV_B[1],
                    //                                   YUV_B[1], YUV_B[1], YUV_B[1], YUV_B[1]);
                    const __m128i V_B = _mm_set_epi16(YUV_B[2], YUV_B[2], YUV_B[2], YUV_B[2],
                                                      YUV_B[2], YUV_B[2], YUV_B[2], YUV_B[2]);
                    _mm_empty();

                    __m128i *dst = (__m128i *)rgb->RH;
                    __m128i *src_y = (__m128i *)yuv->YH;
                    __m128i *src_u = (__m128i *)tmp_u;
                    __m128i *src_v = (__m128i *)tmp_v;
                    for (int i = 0; i < (height * width) / 8; ++i)
                        {
                            tmp_data = _mm_subs_epi16(*src_y, OFFSET_16);  // (Y - 16)
                            tmp = _mm_mulhi_epi16(tmp_data, Y_R);  // R = (Y - 16) * 0.164383
                            *dst = _mm_adds_epi16(tmp, tmp_data);  // R += Y - 16

                            tmp_data = _mm_subs_epi16(*src_u, OFFSET_128);  // (U - 128)
                            tmp = _mm_mulhi_epi16(tmp_data, U_R);           // (U - 128) * 0.017232
                            *dst = _mm_adds_epi16(*dst, tmp);  // R += (U - 128) * 0.017232
                            tmp = _mm_slli_epi16(tmp_data, 1);
                            *dst = _mm_adds_epi16(*dst, tmp);  // R += (U - 128) << 1;

                            dst++;
                            src_y++;
                            src_u++;
                            src_v++;
                        }

                    dst = (__m128i *)rgb->GH;
                    src_y = (__m128i *)yuv->YH;
                    src_u = (__m128i *)tmp_u;
                    src_v = (__m128i *)tmp_v;
                    for (int i = 0; i < (height * width) / 8; ++i)
                        {
                            tmp_data = _mm_subs_epi16(*src_y, OFFSET_16);  // (Y - 16)
                            tmp = _mm_mulhi_epi16(tmp_data, Y_G);  // G = (Y - 16) * 0.164383
                            *dst = _mm_adds_epi16(tmp, tmp_data);  // G += Y - 16

                            tmp_data = _mm_subs_epi16(*src_u, OFFSET_128);  // (U - 128)
                            tmp = _mm_mulhi_epi16(tmp_data, U_G);  // (U - 128) * (-0.391762)
                            *dst = _mm_adds_epi16(*dst, tmp);      // G += (U - 128) * (-0.391762)

                            tmp_data = _mm_subs_epi16(*src_v, OFFSET_128);  // (V - 128)
                            tmp = _mm_mulhi_epi16(tmp_data, V_G);  // (V - 128) * (-0.312968)
                            *dst = _mm_adds_epi16(*dst, tmp);      // G += (V - 128) * (-0.312968)
                            tmp = _mm_srai_epi16(tmp_data, 1);
                            *dst = _mm_subs_epi16(*dst, tmp);  // G -= (V - 128) >> 1;

                            dst++;
                            src_y++;
                            src_u++;
                            src_v++;
                        }

                    dst = (__m128i *)rgb->BH;
                    src_y = (__m128i *)yuv->YH;
                    src_u = (__m128i *)tmp_u;
                    src_v = (__m128i *)tmp_v;
                    for (int i = 0; i < (height * width) / 8; ++i)
                        {
                            tmp_data = _mm_subs_epi16(*src_y, OFFSET_16);  // (Y - 16)
                            tmp = _mm_mulhi_epi16(tmp_data, Y_B);  // B = (Y - 16) * 0.164383
                            *dst = _mm_adds_epi16(tmp, tmp_data);  // B += Y - 16

                            tmp_data = _mm_subs_epi16(*src_v, OFFSET_128);  // (V - 128)
                            tmp = _mm_mulhi_epi16(tmp_data, V_B);           // (V - 128) * 0.096027
                            *dst = _mm_adds_epi16(*dst, tmp);  // B += (V - 128) * 0.096027
                            tmp = _mm_srai_epi16(tmp_data, 1);
                            *dst = _mm_adds_epi16(*dst, tmp);  // G += (V - 128) >> 1;
                            tmp = _mm_slli_epi16(tmp_data, 1);
                            *dst = _mm_adds_epi16(*dst, tmp);  // G += (V - 128) << 1;

                            dst++;
                            src_y++;
                            src_u++;
                            src_v++;
                        }
                    _mm_empty();

                    rgb->desaturation();
                }
                break;
            case MODE_AVX:
                {
                    __m256i tmp, tmp_data, tmp_store, tmp_dst;

                    const __m256i OFFSET_128 =
                        _mm256_set_epi16(128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
                                         128, 128, 128, 128);
                    const __m256i OFFSET_16 = _mm256_set_epi16(16, 16, 16, 16, 16, 16, 16, 16, 16,
                                                               16, 16, 16, 16, 16, 16, 16);
                    const __m256i Y_R =
                        _mm256_set_epi16(YUV_R[0], YUV_R[0], YUV_R[0], YUV_R[0], YUV_R[0], YUV_R[0],
                                         YUV_R[0], YUV_R[0], YUV_R[0], YUV_R[0], YUV_R[0], YUV_R[0],
                                         YUV_R[0], YUV_R[0], YUV_R[0], YUV_R[0]);
                    const __m256i U_R =
                        _mm256_set_epi16(YUV_R[1], YUV_R[1], YUV_R[1], YUV_R[1], YUV_R[1], YUV_R[1],
                                         YUV_R[1], YUV_R[1], YUV_R[1], YUV_R[1], YUV_R[1], YUV_R[1],
                                         YUV_R[1], YUV_R[1], YUV_R[1], YUV_R[1]);
                    // const __m256i V_R =
                    //     _mm256_set_epi16(YUV_R[2], YUV_R[2], YUV_R[2], YUV_R[2], YUV_R[2],
                    //     YUV_R[2],
                    //                      YUV_R[2], YUV_R[2], YUV_R[2], YUV_R[2], YUV_R[2],
                    //                      YUV_R[2], YUV_R[2], YUV_R[2], YUV_R[2], YUV_R[2]);
                    const __m256i Y_G =
                        _mm256_set_epi16(YUV_G[0], YUV_G[0], YUV_G[0], YUV_G[0], YUV_G[0], YUV_G[0],
                                         YUV_G[0], YUV_G[0], YUV_G[0], YUV_G[0], YUV_G[0], YUV_G[0],
                                         YUV_G[0], YUV_G[0], YUV_G[0], YUV_G[0]);
                    const __m256i U_G =
                        _mm256_set_epi16(YUV_G[1], YUV_G[1], YUV_G[1], YUV_G[1], YUV_G[1], YUV_G[1],
                                         YUV_G[1], YUV_G[1], YUV_G[1], YUV_G[1], YUV_G[1], YUV_G[1],
                                         YUV_G[1], YUV_G[1], YUV_G[1], YUV_G[1]);
                    const __m256i V_G =
                        _mm256_set_epi16(YUV_G[2], YUV_G[2], YUV_G[2], YUV_G[2], YUV_G[2], YUV_G[2],
                                         YUV_G[2], YUV_G[2], YUV_G[2], YUV_G[2], YUV_G[2], YUV_G[2],
                                         YUV_G[2], YUV_G[2], YUV_G[2], YUV_G[2]);
                    const __m256i Y_B =
                        _mm256_set_epi16(YUV_B[0], YUV_B[0], YUV_B[0], YUV_B[0], YUV_B[0], YUV_B[0],
                                         YUV_B[0], YUV_B[0], YUV_B[0], YUV_B[0], YUV_B[0], YUV_B[0],
                                         YUV_B[0], YUV_B[0], YUV_B[0], YUV_B[0]);
                    // const __m256i U_B =
                    //     _mm256_set_epi16(YUV_B[1], YUV_B[1], YUV_B[1], YUV_B[1], YUV_B[1],
                    //     YUV_B[1],
                    //                      YUV_B[1], YUV_B[1], YUV_B[1], YUV_B[1], YUV_B[1],
                    //                      YUV_B[1], YUV_B[1], YUV_B[1], YUV_B[1], YUV_B[1]);
                    const __m256i V_B =
                        _mm256_set_epi16(YUV_B[2], YUV_B[2], YUV_B[2], YUV_B[2], YUV_B[2], YUV_B[2],
                                         YUV_B[2], YUV_B[2], YUV_B[2], YUV_B[2], YUV_B[2], YUV_B[2],
                                         YUV_B[2], YUV_B[2], YUV_B[2], YUV_B[2]);
                    _mm_empty();

                    __m256i *dst = (__m256i *)rgb->RH;
                    __m256i *src_y = (__m256i *)yuv->YH;
                    __m256i *src_u = (__m256i *)tmp_u;
                    __m256i *src_v = (__m256i *)tmp_v;
                    for (int i = 0; i < (height * width) / 16; ++i)
                        {
                            tmp_store = _mm256_loadu_si256(src_y);
                            tmp_data = _mm256_subs_epi16(tmp_store, OFFSET_16);  // (Y - 16)
                            tmp = _mm256_mulhi_epi16(tmp_data, Y_R);  // R = (Y - 16) * 0.164383
                            tmp_dst = _mm256_adds_epi16(tmp, tmp_data);  // R += Y - 16

                            tmp_store = _mm256_loadu_si256(src_u);
                            tmp_data = _mm256_subs_epi16(tmp_store, OFFSET_128);  // (U - 128)
                            tmp = _mm256_mulhi_epi16(tmp_data, U_R);    // (U - 128) * 0.017232
                            tmp_dst = _mm256_adds_epi16(tmp_dst, tmp);  // R += (U - 128) * 0.017232

                            tmp = _mm256_slli_epi16(tmp_data, 1);
                            tmp_dst = _mm256_adds_epi16(tmp_dst, tmp);  // R += (U - 128) << 1;
                            _mm256_storeu_si256(dst, tmp_dst);

                            dst++;
                            src_y++;
                            src_u++;
                            src_v++;
                        }

                    dst = (__m256i *)rgb->GH;
                    src_y = (__m256i *)yuv->YH;
                    src_u = (__m256i *)tmp_u;
                    src_v = (__m256i *)tmp_v;
                    for (int i = 0; i < (height * width) / 16; ++i)
                        {
                            tmp_store = _mm256_loadu_si256(src_y);
                            tmp_data = _mm256_subs_epi16(tmp_store, OFFSET_16);  // (Y - 16)
                            tmp = _mm256_mulhi_epi16(tmp_data, Y_G);     // G = (Y - 16) * 0.164383
                            tmp_dst = _mm256_adds_epi16(tmp, tmp_data);  // G += Y - 16

                            tmp_store = _mm256_loadu_si256(src_u);
                            tmp_data = _mm256_subs_epi16(tmp_store, OFFSET_128);  // (U - 128)
                            tmp = _mm256_mulhi_epi16(tmp_data, U_G);  // (U - 128) * (-0.391762)
                            tmp_dst =
                                _mm256_adds_epi16(tmp_dst, tmp);  // G += (U - 128) * (-0.391762)

                            tmp_store = _mm256_loadu_si256(src_v);
                            tmp_data = _mm256_subs_epi16(tmp_store, OFFSET_128);  // (V - 128)
                            tmp = _mm256_mulhi_epi16(tmp_data, V_G);  // (V - 128) * (-0.312968)
                            tmp_dst =
                                _mm256_adds_epi16(tmp_dst, tmp);  // G += (V - 128) * (-0.312968)
                            tmp = _mm256_srai_epi16(tmp_data, 1);
                            tmp_dst = _mm256_subs_epi16(tmp_dst, tmp);  // G -= (V - 128) >> 1;
                            _mm256_storeu_si256(dst, tmp_dst);

                            dst++;
                            src_y++;
                            src_u++;
                            src_v++;
                        }

                    dst = (__m256i *)rgb->BH;
                    src_y = (__m256i *)yuv->YH;
                    src_u = (__m256i *)tmp_u;
                    src_v = (__m256i *)tmp_v;
                    for (int i = 0; i < (height * width) / 16; ++i)
                        {

                            tmp_store = _mm256_loadu_si256(src_y);
                            tmp_data = _mm256_subs_epi16(tmp_store, OFFSET_16);  // (Y - 16)
                            tmp = _mm256_mulhi_epi16(tmp_data, Y_B);     // B = (Y - 16) * 0.164383
                            tmp_dst = _mm256_adds_epi16(tmp, tmp_data);  // B += Y - 16

                            tmp_store = _mm256_loadu_si256(src_v);
                            tmp_data = _mm256_subs_epi16(tmp_store, OFFSET_128);  // (V - 128)
                            tmp = _mm256_mulhi_epi16(tmp_data, V_B);  // (V - 128) * 0.096027

                            tmp_dst = _mm256_adds_epi16(tmp_dst, tmp);  // B += (V - 128) * 0.096027
                            tmp = _mm256_srai_epi16(tmp_data, 1);
                            tmp_dst = _mm256_adds_epi16(tmp_dst, tmp);  // G += (V - 128) >> 1;
                            tmp = _mm256_slli_epi16(tmp_data, 1);
                            tmp_dst = _mm256_adds_epi16(tmp_dst, tmp);  // G += (V - 128) << 1;
                            _mm256_storeu_si256(dst, tmp_dst);

                            dst++;
                            src_y++;
                            src_u++;
                            src_v++;
                        }
                    _mm_empty();

                    rgb->desaturation();
                }
                break;
            default:
                break;
        }
    delete[] tmp_u;
    delete[] tmp_v;
}

void RGB2YUV(RGB *rgb, YUV420 *yuv)
{
    int height = yuv->h;
    int width = yuv->w;
    int uvsize = (height * width) / 4;

    int *tru = new int[uvsize];
    int *trv = new int[uvsize];
    memset(tru, 0, sizeof(int) * uvsize);
    memset(trv, 0, sizeof(int) * uvsize);

    SHWORD RGB_Y[3] = {SHWORD(0.256788 * (1 << 16)), SHWORD(0.004129 * (1 << 16)),
                       SHWORD(0.097906 * (1 << 16))};  // offset: 0 -0.5 0
    SHWORD RGB_U[3] = {SHWORD(0.439216 * (1 << 16)), SHWORD(-0.367788 * (1 << 16)),
                       SHWORD(-0.071427 * (1 << 16))};  // offset: 0 0 0
    SHWORD RGB_V[3] = {SHWORD(-0.148223 * (1 << 16)), SHWORD(-0.290993 * (1 << 16)),
                       SHWORD(0.439216 * (1 << 16))};  // offset: 0 0 0
    SHWORD *tmp_r = new SHWORD[uvsize];
    SHWORD *tmp_g = new SHWORD[uvsize];
    SHWORD *tmp_b = new SHWORD[uvsize];

    for (int i = 0, k = 0; i < height; i += 2)
        {
            for (int j = 0; j < width; j += 2, ++k)
                {
                    tmp_r[k] = rgb->RH[i * width + j];
                    tmp_g[k] = rgb->GH[i * width + j];
                    tmp_b[k] = rgb->BH[i * width + j];
                }
        }

    switch (mode)
        {
            case MODE_NORM:
                {
                    rgb->shword_to_byte();
                    for (int i = 0; i < height; ++i)
                        for (int j = 0; j < width; ++j)
                            {
                                int yi = i * width + j, uvi = (i / 2) * (width / 2) + j / 2;
                                BYTE r = rgb->R[yi], g = rgb->G[yi], b = rgb->B[yi];
                                yuv->Y[yi] = CONVERT_ADJUST(0.257 * r + 0.504 * g + 0.098 * b + 16);
                                tru[uvi] +=
                                    CONVERT_ADJUST(-0.148 * r - 0.291 * g + 0.439 * b + 128);
                                trv[uvi] += CONVERT_ADJUST(0.439 * r - 0.368 * g - 0.071 * b + 128);
                            }
                    for (int i = 0; i < (height * width) / 4; ++i)
                        yuv->U[i] = (BYTE)(tru[i] / 4);
                    for (int i = 0; i < (height * width) / 4; ++i)
                        yuv->V[i] = (BYTE)(trv[i] / 4);
                }
                break;
            case MODE_MMX:
                {
                    __m64 tmp;
                    _mm_empty();

                    const __m64 OFFSET_128 = _mm_set_pi16(128, 128, 128, 128);
                    const __m64 OFFSET_16 = _mm_set_pi16(16, 16, 16, 16);
                    const __m64 R_Y = _mm_set_pi16(RGB_Y[0], RGB_Y[0], RGB_Y[0], RGB_Y[0]);
                    const __m64 G_Y = _mm_set_pi16(RGB_Y[1], RGB_Y[1], RGB_Y[1], RGB_Y[1]);
                    const __m64 B_Y = _mm_set_pi16(RGB_Y[2], RGB_Y[2], RGB_Y[2], RGB_Y[2]);
                    const __m64 R_U = _mm_set_pi16(RGB_U[0], RGB_U[0], RGB_U[0], RGB_U[0]);
                    const __m64 G_U = _mm_set_pi16(RGB_U[1], RGB_U[1], RGB_U[1], RGB_U[1]);
                    const __m64 B_U = _mm_set_pi16(RGB_U[2], RGB_U[2], RGB_U[2], RGB_U[2]);
                    const __m64 R_V = _mm_set_pi16(RGB_V[0], RGB_V[0], RGB_V[0], RGB_V[0]);
                    const __m64 G_V = _mm_set_pi16(RGB_V[1], RGB_V[1], RGB_V[1], RGB_V[1]);
                    const __m64 B_V = _mm_set_pi16(RGB_V[2], RGB_V[2], RGB_V[2], RGB_V[2]);

                    __m64 *dst = (__m64 *)yuv->YH;
                    __m64 *src_r = (__m64 *)rgb->RH;
                    __m64 *src_g = (__m64 *)rgb->GH;
                    __m64 *src_b = (__m64 *)rgb->BH;
                    for (int i = 0; i < (height * width) / 4; i++)
                        {
                            *dst = _m_pmulhw(*src_r, R_Y);  // Y = R * 0.256788

                            tmp = _m_pmulhw(*src_g, G_Y);
                            *dst = _m_paddsw(tmp, *dst);  // Y += G * 0.004129
                            tmp = _m_psrlwi(*src_g, 1);
                            *dst = _m_paddsw(tmp, *dst);  // Y += G >> 1;

                            tmp = _m_pmulhw(*src_b, B_Y);
                            *dst = _m_paddsw(tmp, *dst);  // Y += B * 0.097906

                            *dst = _m_paddsw(*dst, OFFSET_16);  // Y += 16

                            dst++;
                            src_r++;
                            src_g++;
                            src_b++;
                        }

                    dst = (__m64 *)yuv->UH;
                    src_r = (__m64 *)tmp_r;
                    src_g = (__m64 *)tmp_g;
                    src_b = (__m64 *)tmp_b;
                    for (int i = 0; i < height * width / 16; i++)
                        {
                            *dst = _m_pmulhw(*src_r, R_U);  // U = R * 0.439216

                            tmp = _m_pmulhw(*src_g, G_U);
                            *dst = _m_paddsw(tmp, *dst);  // U += G * (-0.367788)

                            tmp = _m_pmulhw(*src_b, B_U);
                            *dst = _m_paddsw(tmp, *dst);  // U += B * (-0.071427)

                            *dst = _m_paddsw(*dst, OFFSET_128);  // U += 128

                            dst++;
                            src_r++;
                            src_g++;
                            src_b++;
                        }

                    dst = (__m64 *)yuv->VH;
                    src_r = (__m64 *)tmp_r;
                    src_g = (__m64 *)tmp_g;
                    src_b = (__m64 *)tmp_b;
                    for (int i = 0; i < (width * height) / 16; i++)
                        {
                            *dst = _m_pmulhw(*src_r, R_V);  // V = R * (-0.148223)

                            tmp = _m_pmulhw(*src_g, G_V);
                            *dst = _m_paddsw(tmp, *dst);  // V += G * (-0.290993)

                            tmp = _m_pmulhw(*src_b, B_V);
                            *dst = _m_paddsw(tmp, *dst);  // V += B * (0.439216)

                            *dst = _m_paddsw(*dst, OFFSET_128);  // V += 128

                            dst++;
                            src_r++;
                            src_g++;
                            src_b++;
                        }
                    _mm_empty();

                    yuv->shword_to_byte();
                }
                break;
            case MODE_SSE2:
                {
                    __m128i tmp;

                    _mm_empty();

                    const __m128i OFFSET_128 =
                        _mm_set_epi16(128, 128, 128, 128, 128, 128, 128, 128);
                    const __m128i OFFSET_16 = _mm_set_epi16(16, 16, 16, 16, 16, 16, 16, 16);
                    const __m128i R_Y = _mm_set_epi16(RGB_Y[0], RGB_Y[0], RGB_Y[0], RGB_Y[0],
                                                      RGB_Y[0], RGB_Y[0], RGB_Y[0], RGB_Y[0]);
                    const __m128i G_Y = _mm_set_epi16(RGB_Y[1], RGB_Y[1], RGB_Y[1], RGB_Y[1],
                                                      RGB_Y[1], RGB_Y[1], RGB_Y[1], RGB_Y[1]);
                    const __m128i B_Y = _mm_set_epi16(RGB_Y[2], RGB_Y[2], RGB_Y[2], RGB_Y[2],
                                                      RGB_Y[2], RGB_Y[2], RGB_Y[2], RGB_Y[2]);
                    const __m128i R_U = _mm_set_epi16(RGB_U[0], RGB_U[0], RGB_U[0], RGB_U[0],
                                                      RGB_U[0], RGB_U[0], RGB_U[0], RGB_U[0]);
                    const __m128i G_U = _mm_set_epi16(RGB_U[1], RGB_U[1], RGB_U[1], RGB_U[1],
                                                      RGB_U[1], RGB_U[1], RGB_U[1], RGB_U[1]);
                    const __m128i B_U = _mm_set_epi16(RGB_U[2], RGB_U[2], RGB_U[2], RGB_U[2],
                                                      RGB_U[2], RGB_U[2], RGB_U[2], RGB_U[2]);
                    const __m128i R_V = _mm_set_epi16(RGB_V[0], RGB_V[0], RGB_V[0], RGB_V[0],
                                                      RGB_V[0], RGB_V[0], RGB_V[0], RGB_V[0]);
                    const __m128i G_V = _mm_set_epi16(RGB_V[1], RGB_V[1], RGB_V[1], RGB_V[1],
                                                      RGB_V[1], RGB_V[1], RGB_V[1], RGB_V[1]);
                    const __m128i B_V = _mm_set_epi16(RGB_V[2], RGB_V[2], RGB_V[2], RGB_V[2],
                                                      RGB_V[2], RGB_V[2], RGB_V[2], RGB_V[2]);

                    __m128i *dst = (__m128i *)yuv->YH;
                    __m128i *src_r = (__m128i *)rgb->RH;
                    __m128i *src_g = (__m128i *)rgb->GH;
                    __m128i *src_b = (__m128i *)rgb->BH;
                    for (int i = 0; i < (width * height) / 8; i++)
                        {
                            *dst = _mm_mulhi_epi16(*src_r, R_Y);  // Y = R * 0.256788

                            tmp = _mm_mulhi_epi16(*src_g, G_Y);
                            *dst = _mm_adds_epi16(tmp, *dst);  // Y += G * 0.004129
                            tmp = _mm_srli_epi16(*src_g, 1);
                            *dst = _mm_adds_epi16(tmp, *dst);  // Y += G >> 1;

                            tmp = _mm_mulhi_epi16(*src_b, B_Y);
                            *dst = _mm_adds_epi16(tmp, *dst);  // Y += B * 0.097906

                            *dst = _mm_adds_epi16(*dst, OFFSET_16);  // Y += 16

                            dst++;
                            src_r++;
                            src_g++;
                            src_b++;
                        }

                    dst = (__m128i *)yuv->UH;
                    src_r = (__m128i *)tmp_r;
                    src_g = (__m128i *)tmp_g;
                    src_b = (__m128i *)tmp_b;
                    for (int i = 0; i < (width * height) / 32; i++)
                        {
                            *dst = _mm_mulhi_epi16(*src_r, R_U);  // U = R * 0.439216

                            tmp = _mm_mulhi_epi16(*src_g, G_U);
                            *dst = _mm_adds_epi16(tmp, *dst);  // U += G * (-0.367788)

                            tmp = _mm_mulhi_epi16(*src_b, B_U);
                            *dst = _mm_adds_epi16(tmp, *dst);  // U += B * (-0.071427)

                            *dst = _mm_adds_epi16(*dst, OFFSET_128);  // U += 128

                            dst++;
                            src_r++;
                            src_g++;
                            src_b++;
                        }

                    dst = (__m128i *)yuv->VH;
                    src_r = (__m128i *)tmp_r;
                    src_g = (__m128i *)tmp_g;
                    src_b = (__m128i *)tmp_b;
                    for (int i = 0; i < (width * height) / 32; i++)
                        {
                            *dst = _mm_mulhi_epi16(*src_r, R_V);  // V = R * (-0.148223)

                            tmp = _mm_mulhi_epi16(*src_g, G_V);
                            *dst = _mm_adds_epi16(tmp, *dst);  // V += G * (-0.290993)

                            tmp = _mm_mulhi_epi16(*src_b, B_V);
                            *dst = _mm_adds_epi16(tmp, *dst);  // V += B * (0.439216)

                            *dst = _mm_adds_epi16(*dst, OFFSET_128);  // V += 128

                            dst++;
                            src_r++;
                            src_g++;
                            src_b++;
                        }

                    _mm_empty();
                    yuv->shword_to_byte();
                }
                break;
            case MODE_AVX:
                {
                    __m256i tmp, tmp_store, tmp_dst;

                    _mm_empty();

                    const __m256i OFFSET_128 =
                        _mm256_set_epi16(128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
                                         128, 128, 128, 128);
                    const __m256i OFFSET_16 = _mm256_set_epi16(16, 16, 16, 16, 16, 16, 16, 16, 16,
                                                               16, 16, 16, 16, 16, 16, 16);
                    const __m256i R_Y =
                        _mm256_set_epi16(RGB_Y[0], RGB_Y[0], RGB_Y[0], RGB_Y[0], RGB_Y[0], RGB_Y[0],
                                         RGB_Y[0], RGB_Y[0], RGB_Y[0], RGB_Y[0], RGB_Y[0], RGB_Y[0],
                                         RGB_Y[0], RGB_Y[0], RGB_Y[0], RGB_Y[0]);
                    const __m256i G_Y =
                        _mm256_set_epi16(RGB_Y[1], RGB_Y[1], RGB_Y[1], RGB_Y[1], RGB_Y[1], RGB_Y[1],
                                         RGB_Y[1], RGB_Y[1], RGB_Y[1], RGB_Y[1], RGB_Y[1], RGB_Y[1],
                                         RGB_Y[1], RGB_Y[1], RGB_Y[1], RGB_Y[1]);
                    const __m256i B_Y =
                        _mm256_set_epi16(RGB_Y[2], RGB_Y[2], RGB_Y[2], RGB_Y[2], RGB_Y[2], RGB_Y[2],
                                         RGB_Y[2], RGB_Y[2], RGB_Y[2], RGB_Y[2], RGB_Y[2], RGB_Y[2],
                                         RGB_Y[2], RGB_Y[2], RGB_Y[2], RGB_Y[2]);
                    const __m256i R_U =
                        _mm256_set_epi16(RGB_U[0], RGB_U[0], RGB_U[0], RGB_U[0], RGB_U[0], RGB_U[0],
                                         RGB_U[0], RGB_U[0], RGB_U[0], RGB_U[0], RGB_U[0], RGB_U[0],
                                         RGB_U[0], RGB_U[0], RGB_U[0], RGB_U[0]);
                    const __m256i G_U =
                        _mm256_set_epi16(RGB_U[1], RGB_U[1], RGB_U[1], RGB_U[1], RGB_U[1], RGB_U[1],
                                         RGB_U[1], RGB_U[1], RGB_U[1], RGB_U[1], RGB_U[1], RGB_U[1],
                                         RGB_U[1], RGB_U[1], RGB_U[1], RGB_U[1]);
                    const __m256i B_U =
                        _mm256_set_epi16(RGB_U[2], RGB_U[2], RGB_U[2], RGB_U[2], RGB_U[2], RGB_U[2],
                                         RGB_U[2], RGB_U[2], RGB_U[2], RGB_U[2], RGB_U[2], RGB_U[2],
                                         RGB_U[2], RGB_U[2], RGB_U[2], RGB_U[2]);
                    const __m256i R_V =
                        _mm256_set_epi16(RGB_V[0], RGB_V[0], RGB_V[0], RGB_V[0], RGB_V[0], RGB_V[0],
                                         RGB_V[0], RGB_V[0], RGB_V[0], RGB_V[0], RGB_V[0], RGB_V[0],
                                         RGB_V[0], RGB_V[0], RGB_V[0], RGB_V[0]);
                    const __m256i G_V =
                        _mm256_set_epi16(RGB_V[1], RGB_V[1], RGB_V[1], RGB_V[1], RGB_V[1], RGB_V[1],
                                         RGB_V[1], RGB_V[1], RGB_V[1], RGB_V[1], RGB_V[1], RGB_V[1],
                                         RGB_V[1], RGB_V[1], RGB_V[1], RGB_V[1]);
                    const __m256i B_V =
                        _mm256_set_epi16(RGB_V[2], RGB_V[2], RGB_V[2], RGB_V[2], RGB_V[2], RGB_V[2],
                                         RGB_V[2], RGB_V[2], RGB_V[2], RGB_V[2], RGB_V[2], RGB_V[2],
                                         RGB_V[2], RGB_V[2], RGB_V[2], RGB_V[2]);

                    __m256i *dst = (__m256i *)yuv->YH;
                    __m256i *src_r = (__m256i *)rgb->RH;
                    __m256i *src_g = (__m256i *)rgb->GH;
                    __m256i *src_b = (__m256i *)rgb->BH;
                    for (int i = 0; i < (width * height) / 16; i++)
                        {
                            tmp_store = _mm256_loadu_si256(src_r);
                            tmp_dst = _mm256_mulhi_epi16(tmp_store, R_Y);  // Y = R * 0.256788

                            tmp_store = _mm256_loadu_si256(src_g);
                            tmp = _mm256_mulhi_epi16(tmp_store, G_Y);
                            tmp_dst = _mm256_adds_epi16(tmp, tmp_dst);  // Y += G * 0.004129
                            tmp = _mm256_srli_epi16(tmp_store, 1);
                            tmp_dst = _mm256_adds_epi16(tmp, tmp_dst);  // Y += G >> 1;

                            tmp_store = _mm256_loadu_si256(src_g);
                            tmp = _mm256_mulhi_epi16(tmp_store, B_Y);
                            tmp_dst = _mm256_adds_epi16(tmp, tmp_dst);  // Y += B * 0.097906

                            tmp_dst = _mm256_adds_epi16(tmp_dst, OFFSET_16);  // Y += 16
                            _mm256_storeu_si256(dst, tmp_dst);

                            dst++;
                            src_r++;
                            src_g++;
                            src_b++;
                        }

                    dst = (__m256i *)yuv->UH;
                    src_r = (__m256i *)tmp_r;
                    src_g = (__m256i *)tmp_g;
                    src_b = (__m256i *)tmp_b;
                    for (int i = 0; i < (width * height) / 64; i++)
                        {
                            tmp_store = _mm256_loadu_si256(src_r);
                            tmp_dst = _mm256_mulhi_epi16(tmp_store, R_U);  // U = R * 0.439216

                            tmp_store = _mm256_loadu_si256(src_g);
                            tmp = _mm256_mulhi_epi16(tmp_store, G_U);
                            tmp_dst = _mm256_adds_epi16(tmp, tmp_dst);  // U += G * (-0.367788)

                            tmp_store = _mm256_loadu_si256(src_b);
                            tmp = _mm256_mulhi_epi16(tmp_store, B_U);
                            tmp_dst = _mm256_adds_epi16(tmp, tmp_dst);  // U += B * (-0.071427)

                            tmp_dst = _mm256_adds_epi16(tmp_dst, OFFSET_128);  // U += 128
                            _mm256_storeu_si256(dst, tmp_dst);

                            dst++;
                            src_r++;
                            src_g++;
                            src_b++;
                        }

                    dst = (__m256i *)yuv->VH;
                    src_r = (__m256i *)tmp_r;
                    src_g = (__m256i *)tmp_g;
                    src_b = (__m256i *)tmp_b;
                    for (int i = 0; i < (width * height) / 64; i++)
                        {
                            tmp_store = _mm256_loadu_si256(src_r);
                            tmp_dst = _mm256_mulhi_epi16(tmp_store, R_V);  // V = R * (-0.148223)

                            tmp_store = _mm256_loadu_si256(src_g);
                            tmp = _mm256_mulhi_epi16(tmp_store, G_V);
                            tmp_dst = _mm256_adds_epi16(tmp, tmp_dst);  // V += G * (-0.290993)

                            tmp_store = _mm256_loadu_si256(src_b);
                            tmp = _mm256_mulhi_epi16(tmp_store, B_V);
                            tmp_dst = _mm256_adds_epi16(tmp, tmp_dst);  // V += B * (0.439216)

                            tmp_dst = _mm256_adds_epi16(tmp_dst, OFFSET_128);  // V += 128
                            _mm256_storeu_si256(dst, tmp_dst);

                            dst++;
                            src_r++;
                            src_g++;
                            src_b++;
                        }
                    _mm_empty();

                    yuv->shword_to_byte();
                }
            default:
                break;
        }
    delete[] tru;
    delete[] trv;
    delete[] tmp_r;
    delete[] tmp_g;
    delete[] tmp_b;
}

void ALPHA_AMALGAMATE(RGB *src, RGB *dst, int _alpha)
{
    int height = src->h;
    int width = src->w;
    int size = height * width;

    switch (mode)
        {
            case MODE_NORM:
                {
                    for (int i = 0; i < size; ++i)
                        dst->RH[i] = (int)(src->RH[i]) * _alpha / 255;
                    for (int i = 0; i < size; ++i)
                        dst->GH[i] = (int)(src->GH[i]) * _alpha / 255;
                    for (int i = 0; i < size; ++i)
                        dst->BH[i] = (int)(src->BH[i]) * _alpha / 255;
                }
                break;
            case MODE_MMX:
                {
                    _mm_empty();
                    __m64 alpha = _mm_set_pi16(_alpha, _alpha, _alpha, _alpha);
                    __m64 tmp;

                    __m64 *dst_ = (__m64 *)dst->RH;
                    __m64 *src_ = (__m64 *)src->RH;
                    for (int i = 0; i < size / 4; i++)
                        {
                            tmp = _m_pmullw(*src_, alpha);
                            tmp = _m_psrlwi(tmp, 8);
                            *dst_ = tmp;
                            dst_++;
                            src_++;
                        }

                    dst_ = (__m64 *)dst->GH;
                    src_ = (__m64 *)src->GH;
                    for (int i = 0; i < size / 4; i++)
                        {
                            tmp = _m_pmullw(*src_, alpha);
                            tmp = _m_psrlwi(tmp, 8);
                            *dst_ = tmp;
                            dst_++;
                            src_++;
                        }

                    dst_ = (__m64 *)dst->BH;
                    src_ = (__m64 *)src->BH;
                    for (int i = 0; i < size / 4; i++)
                        {
                            tmp = _m_pmullw(*src_, alpha);
                            tmp = _m_psrlwi(tmp, 8);
                            *dst_ = tmp;
                            dst_++;
                            src_++;
                        }

                    _mm_empty();
                }
                break;
            case MODE_SSE2:
                {
                    _mm_empty();
                    __m128i alpha = _mm_set_epi16(_alpha, _alpha, _alpha, _alpha, _alpha, _alpha,
                                                  _alpha, _alpha);
                    __m128i tmp;

                    __m128i *dst_ = (__m128i *)dst->RH;
                    __m128i *src_ = (__m128i *)src->RH;
                    for (int i = 0; i < size / 8; i++)
                        {
                            tmp = _mm_mullo_epi16(*src_, alpha);
                            tmp = _mm_srli_epi16(tmp, 8);
                            *dst_ = tmp;
                            dst_++;
                            src_++;
                        }

                    dst_ = (__m128i *)dst->GH;
                    src_ = (__m128i *)src->GH;
                    for (int i = 0; i < size / 8; i++)
                        {
                            tmp = _mm_mullo_epi16(*src_, alpha);
                            tmp = _mm_srli_epi16(tmp, 8);
                            *dst_ = tmp;
                            dst_++;
                            src_++;
                        }

                    dst_ = (__m128i *)dst->BH;
                    src_ = (__m128i *)src->BH;
                    for (int i = 0; i < size / 8; i++)
                        {
                            tmp = _mm_mullo_epi16(*src_, alpha);
                            tmp = _mm_srli_epi16(tmp, 8);
                            *dst_ = tmp;
                            dst_++;
                            src_++;
                        }

                    _mm_empty();
                }
                break;
            case MODE_AVX:
                {
                    _mm_empty();
                    __m256i alpha = _mm256_set_epi16(_alpha, _alpha, _alpha, _alpha, _alpha, _alpha,
                                                     _alpha, _alpha, _alpha, _alpha, _alpha, _alpha,
                                                     _alpha, _alpha, _alpha, _alpha);
                    __m256i tmp, tmp_store;

                    __m256i *dst_ = (__m256i *)dst->RH;
                    __m256i *src_ = (__m256i *)src->RH;
                    for (int i = 0; i < size / 16; i++)
                        {
                            tmp_store = _mm256_loadu_si256(src_);
                            tmp = _mm256_mullo_epi16(tmp_store, alpha);
                            tmp = _mm256_srli_epi16(tmp, 8);
                            _mm256_storeu_si256(dst_, tmp);
                            dst_++;
                            src_++;
                        }

                    dst_ = (__m256i *)dst->GH;
                    src_ = (__m256i *)src->GH;
                    for (int i = 0; i < size / 16; i++)
                        {
                            tmp_store = _mm256_loadu_si256(src_);
                            tmp = _mm256_mullo_epi16(tmp_store, alpha);
                            tmp = _mm256_srli_epi16(tmp, 8);
                            _mm256_storeu_si256(dst_, tmp);
                            dst_++;
                            src_++;
                        }

                    dst_ = (__m256i *)dst->BH;
                    src_ = (__m256i *)src->BH;
                    for (int i = 0; i < size / 16; i++)
                        {
                            tmp_store = _mm256_loadu_si256(src_);
                            tmp = _mm256_mullo_epi16(tmp_store, alpha);
                            tmp = _mm256_srli_epi16(tmp, 8);
                            _mm256_storeu_si256(dst_, tmp);
                            dst_++;
                            src_++;
                        }

                    _mm_empty();
                }
                break;
            default:
                break;
        }
}

void rgb_add(RGB *src1, RGB *src2, RGB *dst)
{
    int size = src1->h * src1->w;

    for (int i = 0; i < size; ++i)
        dst->RH[i] = CONVERT_CUTOFF((int)(src1->RH[i]) + (int)(src2->RH[i]));
    for (int i = 0; i < size; ++i)
        dst->GH[i] = CONVERT_CUTOFF((int)(src1->GH[i]) + (int)(src2->GH[i]));
    for (int i = 0; i < size; ++i)
        dst->BH[i] = CONVERT_CUTOFF((int)(src1->BH[i]) + (int)(src2->BH[i]));
}

void show_all_alpha(YUV420 *yuv_src)
{
    int height = yuv_src->h;
    int width = yuv_src->w;
    int size = height * width;
    int size_yuv = size * 3 / 2;
    RGB *rgb = new RGB(height, width), *rgb_h = new RGB(height, width);
    YUV2RGB(yuv_src, rgb);
    YUV420 *yuv_h = new YUV420(height, width);
    BYTE *block = new BYTE[size_yuv];

    for (int i = 1; i < 256; i += 3)
        {
            ALPHA_AMALGAMATE(rgb, rgb_h, i);
            RGB2YUV(rgb_h, yuv_h);

            yuv_h->save_block(block);
            char namebuf[32];
            sprintf(namebuf, "result/alpha/%d.yuv", i);
            FILE *fyuv_h = fopen(namebuf, "wb+");
            if (fyuv_h == NULL)
                {
                    printf("Open file faild.\n");
                    return;
                }
            fwrite(block, 1, size_yuv, fyuv_h);
            fclose(fyuv_h);
        }
    delete[] block;
    delete yuv_h;
    delete rgb;
    delete rgb_h;
}

void mix(YUV420 *yuv_src1, YUV420 *yuv_src2)
{
    int height = yuv_src1->h;
    int width = yuv_src1->w;
    int size = height * width;
    int size_yuv = size * 3 / 2;
    RGB *rgb_s1 = new RGB(width, height), *rgb_s2 = new RGB(width, height),
        *rgb_d1 = new RGB(width, height), *rgb_d2 = new RGB(width, height),
        *rgb_dst = new RGB(width, height);
    YUV2RGB(yuv_src1, rgb_s1);
    YUV2RGB(yuv_src2, rgb_s2);
    YUV420 *yuv_dst = new YUV420(yuv_src1->h, yuv_src1->w);
    BYTE *block = new BYTE[size_yuv];
    for (int i = 1; i < 256; i += 3)
        {
            ALPHA_AMALGAMATE(rgb_s1, rgb_d1, (255 - i));
            ALPHA_AMALGAMATE(rgb_s2, rgb_d2, i);
            rgb_add(rgb_d1, rgb_d2, rgb_dst);

            RGB2YUV(rgb_dst, yuv_dst);
            yuv_dst->save_block(block);
            char namebuf[32];
            sprintf(namebuf, "result/mix/%d.yuv", i);
            FILE *fyuv_dst = fopen(namebuf, "wb+");
            if (fyuv_dst == NULL)
                {
                    printf("Open file faild.\n");
                    return;
                }
            fwrite(block, 1, size_yuv, fyuv_dst);
            fclose(fyuv_dst);
        }
}

int main()
{
    mode = MODE_AVX;
    // int size = 1920 * 1080;
    int fsize = 1920 * 1080 * 3 / 2;
    BYTE *block_1 = new BYTE[fsize], *block_2 = new BYTE[fsize], *block_h = new BYTE[fsize];
    FILE *f_1 = fopen("demo/dem1.yuv", "rb");
    if (f_1 == NULL)
        {
            printf("Open file faild.\n");
            return 1;
        }
    fread(block_1, 1, fsize, f_1);
    fclose(f_1);
    FILE *f_2 = fopen("demo/dem2.yuv", "rb");
    if (f_2 == NULL)
        {
            printf("Open file faild.\n");
            return 1;
        }
    fread(block_2, 1, fsize, f_2);
    fclose(f_2);
    YUV420 *yuv_1 = new YUV420(1080, 1920), *yuv_2 = new YUV420(1080, 1920),
           *yuv_h = new YUV420(1080, 1920);
    yuv_1->load_from_block(block_1);
    yuv_2->load_from_block(block_2);


    show_all_alpha(yuv_1);
    mix(yuv_1, yuv_2);

    return 0;
}