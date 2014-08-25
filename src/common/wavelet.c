/* Wavelet code */

static inline  __m128
weight_sse(const __m128 *c1, const __m128 *c2, const float param);

/* SSE intrinsics version of dt_fast_expf defined in darktable.h */
static inline __m128
dt_fast_expf_sse(const __m128 x)
{
  const __m128 fone = _mm_set_ps1(0x3f800000u);
  const __m128 femo = _mm_set_ps1(0x00adf880u);
  __m128  f = _mm_add_ps(fone, _mm_mul_ps(x, femo)); // f(n) = i1 + x(n)*(i2-i1)
  __m128i i = _mm_cvtps_epi32(f);                    // i(n) = int(f(n))
  __m128i mask = _mm_srai_epi32(i, 31);              // mask(n) = 0xffffffff if i(n) < 0
  i = _mm_andnot_si128(mask, i);                     // i(n) = 0 if i(n) < 0
  return _mm_castsi128_ps(i);                        // return *(float*)&i
}

// very fast approximation for 2^-x (returns 0 for x > 126)
static inline float
fast_mexp2f(const float x)
{
  const float i1 = (float)0x3f800000u; // 2^0
  const float i2 = (float)0x3f000000u; // 2^-1
  const float k0 = i1 + x * (i2 - i1);
  union floatint_t
  {
    float f;
    uint32_t i;
  } k;
  k.i = k0 >= (float)0x800000u ? k0 : 0;
  return k.f;
}

#define SUM_PIXEL_CONTRIBUTION_COMMON(ii, jj) \
  do { \
    const __m128 f = _mm_set1_ps(filter[(ii)]*filter[(jj)]); \
    const __m128 wp = weight_sse(px, px2, param); \
    const __m128 w = _mm_mul_ps(f, wp); \
    const __m128 pd = _mm_mul_ps(w, *px2); \
    sum = _mm_add_ps(sum, pd); \
    wgt = _mm_add_ps(wgt, w); \
  } while (0)

#define SUM_PIXEL_CONTRIBUTION_WITH_TEST(ii, jj) \
  do { \
    const int iii = (ii)-2; \
    const int jjj = (jj)-2; \
    int x = i + mult*iii; \
    int y = j + mult*jjj; \
    \
    if(x < 0)       x = 0; \
    if(x >= width)  x = width  - 1; \
    if(y < 0)       y = 0; \
    if(y >= height) y = height - 1; \
    \
    px2 = ((__m128 *)in) + x + (size_t)y*width; \
    \
    SUM_PIXEL_CONTRIBUTION_COMMON(ii, jj); \
  } while (0)

#define ROW_PROLOGUE \
  const __m128 *px = ((__m128 *)in) + (size_t)j*width; \
  const __m128 *px2; \
  float *pdetail = detail + (size_t)4*j*width; \
  float *pcoarse = out + (size_t)4*j*width;

#define SUM_PIXEL_PROLOGUE \
  __m128 sum = _mm_setzero_ps(); \
  __m128 wgt = _mm_setzero_ps();

// may replace _mm_div by _mm_mul + _mm_rcp
#define SUM_PIXEL_EPILOGUE \
  sum = _mm_div_ps(sum, wgt); \
  \
  _mm_stream_ps(pdetail, _mm_sub_ps(*px, sum)); \
  _mm_stream_ps(pcoarse, sum); \
  px++; \
  pdetail+=4; \
  pcoarse+=4;

static void
eaw_decompose (float *const out, const float *const in, float *const detail, const int scale,
               const float param, const int32_t width, const int32_t height)
{
  const int mult = 1<<scale;
  static const float filter[5] = {1.0f/16.0f, 4.0f/16.0f, 6.0f/16.0f, 4.0f/16.0f, 1.0f/16.0f};

  /* The first "2*mult" lines use the macro with tests because the 5x5 kernel
   * requires nearest pixel interpolation for at least a pixel in the sum */
#ifdef _OPENMP
  #pragma omp parallel for default(none) schedule(static)
#endif
  for (int j=0; j<2*mult; j++)
  {
    ROW_PROLOGUE

    for(int i=0; i<width; i++)
    {
      SUM_PIXEL_PROLOGUE
      for (int jj=0; jj<5; jj++)
      {
        for (int ii=0; ii<5; ii++)
        {
          SUM_PIXEL_CONTRIBUTION_WITH_TEST(ii, jj);
        }
      }
      SUM_PIXEL_EPILOGUE
    }
  }

#ifdef _OPENMP
  #pragma omp parallel for default(none) schedule(static)
#endif
  for(int j=2*mult; j<height-2*mult; j++)
  {
    ROW_PROLOGUE

    /* The first "2*mult" pixels use the macro with tests because the 5x5 kernel
     * requires nearest pixel interpolation for at least a pixel in the sum */
    for (int i=0; i<2*mult; i++)
    {
      SUM_PIXEL_PROLOGUE
      for (int jj=0; jj<5; jj++)
      {
        for (int ii=0; ii<5; ii++)
        {
          SUM_PIXEL_CONTRIBUTION_WITH_TEST(ii, jj);
        }
      }
      SUM_PIXEL_EPILOGUE
    }

    /* For pixels [2*mult, width-2*mult], we can safely use macro w/o tests
     * to avoid unneeded branching in the inner loops */
    for(int i=2*mult; i<width-2*mult; i++)
    {
      SUM_PIXEL_PROLOGUE
      px2 = ((__m128*)in) + i-2*mult + (size_t)(j-2*mult)*width;
      for (int jj=0; jj<5; jj++)
      {
        for (int ii=0; ii<5; ii++)
        {
          SUM_PIXEL_CONTRIBUTION_COMMON(ii, jj);
          px2 += mult;
        }
        px2 += (width-5)*mult;
      }
      SUM_PIXEL_EPILOGUE
    }

    /* Last two pixels in the row require a slow variant... blablabla */
    for (int i=width-2*mult; i<width; i++)
    {
      SUM_PIXEL_PROLOGUE
      for (int jj=0; jj<5; jj++)
      {
        for (int ii=0; ii<5; ii++)
        {
          SUM_PIXEL_CONTRIBUTION_WITH_TEST(ii, jj);
        }
      }
      SUM_PIXEL_EPILOGUE
    }
  }

  /* The last "2*mult" lines use the macro with tests because the 5x5 kernel
   * requires nearest pixel interpolation for at least a pixel in the sum */
#ifdef _OPENMP
  #pragma omp parallel for default(none) schedule(static)
#endif
  for (int j=height-2*mult; j<height; j++)
  {
    ROW_PROLOGUE

    for(int i=0; i<width; i++)
    {
      SUM_PIXEL_PROLOGUE
      for (int jj=0; jj<5; jj++)
      {
        for (int ii=0; ii<5; ii++)
        {
          SUM_PIXEL_CONTRIBUTION_WITH_TEST(ii, jj);
        }
      }
      SUM_PIXEL_EPILOGUE
    }
  }

  _mm_sfence();
}

#undef SUM_PIXEL_CONTRIBUTION_COMMON
#undef SUM_PIXEL_CONTRIBUTION_WITH_TEST
#undef ROW_PROLOGUE
#undef SUM_PIXEL_PROLOGUE
#undef SUM_PIXEL_EPILOGUE

static void
eaw_synthesize (float *const out, const float *const in, const float *const detail,
                const float *thrsf, const float *boostf, const int32_t width, const int32_t height)
{
  const __m128 threshold = _mm_set_ps(2*thrsf[3], 2*thrsf[2], 2*thrsf[1], 2*thrsf[0]);
  const __m128 boost     = _mm_set_ps(boostf[3], boostf[2], boostf[1], boostf[0]);

  const __m128 threshold2 = _mm_mul_ps(threshold, threshold);

#ifdef _OPENMP
  #pragma omp parallel for default(none) schedule(static)
#endif
  for(int j=0; j<height; j++)
  {
    // TODO: prefetch? _mm_prefetch()
    const __m128 *pin = (__m128 *)in + (size_t)j*width;
    __m128 *pdetail = (__m128 *)detail + (size_t)j*width;
    float *pout = out + (size_t)4*j*width;
    for(int i=0; i<width; i++)
    {
      const __m128 v2 = _mm_mul_ps(*pdetail, *pdetail);
      const __m128 v3 = _mm_mul_ps(v2, *pdetail);
      const __m128 val = _mm_div_ps(v3, _mm_add_ps(v2, threshold2));
      const __m128 amount = _mm_and_ps(_mm_cmpord_ps(val, val), val);
      _mm_stream_ps(pout, _mm_add_ps(*pin, _mm_mul_ps(boost, amount)));
      pdetail ++;
      pin ++;
      pout += 4;
    }
  }
  _mm_sfence();
}
