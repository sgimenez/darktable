/*
    This file is part of darktable,
    copyright (c) 2009--2011 johannes hanika.
    copyright (c) 2014 Stéphane Gimenez

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "bauhaus/bauhaus.h"
#include "common/darktable.h"
#include "common/linsolve.c"
#include "common/maze.h"
#include "develop/imageop.h"
#include "gui/gtk.h"
#include <gtk/gtk.h>
#include <stdlib.h>

#define ETHER_DEBUG 1

#define CA_SHIFT 8 // max allowed CA shift, must be even

const int interpolation_max_width = 8;

// this is the version of the modules parameters,
// and includes version information about compile-time dt
DT_MODULE_INTROSPECTION(1, dt_iop_ether_params_t)

typedef struct dt_iop_ether_params_t
{
  int cacorrect; // rename to ca_type
  int idegree;
  int ipass;
  int ideconv;
} dt_iop_ether_params_t;

typedef struct dt_iop_ether_gui_data_t
{
  GtkWidget *tcombo1;
  GtkWidget *tcombo2;
  GtkWidget *tcombo3;
  GtkWidget *tcombo4;
} dt_iop_ether_gui_data_t;

typedef struct dt_iop_ether_data_t
{
  int cacorrect;
  int idegree;
  int ipass;
  int ideconv;
  int fdegree;
  void *fitdata;
} dt_iop_ether_data_t;

typedef struct dt_iop_ether_global_data_t
{
} dt_iop_ether_global_data_t;

// this returns a translatable name
const char *name()
{
  // make sure you put all your translatable strings into _() !
  return _("etherize");
}

int groups()
{
  return IOP_GROUP_BASIC;
}

int flags()
{
  return IOP_FLAGS_ONE_INSTANCE;
}

static int get_degree()
{
  return dt_conf_get_int("plugins/darkroom/cacorrect/degree");
}

static int get_radial_degree()
{
  return dt_conf_get_int("plugins/darkroom/cacorrect/radial_degree");
}

static int get_visuals()
{
  return dt_conf_get_int("plugins/darkroom/cacorrect/visuals");
}

// not activated but the old dummy value was 50
int legacy_params(dt_iop_module_t *self, const void *const old_params, const int old_version,
                  dt_iop_ether_params_t *new_p, const int new_version)
{
  if(old_version == 1 && new_version == 2)
  {
    new_p->cacorrect = 1;
    new_p->idegree = 2;
    new_p->ipass = 1;
    new_p->ideconv = 2;
    return 0;
  }
  return 1;
}

void modify_roi_in(dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, const dt_iop_roi_t *roi_out,
                   dt_iop_roi_t *roi_in)
{
  dt_iop_ether_data_t *d = piece->data;
  if(1 || !d->fitdata)
  {
    // we want all of it
    roi_in->x = 0;
    roi_in->y = 0;
    roi_in->width = piece->pipe->image.width;
    roi_in->height = piece->pipe->image.height;
  }
  else { // decreases performance since the buffer has to be adapted
    const int margin = CA_SHIFT + interpolation_max_width;
    roi_in->width  = roi_out->width + 2 * margin;
    roi_in->height = roi_out->height + 2 * margin;
    roi_in->x = roi_out->x - margin;
    roi_in->y = roi_out->y - margin;
    if(roi_in->x < 0) { roi_in->width -= roi_in->x; roi_in->x = 0; }
    if(roi_in->y < 0) { roi_in->height -= roi_in->y; roi_in->y = 0; }
    roi_in->width = MIN(piece->pipe->image.width, roi_in->width);
    roi_in->height = MIN(piece->pipe->image.height, roi_in->height);
  }
  roi_in->scale = 1.0f;
}

static void CA_analyse(dt_iop_module_t *const self, dt_dev_pixelpipe_iop_t *const piece,
                       const float *const idata, float *const odata,
                       const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out)
{
  dt_iop_ether_data_t *d = piece->data;

  const int iwidth  = roi_in->width;
  const int iheight = roi_in->height;
  const int isize = MIN(iwidth, iheight);
  const int owidth  = roi_out->width;
  const int oheight = roi_out->height;
  const float oscale = roi_out->scale;
  const int ox = roi_out->x;
  const int oy = roi_out->y;
  const int ix = roi_in->x;
  const int iy = roi_in->y;
  const int ch = piece->colors;

  const uint32_t filters = dt_image_filter(&piece->pipe->image);
  const int ca = (filters & 3) != 1, fb = (filters >> (!ca << 2) & 3) == 2;

  const int sl = CA_SHIFT;
  const int radial = d->cacorrect == 2;
  const int deg = d->fdegree;
  const int degn = lin_size(radial ? 1 : 2, deg);

  const int srad = 22; // size of samples, must be even
  const int subs = srad; // spacing between samples, must be even

  const int TS = (iwidth > 2024 && iheight > 2024) ? 256 : 128;
  const int border = 2 + srad + sl;
  const int ncoeff = radial ? 1 : 2;

  // these small sizes aren't safe to run through the code below.
  if(iheight < 3 * border || iwidth < 3 * border) return;

  int visuals = get_visuals();

#if ETHER_DEBUG
  clock_t time_begin = clock();
  clock_t time_start = time_begin;
  clock_t time_end;
#endif

  // initialize fit coefficients
  float fitmat[2][ncoeff][degn*degn];
  float fitvec[2][ncoeff][degn];
  for(int color = 0; color < 2; color++)
    for(int coeff = 0; coeff < ncoeff; coeff++)
      lin_zero(degn, fitmat[color][coeff], fitvec[color][coeff]);

#if ETHER_DEBUG
  printf("cacorrect: pre-analysis\n");
#endif

  // pre-analysis (todo: allow a second pass to ditch out the outliers)

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int top = 0; top < iheight; top += TS - 2 * border)
    for(int left = 0; left < iwidth; left += TS - 2 * border)
    {
      const int rm = top + TS > iheight ? iheight - top : TS;
      const int cm = left + TS > iwidth ? iwidth - left : TS;

      float rgb[TS][TS] __attribute__((aligned(16)));
      // assume gamma = 2.0 (experimental)
      for(int r = 0; r < rm; r++)
      {
        int c = 0;
        const int inr = top + r;
        const int inc = left + c;
        const float *pi = idata + (size_t)iwidth*inr + inc;
        float *po = &rgb[r][c];
#if __SSE2__
        for(; c + 3 < cm; pi += 4, po += 4, c += 4)
          _mm_store_ps(po, _mm_sqrt_ps(_mm_loadu_ps(pi)));
#endif
        for(; c < cm; pi += 1, po += 1, c += 1)
          *po = sqrt(*pi);
      }

      float la[TS][TS][2];
      // compute derivatives
      for(int r = 2; r < rm - 2; r++)
        for(int c = 2; c < cm - 2; c++)
          {
            la[r][c][0] = (rgb[r + 2][c] - rgb[r - 2][c]) /* /4 */;
            la[r][c][1] = (rgb[r][c + 2] - rgb[r][c - 2]) /* /4 */;
          }

      // thread-local components to be summed up
      float tfitmat[2][ncoeff][degn*degn];
      float tfitvec[2][ncoeff][degn];
      for(int color = 0; color < 2; color++)
        for(int coeff = 0; coeff < ncoeff; coeff++)
          lin_zero(degn, tfitmat[color][coeff], tfitvec[color][coeff]);

      // compute centered offsets for sampling areas
      int rsuboff = (iheight / 4 - border / 2) % subs < subs / 2 ? 0 : subs / 2;
      int csuboff = (iwidth / 4 - border / 2) % subs < subs / 2 ? 0 : subs / 2;
      int roff = (iheight / 4 + rsuboff - (top + border) / 2) % subs;
      int coff = (iwidth / 4 + csuboff - (left + border) / 2) % subs;
      if(roff < 0) roff += subs;
      if(coff < 0) coff += subs;

      for(int r = border + 2 * roff; r < rm - border; r += 2 * subs)
        for(int c = border + 2 * coff; c < cm - border; c += 2 * subs)
        {
          int inr = top + r;
          int inc = left + c;
          if(inr < border || inr > iheight - 1 - border ||
              inc < border || inc > iwidth - 1 - border)
            continue;

          const float y = (2 * (inr + 0.5) - iheight) / isize;
          const float x = (2 * (inc + 0.5) - iwidth) / isize;

          // compute variation map
          float t[2][2*sl+1][2*sl+1];
          for(int color = 0; color < 2; color++)
            for(int i = 0; i < 2 * sl + 1; i++)
              for(int j = 0; j < 2 * sl + 1; j++)
                t[color][i][j] = 0.;

          for(int gi = -srad; gi < srad; gi++)
            for(int gj = -srad + ((gi + ca) & 1); gj < srad; gj += 2)
            {
              const float dgv = la[r + gi][c + gj][0];
              const float dgh = la[r + gi][c + gj][1];
              const int B = (gi + fb) & 1;
              const int R = !B;
              for(int di = -sl + 1; di <= sl; di += 2)
                for(int dj = -sl; dj <= sl; dj += 2)
                  // this test is too costly if the loop is not unrolled
                  /* if(di*di + dj*dj < (int)((sl+0.5)*(sl+0.5))) */
                {
                  float kvR = la[r + gi + di][c + gj + dj][0] - dgv;
                  float khR = la[r + gi + di][c + gj + dj][1] - dgh;
                  float kvB = la[r + gi + dj][c + gj + di][0] - dgv;
                  float khB = la[r + gi + dj][c + gj + di][1] - dgh;
                  t[R][di + sl][dj + sl] += kvR * kvR + khR * khR;
                  t[B][dj + sl][di + sl] += kvB * kvB + khB * khB;
                }
            }

          float bc[2][2] = { };
          float bw[2][2] = { };
          for(int color = 0; color < 2; color++)
          {
            // compute star averages and locate indexed minimum
            int I = 0;
            int J = 0;
            float min = INFINITY;
            for(int di = -sl + 1; di <= sl - 1; di++)
              for(int dj = -sl + 1 + !(di & 1); dj <= sl - 1; dj += 2)
              {
                int d2 = di * di + dj * dj;
                if(d2 < (int)((sl - 0.5) * (sl - 0.5))) // mask
                {
                  float sum =
                    t[color][di + sl + 1][dj + sl] + t[color][di + sl - 1][dj + sl] +
                    t[color][di + sl][dj + sl + 1] + t[color][di + sl][dj + sl - 1];
                  t[color][di + sl][dj + sl] = sum /* /4 */;
                  if(sum < min)
                    if(d2 < (int)((sl - 2.5)*(sl - 2.5))) // mask
                    { min = sum; I = di; J = dj; }
                }
              }

            // find shift estimation and weights
            {
              int i = I + sl;
              int j = J + sl;
              float dv = (t[color][i + 2][j] - t[color][i - 2][j]) /* /4 */;
              float dh = (t[color][i][j + 2] - t[color][i][j - 2]) /* /4 */;
              float d2v = (t[color][i + 2][j] + t[color][i - 2][j]
                           - 2*t[color][i][j])/2 /* /4 */;
              float d2h = (t[color][i][j + 2] + t[color][i][j - 2]
                           - 2*t[color][i][j])/2 /* /4 */;
              float d2m = (t[color][i + 1][j + 1] - t[color][i - 1][j + 1] -
                           t[color][i + 1][j - 1] + t[color][i - 1][j - 1]) /* /4 */;
              float div = 4 * d2v * d2h - d2m * d2m;
              if(div > 0)
              {
                float d2d = d2v - d2h;
                float sqrtD = sqrt(d2d * d2d + d2m * d2m);
                float lb = (d2v + d2h + sqrtD)/2;
                float ls = (d2v + d2h - sqrtD)/2;
                float lb2 = lb * lb;
                float ls2 = ls * ls;
                float pv = -(2 * dv * d2h - dh * d2m) / div;
                float ph = -(2 * dh * d2v - dv * d2m) / div;
                float sv = I + CLAMP(pv, -2.5, +2.5);
                float sh = J + CLAMP(ph, -2.5, +2.5);
                if(radial)
                {
                  float rad = sqrt(x * x + y * y);
                  if(rad != 0.0)
                  {
                    float theta = atan2(d2m, d2d) / 2;
                    float delta = atan2(x, y);
                    float cosd = cos(theta - delta);
                    float sind = sin(theta - delta);
                    bc[color][0] = (y * sv + x * sh) / rad;
                    // weight (deviations are in 1/l^2)
                    bw[color][0] = 1 / (sqrt(cosd * cosd / lb2 + sind * sind / ls2));
                  }
                }
                else
                {
                  /* float sin2th = d2m / sqrtD; */
                  float cos2th = d2d / sqrtD;
                  float sqsinth = (1 - cos2th) / 2;
                  float sqcosth = (1 + cos2th) / 2;
                  bc[color][0] = sv;
                  bc[color][1] = sh;
                  // weights (deviations are in 1/l^2)
                  bw[color][0] = 1/(sqrt(sqcosth / lb2 + sqsinth / ls2));
                  bw[color][1] = 1/(sqrt(sqsinth / lb2 + sqcosth / ls2));
                }
              }
            }
          }

          // show weights
          if(visuals)
          {
            const int rad = 6;
            const int outr = (inr + iy + 0.5) * oscale - 0.5 - oy;
            const int outc = (inc + ix + 0.5) * oscale - 0.5 - ox;
            float *outp = odata + (size_t)owidth*outr + outc;
            if(outr >= rad && outr < oheight - rad &&
                outc >= rad && outc < owidth - rad)
            {
              const float sr =
                sqrtf(bw[0][0] * bw[0][0] + bw[0][1] * bw[0][1]) / 4;
              const float sb =
                sqrtf(bw[1][0] * bw[1][0] + bw[1][1] * bw[1][1]) / 4;
              for(int di = -rad; di <= rad; di++)
                for(int dj = -rad; dj <= rad; dj++)
                {
                  int d2 = di * di + dj * dj;
                  if(d2 < rad * rad)
                  {
                    outp[ch * owidth * di + dj + 0] += fmin(0.5, sr / expf(1.5 * d2));
                    outp[ch * owidth * di + ch * dj + 2] += fmin(0.5, sb / expf(1.5 * d2));
                  }
                }
            }
          }

#if ETHER_DEBUG
          if(y > -0.05 && y < 0.05 && x > -0.05 && x < 0.05)
          {
            printf("x=%+.4lf y=%+.4lf g=%.2lf r  ", x, y, rgb[r][c]);
            printf("%+.5f (%8.4f)  %+.5f (%8.4f)\n",
                   bc[0][0], bw[0][0], bc[0][1], bw[0][1]);
            printf("x=%+.4lf y=%+.4lf g=%.2lf b  ", x, y, rgb[r][c]);
            printf("%+.5f (%8.4f)  %+.5f (%8.4f)\n",
                   bc[1][0], bw[1][0], bc[1][1], bw[1][1]);
          }
#endif

          // fill up parameters for the fitting
          for(int color = 0; color < 2; color++)
            for(int coeff = 0; coeff < ncoeff; coeff++)
              if(radial)
                lin_push1(
                  degn, deg, tfitmat[color][coeff], tfitvec[color][coeff],
                  sqrtf(x * x + y * y), bw[color][coeff], bc[color][coeff]);
              else
                lin_push2(
                  degn, deg, tfitmat[color][coeff], tfitvec[color][coeff],
                  x, y, bw[color][coeff], bc[color][coeff]);
        }

#ifdef _OPENMP
#pragma omp critical
#endif
      {
        for(int color = 0; color < 2; color++)
          for(int coeff = 0; coeff < ncoeff; coeff++)
            lin_add(
              degn,
              tfitmat[color][coeff], tfitvec[color][coeff],
              fitmat[color][coeff], fitvec[color][coeff]);
      }
    }

  // end of pre-analysis

#if ETHER_DEBUG
  time_end = clock();
  printf("time: %f\n", (float)(time_end - time_start) / CLOCKS_PER_SEC);
  time_start = time_end;
  printf("cacorrect: solve\n");
#endif

  // resolve fit coefficients
  float fitsol[2][ncoeff][degn];
  for(int color = 0; color < 2; color++)
    for(int coeff = 0; coeff < ncoeff; coeff++)
    {
      int unk = lin_solve(degn, fitmat[color][coeff], fitvec[color][coeff], fitsol[color][coeff]);
      if(unk)
      {
        printf("cacorrect: unsolved unk=%d color=%d coeff=%d\n", unk, color, coeff);
        return;
      }

      if(radial)
      {
#if ETHER_DEBUG
        printf("cacorrect: radial color=%d offset=%+f\n",
               color, fitsol[color][coeff][0]);
#endif
        fitsol[color][coeff][0] = 0.0;
      }
    }

#if ETHER_DEBUG
  time_end = clock();
  printf("time: %f\n", (float)(time_end - time_start) / CLOCKS_PER_SEC);
  time_start = time_end;
  printf("ether: analyse in %f seconds (cpu time)\n",
         (float)(time_end - time_begin) / CLOCKS_PER_SEC);
#endif

  // store the fit
  d->fitdata = malloc(sizeof(fitsol));
  memcpy(d->fitdata, fitsol, sizeof(fitsol));
}

static void CA_correct(dt_iop_module_t *const self, dt_dev_pixelpipe_iop_t *const piece,
                       const float *const idata, float *const odata,
                       const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out)
{
  dt_iop_ether_data_t *d = piece->data;

  const int iwidth  = roi_in->width;
  const int iheight = roi_in->height;
  const int owidth  = roi_out->width;
  const int oheight = roi_out->height;
  const int isize = MIN(iwidth, iheight);
  const float oscale = roi_out->scale;
  const int ox = roi_out->x;
  const int oy = roi_out->y;
  const int ix = roi_in->x;
  const int iy = roi_in->y;
  const int ch = piece->colors;

#if ETHER_DEBUG
  printf("ether: i.w=%d i.h=%d i.x=%d i.y=%d s=%f\n",
         iwidth, iheight, roi_in->x, roi_in->y, roi_in->scale);
  printf("ether: o.w=%d o.h=%d o.x=%d o.y=%d s=%f\n",
         owidth, oheight, roi_out->x, roi_out->y, roi_out->scale);
#endif

  const uint32_t filters = dt_image_filter(&piece->pipe->image);
  const int ca = (filters & 3) != 1, fb = (filters >> (!ca << 2) & 3) == 2;

  const int sl = CA_SHIFT;
  int generic = d->cacorrect == 1;
  int radial = d->cacorrect == 2;
  const int deg = d->fdegree;
  const int degn = lin_size(radial ? 1 : 2, deg);

  const int TS = (iwidth > 2024 && iheight > 2024) ? 256 : 128;
  const int ncoeff = radial ? 1 : 2;

  int visuals = get_visuals();

  float fitsol[2][ncoeff][degn];
  if(d->cacorrect)
  {
    // do the fitting if not already computed
    if(!d->fitdata)
      CA_analyse(self, piece, idata, odata, roi_in, roi_out);

    // retrieve fit coefficients
    if(d->fitdata)
      memcpy(fitsol, d->fitdata, sizeof(fitsol));
    else
      generic = radial = 0;
  }

#if ETHER_DEBUG
  clock_t time_start = clock();
  clock_t time_end;
#endif

  maze_image_t img;
  img.ch = 3;
  img.data = (float *)idata;
  img.sst = 1;
  img.lst = iwidth;
  img.xmin = 0;
  img.ymin = 0;
  img.xmax = iwidth;
  img.ymax = iheight;
  const int patdata[] =
    {
      0, 1,
      1, 2,
    };
  maze_pattern_t pat;
  pat.x = 2;
  pat.y = 2;
  pat.data = patdata;
  pat.offx = !ca & 1;
  pat.offy = fb & 1;

  maze_image_t buf;
  buf.ch = 3;
  buf.data = malloc(iwidth*iheight*sizeof(float));
  buf.sst = 1;
  buf.lst = iwidth;
  buf.xmin = 0;
  buf.ymin = 0;
  buf.xmax = iwidth;
  buf.ymax = iheight;

  maze_image_t shift;
  shift.ch = 6;
  shift.data = malloc(6*owidth*oheight*sizeof(float));
  shift.sst = 6;
  shift.lst = 6*owidth;
  shift.xmin = 0;
  shift.ymin = 0;
  shift.xmax = owidth;
  shift.ymax = oheight;

  maze_image_t dst;
  dst.ch = ch;
  dst.data = (float *)odata;
  dst.sst = ch;
  dst.lst = ch*owidth;
  dst.xmin = 0;
  dst.ymin = 0;
  dst.xmax = owidth;
  dst.ymax = oheight;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int oby = 0; oby < oheight; oby += TS)
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(int obx = 0; obx < owidth; obx += TS)
    {
      int oex = MIN(owidth, obx + TS);
      int oey = MIN(oheight, oby + TS);

      for(int outr = oby; outr < oey; outr++)
        for(int outc = obx; outc < oex; outc++)
        {
          float bc[3][2] = { { } };
          const float y = (2 * (outr + oy + 0.5) / oscale - iheight) / isize;
          const float x = (2 * (outc + ox + 0.5) / oscale - iwidth) / isize;

          // compute the polynomial
          if(radial)
          {
            const float rad = sqrtf(x * x + y * y);
            if(rad != 0.0)
            {
              float vr = lin_get1(deg, fitsol[0][0], rad);
              float vb = lin_get1(deg, fitsol[1][0], rad);
              bc[0][0] = vr * y / rad;
              bc[0][1] = vr * x / rad;
              bc[2][0] = vb * y / rad;
              bc[2][1] = vb * x / rad;
            }
          }
          else if(generic)
          {
            for(int coeff = 0; coeff < ncoeff; coeff++)
            {
              bc[0][coeff] = lin_get2(deg, fitsol[0][coeff], x, y);
              bc[2][coeff] = lin_get2(deg, fitsol[1][coeff], x, y);
            }
          }

          const float ir = (outr + oy + 0.5) / oscale - iy - 0.5;
          const float ic = (outc + ox + 0.5) / oscale - ix - 0.5;
          float bv[3], bh[3];
          for(int color = 0; color < 3; color++)
          {
            bv[color] = ir + CLAMPS(bc[color][0], -sl, sl);
            bh[color] = ic + CLAMPS(bc[color][1], -sl, sl);
          }
          shift.data[outr * shift.lst + outc * shift.sst + 0] = bh[0];
          shift.data[outr * shift.lst + outc * shift.sst + 1] = bh[1];
          shift.data[outr * shift.lst + outc * shift.sst + 2] = bh[2];
          shift.data[outr * shift.lst + outc * shift.sst + 3] = bv[0];
          shift.data[outr * shift.lst + outc * shift.sst + 4] = bv[1];
          shift.data[outr * shift.lst + outc * shift.sst + 5] = bv[2];

          // shift by interpolation
          float val[3];
          const float r = 2 / M_PI / oscale;
          dt_maze_mosaic_interpolate(&img, &pat, d->idegree, d->ipass, r, bh, bv, val);
          /* dt_maze_mosaic_closest(&img, &pat, bh, bv, val); */
          for (int color = 0; color < 3; color += 2)
            val[color] = fmaxf(0.0, val[color]);

          // show shift norms as isos
          if(visuals)
            for(int color = 0; color < 3; color += 2) // excludes green
            {
              const float norm =
                sqrtf(bv[color] * bv[color] + bh[color] * bh[color]);
              for(float l = 0; l <= sl; l += 0.5)
              {
                const float d = norm - l;
                if((-0.01 < d && d < -0.005) || (0 < d && d < 0.01))
                {
                  val[color] *= d >= 0 ? 1.75 : 1.5;
                  val[color] += d > 0 ? 0.2 : 0.075;
                }
              }
            }

          for(int color = 0; color < 3; color++)
            odata[ch * (size_t)owidth * outr + ch * outc + color] = val[color];
        }
    }

  printf("ether: deconvolve\n");

  float r = 0.90;
  float rsrc[3] = { 2.0f*r, sqrt(2.0f)*r, 2.0f*r };
  float rdst[3] = { 1.0f*r / oscale, 1.0f*r / oscale, 1.0f*r / oscale };
  for(int k = 0; k < d->ideconv; k++)
    dt_maze_mosaic_deconvolve(&img, &pat, &buf, &shift, &dst, rsrc, rdst);

  free(buf.data);
  free(shift.data);

#if ETHER_DEBUG
  time_end = clock();
  printf("ether: correct in %f seconds (cpu time)\n",
         (float)(time_end - time_start) / CLOCKS_PER_SEC);
#endif

}

void process(dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece,
             const float *const in, float *const out,
             const dt_iop_roi_t *roi_in, const dt_iop_roi_t *roi_out)
{
  CA_correct(self, piece, in, out, roi_in, roi_out);
}

void reload_defaults(dt_iop_module_t *module)
{
  // init defaults:
  dt_iop_ether_params_t p;
  p.cacorrect = 1;
  p.idegree = 2;
  p.ipass = 1;
  p.ideconv = 2;
  memcpy(module->params, &p, sizeof(dt_iop_ether_params_t));
  memcpy(module->default_params, &p, sizeof(dt_iop_ether_params_t));

  module->hide_enable_button = 1;
  // (module->dev->image_storage.filters == 9u);
}

void init(dt_iop_module_t *module)
{
  module->params_size = sizeof(dt_iop_ether_params_t);
  module->params = malloc(sizeof(dt_iop_ether_params_t));
  module->default_params = malloc(sizeof(dt_iop_ether_params_t));

  module->default_enabled = 1;

  // we come instead of demosaicing.
  module->priority = 140; // module order created by iop_dependencies.py, do not edit!
}

void cleanup(dt_iop_module_t *module)
{
  free(module->params);
}

void commit_params(dt_iop_module_t *self, dt_iop_ether_params_t *p,
                   dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  dt_iop_ether_data_t *d = piece->data;

  // non raw and preview pipe do not have mosaiced data:
  if(!(pipe->image.flags & DT_IMAGE_RAW)
     || dt_dev_pixelpipe_uses_downsampled_input(pipe))
  {
    piece->enabled = 0;
    return;
  }

  // drop the current fit
  if(d->fitdata)
  {
    free(d->fitdata);
    d->fitdata = NULL;
  }

  // set up parameters
  d->cacorrect = p->cacorrect;
  d->idegree = p->idegree;
  d->ipass = p->ipass;
  d->ideconv = p->ideconv;
  if (p->cacorrect == 1)
  {
    int rdegree = get_radial_degree();
    if (rdegree <= 0) rdegree = 4; // default value
    d->fdegree = rdegree; // order of the radial polynomial fit
  }
  else
  {
    int gdegree = get_degree();
    if (gdegree < 0) gdegree = 4; // default value
    d->fdegree = gdegree; // order of the 2-dimentional polynomial fit
  }
}

void init_pipe(dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  piece->data = malloc(sizeof(dt_iop_ether_data_t));
  dt_iop_ether_data_t *d = piece->data;
  d->fitdata = NULL;
  //printf("ether: init_pipe\n");
}

void cleanup_pipe(dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  dt_iop_ether_data_t *d = piece->data;
  if(d->fitdata)
    free(d->fitdata);
  free(piece->data);
  piece->data = NULL;
  //printf("ether: cleanup_pipe\n");
}

void gui_update(dt_iop_module_t *self)
{
  dt_iop_ether_gui_data_t *g = self->gui_data;
  dt_iop_ether_params_t *p = (void *)self->params;
  dt_bauhaus_combobox_set(g->tcombo1, p->cacorrect);
  dt_bauhaus_combobox_set(g->tcombo2, p->idegree);
  dt_bauhaus_combobox_set(g->tcombo3, p->ipass);
  dt_bauhaus_combobox_set(g->tcombo4, p->ideconv);
}

static void cacorrect_changed(GtkWidget *widget, dt_iop_module_t *self)
{
  dt_iop_ether_gui_data_t *g = self->gui_data;
  dt_iop_ether_params_t *p = (void *)self->params;
  if(self->dt->gui->reset) return;
  p->cacorrect = dt_bauhaus_combobox_get(g->tcombo1);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void idegree_changed(GtkWidget *widget, dt_iop_module_t *self)
{
  dt_iop_ether_gui_data_t *g = self->gui_data;
  dt_iop_ether_params_t *p = (void *)self->params;
  if(self->dt->gui->reset) return;
  p->idegree = dt_bauhaus_combobox_get(g->tcombo2);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void ipass_changed(GtkWidget *widget, dt_iop_module_t *self)
{
  dt_iop_ether_gui_data_t *g = self->gui_data;
  dt_iop_ether_params_t *p = (void *)self->params;
  if(self->dt->gui->reset) return;
  p->ipass = dt_bauhaus_combobox_get(g->tcombo3);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void ideconv_changed(GtkWidget *widget, dt_iop_module_t *self)
{
  dt_iop_ether_gui_data_t *g = self->gui_data;
  dt_iop_ether_params_t *p = (void *)self->params;
  if(self->dt->gui->reset) return;
  p->ideconv = dt_bauhaus_combobox_get(g->tcombo4);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

void gui_init(dt_iop_module_t *self)
{
  self->gui_data = malloc(sizeof(dt_iop_ether_gui_data_t));
  dt_iop_ether_gui_data_t *g = self->gui_data;

  self->widget = gtk_box_new(GTK_ORIENTATION_VERTICAL, DT_BAUHAUS_SPACE);

  g->tcombo1 = dt_bauhaus_combobox_new(self);
  dt_bauhaus_widget_set_label(g->tcombo1, NULL, _("chromatic correction"));
  dt_bauhaus_combobox_clear(g->tcombo1);
  dt_bauhaus_combobox_add(g->tcombo1, _("none"));
  dt_bauhaus_combobox_add(g->tcombo1, _("generic"));
  dt_bauhaus_combobox_add(g->tcombo1, _("radial"));
  g_signal_connect(G_OBJECT(g->tcombo1), "value-changed", G_CALLBACK(cacorrect_changed), (gpointer)self);
  gtk_box_pack_start(GTK_BOX(self->widget), g->tcombo1, TRUE, TRUE, 0);

  g->tcombo2 = dt_bauhaus_combobox_new(self);
  dt_bauhaus_widget_set_label(g->tcombo2, NULL, _("interpolation degree"));
  dt_bauhaus_combobox_clear(g->tcombo2);
  dt_bauhaus_combobox_add(g->tcombo2, _("0"));
  dt_bauhaus_combobox_add(g->tcombo2, _("1"));
  dt_bauhaus_combobox_add(g->tcombo2, _("2"));
  dt_bauhaus_combobox_add(g->tcombo2, _("3"));
  dt_bauhaus_combobox_add(g->tcombo2, _("4"));
  g_signal_connect(G_OBJECT(g->tcombo2), "value-changed", G_CALLBACK(idegree_changed), (gpointer)self);
  gtk_box_pack_start(GTK_BOX(self->widget), g->tcombo2, TRUE, TRUE, 0);

  g->tcombo3 = dt_bauhaus_combobox_new(self);
  dt_bauhaus_widget_set_label(g->tcombo3, NULL, _("weighting iterations"));
  dt_bauhaus_combobox_clear(g->tcombo3);
  dt_bauhaus_combobox_add(g->tcombo3, _("0"));
  dt_bauhaus_combobox_add(g->tcombo3, _("1"));
  dt_bauhaus_combobox_add(g->tcombo3, _("2"));
  dt_bauhaus_combobox_add(g->tcombo3, _("3"));
  g_signal_connect(G_OBJECT(g->tcombo3), "value-changed", G_CALLBACK(ipass_changed), (gpointer)self);
  gtk_box_pack_start(GTK_BOX(self->widget), g->tcombo3, TRUE, TRUE, 0);

  g->tcombo4 = dt_bauhaus_combobox_new(self);
  dt_bauhaus_widget_set_label(g->tcombo4, NULL, _("deconv iterations"));
  dt_bauhaus_combobox_clear(g->tcombo4);
  dt_bauhaus_combobox_add(g->tcombo4, _("0"));
  dt_bauhaus_combobox_add(g->tcombo4, _("1"));
  dt_bauhaus_combobox_add(g->tcombo4, _("2"));
  dt_bauhaus_combobox_add(g->tcombo4, _("3"));
  g_signal_connect(G_OBJECT(g->tcombo4), "value-changed", G_CALLBACK(ideconv_changed), (gpointer)self);
  gtk_box_pack_start(GTK_BOX(self->widget), g->tcombo4, TRUE, TRUE, 0);

}

void gui_cleanup(dt_iop_module_t *self)
{
  free(self->gui_data);
}

// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.sh
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-space on;
