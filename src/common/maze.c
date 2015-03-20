/* --------------------------------------------------------------------------
    This file is part of darktable,
    copyright (c) 2014 Stéphane Gimenez <dev@gim.name>

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
* ------------------------------------------------------------------------*/

#include "common/maze.h"

#include "common/darktable.h"
#include "common/linsolve.c"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>

static inline float envelop(float x, float r)
{
  float t = 1 - x*x/(r*r);
  return t*t;
}

void
dt_maze_mosaic_closest(
  const maze_image_t *img,
  const maze_pattern_t *pat,
  const float *x,
  const float *y,
  float *val)
{
  for (int ch = 0; ch < img->ch; ch++)
  {
    const int bx = MAX(floor(x[ch] - pat->x) + 1, img->xmin);
    const int by = MAX(floor(y[ch] - pat->y) + 1, img->ymin);
    const int ex = MIN(ceil(x[ch] + pat->x), img->xmax);
    const int ey = MIN(ceil(y[ch] + pat->y), img->ymax);
    for (int sx = bx; sx < ex; sx++)
      for (int sy = by; sy < ey; sy++)
      {
        float d2 = pat->x * pat->x + pat->y * pat->y + 1;
        int i = pat->x*((pat->offy+sy)%pat->y)+(pat->offx+sx)%pat->x;
        if (pat->data[i] == ch)
        {
          float dx = x[ch] - sx;
          float dy = y[ch] - sy;
          float dd2 = dx * dx + dy * dy;
          if (dd2 < d2)
          {
            d2 = dd2;
            val[ch] = img->data[img->lst * sy + img->sst * sx];
          }
        }
      }
  }
}

// todo: integrate this to interpolate
void
dt_maze_mosaic_downsample(
  const maze_image_t *img,
  const maze_pattern_t *pat,
  const float r,
  const float *x,
  const float *y,
  float *val)
{
  const float ir = r; // todo: compensate smoothing
  const float ir2 = ir * ir;
  for (int ch = 0; ch < img->ch; ch++)
  {
    float v = 0.0;
    float n = 0.0;
    const int bx = MAX(floor(x[ch] - ir) + 1, img->xmin);
    const int by = MAX(floor(y[ch] - ir) + 1, img->ymin);
    const int ex = MIN(ceil(x[ch] + ir), img->xmax);
    const int ey = MIN(ceil(y[ch] + ir), img->ymax);
    for (int sx = bx; sx < ex; sx++)
      for (int sy = by; sy < ey; sy++)
      {
        int i = pat->x*((pat->offy+sy)%pat->y)+(pat->offx+sx)%pat->x;
        if (pat->data[i] == ch)
        {
          float dx = x[ch] - sx;
          float dy = y[ch] - sy;
          float dd = dx*dx + dy*dy;
          if (dd >= ir2) continue;
          float k = envelop(dd, ir2);
          v += k*img->data[img->lst*sy + img->sst*sx];
          n += k;
        }
      }
    val[ch] = v/n;
  }
}

float dt_maze_p(float rsrc, float rdst, float d2)
{
  float r = rsrc + rdst;
  float r2 = r * r;
  return fmax(1.0f - d2 / r2, 0.0f);
}

void
dt_maze_mosaic_deconvolve(
  const maze_image_t *src,
  const maze_pattern_t *pat,
  const maze_image_t *buf,
  const maze_image_t *shift,
  const maze_image_t *dst,
  const float *const rsrc,
  const float *const rdst)
{
  float r[3];
  for(int l = 0; l < src->ch; l++)
    r[l] = rsrc[l] + rdst[l];

  for(int ix = src->xmin; ix < src->xmax; ix++)
    for(int iy = src->ymin; iy < src->ymax; iy++)
      buf->data[iy * buf->lst + ix * buf->sst] = 0.0f;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int jx = dst->xmin; jx < dst->xmax; jx++)
    for(int jy = dst->ymin; jy < dst->ymax; jy++)
    {
      float *tx = shift->data + jy * shift->lst + jx * shift->sst + 0;
      float *ty = shift->data + jy * shift->lst + jx * shift->sst + src->ch;
      for(int px = 0; px < pat->x; px++)
        for(int py = 0; py < pat->y; py++)
        {
          int l = pat->data[py*pat->x + px];
          int xb = MAX(src->xmin, tx[l] - r[l]);
          int xe = MIN(src->xmax, tx[l] + r[l] + 1);
          int yb = MAX(src->ymin, ty[l] - r[l]);
          int ye = MIN(src->ymax, ty[l] + r[l] + 1);
          int xo = (px - xb - pat->offx)%pat->x;
          int yo = (py - yb - pat->offy)%pat->y;
          if (xo < 0) xo += pat->x;
          if (yo < 0) yo += pat->y;
          for(int ix = xb + xo; ix < xe; ix += pat->x)
            for(int iy = yb + yo; iy < ye; iy += pat->y)
            {
              float dx = ix - tx[l];
              float dy = iy - ty[l];
              float d2 = dx * dx + dy * dy;
              buf->data[iy * buf->lst + ix * buf->sst] +=
                dt_maze_p(rsrc[l], rdst[l], d2) *
                dst->data[jy * dst->lst + jx * dst->sst + l];
            }
        }
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int jx = dst->xmin; jx < dst->xmax; jx++)
    for(int jy = dst->ymin; jy < dst->ymax; jy++)
    {
      float *tx = shift->data + jy * shift->lst + jx * shift->sst + 0;
      float *ty = shift->data + jy * shift->lst + jx * shift->sst + src->ch;
      float k[src->ch];
      for(int l = 0; l < src->ch; l++)
        k[l] = 0.0f;
      for(int px = 0; px < pat->x; px++)
        for(int py = 0; py < pat->y; py++)
        {
          int l = pat->data[py*pat->x + px];
          int xb = MAX(src->xmin, tx[l] - r[l]);
          int xe = MIN(src->xmax, tx[l] + r[l] + 1);
          int yb = MAX(src->ymin, ty[l] - r[l]);
          int ye = MIN(src->ymax, ty[l] + r[l] + 1);
          int xo = (px - xb - pat->offx)%pat->x;
          int yo = (py - yb - pat->offy)%pat->y;
          if (xo < 0) xo += pat->x;
          if (yo < 0) yo += pat->y;
          for(int ix = xb + xo; ix < xe; ix += pat->x)
            for(int iy = yb + yo; iy < ye; iy += pat->y)
            {
              float dx = ix - tx[l];
              float dy = iy - ty[l];
              float d2 = dx * dx + dy * dy;
              k[l] += src->data[iy * src->lst + ix * src->sst]
                / buf->data[iy * buf->lst + ix * buf->sst]
                * dt_maze_p(rsrc[l], rdst[l], d2);
            }
        }
      /* k[0]*=4; */
      /* k[1]*=2; */
      /* k[2]*=4; */
      float ka = (k[0]*4 +  k[1]*2 + k[2]*4)/3;
      for(int l = 0; l < src->ch; l++)
      {
        dst->data[jy * dst->lst + jx * dst->sst + l] *= ka;//[l];
        /* printf("%f %f %f\n",k[0],k[1],k[2]); */
      }
    }
}

void
dt_maze_mosaic_interpolate(
  const maze_image_t *img,
  const maze_pattern_t *pat,
  const int deg,
  const int iterations,
  const float r,
  const float *x,
  const float *y,
  float *val)
{
  const float margin = 8.0;

  const int degn = lin_size(2, deg) + img->ch - 1;

  float er[img->ch], mr[img->ch], ir[img->ch];
  for (int ch = 0; ch < img->ch; ch++)
  {
    er[ch] = ch == 1 ? MAX(1.44, r) : MAX(2.0, r);
    mr[ch] = MAX(4.0, deg); // margin radius
    ir[ch] = MAX(mr[ch], r); // input radius
  }

  for (int ch = 0; ch < img->ch; ch++)
    val[ch] = 0.0;

  // check boundaries
  int outside = 0;
  for (int ch = 0; ch < img->ch; ch++)
  {
    if (x[ch] < img->xmin - margin ||
        y[ch] < img->ymin - margin ||
        x[ch] > img->xmax - 1 + margin ||
        y[ch] > img->ymax - 1 + margin)
      return;
    if (x[ch] >= img->xmin + mr[ch] &&
        y[ch] >= img->ymin + mr[ch] &&
        x[ch] <= img->xmax - 1 - mr[ch] &&
        y[ch] <= img->ymax - 1 - mr[ch])
      break;
    float d = MAX(
      MAX(img->xmin + mr[ch] - x[ch], x[ch] - img->xmax + mr[ch] + 1),
      MAX(img->ymin + mr[ch] - y[ch], y[ch] - img->ymax + mr[ch] + 1));
    ir[ch] = MAX(ir[ch], 1.5*d);
    outside = 1;
  }

  // retreive the data
  int vcount = 0;
  int vcountmax = 1024;
  float sk[img->ch];
  float sv[img->ch];
  int vc[vcountmax];
  float vk[vcountmax], ve[vcountmax];
  float vx[vcountmax], vy[vcountmax];
  float vv[vcountmax];
  for (int ch = 0; ch < img->ch; ch++)
  {
    sk[ch] = 0.0;
    sv[ch] = 0.0;
    const float ir2 = ir[ch]*ir[ch];
    const float er2 = er[ch]*er[ch];
    const int bx = MAX(floor(x[ch] - ir[ch]) + 1, img->xmin);
    const int by = MAX(floor(y[ch] - ir[ch]) + 1, img->ymin);
    const int ex = MIN(ceil(x[ch] + ir[ch]), img->xmax);
    const int ey = MIN(ceil(y[ch] + ir[ch]), img->ymax);
    for (int sx = bx; sx < ex; sx++)
      for (int sy = by; sy < ey; sy++)
      {
        int i = pat->x*((pat->offy+sy)%pat->y)+(pat->offx+sx)%pat->x;
        if (pat->data[i] == ch)
        {
          float dx = x[ch] - sx;
          float dy = y[ch] - sy;
          float dd = dx*dx + dy*dy;
          if (dd >= ir2) continue;
          int c = vcount++;
          vc[c] = ch;
          vx[c] = dx;
          vy[c] = dy;
          vv[c] = img->data[img->lst*sy + img->sst*sx];
          vk[c] = envelop(dd, ir2) / (1 + dd/er2);
          ve[c] = 50*dd/er2;
          sk[ch] += vk[c];
          sv[ch] += vk[c]*vv[c];
        }
      }
  }
  assert(vcount < vcountmax);

  int downscale = 0; // debug: should be 1
  for (int ch = 0; ch < img->ch; ch++)
    downscale &= sk[ch] > degn;
  if (outside || downscale)
  {
    for (int ch = 0; ch < img->ch; ch++)
      val[ch] = sv[ch] / sk[ch];
    return;
  }

  float st[degn];
  float vt[degn];
  float mt[degn*degn];

  lin_zero(degn, mt, vt);
  for (int c = 0; c < vcount; c++)
    lin_push2m(degn, img->ch, deg, mt, vt, vc[c], vx[c], vy[c], vk[c], vv[c]);
  int unk = lin_solve(degn, mt, vt, st);
  if (unk != 0)
  {
    printf("mosaic interpolate: unsolved unk=%d count=%d\n", unk, vcount);
    lin_print(degn, mt, vt);
    return;
  }

  // adjust the weights
  for (int i = 0; i < iterations; i++)
  {
    lin_zero(degn, mt, vt);
    for (int c = 0; c < vcount; c++)
    {
      float av = lin_get2m(deg, img->ch, st, vc[c], vx[c], vy[c]);
      float dv = vv[c] - av;
      float dn = 1; // todo: relativize
      float d = dv / dn;
      float e = 1 + ve[c]*d*d;
      lin_push2m(degn, img->ch, deg, mt, vt,
                 vc[c], vx[c], vy[c], vk[c]/e, vv[c]);
    }
    if (lin_solve(degn, mt, vt, st) != 0) break;
  }

  for (int ch = 0; ch < img->ch; ch++)
    val[ch] = lin_get2m_mean(deg, img->ch, st, ch, 0.5*r); // anti-aliasing
}

float
dt_maze_downsample(
  const maze_image_t *img,
  const float r,
  const int ch,
  const float x,
  const float y,
  float val)
{
  const float ir = r; // todo: compensate smoothing
  const float ir2 = ir * ir;
  const int bx = MAX(floor(x - ir) + 1, img->xmin);
  const int by = MAX(floor(y - ir) + 1, img->ymin);
  const int ex = MIN(ceil(x + ir), img->xmax);
  const int ey = MIN(ceil(y + ir), img->ymax);
  float v = 0.0;
  float n = 0.0;
  for (int sx = bx; sx < ex; sx++)
    for (int sy = by; sy < ey; sy++)
    {
      float dx = x - sx;
      float dy = y - sy;
      float dd = dx*dx + dy*dy;
      if (dd >= ir2) continue;
      float k = envelop(dd, ir2);
      v += k*img->data[img->lst*sy + img->sst*sx + ch];
      n += k;
    }
  return v/n;
}

void
dt_maze_interpolate(
  const maze_image_t *img,
  const int deg,
  const int iterations,
  const float margin,
  const float r,
  const float *x,
  const float *y,
  float *val)
{
  const int degn = lin_size(2, deg);

  float er = MAX(1.0, r); // focus radius
  float mr = MAX(2.0, deg);
  float ir = MAX(mr, r); // actual sampling radius

  for (int ch = 0; ch < img->ch; ch++)
    val[ch] = 0.0;

  // check boundaries
  int outside = 0;
  for (int ch = 0; ch < img->ch; ch++)
  {
    if (x[ch] < img->xmin - margin ||
        y[ch] < img->ymin - margin ||
        x[ch] > img->xmax - 1 + margin ||
        y[ch] > img->ymax - 1 + margin)
      return;
    if (x[ch] >= img->xmin + mr &&
        y[ch] >= img->ymin + mr &&
        x[ch] <= img->xmax - 1 - mr &&
        y[ch] <= img->ymax - 1 - mr)
      break;
    float d = MAX(
      MAX(img->xmin + mr - x[ch], x[ch] - img->xmax + mr + 1),
      MAX(img->ymin + mr - y[ch], y[ch] - img->ymax + mr + 1));
    ir = MAX(ir, 1.5*d);
    outside = 1;
  }

  // retreive the data
  int vcount[img->ch];
  int vcountmax = 256;
  float sk[img->ch];
  float sv[img->ch];
  float vk[img->ch][vcountmax], ve[img->ch][vcountmax];
  float vx[img->ch][vcountmax], vy[img->ch][vcountmax];
  float vv[img->ch][vcountmax];
  for (int ch = 0; ch < img->ch; ch++)
  {
    sk[ch] = 0.0;
    sv[ch] = 0.0;
    vcount[ch] = 0;
    float ir2 = ir*ir;
    float er2 = er*er;
    const int bx = MAX(floor(x[ch] - ir) + 1, img->xmin);
    const int by = MAX(floor(y[ch] - ir) + 1, img->ymin);
    const int ex = MIN(ceil(x[ch] + ir), img->xmax);
    const int ey = MIN(ceil(y[ch] + ir), img->ymax);
    for (int sx = bx; sx < ex; sx++)
      for (int sy = by; sy < ey; sy++)
      {
        float dx = x[ch] - sx;
        float dy = y[ch] - sy;
        float dd = dx*dx + dy*dy;
        if (dd >= ir2) continue;
        int c = vcount[ch]++;
        vx[ch][c] = dx;
        vy[ch][c] = dy;
        vv[ch][c] = img->data[img->lst*sy + img->sst*sx + ch];
        vk[ch][c] = envelop(dd, ir2) / (1 + dd*dd/(er2*er2));
        ve[ch][c] = 50*dd/er2;
        sk[ch] += vk[ch][c];
        sv[ch] += vk[ch][c]*vv[ch][c];
      }
    assert(vcount[ch] < vcountmax);
  }

  int downscale = 0; // debug: should be 1
  for (int ch = 0; ch < img->ch; ch++)
    downscale &= sk[ch] > deg*deg;
  if (outside || downscale)
  {
    for (int ch = 0; ch < img->ch; ch++)
      val[ch] = sv[ch] / sk[ch];
    return;
  }

  for (int ch = 0; ch < img->ch; ch++)
  {
    float st[degn];
    float vt[degn];
    float mt[degn*degn];

    lin_zero(degn, mt, vt);
    for (int c = 0; c < vcount[ch]; c++)
      lin_push2(degn, deg, mt, vt, vx[ch][c], vy[ch][c], vk[ch][c], vv[ch][c]);
    int unk = lin_solve(degn, mt, vt, st);
    if (unk != 0)
    {
      printf("interpolate: unsolved unk=%d count=%d\n", unk, vcount[ch]);
      lin_print(degn, mt, vt);
      return;
    }

    // adjust the weights
    for (int i = 0; i < iterations; i++)
    {
      lin_zero(degn, mt, vt);
      for (int c = 0; c < vcount[ch]; c++)
      {
        float av = lin_get2(deg, st, vx[ch][c], vy[ch][c]);
        float dv = vv[ch][c] - av;
        float dn = 1; // todo: relativize
        float d = dv / dn;
        float e = 1 + ve[ch][c]*d*d;
        lin_push2(degn, deg, mt, vt,
                  vx[ch][c], vy[ch][c], vk[ch][c]/e, vv[ch][c]);
      }
      if (lin_solve(degn, mt, vt, st) != 0) break;
    }

    val[ch] = lin_get2_mean(deg, st, 0.5*r); // anti-aliasing
  }
}

