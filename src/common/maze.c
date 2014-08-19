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
  return (1 - x/r)*(1 - x/r);
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


void
dt_maze_mosaic_interpolate(
  const maze_image_t *img,
  const maze_pattern_t *pat,
  const int deg,
  const float r,
  const float *x,
  const float *y,
  float *val)
{
  const float margin = 8.0;
  const int weighting_iterations = 1;

  const int degn = lin_size(2, deg) + img->ch - 1;

  float er = MAX(1.0, r); // focus radius
  float r0[img->ch], ir[img->ch];
  for (int ch = 0; ch < img->ch; ch++)
  {
#if 0
    float adj = pat->quant[ch];
#else
    float adj = 1.0/img->ch;
#endif
    r0[ch] = sqrtf(2*degn*adj/M_PI + img->ch + 2); // min sampling radius
    ir[ch] = MAX(r0[ch], r); // actual sampling radius
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
    if (x[ch] >= img->xmin + r0[ch] &&
        y[ch] >= img->ymin + r0[ch] &&
        x[ch] <= img->xmax - 1 - r0[ch] &&
        y[ch] <= img->ymax - 1 - r0[ch])
      break;
    float d = MAX(
      MAX(img->xmin + r0[ch] - x[ch], x[ch] - img->xmax + r0[ch] + 1),
      MAX(img->ymin + r0[ch] - y[ch], y[ch] - img->ymax + r0[ch] + 1));
    ir[ch] = MAX(ir[ch], 1.5*d);
    outside = 1;
  }
  (void)outside;

  // retreive the data
  int vcount = 0;
  int vcountmax = 1024;
  int vc[vcountmax];
  float vk[vcountmax], ve[vcountmax];
  float vx[vcountmax], vy[vcountmax];
  float vv[vcountmax];
  for (int ch = 0; ch < img->ch; ch++)
  {
    const float ir2 = ir[ch]*ir[ch];
    const float er2 = er*er;
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
          ve[c] = 5000*dd/er2;
        }
      }
  }

  assert(vcount < vcountmax);
  if (vcount <= degn)
  {
    printf("mosaic interpolate: insufficient data count=%d\n", vcount);
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
  for (int i = 0; i < weighting_iterations; i++)
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
  const float margin,
  const float r,
  const float *x,
  const float *y,
  float *val)
{
  const int weighting_iterations = 1;

  const int degn = lin_size(2, deg);

  float er = MAX(1.0, r); // focus radius
  float r0 = sqrtf(1.75*degn/M_PI + 2); // min sampling radius
  float ir = MAX(r0, r); // actual sampling radius

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
    if (x[ch] >= img->xmin + r0 &&
        y[ch] >= img->ymin + r0 &&
        x[ch] <= img->xmax - 1 - r0 &&
        y[ch] <= img->ymax - 1 - r0)
      break;
    float d = MAX(
      MAX(img->xmin + r0 - x[ch], x[ch] - img->xmax + r0 + 1),
      MAX(img->ymin + r0 - y[ch], y[ch] - img->ymax + r0 + 1));
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
    vcount[ch] = 0;
    sk[ch] = 0.0;
    sv[ch] = 0.0;
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
        vk[ch][c] = envelop(dd, ir2) / (1 + dd/er2);
        ve[ch][c] = 5000*dd/er2;
        sk[ch] += vk[ch][c];
        sv[ch] += vk[ch][c]*vv[ch][c];
      }

      assert(vcount[ch] < vcountmax);
  }

  for (int ch = 0; ch < img->ch; ch++)
    if (sk[ch] < 1.0)
      return; // outside data range

  for (int ch = 0; ch < img->ch; ch++)
  {
    if (outside || sk[ch] > 3*degn)
    {
      if (!outside)
        printf("interpolate: downsampling\n");
      val[ch] = sv[ch] / sk[ch];
      continue;
    }

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
    for (int i = 0; i < weighting_iterations; i++)
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

