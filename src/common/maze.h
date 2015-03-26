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

#ifndef MAZE_H
#define MAZE_H

typedef struct {
  float *data;
  int ch;
  int sst;
  int lst;
  int xmin;
  int ymin;
  int xmax;
  int ymax;
} maze_image_t;

typedef struct {
  int x;
  int y;
  int ch;
  const int *data;
  int offx;
  int offy;
} maze_pattern_t;

typedef struct {
  float *data;
  int width;
  int height;
  float xmin;
  float ymin;
  float xmax;
  float ymax;
} maze_trans_t;

void
dt_maze_interpolate(
  const maze_image_t *img,
  const int degree,
  const int iterations,
  const float margin,
  const float r,
  const float *x,
  const float *y,
  float *val);

void
dt_maze_mosaic_closest(
  const maze_image_t *img,
  const maze_pattern_t *pat,
  const float *x,
  const float *y,
  float *val);

void
dt_maze_mosaic_interpolate(
  const maze_image_t *img,
  const maze_pattern_t *pat,
  const int degree,
  const int iterations,
  const float r,
  const float *x,
  const float *y,
  float *val);

void
dt_maze_mosaic_downsample(
  const maze_image_t *img,
  const maze_pattern_t *pat,
  const float r,
  const float *x,
  const float *y,
  float *val);

void
dt_maze_mosaic_deconvolve(
  const maze_image_t *src,
  const maze_pattern_t *pat,
  const maze_image_t *buf,
  const maze_image_t *shift,
  const maze_image_t *dst,
  const maze_trans_t **tr,
  const float chroma);

void
dt_maze_deconvolve(
  const maze_image_t *src,
  const maze_image_t *buf,
  const maze_image_t *shift,
  const maze_image_t *dst,
  const maze_trans_t **tr);

void
dt_maze_trans_build(
  maze_trans_t *t,
  const maze_pattern_t *pat,
  float rs,
  float rd);

#endif
