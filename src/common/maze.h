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
  const int *data;
  const float *quant;
  int offx;
  int offy;
} maze_pattern_t;

void
dt_maze_interpolate(
  const maze_image_t *img,
  const int degree,
  const float margin,
  const float r,
  const float *x,
  const float *y,
  float *val);

void
dt_maze_mosaic_interpolate(
  const maze_image_t *img,
  const maze_pattern_t *pat,
  const int degree,
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

#endif
