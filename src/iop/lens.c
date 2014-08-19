/*
    This file is part of darktable,
    copyright (c) 2009--2012 johannes hanika.
    copyright (c) 2014-2015 LebedevRI.
    copyright (c) 2014-2015 Stéphane Gimenez

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
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <gtk/gtk.h>
#include <inttypes.h>
#include <ctype.h>
#include <lensfun/lensfun.h>
#include "develop/develop.h"
#include "develop/imageop.h"
#include "develop/tiling.h"
#include "common/opencl.h"
#include "common/interpolation.h"
#include "common/maze.h"
#include "control/control.h"
#include "dtgtk/button.h"
#include "dtgtk/resetlabel.h"
#include "bauhaus/bauhaus.h"
#include "gui/accelerators.h"
#include "gui/gtk.h"
#include "gui/draw.h"

#define LENS_DEBUG 1

#if LF_VERSION < ((0 << 24) | (2 << 16) | (9 << 8) | 0)
#define LF_SEARCH_SORT_AND_UNIQUIFY 2
#endif

#ifndef __GNUC_PREREQ
// on OSX, gcc-4.6 and clang chokes if this is not here.
#if defined __GNUC__ && defined __GNUC_MINOR__
#define __GNUC_PREREQ(maj, min) ((__GNUC__ << 16) + __GNUC_MINOR__ >= ((maj) << 16) + (min))
#else
#define __GNUC_PREREQ(maj, min) 0
#endif
#endif

DT_MODULE_INTROSPECTION(5, dt_iop_lensfun_params_t)

#define LENS_MOZ 0

typedef enum dt_iop_lensfun_modflag_t
{
  LENSFUN_MODFLAG_NONE = 0,
  LENSFUN_MODFLAG_ALL = LF_MODIFY_DISTORTION | LF_MODIFY_TCA | LF_MODIFY_VIGNETTING,
  LENSFUN_MODFLAG_DIST_TCA = LF_MODIFY_DISTORTION | LF_MODIFY_TCA,
  LENSFUN_MODFLAG_DIST_VIGN = LF_MODIFY_DISTORTION | LF_MODIFY_VIGNETTING,
  LENSFUN_MODFLAG_TCA_VIGN = LF_MODIFY_TCA | LF_MODIFY_VIGNETTING,
  LENSFUN_MODFLAG_DIST = LF_MODIFY_DISTORTION,
  LENSFUN_MODFLAG_TCA = LF_MODIFY_TCA,
  LENSFUN_MODFLAG_VIGN = LF_MODIFY_VIGNETTING,
  LENSFUN_MODFLAG_MASK = LF_MODIFY_DISTORTION | LF_MODIFY_TCA | LF_MODIFY_VIGNETTING
} dt_iop_lensfun_modflag_t;

typedef struct dt_iop_lensfun_modifier_t
{
  char name[40];
  int pos; // position in combo box
  int modflag;
} dt_iop_lensfun_modifier_t;

typedef struct dt_iop_lensfun_params_t
{
  int modify_flags;
  int inverse;
  float scale;
  float crop;
  float focal;
  float aperture;
  float distance;
  lfLensType target_geom;
  char camera[128];
  char lens[128];
  int tca_override;
  float tca_r, tca_b;
  int modified; // did user changed anything from automatically detected?
} dt_iop_lensfun_params_t;

typedef struct dt_iop_lensfun_gui_data_t
{
  const lfCamera *camera;
  GtkWidget *lens_param_box;
  GtkWidget *detection_warning;
  GtkWidget *cbe[3];
  GtkButton *camera_model;
  GtkMenu *camera_menu;
  GtkButton *lens_model;
  GtkMenu *lens_menu;
  GtkWidget *modflags, *target_geom, *reverse, *tca_r, *tca_b, *scale;
  GtkWidget *find_lens_button;
  GtkWidget *find_camera_button;
  GList *modifiers;
  GtkLabel *message;
  int corrections_done;
  dt_pthread_mutex_t lock;
} dt_iop_lensfun_gui_data_t;

typedef struct dt_iop_lensfun_global_data_t
{
  lfDatabase *db;
  int kernel_lens_distort_bilinear;
  int kernel_lens_distort_bicubic;
  int kernel_lens_distort_lanczos2;
  int kernel_lens_distort_lanczos3;
  int kernel_lens_vignette;
} dt_iop_lensfun_global_data_t;

typedef struct dt_iop_lensfun_data_t
{
  lfLens *lens;
  int modify_flags;
  int inverse;
  float scale;
  float crop;
  float focal;
  float aperture;
  float distance;
  lfLensType target_geom;
  gboolean do_nan_checks;
} dt_iop_lensfun_data_t;

const char *name()
{
  return _("lens correction");
}

int groups()
{
  return IOP_GROUP_CORRECT;
}

int operation_tags()
{
  return IOP_TAG_DISTORT;
}

int flags()
{
  return IOP_FLAGS_ALLOW_TILING | IOP_FLAGS_TILING_FULL_ROI |
    IOP_FLAGS_ONE_INSTANCE;
}

int
output_bpp(dt_iop_module_t *module, dt_dev_pixelpipe_t *pipe,
           dt_dev_pixelpipe_iop_t *piece)
{
  const int moz = LENS_MOZ && (piece->pipe->image.flags & DT_IMAGE_RAW) &&
    !dt_dev_pixelpipe_uses_downsampled_input(piece->pipe);
  return (moz ? 1 : piece->colors)*sizeof(float);
}

void init_key_accels(dt_iop_module_so_t *self)
{
  dt_accel_register_slider_iop(self, FALSE, NC_("accel", "scale"));
  dt_accel_register_slider_iop(self, FALSE, NC_("accel", "TCA R"));
  dt_accel_register_slider_iop(self, FALSE, NC_("accel", "TCA B"));

  dt_accel_register_iop(self, FALSE, NC_("accel", "find camera"), 0, 0);
  dt_accel_register_iop(self, FALSE, NC_("accel", "find lens"), 0, 0);
  dt_accel_register_iop(self, FALSE, NC_("accel", "auto scale"), 0, 0);
  dt_accel_register_iop(self, FALSE, NC_("accel", "camera model"), 0, 0);
  dt_accel_register_iop(self, FALSE, NC_("accel", "lens model"), 0, 0);
  dt_accel_register_iop(self, FALSE, NC_("accel", "select corrections"), 0, 0);
}

void connect_key_accels(dt_iop_module_t *self)
{
  dt_iop_lensfun_gui_data_t *g = (dt_iop_lensfun_gui_data_t *)self->gui_data;

  dt_accel_connect_button_iop(self, "find lens", GTK_WIDGET(g->find_lens_button));
  dt_accel_connect_button_iop(self, "lens model", GTK_WIDGET(g->lens_model));
  dt_accel_connect_button_iop(self, "camera model", GTK_WIDGET(g->camera_model));
  dt_accel_connect_button_iop(self, "find camera", GTK_WIDGET(g->find_camera_button));
  dt_accel_connect_button_iop(self, "select corrections", GTK_WIDGET(g->modflags));

  dt_accel_connect_slider_iop(self, "scale", GTK_WIDGET(g->scale));
  dt_accel_connect_slider_iop(self, "tca R", GTK_WIDGET(g->tca_r));
  dt_accel_connect_slider_iop(self, "tca B", GTK_WIDGET(g->tca_b));
}

int legacy_params(dt_iop_module_t *self, const void *const old_params, const int old_version,
                  void *new_params, const int new_version)
{
  if(old_version == 2 && new_version == 5)
  {
    // legacy params of version 2; version 1 comes from ancient times and seems to be forgotten by now
    typedef struct dt_iop_lensfun_params_v2_t
    {
      int modify_flags;
      int inverse;
      float scale;
      float crop;
      float focal;
      float aperture;
      float distance;
      lfLensType target_geom;
      char camera[52];
      char lens[52];
      int tca_override;
      float tca_r, tca_b;
    } dt_iop_lensfun_params_v2_t;

    const dt_iop_lensfun_params_v2_t *o = old_params;
    dt_iop_lensfun_params_t *n = new_params;
    dt_iop_lensfun_params_t *d = self->default_params;

    *n = *d; // start with a fresh copy of default parameters

    n->modify_flags = o->modify_flags;
    n->inverse = o->inverse;
    n->scale = o->scale;
    n->crop = o->crop;
    n->focal = o->focal;
    n->aperture = o->aperture;
    n->distance = o->distance;
    n->target_geom = o->target_geom;
    n->tca_override = o->tca_override;
    strncpy(n->camera, o->camera, sizeof(n->camera));
    strncpy(n->lens, o->lens, sizeof(n->lens));
    n->modified = 1;

    // old versions had R and B swapped
    n->tca_r = o->tca_b;
    n->tca_b = o->tca_r;

    return 0;
  }
  if(old_version == 3 && new_version == 5)
  {
    typedef struct dt_iop_lensfun_params_v3_t
    {
      int modify_flags;
      int inverse;
      float scale;
      float crop;
      float focal;
      float aperture;
      float distance;
      lfLensType target_geom;
      char camera[128];
      char lens[128];
      int tca_override;
      float tca_r, tca_b;
    } dt_iop_lensfun_params_v3_t;

    const dt_iop_lensfun_params_v3_t *o = old_params;
    dt_iop_lensfun_params_t *n = new_params;
    dt_iop_lensfun_params_t *d = self->default_params;

    *n = *d; // start with a fresh copy of default parameters

    memcpy(n, o, sizeof(dt_iop_lensfun_params_t) - sizeof(int));
    n->modified = 0;

    // old versions had R and B swapped
    n->tca_r = o->tca_b;
    n->tca_b = o->tca_r;

    return 0;
  }

  if(old_version == 4 && new_version == 5)
  {
    typedef struct dt_iop_lensfun_params_v4_t
    {
      int modify_flags;
      int inverse;
      float scale;
      float crop;
      float focal;
      float aperture;
      float distance;
      lfLensType target_geom;
      char camera[128];
      char lens[128];
      int tca_override;
      float tca_r, tca_b;
      int modified;
    } dt_iop_lensfun_params_v4_t;

    const dt_iop_lensfun_params_v4_t *o = old_params;
    dt_iop_lensfun_params_t *n = new_params;
    dt_iop_lensfun_params_t *d = self->default_params;

    *n = *d; // start with a fresh copy of default parameters

    memcpy(n, o, sizeof(dt_iop_lensfun_params_t));

    // old versions had R and B swapped
    n->tca_r = o->tca_b;
    n->tca_b = o->tca_r;

    return 0;
  }

  return 1;
}

static void get_sanitized_lens(char *dest, int length, const char *orig)
{
  char *found_or = strstr(orig, " or ");
  char *found_parenthesis = strstr(orig, " (");
  int l = length;
  if(found_or || found_parenthesis)
  {
    size_t pos_or = (size_t)(found_or - orig);
    size_t pos_parenthesis = (size_t)(found_parenthesis - orig);
    size_t pos = pos_or < pos_parenthesis ? pos_or : pos_parenthesis;
    if(pos > 0 && pos < length)
      l = pos;
  }
  memset(dest, 0, length);
  strncpy(dest, orig, l-1); // ensures that it's null-terminated
}

static int mod_dist(int modflags)
{
  return modflags & (LF_MODIFY_TCA | LF_MODIFY_DISTORTION | LF_MODIFY_GEOMETRY | LF_MODIFY_SCALE);
}

static int mod_vign(int modflags)
{
  return modflags & (LF_MODIFY_VIGNETTING);
}

void process(dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece,
             float *idata, float *odata,
             const dt_iop_roi_t *roi_in, const dt_iop_roi_t *roi_out)
{
  const dt_iop_lensfun_data_t *const d = piece->data;
  dt_iop_lensfun_gui_data_t *g = self->gui_data;

  const int moz = LENS_MOZ && (piece->pipe->image.flags & DT_IMAGE_RAW) &&
    !dt_dev_pixelpipe_uses_downsampled_input(piece->pipe);

  const unsigned int ch = moz ? 1 : piece->colors;
  const unsigned int pixelformat =
    moz ? LF_CR_1(INTENSITY) :
    ch == 3 ? LF_CR_3(RED, GREEN, BLUE) : LF_CR_4(RED, GREEN, BLUE, UNKNOWN);

  const int iwidth  = roi_in->width;
  const int iheight = roi_in->height;
  const float iscale = roi_in->scale;
  const int owidth  = roi_out->width;
  const int oheight = roi_out->height;
  const float oscale = roi_out->scale;

  const float orig_w = iscale * piece->iwidth;
  const float orig_h = iscale * piece->iheight;

#if LENS_DEBUG
  printf("lens: i.w=%d i.h=%d i.x=%d i.y=%d s=%f\n",
         iwidth, iheight, roi_in->x, roi_in->y, roi_in->scale);
  printf("lens: o.w=%d o.h=%d o.x=%d o.y=%d s=%f\n",
         owidth, oheight, roi_out->x, roi_out->y, roi_out->scale);
  printf("lens: orig_w=%f orig_h=%f moz=%d\n", orig_w, orig_h, moz);
#endif

  dt_pthread_mutex_lock(&darktable.plugin_threadsafe);
  lfModifier *modifier = lf_modifier_new(d->lens, d->crop, orig_w, orig_h);
  int modflags =
    lf_modifier_initialize(
      modifier, d->lens, LF_PF_F32,
      d->focal, d->aperture, d->distance, d->scale,
      d->target_geom, d->modify_flags, d->inverse);
  dt_pthread_mutex_unlock(&darktable.plugin_threadsafe);

  const int mod_d = mod_dist(modflags);
  const int mod_v = mod_vign(modflags);

  const uint32_t filters = dt_image_filter(&piece->pipe->image);
  const int ca = (filters & 3) != 1, fb = (filters >> (!ca << 2) & 3) == 2;

  // acquire temp memory for image buffer
  float *tmpbuf = dt_alloc_align(16, ch * iwidth * iheight * sizeof(float));
  memcpy(tmpbuf, idata, ch * iwidth * iheight * sizeof(float));

#if LENS_DEBUG
  printf("lens: start\n");
#endif

  if(mod_v && !d->inverse)
  {
#if LENS_DEBUG
    printf("lens: vignette\n");
#endif
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(int y = 0; y < iheight; y++)
    {
      lf_modifier_apply_color_modification(
        modifier, tmpbuf + ch*iwidth*y,
        roi_in->x, roi_in->y + y, iwidth, 1, pixelformat, iwidth);
    }
  }

  if(mod_d)
  {
    maze_image_t img;
    img.ch = moz ? 3 : ch;
    img.data = tmpbuf;
    img.sst = ch;
    img.lst = ch*iwidth;
    img.xmin = 0;
    img.ymin = 0;
    img.xmax = iwidth;
    img.ymax = iheight;
    maze_pattern_t pat;
    if (moz)
    {
      const int patdata[] =
        {
          0, 1,
          1, 2,
        };
      const float patquant[] = { 2.0/5.0, 1.0/5.0, 2.0/5.0 };
      pat.x = 2;
      pat.y = 2;
      pat.data = patdata;
      pat.quant = patquant;
      pat.offx = !ca & 1;
      pat.offy = fb & 1;
    }

    // acquire temp memory for image buffer
    const size_t dbuf_req = owidth * 2 * 3 * sizeof(float);
    void *dbuf = dt_alloc_align(16, dbuf_req * dt_get_num_threads());

#if LENS_DEBUG
    printf("lens: distortion\n");
#endif
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(int y = 0; y < oheight; y++)
    {
      float *pi = (float *)(dbuf + dbuf_req * dt_get_thread_num());
      lf_modifier_apply_subpixel_geometry_distortion(
        modifier, roi_out->x, roi_out->y + y, owidth, 1, pi);

      if(moz)
        for(int x = 0; x < owidth; x++, pi += 6)
        {
          const int color = (x + y + ca) & 1 ? 1 : (x + y + fb) & 1 ? 2 : 0;
          float px[3], py[3];
          for (int i = 0; i < 3; i++)
          {
            px[i] = pi[2 * i + 0] - roi_in->x;
            py[i] = pi[2 * i + 1] - roi_in->y;
          }
          const float r = 2 / M_PI * iscale / oscale;
          const float rmax = 3.0;
          float val[3];
          if (r > rmax)
            dt_maze_mosaic_downsample(&img, &pat, r, px, py, val);
          else
            dt_maze_mosaic_interpolate(&img, &pat, 2, r, px, py, val);
          odata[ch * owidth * y + ch * x] = fmaxf(0.0, val[color]);
        }
      else // demosaicaized
        for (int x = 0; x < owidth; x++, pi += 6)
        {
          float px[3], py[3];
          for (int i = 0; i < ch; i++)
          {
            const int index = i > 2 ? 1 : i;
            px[i] = pi[2 * index + 0] - roi_in->x;
            py[i] = pi[2 * index + 1] - roi_in->y;
          }
          float val[3];
          const int degree = 2;
          const float margin = 10*iscale;
          const float r = 2/M_PI*iscale/oscale;
          dt_maze_interpolate(&img, degree, margin, r, px, py, val);
          for (int i = 0; i < ch; i++)
            odata[ch*owidth*y+ch*x+i] = fmaxf(0.0, val[i]);
        }
    }

    dt_free_align(dbuf);
  }
  else
    memcpy(odata, tmpbuf, ch * iwidth * iheight * sizeof(float));

  dt_free_align(tmpbuf);

  if(mod_v && d->inverse)
  {
#if LENS_DEBUG
    printf("lens: vignette\n");
#endif
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(int y = 0; y < oheight; y++)
    {
      lf_modifier_apply_color_modification(
        modifier, odata + ch * owidth * y, roi_out->x, roi_out->y + y,
        owidth, 1, pixelformat, owidth);  // fixme
    }
  }

#if LENS_DEBUG
  printf("lens: end\n");
#endif

  lf_modifier_destroy(modifier);

  if(self->dev->gui_attached && g)
  {
    dt_pthread_mutex_lock(&g->lock);
    g->corrections_done = (modflags & LENSFUN_MODFLAG_MASK);
    dt_pthread_mutex_unlock(&g->lock);
  }

}

void tiling_callback(dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece,
                     const dt_iop_roi_t *roi_in, const dt_iop_roi_t *roi_out,
                     dt_develop_tiling_t *tiling)
{
  const int moz = LENS_MOZ && (piece->pipe->image.flags & DT_IMAGE_RAW) &&
    !dt_dev_pixelpipe_uses_downsampled_input(piece->pipe);

  tiling->factor = 3.5f;  // in + out + tmpbuf
  tiling->maxbuf = 1.5f;
  tiling->overhead = 0;
  tiling->overlap = 0;
  tiling->xalign = moz ? 2 : 1;
  tiling->yalign = moz ? 2 : 1;
  return;
}

int
distort_transform(dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece,
                  float *points, size_t points_count)
{
  dt_iop_lensfun_data_t *d = piece->data;

  if(!d->lens->Maker || d->crop <= 0.0f) return 0;

  const float orig_w = piece->buf_in.width, orig_h = piece->buf_in.height;
  lfModifier *modifier = lf_modifier_new(d->lens, d->crop, orig_w, orig_h);

  int modflags = lf_modifier_initialize(modifier, d->lens, LF_PF_F32, d->focal, d->aperture, d->distance,
                                        d->scale, d->target_geom, d->modify_flags, !d->inverse);
  float *buf = malloc(2 * 3 * sizeof(float));

  for(size_t i = 0; i < points_count * 2; i += 2)
  {
    if(mod_dist(modflags))
    {
      lf_modifier_apply_subpixel_geometry_distortion(modifier, points[i], points[i + 1], 1, 1, buf);
      points[i] = buf[0];
      points[i + 1] = buf[3];
    }
  }
  free(buf);
  lf_modifier_destroy(modifier);

  return 1;
}

int distort_backtransform(dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece,
                          float *points, size_t points_count)
{
  dt_iop_lensfun_data_t *d = (dt_iop_lensfun_data_t *)piece->data;
  if(!d->lens->Maker || d->crop <= 0.0f) return 0;

  const float orig_w = piece->buf_in.width, orig_h = piece->buf_in.height;
  lfModifier *modifier = lf_modifier_new(d->lens, d->crop, orig_w, orig_h);

  int modflags = lf_modifier_initialize(modifier, d->lens, LF_PF_F32, d->focal, d->aperture, d->distance,
                                        d->scale, d->target_geom, d->modify_flags, d->inverse);
  float *buf = malloc(2 * 3 * sizeof(float));

  for(size_t i = 0; i < points_count * 2; i += 2)
  {
    if(mod_dist(modflags))
    {
      lf_modifier_apply_subpixel_geometry_distortion(modifier, points[i], points[i + 1], 1, 1, buf);
      points[i] = buf[0];
      points[i + 1] = buf[3];
    }
  }
  free(buf);
  lf_modifier_destroy(modifier);
  return 1;
}

void modify_roi_out(dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece,
                    dt_iop_roi_t *roi_out, const dt_iop_roi_t *roi_in)
{
  *roi_out = *roi_in;
}

void modify_roi_in(dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece,
                   const dt_iop_roi_t *roi_out, dt_iop_roi_t *roi_in)
{
  *roi_in = *roi_out;

  // inverse transform with given params

  dt_iop_lensfun_data_t *d = piece->data;

  const float orig_w = roi_in->scale*piece->iwidth;
  const float orig_h = roi_in->scale*piece->iheight;

  const int moz = LENS_MOZ && (piece->pipe->image.flags & DT_IMAGE_RAW) &&
    !dt_dev_pixelpipe_uses_downsampled_input(piece->pipe);

  const int imargin = moz ? 12 : 8;

  lfModifier *modifier = lf_modifier_new(d->lens, d->crop, orig_w, orig_h);

  float xm = FLT_MAX, xM = -FLT_MAX, ym = FLT_MAX, yM = -FLT_MAX;

  int modflags = lf_modifier_initialize(modifier, d->lens, LF_PF_F32, d->focal, d->aperture, d->distance,
                                        d->scale, d->target_geom, d->modify_flags, d->inverse);
  // todo: extend input area

  if(mod_dist(modflags))
  {
    // acquire temp memory for distorted pixel coords
    const size_t dbuf_len = sizeof(float)*roi_in->width*2*3;
    float *dbuf = dt_alloc_align(16, dbuf_len);
    for(int y = 0; y < roi_out->height; y++)
    {
      lf_modifier_apply_subpixel_geometry_distortion(
        modifier, roi_out->x, roi_out->y+y, roi_out->width, 1, dbuf);
      const float *pi = dbuf;
      // reverse transform the global coords from lf to our buffer
      for(int x = 0; x < roi_out->width; x++)
      {
        for(int c = 0; c < 3; c++, pi += 2)
        {
          xm = MIN(xm, pi[0]);
          xM = MAX(xM, pi[0]);
          ym = MIN(ym, pi[1]);
          yM = MAX(yM, pi[1]);
        }
      }
    }
    dt_free_align(dbuf);

    roi_in->x = floor(fmaxf(0.0f, xm - imargin));
    roi_in->y = floor(fmaxf(0.0f, ym - imargin));
    roi_in->width = ceil(fminf(orig_w, xM + imargin)) - roi_in->x;
    roi_in->height = ceil(fminf(orig_h, yM + imargin)) - roi_in->y;
  }
  lf_modifier_destroy(modifier);
}

void auto_fill(dt_iop_module_t *self, dt_iop_lensfun_params_t *p,
               dt_iop_lensfun_data_t *d)
{
  dt_iop_lensfun_global_data_t *gd = self->data;
  dt_iop_lensfun_gui_data_t *g = self->gui_data;
  lfDatabase *dt_iop_lensfun_db = gd->db;
  const dt_image_t *img = &self->dev->image_storage;

  char s_camera[52];
  if(p->camera[0])
    g_strlcpy(s_camera, p->camera, sizeof(s_camera));
  else
    g_strlcpy(s_camera, img->exif_model, sizeof(s_camera));

  char s_lens[52];
  if(p->lens[0])
    g_strlcpy(s_lens, p->lens, sizeof(s_lens));
  else
    get_sanitized_lens(s_lens, sizeof(s_lens), img->exif_lens);

  d->crop = p->crop == -1 ? img->exif_crop : p->crop;
  d->aperture = p->aperture == -1 ? img->exif_aperture : p->aperture;
  d->focal = p->focal == -1 ? img->exif_focal_length : p->focal;
  d->distance = p->distance == -1 ?
    //if we did not find focus_distance in EXIF, lets default to 4.0
    (img->exif_focus_distance == 0.0f ? 4.0f : img->exif_focus_distance) :
    p->distance;
  d->modify_flags = p->modify_flags;
  d->inverse      = p->inverse;
  d->target_geom  = p->target_geom;

  // get camera
  const lfCamera *camera = NULL;
  if(s_camera[0])
  {
    const lfCamera **cam = NULL;
    dt_pthread_mutex_lock(&darktable.plugin_threadsafe);
    cam = lf_db_find_cameras_ext(dt_iop_lensfun_db, NULL, s_camera, 0);
    if(cam)
    {
      camera = cam[0];
      if(p->crop == -1)
        d->crop = cam[0]->CropFactor;
    }
    lf_free(cam);
    dt_pthread_mutex_unlock(&darktable.plugin_threadsafe);
  }

  if(g)
    g->camera = camera;

  // build lens
  if(d->lens)
    lf_lens_destroy(d->lens);
  d->lens = lf_lens_new();
  if(s_lens[0])
  {
    dt_pthread_mutex_lock(&darktable.plugin_threadsafe);
    const lfLens **lens = lf_db_find_lenses_hd(
      dt_iop_lensfun_db, camera, NULL, s_lens, LF_SEARCH_SORT_AND_UNIQUIFY);
    dt_pthread_mutex_unlock(&darktable.plugin_threadsafe);
    if(lens)
    {
      lf_lens_copy(d->lens, lens[0]);
      if(p->tca_override)
      {
        // add manual d->lens stuff:
        lfLensCalibTCA tca = { 0 };
        tca.Focal = 0;
        tca.Model = LF_TCA_MODEL_LINEAR;
        tca.Terms[0] = p->tca_r;
        tca.Terms[1] = p->tca_b;
        if(d->lens->CalibTCA)
          while(d->lens->CalibTCA[0]) lf_lens_remove_calib_tca(d->lens, 0);
        lf_lens_add_calib_tca(d->lens, &tca);
      }
      lf_free(lens);
    }
  }

  if(p->scale != -1)
    d->scale = p->scale;
  else
  {
    float scale = 1.0;
    if(d->lens)
    {
      dt_pthread_mutex_lock(&darktable.plugin_threadsafe);
      // create dummy modifier

      // FIXME: get those from rawprepare IOP somehow !!!
      const int iwd = img->width - img->crop_x - img->crop_width,
                iht = img->height - img->crop_y - img->crop_height;

      lfModifier *modifier = lf_modifier_new(d->lens, d->crop, iwd, iht);
      (void)lf_modifier_initialize(
        modifier, d->lens, LF_PF_F32,
        d->focal, d->aperture,
        d->distance, 1.0f,
        d->target_geom, d->modify_flags, d->inverse);
      scale = lf_modifier_get_auto_scale(modifier, d->inverse);
      lf_modifier_destroy(modifier);
      dt_pthread_mutex_unlock(&darktable.plugin_threadsafe);
    }
    d->scale = scale;
  }
}

void commit_params(dt_iop_module_t *self, dt_iop_lensfun_params_t *p,
                   dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  dt_iop_lensfun_data_t *d = piece->data;
#ifdef HAVE_GEGL
  // pull in new params to gegl
#error "lensfun needs to be ported to GEGL!"
#else
  auto_fill(self, p, d);
#endif

  if(!d->lens->Maker || d->crop <= 0.0f)
    piece->enabled = 0;
}

void init_pipe(dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
#ifdef HAVE_GEGL
#error "lensfun needs to be ported to GEGL!"
#else
  piece->data = malloc(sizeof(dt_iop_lensfun_data_t));
  dt_iop_lensfun_data_t *d = piece->data;
  d->lens = NULL;
  self->commit_params(self, self->default_params, pipe, piece);
#endif
}

void cleanup_pipe(dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
#ifdef HAVE_GEGL
#error "lensfun needs to be ported to GEGL!"
#else
  dt_iop_lensfun_data_t *d = piece->data;
  if(d->lens)
    lf_lens_destroy(d->lens);
  free(piece->data);
  piece->data = NULL;
#endif
}

void init_global(dt_iop_module_so_t *module)
{
  const int program = 2; // basic.cl, from programs.conf
  dt_iop_lensfun_global_data_t *gd
      = (dt_iop_lensfun_global_data_t *)calloc(1, sizeof(dt_iop_lensfun_global_data_t));
  module->data = gd;
  gd->kernel_lens_distort_bilinear =
    dt_opencl_create_kernel(program, "lens_distort_bilinear");
  gd->kernel_lens_distort_bicubic =
    dt_opencl_create_kernel(program, "lens_distort_bicubic");
  gd->kernel_lens_distort_lanczos2 =
    dt_opencl_create_kernel(program, "lens_distort_lanczos2");
  gd->kernel_lens_distort_lanczos3 =
    dt_opencl_create_kernel(program, "lens_distort_lanczos3");
  gd->kernel_lens_vignette =
    dt_opencl_create_kernel(program, "lens_vignette");

  lfDatabase *dt_iop_lensfun_db = lf_db_new();
  gd->db = (void *)dt_iop_lensfun_db;
#if defined(__MACH__) || defined(__APPLE__)
#else
  if(lf_db_load(dt_iop_lensfun_db) != LF_NO_ERROR)
#endif
  {
    char path[PATH_MAX] = { 0 };
    dt_loc_get_datadir(path, sizeof(path));
    char *c = path + strlen(path);
    for(; c > path && *c != '/'; c--)
      ;
    sprintf(c, "/lensfun");
    dt_iop_lensfun_db->HomeDataDir = g_strdup(path);
    if(lf_db_load(dt_iop_lensfun_db) != LF_NO_ERROR)
      fprintf(stderr, "[iop_lens]: could not load lensfun database in `%s'!\n", path);
  }
}

void reload_defaults(dt_iop_module_t *module)
{
  dt_iop_lensfun_params_t tmp;
  memset(tmp.lens, 0, sizeof(tmp.lens)); // auto
  memset(tmp.camera, 0, sizeof(tmp.camera)); // auto
  tmp.crop     = -1; // auto
  tmp.aperture = -1;
  tmp.focal    = -1;
  tmp.scale    = -1;
  tmp.inverse  = 0;
  tmp.modify_flags =
    /* LF_MODIFY_TCA | */ LF_MODIFY_VIGNETTING | LF_MODIFY_DISTORTION | LF_MODIFY_GEOMETRY | LF_MODIFY_SCALE;
  tmp.distance = -1; //auto
  tmp.target_geom = LF_RECTILINEAR;
  tmp.tca_override = 0;
  tmp.tca_r = 1.0;
  tmp.tca_b = 1.0;

  memcpy(module->params, &tmp, sizeof(dt_iop_lensfun_params_t));
  memcpy(module->default_params, &tmp, sizeof(dt_iop_lensfun_params_t));
  module->default_enabled = 0;
}

void init(dt_iop_module_t *module)
{
  module->params = malloc(sizeof(dt_iop_lensfun_params_t));
  module->default_params = malloc(sizeof(dt_iop_lensfun_params_t));
  module->default_enabled = 0;
  module->params_size = sizeof(dt_iop_lensfun_params_t);
  module->gui_data = NULL;
#if LENS_MOZ
  module->priority = 130; // changed
#else
  module->priority = 170; // changed
#endif
}

void cleanup(dt_iop_module_t *module)
{
  free(module->gui_data);
  module->gui_data = NULL;
  free(module->params);
  module->params = NULL;
}

void cleanup_global(dt_iop_module_so_t *module)
{
  dt_iop_lensfun_global_data_t *gd = module->data;
  lfDatabase *dt_iop_lensfun_db = gd->db;
  lf_db_destroy(dt_iop_lensfun_db);

  dt_opencl_free_kernel(gd->kernel_lens_distort_bilinear);
  dt_opencl_free_kernel(gd->kernel_lens_distort_bicubic);
  dt_opencl_free_kernel(gd->kernel_lens_distort_lanczos2);
  dt_opencl_free_kernel(gd->kernel_lens_distort_lanczos3);
  dt_opencl_free_kernel(gd->kernel_lens_vignette);
  free(module->data);
  module->data = NULL;
}


/// ############################################################
/// gui stuff: inspired by ufraws lensfun tab:


/* simple function to compute the floating-point precision
   which is enough for "normal use". The criteria is to have
   about 3 leading digits after the initial zeros.  */
static int precision(double x, double adj)
{
  x *= adj;

  if(x == 0) return 1;
  if(x < 1.0)
    if(x < 0.1)
      if(x < 0.01)
        return 5;
      else
        return 4;
    else
      return 3;
  else if(x < 100.0)
    if(x < 10.0)
      return 2;
    else
      return 1;
  else
    return 0;
}

/* -- ufraw ptr array functions -- */

static int ptr_array_insert_sorted(GPtrArray *array, const void *item, GCompareFunc compare)
{
  int length = array->len;
  g_ptr_array_set_size(array, length + 1);
  const void **root = (const void **)array->pdata;

  int m = 0, l = 0, r = length - 1;

  // Skip trailing NULL, if any
  if(l <= r && !root[r]) r--;

  while(l <= r)
  {
    m = (l + r) / 2;
    int cmp = compare(root[m], item);

    if(cmp == 0)
    {
      ++m;
      goto done;
    }
    else if(cmp < 0)
      l = m + 1;
    else
      r = m - 1;
  }
  if(r == m) m++;

done:
  memmove(root + m + 1, root + m, (length - m) * sizeof(void *));
  root[m] = item;
  return m;
}

static int ptr_array_find_sorted(const GPtrArray *array, const void *item, GCompareFunc compare)
{
  int length = array->len;
  void **root = array->pdata;

  int l = 0, r = length - 1;
  int m = 0, cmp = 0;

  if(!length) return -1;

  // Skip trailing NULL, if any
  if(!root[r]) r--;

  while(l <= r)
  {
    m = (l + r) / 2;
    cmp = compare(root[m], item);

    if(cmp == 0)
      return m;
    else if(cmp < 0)
      l = m + 1;
    else
      r = m - 1;
  }

  return -1;
}

static void ptr_array_insert_index(GPtrArray *array, const void *item, int index)
{
  const void **root;
  int length = array->len;
  g_ptr_array_set_size(array, length + 1);
  root = (const void **)array->pdata;
  memmove(root + index + 1, root + index, (length - index) * sizeof(void *));
  root[index] = item;
}

/* -- end ufraw ptr array functions -- */

/* -- camera -- */

static void camera_set(dt_iop_module_t *self, const lfCamera *cam)
{
  dt_iop_lensfun_gui_data_t *g = self->gui_data;
  dt_iop_lensfun_params_t *p = (dt_iop_lensfun_params_t *)self->params;
  gchar *fm;
  const char *maker, *model, *variant;
  char _variant[100];

  if(!cam)
  {
    gtk_button_set_label(GTK_BUTTON(g->camera_model), "");
    gtk_label_set_ellipsize(GTK_LABEL(gtk_bin_get_child(GTK_BIN(g->camera_model))), PANGO_ELLIPSIZE_END);
    g_object_set(G_OBJECT(g->camera_model), "tooltip-text", "", (char *)NULL);
    return;
  }

  memset(p->camera, 0, sizeof(p->camera));
  g_strlcpy(p->camera, cam->Model, sizeof(p->camera));
  p->crop = cam->CropFactor;
  g->camera = cam;

  maker = lf_mlstr_get(cam->Maker);
  model = lf_mlstr_get(cam->Model);
  variant = lf_mlstr_get(cam->Variant);

  if(model)
  {
    if(maker)
      fm = g_strdup_printf("%s, %s", maker, model);
    else
      fm = g_strdup_printf("%s", model);
    gtk_button_set_label(GTK_BUTTON(g->camera_model), fm);
    gtk_label_set_ellipsize(GTK_LABEL(gtk_bin_get_child(GTK_BIN(g->camera_model))), PANGO_ELLIPSIZE_END);
    g_free(fm);
  }

  if(variant)
    snprintf(_variant, sizeof(_variant), " (%s)", variant);
  else
    _variant[0] = 0;

  fm = g_strdup_printf(_("maker: %s\n"
                         "model: %s%s\n"
                         "mount: %s\n"
                         "crop factor: %.1f"),
                       maker, model, _variant, cam->Mount, cam->CropFactor);
  g_object_set(G_OBJECT(g->camera_model), "tooltip-text", fm, (char *)NULL);
  g_free(fm);
}

static void camera_menu_select(GtkMenuItem *menuitem, gpointer user_data)
{
  dt_iop_module_t *self = user_data;
  camera_set(self, (lfCamera *)g_object_get_data(G_OBJECT(menuitem), "lfCamera"));
  if(!darktable.gui->reset)
    dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void camera_menu_fill(dt_iop_module_t *self, const lfCamera *const *camlist)
{
  dt_iop_lensfun_gui_data_t *g = self->gui_data;
  unsigned i;
  GPtrArray *makers, *submenus;

  if(g->camera_menu)
  {
    gtk_widget_destroy(GTK_WIDGET(g->camera_menu));
    g->camera_menu = NULL;
  }

  /* Count all existing camera makers and create a sorted list */
  makers = g_ptr_array_new();
  submenus = g_ptr_array_new();
  for(i = 0; camlist[i]; i++)
  {
    GtkWidget *submenu, *item;
    const char *m = lf_mlstr_get(camlist[i]->Maker);
    int idx = ptr_array_find_sorted(makers, m, (GCompareFunc)g_utf8_collate);
    if(idx < 0)
    {
      /* No such maker yet, insert it into the array */
      idx = ptr_array_insert_sorted(makers, m, (GCompareFunc)g_utf8_collate);
      /* Create a submenu for cameras by this maker */
      submenu = gtk_menu_new();
      ptr_array_insert_index(submenus, submenu, idx);
    }

    submenu = g_ptr_array_index(submenus, idx);
    /* Append current camera name to the submenu */
    m = lf_mlstr_get(camlist[i]->Model);
    if(!camlist[i]->Variant)
      item = gtk_menu_item_new_with_label(m);
    else
    {
      gchar *fm = g_strdup_printf("%s (%s)", m, camlist[i]->Variant);
      item = gtk_menu_item_new_with_label(fm);
      g_free(fm);
    }
    gtk_widget_show(item);
    g_object_set_data(G_OBJECT(item), "lfCamera", (void *)camlist[i]);
    g_signal_connect(G_OBJECT(item), "activate", G_CALLBACK(camera_menu_select), self);
    gtk_menu_shell_append(GTK_MENU_SHELL(submenu), item);
  }

  g->camera_menu = GTK_MENU(gtk_menu_new());
  for(i = 0; i < makers->len; i++)
  {
    GtkWidget *item = gtk_menu_item_new_with_label(g_ptr_array_index(makers, i));
    gtk_widget_show(item);
    gtk_menu_shell_append(GTK_MENU_SHELL(g->camera_menu), item);
    gtk_menu_item_set_submenu(GTK_MENU_ITEM(item), (GtkWidget *)g_ptr_array_index(submenus, i));
  }

  g_ptr_array_free(submenus, TRUE);
  g_ptr_array_free(makers, TRUE);
}

static void parse_maker_model(const char *txt, char *make, size_t sz_make, char *model, size_t sz_model)
{
  const gchar *sep;

  while(txt[0] && isspace(txt[0]))
    txt++;
  sep = strchr(txt, ',');
  if(sep)
  {
    size_t len = sep - txt;
    if(len > sz_make - 1)
      len = sz_make - 1;
    memcpy(make, txt, len);
    make[len] = 0;

    while(*++sep && isspace(sep[0]))
      ;
    len = strlen(sep);
    if(len > sz_model - 1)
      len = sz_model - 1;
    memcpy(model, sep, len);
    model[len] = 0;
  }
  else
  {
    size_t len = strlen(txt);
    if(len > sz_model - 1) len = sz_model - 1;
    memcpy(model, txt, len);
    model[len] = 0;
    make[0] = 0;
  }
}

static void camera_menusearch_clicked(GtkWidget *button, gpointer user_data)
{
  dt_iop_module_t *self = user_data;
  dt_iop_lensfun_global_data_t *gd = self->data;
  lfDatabase *dt_iop_lensfun_db = gd->db;
  dt_iop_lensfun_gui_data_t *g = self->gui_data;

  (void)button;

  const lfCamera *const *camlist;
  dt_pthread_mutex_lock(&darktable.plugin_threadsafe);
  camlist = lf_db_get_cameras(dt_iop_lensfun_db);
  dt_pthread_mutex_unlock(&darktable.plugin_threadsafe);
  if(!camlist) return;
  camera_menu_fill(self, camlist);
  gtk_menu_popup(GTK_MENU(g->camera_menu), NULL, NULL, NULL, NULL, 0, gtk_get_current_event_time());
}

static void camera_autosearch_clicked(GtkWidget *button, gpointer user_data)
{
  dt_iop_module_t *self = user_data;
  dt_iop_lensfun_global_data_t *gd = self->data;
  lfDatabase *dt_iop_lensfun_db = gd->db;
  dt_iop_lensfun_gui_data_t *g = self->gui_data;
  char make[200], model[200];
  const gchar *txt = ((dt_iop_lensfun_params_t*)self->default_params)->camera;

  (void)button;

  if(txt[0] == '\0')
  {
    const lfCamera *const *camlist;
    dt_pthread_mutex_lock(&darktable.plugin_threadsafe);
    camlist = lf_db_get_cameras(dt_iop_lensfun_db);
    dt_pthread_mutex_unlock(&darktable.plugin_threadsafe);
    if(!camlist) return;
    camera_menu_fill(self, camlist);
  }
  else
  {
    parse_maker_model(txt, make, sizeof(make), model, sizeof(model));
    dt_pthread_mutex_lock(&darktable.plugin_threadsafe);
    const lfCamera **camlist = lf_db_find_cameras_ext(dt_iop_lensfun_db, make, model, 0);
    dt_pthread_mutex_unlock(&darktable.plugin_threadsafe);
    if(!camlist) return;
    camera_menu_fill(self, camlist);
    lf_free(camlist);
  }

  gtk_menu_popup(GTK_MENU(g->camera_menu), NULL, NULL, NULL, NULL, 0, gtk_get_current_event_time());
}

/* -- end camera -- */

static void lens_comboentry_focal_update(GtkWidget *widget, dt_iop_module_t *self)
{
  dt_iop_lensfun_params_t *p = (dt_iop_lensfun_params_t *)self->params;
  const char *text = dt_bauhaus_combobox_get_text(widget);
  if(dt_bauhaus_combobox_get(widget) == 0)
    p->focal = -1;
  else
    if(text)
      (void)sscanf(text, "%f", &p->focal);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void lens_comboentry_aperture_update(GtkWidget *widget, dt_iop_module_t *self)
{
  dt_iop_lensfun_params_t *p = (dt_iop_lensfun_params_t *)self->params;
  const char *text = dt_bauhaus_combobox_get_text(widget);
  if(dt_bauhaus_combobox_get(widget) == 0)
    p->aperture = -1;
  else
    if(text)
      (void)sscanf(text, "%f", &p->aperture);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void lens_comboentry_distance_update(GtkWidget *widget, dt_iop_module_t *self)
{
  dt_iop_lensfun_params_t *p = (dt_iop_lensfun_params_t *)self->params;
  const char *text = dt_bauhaus_combobox_get_text(widget);
  if(dt_bauhaus_combobox_get(widget) == 0)
    p->distance = -1;
  else
    if(text)
      (void)sscanf(text, "%f", &p->distance);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void delete_children(GtkWidget *widget, gpointer data)
{
  (void)data;
  gtk_widget_destroy(widget);
}

static void lens_set(dt_iop_module_t *self, const lfLens *lens)
{
  dt_iop_lensfun_gui_data_t *g = self->gui_data;
  gchar *fm;
  const char *maker, *model;
  unsigned i;

  gdouble focal_values[]
      = { -INFINITY, 4.5, 8,   10,  12,  14,  15,  16,  17,  18,  20,  24,  28,   30,      31,  35,
          38,        40,  43,  45,  50,  55,  60,  70,  75,  77,  80,  85,  90,   100,     105, 110,
          120,       135, 150, 200, 210, 240, 250, 300, 400, 500, 600, 800, 1000, INFINITY };
  gdouble aperture_values[]
      = { -INFINITY, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.8, 2,  2.2, 2.5, 2.8, 3.2, 3.4, 4,  4.5, 5.0,
          5.6,       6.3, 7.1, 8,   9, 10,  11,  13,  14,  16, 18,  20,  22,  25,  29,  32, 38,  INFINITY };

  if(!lens)
  {
    gtk_widget_set_sensitive(GTK_WIDGET(g->modflags), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(g->target_geom), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(g->scale), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(g->reverse), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(g->tca_r), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(g->tca_b), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(g->message), FALSE);

    gtk_container_foreach(GTK_CONTAINER(g->detection_warning), delete_children, NULL);

    GtkLabel *label;

    label = GTK_LABEL(
      gtk_label_new(_("camera/lens not found - please select manually")));

    g_object_set(G_OBJECT(label), "tooltip-text", _("try to locate your camera/lens in the above two menus"),
                 (char *)NULL);

    gtk_box_pack_start(GTK_BOX(g->detection_warning), GTK_WIDGET(label),
                       FALSE, FALSE, 0);

    gtk_widget_hide(g->lens_param_box);
    gtk_widget_show_all(g->detection_warning);
    return;
  }
  else
  {
    gtk_widget_set_sensitive(GTK_WIDGET(g->modflags), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(g->target_geom), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(g->scale), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(g->reverse), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(g->tca_r), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(g->tca_b), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(g->message), TRUE);
  }

  maker = lf_mlstr_get(lens->Maker);
  model = lf_mlstr_get(lens->Model);

  if(model)
  {
    if(maker)
      fm = g_strdup_printf("%s, %s", maker, model);
    else
      fm = g_strdup_printf("%s", model);
      gtk_button_set_label(GTK_BUTTON(g->lens_model), fm);
    gtk_label_set_ellipsize(GTK_LABEL(gtk_bin_get_child(GTK_BIN(g->lens_model))), PANGO_ELLIPSIZE_END);
    g_free(fm);
  }

  char focal[100], aperture[100], mounts[200];

  if(lens->MinFocal < lens->MaxFocal)
    snprintf(focal, sizeof(focal), "%g-%gmm", lens->MinFocal, lens->MaxFocal);
  else
    snprintf(focal, sizeof(focal), "%gmm", lens->MinFocal);
  if(lens->MinAperture < lens->MaxAperture)
    snprintf(aperture, sizeof(aperture), "%g-%g", lens->MinAperture, lens->MaxAperture);
  else
    snprintf(aperture, sizeof(aperture), "%g", lens->MinAperture);

  mounts[0] = 0;
  if(lens->Mounts)
    for(i = 0; lens->Mounts[i]; i++)
    {
      if(i > 0)
        g_strlcat(mounts, ", ", sizeof(mounts));
      g_strlcat(mounts, lens->Mounts[i], sizeof(mounts));
    }

  fm = g_strdup_printf(_("maker: %s\n"
                         "model: %s\n"
                         "focal range: %s\n"
                         "aperture: %s\n"
                         "crop factor: %.1f\n"
                         "type: %s\n"
                         "mounts: %s"),
                       maker ? maker : "?", model ? model : "?", focal, aperture, lens->CropFactor,
                       lf_get_lens_type_desc(lens->Type, NULL), mounts);
  g_object_set(G_OBJECT(g->lens_model), "tooltip-text", fm, (char *)NULL);
  g_free(fm);

  /* Create the focal/aperture/distance combo boxes */
  gtk_container_foreach(GTK_CONTAINER(g->lens_param_box), delete_children, NULL);

  int ffi = 1, fli = -1;
  for(i = 1; i < sizeof(focal_values) / sizeof(gdouble) - 1; i++)
  {
    if(focal_values[i] < lens->MinFocal) ffi = i + 1;
    if(focal_values[i] > lens->MaxFocal && fli == -1) fli = i;
  }
  if(focal_values[ffi] > lens->MinFocal)
  {
    focal_values[ffi - 1] = lens->MinFocal;
    ffi--;
  }
  if(lens->MaxFocal == 0 || fli < 0) fli = sizeof(focal_values) / sizeof(gdouble) - 2;
  if(focal_values[fli + 1] < lens->MaxFocal)
  {
    focal_values[fli + 1] = lens->MaxFocal;
    ffi++;
  }
  if(fli < ffi) fli = ffi + 1;


  GtkWidget *w;
  char txt[30];

  // focal length
  w = dt_bauhaus_combobox_new(self);
  dt_bauhaus_widget_set_label(w, NULL, _("mm"));
  g_object_set(G_OBJECT(w), "tooltip-text", _("focal length (mm)"), (char *)NULL);
  snprintf(txt, sizeof(txt), "auto");
  dt_bauhaus_combobox_add(w, txt);
  for(int k = 0; k < fli - ffi; k++)
  {
    snprintf(txt, sizeof(txt), "%.*f", precision(focal_values[ffi+k], 10.0), focal_values[ffi+k]);
    dt_bauhaus_combobox_add(w, txt);
  }
  g_signal_connect(G_OBJECT(w), "value-changed", G_CALLBACK(lens_comboentry_focal_update), self);
  gtk_box_pack_start(GTK_BOX(g->lens_param_box), w, TRUE, TRUE, 0);
  dt_bauhaus_combobox_set_editable(w, 1);
  g->cbe[0] = w;


  // f-stop
  ffi = 1, fli = sizeof(aperture_values) / sizeof(gdouble) - 1;
  for(i = 1; i < sizeof(aperture_values) / sizeof(gdouble) - 1; i++)
    if(aperture_values[i] < lens->MinAperture) ffi = i + 1;
  if(aperture_values[ffi] > lens->MinAperture)
  {
    aperture_values[ffi - 1] = lens->MinAperture;
    ffi--;
  }

  w = dt_bauhaus_combobox_new(self);
  dt_bauhaus_widget_set_label(w, NULL, _("f/"));
  g_object_set(G_OBJECT(w), "tooltip-text", _("f-number (aperture)"), (char *)NULL);
  snprintf(txt, sizeof(txt), "auto");
  dt_bauhaus_combobox_add(w, txt);
  for(size_t k = 0; k < fli - ffi; k++)
  {
    snprintf(txt, sizeof(txt), "%.*f", precision(aperture_values[ffi + k], 10.0), aperture_values[ffi + k]);
    dt_bauhaus_combobox_add(w, txt);
  }
  g_signal_connect(G_OBJECT(w), "value-changed", G_CALLBACK(lens_comboentry_aperture_update), self);
  gtk_box_pack_start(GTK_BOX(g->lens_param_box), w, TRUE, TRUE, 0);
  dt_bauhaus_combobox_set_editable(w, 1);
  g->cbe[1] = w;

  w = dt_bauhaus_combobox_new(self);
  dt_bauhaus_widget_set_label(w, NULL, _("d"));
  g_object_set(G_OBJECT(w), "tooltip-text", _("distance to subject"), (char *)NULL);
  snprintf(txt, sizeof(txt), "auto");
  dt_bauhaus_combobox_add(w, txt);
  float val = 0.25f;
  for(int k = 0; k < 25; k++)
  {
    if(val > 1000.0f) val = 1000.0f;
    snprintf(txt, sizeof(txt), "%.*f", precision(val, 10.0), val);
    dt_bauhaus_combobox_add(w, txt);
    if(val >= 1000.0f) break;
    val *= sqrtf(2.0f);
  }
  g_signal_connect(G_OBJECT(w), "value-changed", G_CALLBACK(lens_comboentry_distance_update), self);
  gtk_box_pack_start(GTK_BOX(g->lens_param_box), w, TRUE, TRUE, 0);
  dt_bauhaus_combobox_set_editable(w, 1);
  g->cbe[2] = w;

  gtk_widget_hide(g->detection_warning);
  gtk_widget_show_all(g->lens_param_box);
}

static void lens_menu_select(GtkMenuItem *menuitem, gpointer user_data)
{
  dt_iop_module_t *self = user_data;
  dt_iop_lensfun_params_t *p = (void *)self->params;
  lfLens *lens = g_object_get_data(G_OBJECT(menuitem), "lfLens");
  lens_set(self, lens);
  memset(p->lens, 0, sizeof(p->lens));
  g_strlcpy(p->lens, lens->Model, sizeof(p->lens));
  if(!darktable.gui->reset)
    dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void lens_menu_fill(dt_iop_module_t *self, const lfLens *const *lenslist)
{
  dt_iop_lensfun_gui_data_t *g = self->gui_data;
  unsigned i;
  GPtrArray *makers, *submenus;

  if(g->lens_menu)
  {
    gtk_widget_destroy(GTK_WIDGET(g->lens_menu));
    g->lens_menu = NULL;
  }

  /* Count all existing lens makers and create a sorted list */
  makers = g_ptr_array_new();
  submenus = g_ptr_array_new();
  for(i = 0; lenslist[i]; i++)
  {
    GtkWidget *submenu, *item;
    const char *m = lf_mlstr_get(lenslist[i]->Maker);
    int idx = ptr_array_find_sorted(makers, m, (GCompareFunc)g_utf8_collate);
    if(idx < 0)
    {
      /* No such maker yet, insert it into the array */
      idx = ptr_array_insert_sorted(makers, m, (GCompareFunc)g_utf8_collate);
      /* Create a submenu for lenses by this maker */
      submenu = gtk_menu_new();
      ptr_array_insert_index(submenus, submenu, idx);
    }

    submenu = g_ptr_array_index(submenus, idx);
    /* Append current lens name to the submenu */
    item = gtk_menu_item_new_with_label(lf_mlstr_get(lenslist[i]->Model));
    gtk_widget_show(item);
    g_object_set_data(G_OBJECT(item), "lfLens", (void *)lenslist[i]);
    g_signal_connect(G_OBJECT(item), "activate", G_CALLBACK(lens_menu_select), self);
    gtk_menu_shell_append(GTK_MENU_SHELL(submenu), item);
  }

  g->lens_menu = GTK_MENU(gtk_menu_new());
  for(i = 0; i < makers->len; i++)
  {
    GtkWidget *item = gtk_menu_item_new_with_label(g_ptr_array_index(makers, i));
    gtk_widget_show(item);
    gtk_menu_shell_append(GTK_MENU_SHELL(g->lens_menu), item);
    gtk_menu_item_set_submenu(GTK_MENU_ITEM(item), (GtkWidget *)g_ptr_array_index(submenus, i));
  }

  g_ptr_array_free(submenus, TRUE);
  g_ptr_array_free(makers, TRUE);
}

static void lens_menusearch_clicked(GtkWidget *button, gpointer user_data)
{
  dt_iop_module_t *self = user_data;
  dt_iop_lensfun_global_data_t *gd = self->data;
  lfDatabase *dt_iop_lensfun_db = gd->db;
  dt_iop_lensfun_gui_data_t *g = self->gui_data;
  const lfLens **lenslist;

  (void)button;

  dt_pthread_mutex_lock(&darktable.plugin_threadsafe);
  lenslist = lf_db_find_lenses_hd(dt_iop_lensfun_db, g->camera, NULL, NULL, LF_SEARCH_SORT_AND_UNIQUIFY);
  dt_pthread_mutex_unlock(&darktable.plugin_threadsafe);
  if(!lenslist) return;
  lens_menu_fill(self, lenslist);
  lf_free(lenslist);

  gtk_menu_popup(GTK_MENU(g->lens_menu), NULL, NULL, NULL, NULL, 0, gtk_get_current_event_time());
}

static void lens_autosearch_clicked(GtkWidget *button, gpointer user_data)
{
  dt_iop_module_t *self = user_data;
  dt_iop_lensfun_global_data_t *gd = self->data;
  lfDatabase *dt_iop_lensfun_db = gd->db;
  dt_iop_lensfun_gui_data_t *g = self->gui_data;
  const lfLens **lenslist;
  char make[200], model[200];
  const gchar *txt = ((dt_iop_lensfun_params_t *)self->default_params)->lens;

  (void)button;

  parse_maker_model(txt, make, sizeof(make), model, sizeof(model));
  dt_pthread_mutex_lock(&darktable.plugin_threadsafe);
  lenslist = lf_db_find_lenses_hd(dt_iop_lensfun_db, g->camera, make[0] ? make : NULL,
                                  model[0] ? model : NULL, LF_SEARCH_SORT_AND_UNIQUIFY);
  dt_pthread_mutex_unlock(&darktable.plugin_threadsafe);
  if(!lenslist) return;
  lens_menu_fill(self, lenslist);
  lf_free(lenslist);

  gtk_menu_popup(GTK_MENU(g->lens_menu), NULL, NULL, NULL, NULL, 0, gtk_get_current_event_time());
}

/* -- end lens -- */

static void target_geometry_changed(GtkWidget *widget, gpointer user_data)
{
  dt_iop_module_t *self = user_data;
  dt_iop_lensfun_params_t *p = (dt_iop_lensfun_params_t *)self->params;

  int pos = dt_bauhaus_combobox_get(widget);
  p->target_geom = pos + LF_UNKNOWN + 1;
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void modflags_changed(GtkWidget *widget, gpointer user_data)
{
  dt_iop_module_t *self = user_data;
  if(self->dt->gui->reset) return;
  dt_iop_lensfun_params_t *p = (dt_iop_lensfun_params_t *)self->params;
  dt_iop_lensfun_gui_data_t *g = self->gui_data;
  int pos = dt_bauhaus_combobox_get(widget);
  GList *modifiers = g->modifiers;
  while(modifiers)
  {
    // could use g_list_nth. this seems safer?
    dt_iop_lensfun_modifier_t *mm = modifiers->data;
    if(mm->pos == pos)
    {
      p->modify_flags = (p->modify_flags & ~LENSFUN_MODFLAG_MASK) | mm->modflag;
      dt_dev_add_history_item(darktable.develop, self, TRUE);
      break;
    }
    modifiers = g_list_next(modifiers);
  }
}

static void reverse_toggled(GtkWidget *widget, gpointer user_data)
{
  dt_iop_module_t *self = user_data;
  dt_iop_lensfun_params_t *p = (dt_iop_lensfun_params_t *)self->params;
  p->inverse = dt_bauhaus_combobox_get(widget);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void tca_changed(GtkWidget *slider, dt_iop_module_t *self)
{
  dt_iop_lensfun_params_t *p = (dt_iop_lensfun_params_t *)self->params;
  dt_iop_lensfun_gui_data_t *g = (dt_iop_lensfun_gui_data_t *)self->gui_data;
  const float val = dt_bauhaus_slider_get(slider);
  if(slider == g->tca_r)
    p->tca_r = val;
  else
    p->tca_b = val;
  if(p->tca_r != 1.0 || p->tca_b != 1.0) p->tca_override = 1;
  p->modified = 1;
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void scale_changed(GtkWidget *slider, gpointer user_data)
{
  dt_iop_module_t *self = user_data;
  dt_iop_lensfun_params_t *p = (void *)self->params;
  p->scale = dt_bauhaus_slider_get(slider);
  dt_bauhaus_widget_set_label(slider, NULL, _("scale"));
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void autoscale_pressed(GtkWidget *button, gpointer user_data)
{
  dt_iop_module_t *self = user_data;
  dt_iop_lensfun_gui_data_t *g = self->gui_data;
  dt_iop_lensfun_params_t *p = (void *)self->params;
  p->scale = -1;
  dt_iop_lensfun_data_t d;
  d.lens = NULL;
  auto_fill(self, p, &d);
  lf_lens_destroy(d.lens);
  dt_bauhaus_slider_set(g->scale, d.scale);
  dt_bauhaus_widget_set_label(g->scale, NULL, _("scale (auto)"));
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static gboolean draw(GtkWidget *widget, cairo_t *cr, dt_iop_module_t *self)
{
  dt_iop_lensfun_gui_data_t *g = self->gui_data;
  if(darktable.gui->reset) return FALSE;

  dt_pthread_mutex_lock(&g->lock);
  const int corrections_done = g->corrections_done;
  g->corrections_done = -1;
  dt_pthread_mutex_unlock(&g->lock);

  if(corrections_done == -1) return FALSE;

  char *message = "";
  GList *modifiers = g->modifiers;
  while(modifiers)
  {
    // could use g_list_nth. this seems safer?
    dt_iop_lensfun_modifier_t *mm = modifiers->data;
    if(mm->modflag == corrections_done)
    {
      message = mm->name;
      break;
    }
    modifiers = g_list_next(modifiers);
  }

  darktable.gui->reset = 1;
  gtk_label_set_text(g->message, message);
  darktable.gui->reset = 0;

  return FALSE;
}


void gui_init(struct dt_iop_module_t *self)
{
  self->gui_data = malloc(sizeof(dt_iop_lensfun_gui_data_t));
  dt_iop_lensfun_gui_data_t *g = self->gui_data;
  dt_iop_lensfun_params_t *p = (dt_iop_lensfun_params_t *)self->params;

  dt_pthread_mutex_init(&g->lock, NULL);

  g->camera = NULL;
  g->camera_menu = NULL;
  g->lens_menu = NULL;
  g->modifiers = NULL;

  dt_pthread_mutex_lock(&g->lock);
  g->corrections_done = -1;
  dt_pthread_mutex_unlock(&g->lock);

  // initialize modflags options
  int pos = -1;
  dt_iop_lensfun_modifier_t *modifier;
  modifier = g_malloc0(sizeof(dt_iop_lensfun_modifier_t));
  dt_utf8_strlcpy(modifier->name, _("none"), sizeof(modifier->name));
  g->modifiers = g_list_append(g->modifiers, modifier);
  modifier->modflag = LENSFUN_MODFLAG_NONE;
  modifier->pos = ++pos;

  modifier = g_malloc0(sizeof(dt_iop_lensfun_modifier_t));
  dt_utf8_strlcpy(modifier->name, _("all"), sizeof(modifier->name));
  g->modifiers = g_list_append(g->modifiers, modifier);
  modifier->modflag = LENSFUN_MODFLAG_ALL;
  modifier->pos = ++pos;

  modifier = g_malloc0(sizeof(dt_iop_lensfun_modifier_t));
  dt_utf8_strlcpy(modifier->name, _("distortion & TCA"), sizeof(modifier->name));
  g->modifiers = g_list_append(g->modifiers, modifier);
  modifier->modflag = LENSFUN_MODFLAG_DIST_TCA;
  modifier->pos = ++pos;

  modifier = g_malloc0(sizeof(dt_iop_lensfun_modifier_t));
  dt_utf8_strlcpy(modifier->name, _("distortion & vignetting"), sizeof(modifier->name));
  g->modifiers = g_list_append(g->modifiers, modifier);
  modifier->modflag = LENSFUN_MODFLAG_DIST_VIGN;
  modifier->pos = ++pos;

  modifier = g_malloc0(sizeof(dt_iop_lensfun_modifier_t));
  dt_utf8_strlcpy(modifier->name, _("TCA & vignetting"), sizeof(modifier->name));
  g->modifiers = g_list_append(g->modifiers, modifier);
  modifier->modflag = LENSFUN_MODFLAG_TCA_VIGN;
  modifier->pos = ++pos;

  modifier = g_malloc0(sizeof(dt_iop_lensfun_modifier_t));
  dt_utf8_strlcpy(modifier->name, _("only distortion"), sizeof(modifier->name));
  g->modifiers = g_list_append(g->modifiers, modifier);
  modifier->modflag = LENSFUN_MODFLAG_DIST;
  modifier->pos = ++pos;

  modifier = g_malloc0(sizeof(dt_iop_lensfun_modifier_t));
  dt_utf8_strlcpy(modifier->name, _("only TCA"), sizeof(modifier->name));
  g->modifiers = g_list_append(g->modifiers, modifier);
  modifier->modflag = LENSFUN_MODFLAG_TCA;
  modifier->pos = ++pos;

  modifier = g_malloc0(sizeof(dt_iop_lensfun_modifier_t));
  dt_utf8_strlcpy(modifier->name, _("only vignetting"), sizeof(modifier->name));
  g->modifiers = g_list_append(g->modifiers, modifier);
  modifier->modflag = LENSFUN_MODFLAG_VIGN;
  modifier->pos = ++pos;

  GtkWidget *button;

  self->widget = gtk_box_new(GTK_ORIENTATION_VERTICAL, DT_BAUHAUS_SPACE);

  g_signal_connect(G_OBJECT(self->widget), "draw", G_CALLBACK(draw), self);

  // camera selector
  GtkWidget *hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  g->camera_model = GTK_BUTTON(gtk_button_new_with_label(self->dev->image_storage.exif_model));
  dt_gui_key_accel_block_on_focus_connect(GTK_WIDGET(g->camera_model));
  gtk_label_set_ellipsize(GTK_LABEL(gtk_bin_get_child(GTK_BIN(g->camera_model))), PANGO_ELLIPSIZE_END);
  g_signal_connect(G_OBJECT(g->camera_model), "clicked", G_CALLBACK(camera_menusearch_clicked), self);
  gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(g->camera_model), TRUE, TRUE, 0);
  button = dtgtk_button_new(dtgtk_cairo_paint_solid_triangle, CPF_STYLE_FLAT | CPF_DIRECTION_DOWN);
  g->find_camera_button = button;
  gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 0);
  g_object_set(G_OBJECT(button), "tooltip-text", _("find camera"), (char *)NULL);
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(camera_autosearch_clicked), self);
  gtk_box_pack_start(GTK_BOX(self->widget), hbox, TRUE, TRUE, 0);

  // lens selector
  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  g->lens_model = GTK_BUTTON(gtk_button_new_with_label(self->dev->image_storage.exif_lens));
  dt_gui_key_accel_block_on_focus_connect(GTK_WIDGET(g->lens_model));
  gtk_label_set_ellipsize(GTK_LABEL(gtk_bin_get_child(GTK_BIN(g->lens_model))), PANGO_ELLIPSIZE_END);
  g_signal_connect(G_OBJECT(g->lens_model), "clicked", G_CALLBACK(lens_menusearch_clicked), self);
  gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(g->lens_model), TRUE, TRUE, 0);
  button = dtgtk_button_new(dtgtk_cairo_paint_solid_triangle, CPF_STYLE_FLAT | CPF_DIRECTION_DOWN);
  g->find_lens_button = GTK_WIDGET(button);
  gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 0);
  g_object_set(G_OBJECT(button), "tooltip-text", _("find lens"), (char *)NULL);
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(lens_autosearch_clicked), self);
  gtk_box_pack_start(GTK_BOX(self->widget), hbox, TRUE, TRUE, 0);


  // lens properties
  g->lens_param_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_pack_start(GTK_BOX(self->widget), g->lens_param_box, TRUE, TRUE, 0);

  // camera/lens not detected warning box
  g->detection_warning = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_pack_start(GTK_BOX(self->widget), g->detection_warning, TRUE, TRUE, 0);

  // selector for correction type (modflags): one or more out of
  // distortion, TCA, vignetting
  g->modflags = dt_bauhaus_combobox_new(self);
  dt_bauhaus_widget_set_label(g->modflags, NULL, _("corrections"));
  gtk_box_pack_start(GTK_BOX(self->widget), g->modflags, TRUE, TRUE, 0);
  g_object_set(G_OBJECT(g->modflags), "tooltip-text", _("which corrections to apply"), (char *)NULL);
  GList *l = g->modifiers;
  while(l)
  {
    dt_iop_lensfun_modifier_t *modifier = l->data;
    dt_bauhaus_combobox_add(g->modflags, modifier->name);
    l = g_list_next(l);
  }
  dt_bauhaus_combobox_set(g->modflags, 0);
  g_signal_connect(G_OBJECT(g->modflags), "value-changed", G_CALLBACK(modflags_changed), (gpointer)self);

  // target geometry
  g->target_geom = dt_bauhaus_combobox_new(self);
  dt_bauhaus_widget_set_label(g->target_geom, NULL, _("geometry"));
  gtk_box_pack_start(GTK_BOX(self->widget), g->target_geom, TRUE, TRUE, 0);
  g_object_set(G_OBJECT(g->target_geom), "tooltip-text", _("target geometry"), (char *)NULL);
  dt_bauhaus_combobox_add(g->target_geom, _("rectilinear"));
  dt_bauhaus_combobox_add(g->target_geom, _("fish-eye"));
  dt_bauhaus_combobox_add(g->target_geom, _("panoramic"));
  dt_bauhaus_combobox_add(g->target_geom, _("equirectangular"));
#if LF_VERSION >= ((0 << 24) | (2 << 16) | (6 << 8) | 0)
  dt_bauhaus_combobox_add(g->target_geom, _("orthographic"));
  dt_bauhaus_combobox_add(g->target_geom, _("stereographic"));
  dt_bauhaus_combobox_add(g->target_geom, _("equisolid angle"));
  dt_bauhaus_combobox_add(g->target_geom, _("thoby fish-eye"));
#endif
  g_signal_connect(G_OBJECT(g->target_geom), "value-changed", G_CALLBACK(target_geometry_changed),
                   (gpointer)self);

  // scale
  g->scale = dt_bauhaus_slider_new_with_range(self, 0.1, 2.0, 0.005, 1.0, 3);
  g_object_set(G_OBJECT(g->scale), "tooltip-text", _("auto scale"), (char *)NULL);
  dt_bauhaus_widget_set_label(g->scale, NULL, _("scale"));
  g_signal_connect(G_OBJECT(g->scale), "value-changed", G_CALLBACK(scale_changed), self);
  g_signal_connect(G_OBJECT(g->scale), "quad-pressed", G_CALLBACK(autoscale_pressed), self);
  gtk_box_pack_start(GTK_BOX(self->widget), g->scale, TRUE, TRUE, 0);
  dt_bauhaus_widget_set_quad_paint(g->scale, dtgtk_cairo_paint_refresh, 0);

  // reverse direction
  g->reverse = dt_bauhaus_combobox_new(self);
  dt_bauhaus_widget_set_label(g->reverse, NULL, _("mode"));
  gtk_box_pack_start(GTK_BOX(self->widget), g->reverse, TRUE, TRUE, 0);
  g_object_set(G_OBJECT(g->reverse), "tooltip-text", _("correct distortions or apply them"), (char *)NULL);
  dt_bauhaus_combobox_add(g->reverse, _("correct"));
  dt_bauhaus_combobox_add(g->reverse, _("distort"));
  g_signal_connect(G_OBJECT(g->reverse), "value-changed", G_CALLBACK(reverse_toggled), (gpointer)self);

  // override linear tca (if not 1.0):
  g->tca_r = dt_bauhaus_slider_new_with_range(self, 0.99, 1.01, 0.0001, p->tca_r, 5);
  g_object_set(G_OBJECT(g->tca_r), "tooltip-text", _("Transversal Chromatic Aberration red"), (char *)NULL);
  dt_bauhaus_widget_set_label(g->tca_r, NULL, _("TCA red"));
  g_signal_connect(G_OBJECT(g->tca_r), "value-changed", G_CALLBACK(tca_changed), self);
  gtk_box_pack_start(GTK_BOX(self->widget), g->tca_r, TRUE, TRUE, 0);

  g->tca_b = dt_bauhaus_slider_new_with_range(self, 0.99, 1.01, 0.0001, p->tca_b, 5);
  g_object_set(G_OBJECT(g->tca_b), "tooltip-text", _("Transversal Chromatic Aberration blue"), (char *)NULL);
  dt_bauhaus_widget_set_label(g->tca_b, NULL, _("TCA blue"));
  g_signal_connect(G_OBJECT(g->tca_b), "value-changed", G_CALLBACK(tca_changed), self);
  gtk_box_pack_start(GTK_BOX(self->widget), g->tca_b, TRUE, TRUE, 0);

  // message box to inform user what corrections have been done
  // this is useful as depending on lensfuns
  // profile only some of the lens flaws can be corrected
  GtkBox *hbox1 = GTK_BOX(gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0));
  GtkLabel *label = GTK_LABEL(gtk_label_new(_("corrections done: ")));
  g_object_set(G_OBJECT(label), "tooltip-text", _("which corrections have actually been done"), (char *)NULL);
  gtk_box_pack_start(GTK_BOX(hbox1), GTK_WIDGET(label), FALSE, FALSE, 0);
  g->message = GTK_LABEL(gtk_label_new("")); // This gets filled in by process
  gtk_box_pack_start(GTK_BOX(hbox1), GTK_WIDGET(g->message), FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(self->widget), GTK_WIDGET(hbox1), TRUE, TRUE, 0);
}

void gui_update(struct dt_iop_module_t *self)
{
  // let gui elements reflect params
  dt_iop_lensfun_gui_data_t *g = self->gui_data;
  dt_iop_lensfun_params_t *p = (dt_iop_lensfun_params_t *)self->params;
  dt_iop_lensfun_global_data_t *gd = self->data;
  lfDatabase *dt_iop_lensfun_db = gd->db;
  dt_iop_lensfun_data_t d;
  const dt_image_t *img = &self->dev->image_storage;

  d.lens = NULL;
  auto_fill(self, p, &d);
  lf_lens_destroy(d.lens);

  gtk_button_set_label(g->camera_model, p->camera);
  gtk_button_set_label(g->lens_model, p->lens);

  gtk_label_set_ellipsize(GTK_LABEL(gtk_bin_get_child(GTK_BIN(g->camera_model))), PANGO_ELLIPSIZE_END);
  gtk_label_set_ellipsize(GTK_LABEL(gtk_bin_get_child(GTK_BIN(g->lens_model))), PANGO_ELLIPSIZE_END);
  g_object_set(G_OBJECT(g->camera_model), "tooltip-text", "", (char *)NULL);
  g_object_set(G_OBJECT(g->lens_model), "tooltip-text", "", (char *)NULL);

  dt_pthread_mutex_lock(&g->lock);
  g->corrections_done = -1;
  dt_pthread_mutex_unlock(&g->lock);
  gtk_label_set_text(g->message, "");

  int modflag = p->modify_flags & LENSFUN_MODFLAG_MASK;
  GList *modifiers = g->modifiers;
  while(modifiers)
  {
    // could use g_list_nth. this seems safer?
    dt_iop_lensfun_modifier_t *mm = modifiers->data;
    if(mm->modflag == modflag)
    {
      dt_bauhaus_combobox_set(g->modflags, mm->pos);
      break;
    }
    modifiers = g_list_next(modifiers);
  }

  dt_bauhaus_combobox_set(g->target_geom, p->target_geom - LF_UNKNOWN - 1);
  dt_bauhaus_combobox_set(g->reverse, p->inverse);
  dt_bauhaus_slider_set(g->tca_r, p->tca_r);
  dt_bauhaus_slider_set(g->tca_b, p->tca_b);

  if(g->camera)
  {
    char make[200], model[200];
    if(p->lens[0])
    {
      const gchar *txt = gtk_button_get_label(GTK_BUTTON(g->lens_model));
      parse_maker_model(txt, make, sizeof(make), model, sizeof(model));
    }
    else
    {
      char txt[52];
      get_sanitized_lens(txt, sizeof(txt), img->exif_lens);
      parse_maker_model(txt, make, sizeof(make), model, sizeof(model));
    }
    dt_pthread_mutex_lock(&darktable.plugin_threadsafe);
    const lfLens **lenslist = lf_db_find_lenses_hd(
      dt_iop_lensfun_db, g->camera, make[0] ? make : NULL,
      model[0] ? model : NULL, LF_SEARCH_SORT_AND_UNIQUIFY);
    if(lenslist)
      lens_set(self, lenslist[0]);
    else
      lens_set(self, NULL);
    lf_free(lenslist);
    dt_pthread_mutex_unlock(&darktable.plugin_threadsafe);
  }
  else
  {
    dt_pthread_mutex_lock(&darktable.plugin_threadsafe);
    lens_set(self, NULL);
    dt_pthread_mutex_unlock(&darktable.plugin_threadsafe);
  }

  if(!p->camera[0])
    gtk_button_set_label(g->camera_model, "Camera from EXIF Data");
  if(!p->lens[0])
    gtk_button_set_label(g->lens_model, "Lens from EXIF Data");

  dt_bauhaus_slider_set(g->scale, d.scale);
  if(p->scale == -1)
    dt_bauhaus_widget_set_label(g->scale, NULL, _("scale (auto)"));
  else
    dt_bauhaus_widget_set_label(g->scale, NULL, _("scale"));
}

void gui_cleanup(struct dt_iop_module_t *self)
{
  dt_iop_lensfun_gui_data_t *g = self->gui_data;
  dt_gui_key_accel_block_on_focus_disconnect(GTK_WIDGET(g->lens_model));
  dt_gui_key_accel_block_on_focus_disconnect(GTK_WIDGET(g->camera_model));
  while(g->modifiers)
  {
    g_free(g->modifiers->data);
    g->modifiers = g_list_delete_link(g->modifiers, g->modifiers);
  }

  dt_pthread_mutex_destroy(&g->lock);

  free(self->gui_data);
  self->gui_data = NULL;
}

// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.sh
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-space on;
