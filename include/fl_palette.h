// Latest modification 21/03/2008
// ========================================================================
// version: 0.1
// name:    fl_pallete.h
//
// Copyrigth 2008 by Edmanuel Torres A. (eetorres@gmail.com)
//
// This widget is  free  software;  you  can  redistribute  it and/or
// modify it  under  the  terms  of  the   GNU Library General Public
// License  as  published  by  the  Free  Software Foundation; either
// version 2 of the License,  or  (at your option)  any later version
// or much better FLTK license, which allow you static linking.
//
// This  library  is  distributed  in the hope that it will be useful,
// but  WITHOUT ANY WARRANTY;  without  even  the  implied warranty of
// MERCHANTABILITY  or FITNESS FOR A PARTICULAR PURPOSE.   See the GNU
// Library General Public License for more details.
//
// You should have  received a copy  of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA.
//
// Sent me suggestions, modifications or bugs. Don't hesitate to contact
// me for any question, I will be very grateful with your feedbacks.
// Get the last version at http://eetorres.googlepages.com/fltk_en
//
//========================================================================

#ifndef _Fl_Palette_H_
#define _Fl_Palette_H_

#include <msmvtl/tmatrix.h>
#include <msmvtl/tmmath.h>

typedef struct {
  double r,g,b;
} gm_rgb;

class fl_palette{

public:
    //
    unsigned int xcells, ycells, _lvls, _ipalette;
    double xdelta1, ydelta1, xdelta2, ydelta2;
    double y0, y1, x0, x1, x2, x3, x4, x5;
    double w0, dz, zmin, zmax;
    gm_rgb rgb;    
    gm_rgb cf1, cf2;

    std::vector<gm_rgb> _color_palette;

    fl_palette();
    ~fl_palette(){};
    // set
    void set_color(unsigned int);
    // get
    gm_rgb get_color(double);
    void initialize(double mn, double mx, unsigned int l){zmin=mn; zmax=mx; _lvls=l;};

//private:
    //
    gm_rgb  hsv_palette_(double val);
    gm_rgb  rgb_palette_(double val);
    gm_rgb  linear_palette_(double val);
    gm_rgb  earth_palette_(double val);
    gm_rgb  terrain_palette_(double val);
    //
    gm_rgb  palette_selection_(double val);
    double color_interpolation_(double,double,double);
    //
    // Color functions
    void update_palette_(void);
    void set(double w);
};

#endif

////

