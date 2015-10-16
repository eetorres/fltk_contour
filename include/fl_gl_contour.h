// Latest modification 19/03/2008
// ========================================================================
// version: 1.9
// name:    fl_gl_contour.h
// note:    this widget is a c++ class based on the Fast Light Tool Kit 
//          (FLTK) - www.fltk.org
//
// Copyrigth 2002 by Edmanuel Torres A. (eetorres@gmail.com)
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


#ifndef _Fl_2D_Gl_Contour_H_
#define _Fl_Gl_Contour_H_

#include <FL/Fl.H>
#include <geometry.h>
#include <fl_contour.h>

#define HAVE_GL 1

#if HAVE_GL
    #include <FL/Fl_Gl_Window.H>
    #include <FL/gl.h>
#else
    #include <FL/Fl_Box.H>
#endif /* HAVE_GL */

#include <stdlib.h>

#if HAVE_GL
class fl_gl_contour : public Fl_Gl_Window , public fl_contour{
#else
class fl_gl_contour : public Fl_Box {
#endif /* HAVE_GL */

public:
    
    fl_gl_contour(int x,int y,int w,int h,const char *l=0);
    fl_gl_contour(int,int,int,int,TMatrix<gm_real>&,const char* l=0);

    void graph_cb(void);
    void graph_pw(void);
    void drawData(void);
    void drawBox(void);
    void drawPalette(void);
    void drawGridPoints(void);
    void drawobjects(void);
//    void redraw(void){draw();};
    void graph_2d(bool _b){ 	_if_2d=_b;};
    void graph_3d(bool _b){ 	_if_3d=_b;};
    void graph_cut(bool _b){ 	_draw_cut = _b;};
    void graph_dat(bool _b){ 	_draw_dat = _b;};
    void graph_mpt(bool _b){ 	_draw_mpt = _b;};
    void graph_box(bool _b){ 	_draw_box = _b;};
    void graph_pal(bool _b){ 	_draw_pal = _b;};
    void graph_dpt(bool _b){ 	_depth = _b;};
    
    void bgColor(float r,float g,float b);
    void fgColor(float r,float g,float b);
    // Scale function
    void scale(float _f){	_scl = _f;};
    void zoom(float _f){	_zoom = _f;};
    void xzoom(float _f){	xZoom = _f;};
    void yzoom(float _f){	yZoom = _f;};
    // Set the rotation about the vertical (y ) axis
    void v_angle(float angle){	vAng=angle;};
    // Return the rotation about the vertical (y ) axis
    float v_angle(){	return vAng;};
    // Set the rotation about the horizontal (x ) axis
    void h_angle(float angle){	hAng=angle;};
    // the rotation about the horizontal (x ) axis
    float h_angle(){return hAng;};
    void panx(float _f){	xshift=_f;};
    void pany(float _f){	yshift=_f;};
    void actualize(void){	_mth_act=3;};
    //
    void intp_method(gm_int _i){  intp_meth = _i;};
    void graph_method(gm_int _i){ graph_meth = _i;};
    void graph_3d(gm_int _i){ _gph_3d = _i;};
    void graph_2d(gm_int _i){ _gph_2d = _i;};

#if HAVE_GL
    void drawMesh(gm_bool);
    void drawCntr(gm_bool);
    void drawMap(gm_bool);
    void drawCntMap(gm_bool);
#else
    void drawMesh() {	printf("Is not posible start any OpenGL device interface\n");}
    void drawCntr() {	printf("Is not posible start any OpenGL device interface\n");}
    void drawMap() {	printf("Is not posible start any OpenGL device interface\n");}
#endif /* HAVE_GL */
#if HAVE_GL
    void draw();    
#endif /* HAVE_GL */

protected:
    int handle(int);

private:
    gm_real _xi, _xf, _yi, _yf;
    gm_int graph_meth, intp_meth, _gph_2d, _gph_3d;
    unsigned int _mth_act;
    
    bool _depth, gm_light, _if_2d, _if_3d;
    bool is_n_data, _draw_mpt, _draw_cut;
    bool  _draw_dat, _draw_box, _draw_pal;            //
        
    float vAng,hAng, y_off, _scl;
    float xZoom,yZoom, _zoom;
    float xshift,yshift,zshift;
    float fgred, fggreen, fgblue;
    float bgred, bggreen, bgblue;
    float boxv0[3];
    float boxv1[3];
    ////
    void init();
    void drawCutMesh(void);
    
};

#endif

