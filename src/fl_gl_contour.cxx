// Latest modification 19/03/2008
// ========================================================================
// version: 1.9
// name:    fl_gl_contour.cxx
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


#include <fl_gl_contour.h>

#define ALPHA 0.5
#define HAVE_GL 1

#define _SHOW_LINES_ 1
//#define _SHOW_MESSAGES 1

#if HAVE_GL
fl_gl_contour::fl_gl_contour(int x,int y,int w,int h,const char *l)
            : Fl_Gl_Window(x,y,w,h,l)
#else
fl_gl_contour::fl_gl_contour(int x,int y,int w,int h,const char *l) : Fl_Box(x,y,w,h,l)
#endif /* HAVE_GL */
{
    init();    
}

#if HAVE_GL
fl_gl_contour::fl_gl_contour(int x,int y,int w,int h,TMatrix<gm_real>& mtx,const char *l)
            : Fl_Gl_Window(x,y,w,h,l)
#else
fl_gl_contour::fl_gl_contour(int x,int y,int w,int h,TMatrix<gm_real>& mtx,const char *l) : Fl_Box(x,y,w,h,l)
#endif /* HAVE_GL */
{
    init();
    set_input_data(mtx);
#if !HAVE_GL
    label("OpenGL is required for this demo to operate.");
    align(FL_ALIGN_WRAP | FL_ALIGN_INSIDE);
#endif /* !HAVE_GL */
}

void fl_gl_contour::init() {
    intp_meth    = 0;
    graph_meth   = 0;
    _gph_3d   = 0;
    _gph_2d   = 0;
    _choose      = 0;
    _mth_act     = 3;

    xcells   = 40;
    ycells   = 40;
    _lvls    = 40;
    

    vAng    = 95.0;//Z
    hAng    = 220.0;//XY
    
    xZoom   = 1.0;
    yZoom   = 1.0;
  
    _scl     = 1;
    _zoom    = 0.8;
    xshift   = 0.0;
    yshift   = 0.0;
    zshift   = 0.0;
    y_off    = 1.0;

    // Foreground colors
    fgred     = 0.5;
    fggreen   = 1.0;
    fgblue    = 0.5;
    // Background colors
    bgred     = 0.0;
    bggreen   = 0.0;
    bgblue    = 0.0;
  
    _depth     = true;
    is_n_data    = false;
    _draw_mpt    = false;
    _draw_dat    = false;
    _draw_box    = false;
    _draw_pal    = false;
    _draw_cut    = false;
    gm_light     = true;
    //
    _if_2d   = true;
    _if_3d   = true;
    
    set_max_vertex(500);
    set_max_nearest_neighbour(5);
    
#if !HAVE_GL
    label("WARNING: OpenGL is required for this demo to operate.");
    align(FL_ALIGN_WRAP | FL_ALIGN_INSIDE);
#endif /* !HAVE_GL */
}


////////////////////////////OPENGL////////////////////////////////

#if HAVE_GL

void fl_gl_contour::drawMesh(gm_bool _sw) {
    unsigned int l;
    gm_rgb _rgb;
    glColor3f(0.8,0.1,0.1);
    /*// Delaunay
    for(unsigned int k=0;k<num_tri;k++){
  glBegin(GL_LINE_STRIP);
  p1 = delaunay_vertex[vtri[k].x()];
  p2 = delaunay_vertex[vtri[k].y()];
  p3 = delaunay_vertex[vtri[k].z()];
  //glVertex3f( p1.x(), p1.y(), p1.z());
  //glVertex3f( p2.x(), p2.y(), p2.z());
  //glVertex3f( p3.x(), p3.y(), p3.z());
  //
  glVertex3f( p1.x(), p1.y(), 0.0);
  glVertex3f( p2.x(), p2.y(), 0.0);
  glVertex3f( p3.x(), p3.y(), 0.0);
  glVertex3f( p1.x(), p1.y(), 0.0);
  glEnd();
    }*/
    glLineWidth(1.0);
    for(l=0; l<triangle_color_mesh.size(); l++){
  if(_sw){
      glBegin(GL_LINE_STRIP);
      _rgb = _cpalette.get_color(triangle_color_mesh[l][0].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][0].x(), triangle_color_mesh[l][0].y(), _scl*triangle_color_mesh[l][0].z());
      _rgb = _cpalette.get_color(triangle_color_mesh[l][1].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][1].x(), triangle_color_mesh[l][1].y(), _scl*triangle_color_mesh[l][1].z());
      _rgb = _cpalette.get_color(triangle_color_mesh[l][2].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][2].x(), triangle_color_mesh[l][2].y(), _scl*triangle_color_mesh[l][2].z());
      _rgb = _cpalette.get_color(triangle_color_mesh[l][0].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][0].x(), triangle_color_mesh[l][0].y(), _scl*triangle_color_mesh[l][0].z());
      glEnd();
  }else{
      glBegin(GL_LINE_STRIP);
      _rgb = _cpalette.get_color(triangle_color_mesh[l][0].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][0].x(), triangle_color_mesh[l][0].y(), -0.5);
      _rgb = _cpalette.get_color(triangle_color_mesh[l][1].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][1].x(), triangle_color_mesh[l][1].y(), -0.5);
      _rgb = _cpalette.get_color(triangle_color_mesh[l][2].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][2].x(), triangle_color_mesh[l][2].y(), -0.5);
      _rgb = _cpalette.get_color(triangle_color_mesh[l][0].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][0].x(), triangle_color_mesh[l][0].y(), -0.5);
      glEnd();
  }
    }
    //drawobjects();
#if _SHOW_MESSAGES
    printf("drawMesh();\n");
#endif /* !HAVE_GL */
}// end drawMesh

void fl_gl_contour::drawMap(gm_bool _sw) {
    unsigned int l;
    gm_rgb _rgb;
    //
    //if(_draw_dat) drawData();
    //if(_draw_box) drawBox();
    glBegin(GL_TRIANGLES);
    for(l=0; l<triangle_color_mesh.size(); l++){
  if(_sw){
      _rgb = _cpalette.get_color(triangle_color_mesh[l][0].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][0].x(), triangle_color_mesh[l][0].y(), _scl*triangle_color_mesh[l][0].z());
  
      _rgb = _cpalette.get_color(triangle_color_mesh[l][1].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][1].x(), triangle_color_mesh[l][1].y(), _scl*triangle_color_mesh[l][1].z());
  
      _rgb = _cpalette.get_color(triangle_color_mesh[l][2].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][2].x(), triangle_color_mesh[l][2].y(), _scl*triangle_color_mesh[l][2].z());
  }else{
      _rgb = _cpalette.get_color(triangle_color_mesh[l][0].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][0].x(), triangle_color_mesh[l][0].y(), -0.5);
  
      _rgb = _cpalette.get_color(triangle_color_mesh[l][1].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][1].x(), triangle_color_mesh[l][1].y(), -0.5);
  
      _rgb = _cpalette.get_color(triangle_color_mesh[l][2].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][2].x(), triangle_color_mesh[l][2].y(), -0.5);
  }
    }
    glEnd();
    //drawobjects();    
#if _SHOW_MESSAGES
    printf("drawMesh();\n");
#endif /* !HAVE_GL */
}// end drawMap

void fl_gl_contour::drawCntr(gm_bool _sw) {
    gm_rgb _rgb;
    glLineWidth(1.0);
    glBegin(GL_LINES);
    for(unsigned int l=0; l<line_3d_contour.size(); l++){
  _rgb = _cpalette.get_color(line_2d_contour[l]-zmin);
  glColor3f(_rgb.r,_rgb.g,_rgb.b);
  if(_sw){
      glVertex3f(line_3d_contour[l][0].x(),line_3d_contour[l][0].y(),_scl*line_3d_contour[l][0].z());
      glVertex3f(line_3d_contour[l][1].x(),line_3d_contour[l][1].y(),_scl*line_3d_contour[l][1].z());
  }else{
      glVertex3f(line_3d_contour[l][0].x(),line_3d_contour[l][0].y(),-0.5);
      glVertex3f(line_3d_contour[l][1].x(),line_3d_contour[l][1].y(),-0.5);
  }
    }
    glEnd();
    //drawobjects();
#if _SHOW_MESSAGES
    printf("drawCntr();\n");
#endif /* !HAVE_GL */
}// end drawMesh

/*******************************************/
void fl_gl_contour::drawCntMap(gm_bool _sw) {
    unsigned int l;
    gm_rgb _rgb;
    //printf("Number of triangles = %i\n",triangle_color_mesh.size());
    glBegin(GL_TRIANGLES);
    if(_sw){
  for(l=0; l<triangle_color_mesh.size(); l++){
      _rgb = _cpalette.get_color(triangle_color_mesh[l][0].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][0].x(), triangle_color_mesh[l][0].y(), _scl*triangle_color_mesh[l][0].z());
      _rgb = _cpalette.get_color(triangle_color_mesh[l][1].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][1].x(), triangle_color_mesh[l][1].y(), _scl*triangle_color_mesh[l][1].z());
      _rgb = _cpalette.get_color(triangle_color_mesh[l][2].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][2].x(), triangle_color_mesh[l][2].y(), _scl*triangle_color_mesh[l][2].z());
  }
    }else{
  for(l=0; l<triangle_color_mesh.size(); l++){
      _rgb = _cpalette.get_color(triangle_color_mesh[l][0].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][0].x(), triangle_color_mesh[l][0].y(), -0.5);
      _rgb = _cpalette.get_color(triangle_color_mesh[l][1].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][1].x(), triangle_color_mesh[l][1].y(), -0.5);
      _rgb = _cpalette.get_color(triangle_color_mesh[l][2].z()-zmin);
      glColor3f(_rgb.r,_rgb.g,_rgb.b);
      glVertex3f( triangle_color_mesh[l][2].x(), triangle_color_mesh[l][2].y(), -0.5);
  }
    }
    glEnd();
    
    glColor3f(0.5,0.5,0.5);
    glLineWidth(2.0);
    glBegin(GL_LINES);
    for(unsigned int l=0; l<line_3d_contour.size(); l++){
  if(_sw){
      glVertex3f(line_3d_contour[l][0].x(),line_3d_contour[l][0].y(),_scl*line_3d_contour[l][0].z());
      glVertex3f(line_3d_contour[l][1].x(),line_3d_contour[l][1].y(),_scl*line_3d_contour[l][1].z());
  }else{
      glVertex3f(line_3d_contour[l][0].x(),line_3d_contour[l][0].y(),-0.5);
      glVertex3f(line_3d_contour[l][1].x(),line_3d_contour[l][1].y(),-0.5);
  }
    }
    glEnd();
    //drawobjects();
#if _SHOW_MESSAGES
    printf("drawMesh();\n");
#endif /* !HAVE_GL */
}// end drawMap

void fl_gl_contour::drawData(void){
    //gm_rgb _rgb;
    glPointSize(1.5);
    glBegin(GL_POINTS);
    glColor3f(1.0,0.0,0.0);
    for(unsigned int l=1; l<delaunay_vertex.size()-3; l++){
  //_rgb = _cpalette.get_color(delaunay_vertex[l].z());
        //glColor3f(_rgb.r,_rgb.g,_rgb.b);
        glVertex3f(delaunay_vertex[l].x(),delaunay_vertex[l].y(),_scl* (delaunay_vertex[l].z()-0.5) );
    }
    glEnd();
}

void fl_gl_contour::drawGridPoints(void){
    glColor3f(0.0,1.0,0.0);
    glPointSize(1.0);
    glBegin(GL_POINTS);
    for(unsigned int l=0; l<square_grid.size(); l++)
  for(unsigned int l2=0; l2<4; l2++){
      glVertex3f(square_grid[l][l2].x(),square_grid[l][l2].y(), _scl*(square_grid[l][l2].z()));
  }
    glEnd();
}

void fl_gl_contour::drawBox(void){
    glLineWidth(1.0);
    glBegin(GL_LINES);
    glColor3f(0.0,0.0,0.5);
    // top
    glVertex3f(-0.5,-0.5,0.5);
    glVertex3f( 0.5,-0.5,0.5);
    glVertex3f(-0.5,-0.5,0.5);
    glVertex3f(-0.5, 0.5,0.5);
    glVertex3f(-0.5,0.5,0.5);
    glVertex3f( 0.5,0.5,0.5);
    glVertex3f(0.5,-0.5,0.5);
    glVertex3f(0.5, 0.5,0.5);
    // bottom
    glVertex3f(-0.5,-0.5,-0.5);
    glVertex3f( 0.5,-0.5,-0.5);
    glVertex3f(-0.5,-0.5,-0.5);
    glVertex3f(-0.5, 0.5,-0.5);
    glVertex3f(-0.5,0.5,-0.5);
    glVertex3f( 0.5,0.5,-0.5);
    glVertex3f(0.5,-0.5,-0.5);
    glVertex3f(0.5, 0.5,-0.5);
    //vertex
    glVertex3f(-0.5,-0.5,-0.5);
    glVertex3f(-0.5,-0.5, 0.5);
    glVertex3f(0.5,0.5,-0.5);
    glVertex3f(0.5,0.5, 0.5);
    glVertex3f(-0.5,0.5,-0.5);
    glVertex3f(-0.5,0.5, 0.5);
    glVertex3f(0.5,-0.5,-0.5);
    glVertex3f(0.5,-0.5, 0.5);
    glEnd();
    // (1,1,1)
    glPointSize(5.0);
    glBegin(GL_POINTS);
    glColor3f(1.0,0.0,0.0); 
    glVertex3f(0.5,0.5,0.5); 
    glEnd();
}

void fl_gl_contour::drawCutMesh(void){
    glLineWidth(2.0);
    glColor3f(1.0,0.0,0.0);

    TMatrix<gm_real> _cut_box = get_normalize_submesh_limits();
    _xi = _cut_box[0][0];
    _xf = _cut_box[0][1];
    _yi = _cut_box[1][0];
    _yf = _cut_box[1][1];
    // Bottom
    glBegin(GL_LINE_LOOP);
    glVertex3f(_xi,_yi,-0.5);
    glVertex3f(_xi,_yf,-0.5);
    glVertex3f(_xf,_yf,-0.5);
    glVertex3f(_xf,_yi,-0.5);
    glEnd();
}

void fl_gl_contour::drawPalette(void){
    gm_rgb _rgb;
    double _ds;
    double _hor=0.45, _h=0.05;
    double _ver=-0.4, _w=0.8;
    
    glLoadIdentity();
    glPushMatrix();
    glLineWidth(1.0);
    _ds = _w/_lvls;
    glBegin(GL_QUADS);
    for(unsigned int _n=0;_n<_lvls;_n++){
  _rgb = _cpalette.get_color(_n*dz);
  glColor3f(_rgb.r,_rgb.g,_rgb.b);
  glVertex3f(_ver+_n*_ds, _hor   , 0.8);
  glVertex3f(_ver+_w    , _hor   , 0.8);
  glVertex3f(_ver+_w    , _hor+_h, 0.8);
  glVertex3f(_ver+_n*_ds, _hor+_h, 0.8);
    }
    glEnd();

    glBegin(GL_LINES);
    glColor3f(0.5,0.5,0.0);
    glVertex3f(_ver,_hor+_h,0.8);
    glVertex3f(_ver+_w, _hor+_h, 0.8);
    glVertex3f(_ver+_w, _hor   , 0.8);
    glVertex3f(_ver+_w, _hor+_h, 0.8);
    glVertex3f(_ver   , _hor   , 0.8);
    glVertex3f(_ver   , _hor+_h, 0.8);
    glVertex3f(_ver   , _hor   , 0.8);
    glVertex3f(_ver+_w, _hor   , 0.8);
    glEnd();
    glPopMatrix();
}

void fl_gl_contour::drawobjects(void){
    if(_draw_cut) drawCutMesh();
    if(_draw_dat) drawData();
    if(_draw_box) drawBox();
    if(_draw_mpt) drawGridPoints();
    if(_draw_pal) drawPalette();
}

void fl_gl_contour::draw() {
    if (!valid()) {
        glLoadIdentity();
      glClearColor(bgred,bggreen,bgblue,1.0);
  //if(_depth){
      glEnable(GL_DEPTH_TEST);    // Enable Depth testing
      glDepthFunc(GL_LEQUAL);
      glShadeModel(GL_SMOOTH);             // Use smooth shading
      glHint(GL_SHADE_MODEL,GL_NICEST); // Set the smooth shaiding to the best we can have
  //}
  glViewport (0, 0, w(), h());
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity();
  if (w() <= h())
          glOrtho (-0.6, 0.6, -0.6*(GLfloat)h()/(GLfloat)w(), 0.6*(GLfloat)h()/(GLfloat)w(), -10.0, 10.0);
  else
          glOrtho (-0.6*(GLfloat)w()/(GLfloat)h(), 0.6*(GLfloat)w()/(GLfloat)h(), -0.6, 0.6, -10.0, 10.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
    }
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPushMatrix();
    //
    glTranslatef(xshift, yshift, 0);
    glRotatef(hAng,0,1,0); 
    glRotatef(vAng,-1,0,0);
    glScalef(_zoom,_zoom,_zoom);

    if(is_n_data){
  if(_if_3d)
  switch(_gph_3d){
      case 0:
    drawMesh(true);
      break;
      case 1:
    drawCntr(true);
      break;
      case 2:
    drawMap(true);
      break;
      default:
    drawCntMap(true);
      break;
  }
  if(_if_2d)
  switch(_gph_2d){
      case 0:
    drawMesh(false);
      break;
      case 1:
    drawCntr(false);
      break;
      case 2:
    drawMap(false);
      break;
      default:
    drawCntMap(false);
      break;
  }
  drawobjects();
    }
  
    //drawMesh();
    glPopMatrix();
    glFlush();
    redraw();
}

#endif /* HAVE_GL */

////////////////////////////HANDLE EVENTS////////////////////////////////

int fl_gl_contour::handle(int event) {
    static int last_x;
    static int last_y;
    int delta_x, delta_y;
  //... position in Fl::event_x() and Fl::event_y()
  // get current mouse position and process event
    int x = Fl::event_x();
    int y = Fl::event_y();
  
  switch(event) {
  case FL_PUSH:
    //... mouse down event ...
    // save mouse position to track drag events
    last_x = x;
    last_y = y;
    return 1;
  case FL_DRAG:
    delta_x = x - last_x;
    delta_y = y - last_y;
    last_x = x;
    last_y = y;
    hAng += 0.2*delta_x;
    vAng += 0.2*delta_y;
    redraw();
    return 1;
    /*case FL_RELEASE:   
    ... mouse up event ...
    return 1;
    case FL_FOCUS :
    case FL_UNFOCUS :
    ... Return 1 if you want keyboard events, 0 otherwise
    return 1;
    case FL_KEYBOARD:
    ... keypress, key is in Fl::event_key(), ascii in Fl::event_text()
    ... Return 1 if you understand/use the keyboard event, 0 otherwise...
    return 1;
    case FL_SHORTCUT:
    ... shortcut, key is in Fl::event_key(), ascii in Fl::event_text()
    ... Return 1 if you understand/use the shortcut event, 0 otherwise...
    return 1;*/
  default:
    // pass other events to the base class...
    return Fl_Gl_Window::handle(event);
  }
}

//////////////////////////////UTILS///////////////////////////////

void fl_gl_contour::bgColor(float r,float g,float b){
    bgred=r; bggreen=g; bgblue=b;
}

void fl_gl_contour::fgColor(float r,float g,float b){
    fgred=r; fggreen=g; fgblue=b;
}

void fl_gl_contour::graph_cb(void){
    initialize_data();
    switch(intp_meth) {
  case 0:
      if(_mth_act != 0){
    inverse_distance_interp();
    _mth_act = 0;
      }
  break;
  case 1:
      if(_mth_act != 1){
    nearest_neighbor_interp();
    _mth_act = 1;
      }
  break;
  default:
      if(_mth_act != 2){
    lineal_interpolation();
    _mth_act = 2;
      }
  break;
    }
    //
    switch(_gph_3d) {
  case 0:
      eval_color_map();
  break;
  case 1:
          eval_contour_map(GM_3D);
  break;
  case 2:
      eval_color_map();
  break;
  default:
      eval_contour_map(GM_3D);
      eval_color_map();
  break;
    }
    if(_gph_2d != _gph_3d) 
    switch(_gph_2d) {
  case 0:
      eval_color_map();
  break;
  case 1:
          eval_contour_map(GM_3D);
  break;
  case 2:
      eval_color_map();
  break;
  default:
      eval_contour_map(GM_3D);
      eval_color_map();
  break;
    }
    
    redraw();
    is_n_data = 1;
#if _SHOW_MESSAGES
    printf("graph_cb();\n");
#endif /* !HAVE_GL */
}

/////




////

