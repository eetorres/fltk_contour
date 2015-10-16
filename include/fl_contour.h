// Latest modification 20/03/2008
// ========================================================================
// version: 2.1
// name:    fl_contour.h
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


#ifndef _Fl_Contour_H_
#define _Fl_Contour_H_

#include <geometry.h>
#include <fl_palette.h>
#include <msmvtl/tmatrix.h>
#include <msmvtl/tmmath.h>
#include <msmvtl/tmmanip.h>

extern const gm_uint typet[3][3][3];

class fl_contour{
	
public:
    //
    Geom_3D_Mesh<real>      triangle_color_mesh;
    Geom_Edges<gm_real,0,3> delaunay_vertex;
    Geom_Group<gm_real,4,3> square_grid;
    Geom_Contour<gm_real,3> line_3d_contour;
    vector<gm_real>         line_2d_contour;
    //
    gm_uint xcells, ycells, _lvls;
    gm_real dz, zmin;
    gm_short _choose;
    //
    fl_contour();
    virtual ~fl_contour();
    // Initialization functions
    void initialize_data(void);
    // Interpolation functions
    void inverse_distance_interp(void);
    void nearest_neighbor_interp(void);
    void lineal_interpolation(void);
    // eval
    void eval_contour_map(bool);
    void eval_color_map(void);
    // set
    void set_input_data(TMatrix<gm_real>& mtdata);
    void set_max_vertex(gm_uint u){ EMAX=u;};
    void set_remove_duplicated(bool b){ _delete_duplicate	= b;};
    void set_average_duplicated(bool b){ _average_duplicate 	= b;};
    void set_duplicated_tolerance(gm_real r){ _duplicate_tolerance = r*r;};
    void set_max_nearest_neighbour(gm_uint u){ _near_points = u;};
    void set_number_x_grid(gm_uint u){ xcells 	= u;};
    void set_number_y_grid(gm_uint u){ ycells 	= u;};
    void set_number_z_grid(gm_uint u){ _lvls 	= u;};
    void set_submesh_limits(int,int,real);
    void set_palette(gm_uint);
    // get
    gm_uint get_max_vertex(void){ return EMAX;};
    gm_real get_grid_limits(int _i, int _j){return _m_grid_limits[_i][_j];};
    gm_real get_x_min(void){return _m_grid_limits[0][0];};
    gm_real get_x_max(void){return _m_grid_limits[0][1];};
    gm_real get_y_min(void){return _m_grid_limits[1][0];};
    gm_real get_y_max(void){return _m_grid_limits[1][1];};
    gm_real get_z_min(void){return _m_grid_limits[2][0];};
    gm_real get_z_max(void){return _m_grid_limits[2][1];};
    gm_real get_submesh_limits(int,int);
    TMatrix<gm_real> get_normalize_submesh_limits(void);
    //
    void normalize_submesh_limits(void);
    void save_xyz_submesh(char*);

    fl_palette _cpalette;

private:

    virtual void draw(){};
    int handle(int event);
    //
    TMatrix<gm_real> _m_input_data, _m_grid_limits;
    TMatrix<gm_real> _cut_lim, _nor_lim;
    vector<gm_real>  _type_triangle_list;
    Geom_3D_Mesh<gm_real> _triangle_mesh;
    //
    gm_int CONTOUR_INIT_X, CONTOUR_INIT_Y;
    gm_int CONTOUR_INIT_WIDTH, CONTOUR_INIT_HEIGHT;
    gm_int POSY, POSX;
    gm_uint EMAX;
    //
    Geom_Edges<gm_uint,0,3>   _vtriangle;
    Geom_Edges<gm_real,4,3>   _vsquare;
    //
    Geom_Vertex<gm_real,2>    _cc;
    Geom_Vertex<gm_real,3>    _p1, _p2, _p3;
    Geom_Edges<gm_real,2,2>   _line;
    Geom_3D_Triangle<gm_real> _Trg;
    //
    gm_uint _near_points, _number_vertices;
    gm_uint trimax, l[3], _view_mode, _view_method, _number_triangles;
    gm_real xmin, xmax, ymin, ymax, zmax;
    gm_real _delta_x, _delta_y, dq, sin1, sin2;
    gm_real _duplicate_tolerance, _ss, _radio, zscale;
    //
    bool _found[4], _delete_duplicate, _average_duplicate;
    //
    void initialize_levels_();
    void compute_triangulation_(void);
    void eval_big_triangle_();
    void delaunay_xz_(void);
    //
    inline void line_interpolation_(gm_real, gm_uint);
    inline Geom_Vertex<gm_real,2> linear_interpolation_(gm_int,gm_int,gm_real,gm_real);
    //
    inline gm_uint triangle_type_(gm_real);
    inline gm_uint is_triange_in_triangle_(void);
    inline bool is_point_in_triangle_(gm_real,gm_real); 
    inline gm_int  which_side_(gm_real,gm_real,gm_real,gm_real,gm_real,gm_real);
    //
    inline bool eval_cell_curvature_(Geom_Edges<gm_real, 4, 3>&);
    inline bool eval_circumcircle_(Geom_Triangle<gm_real>& o,Geom_Vertex<gm_real,3>&);
    inline void eval_centroide_(Geom_Triangle<gm_real>&);
    inline gm_real eval_square_radius_(Geom_Triangle<gm_real>&);
    inline gm_real eval_plane_point_(gm_real,gm_real);
    //
    void add_3d_triangle_(void);
    void add_3d_triangle_mesh_(void);
    void add_2d_line_(Geom_Edges<gm_real, 2, 2>&, gm_real);
    void add_3d_line_(Geom_Edges<gm_real, 2, 2>&, gm_real);
    // Color functions
    void update_palette_(void);
    void set_view_mode_(const int mode);
    void set_method_mode_(const int mode);
    void set_z_levels_(const gm_uint);
    void set_2d_triangle_(Geom_Edges<gm_real, 4, 3>&, int);
};

#endif

////

