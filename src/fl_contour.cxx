// Latest modification 20/03/2008
// ========================================================================
// version: 2.1
// name:    gl_contour.cxx
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

#include <fl_contour.h>
//#define _SHOW_MESSAGE_ 1
//#define _DEBUG 1

const gm_uint typet[3][3][3] = {{{0,0,8},{0,2,6},{9,5,7}},{{0,3,4},{1,0,1},{4,3,0}},{{7,5,9},{6,3,0},{8,0,0}}};
const gm_uint sidet[4][3] = {{0,2,1},{0,3,2},{3,2,1},{1,0,3}};

fl_contour::fl_contour(){
    _delete_duplicate = false;
    _average_duplicate = true;
    _duplicate_tolerance = 1e-9;
    set_palette(0);
}

fl_contour::~fl_contour(){
}

void fl_contour::set_input_data(TMatrix<gm_real>& mdt){
    //TMatrix<real> _mmt;
    real _ss1, _ss2;
    //
    _m_grid_limits = minterval(mdt); // getting scaling parameters
    _cut_lim = _m_grid_limits;
    //nor_lim();
    //
    _ss1 = (_m_grid_limits[0][1] - _m_grid_limits[0][0]);
    _ss2 = (_m_grid_limits[1][1] - _m_grid_limits[1][0]);
    _ss = max(_ss1,_ss2);
    //cout<<"Scaling parameters="<<_m_grid_limits;
    //cout<<"INPUT="<<mdt;
    _m_input_data = mnormxy(mdt,_m_grid_limits,_ss);
    //cout<<"DATA="<<_m_input_data;
    //_mmt = mscalexy(_m_input_data,_m_grid_limits,_ss);
    //cout<<"SCALE="<<_mmt;
    initialize_data();
#ifdef _SHOW_MESSAGE_
    for(gm_uint i=1;i<_m_input_data.cols();i++){
        printf("(%f,%f,%f)\n",_m_input_data[i][0],_m_input_data[i][1],_m_input_data[i][2]);
    }
    printf("input_data\n");
#endif
    normalize_submesh_limits();
}

void fl_contour::initialize_data(void){
    compute_triangulation_();
    initialize_levels_();
#ifdef _SHOW_MESSAGE_
    printf("initialize_data\n");
#endif
}

void fl_contour::compute_triangulation_(void){
    gm_uint i, j;
    Geom_Vertex<gm_real,3> t3p;
    gm_int duplicate = 0;
    gm_uint _duplicate_value;

    xmin = _m_input_data.get_col(0).get_min();
    xmax = _m_input_data.get_col(0).get_max();
    ymin = _m_input_data.get_col(1).get_min();
    ymax = _m_input_data.get_col(1).get_max();
#ifdef _SHOW_MESSAGE_
    printf("size = %i - cols = %i\n",_m_input_data.rows(),_m_input_data.cols());
    printf("xmin = %f\n",xmin);
    printf("xmax = %f\n",xmax);
    printf("ymin = %f\n",ymin);
    printf("ymax = %f\n",ymax);
#endif
  _number_vertices=_m_input_data.rows();
  delaunay_vertex.resize(_number_vertices+4);
  trimax = 20 * (_number_vertices+10);
  _vtriangle.resize(trimax);
  for(i=0;i<_number_vertices;i++){
    t3p.x()=_m_input_data[i][0];
    t3p.y()=_m_input_data[i][1];
    t3p.z()=_m_input_data[i][2];
    delaunay_vertex.addlnpnt(t3p,i+1);
  }
  // sort the samples in x coordinate
  for(i=1;i<_number_vertices;i++){
    for(j=i+1;j<=_number_vertices;j++){
      if(delaunay_vertex[j].x() < delaunay_vertex[i].x()){
        t3p = delaunay_vertex[i];
        delaunay_vertex[i] = delaunay_vertex[j];
        delaunay_vertex[j] = t3p;
      }
    }
  }
  if(_delete_duplicate || _average_duplicate){
    // check for coincident values and delete them
    real _x, _y, _rr;
    for(i=1;i<_number_vertices;i++){
      _duplicate_value = 1;
      for(j=i+1;j<=_number_vertices;j++){
         _x = (delaunay_vertex[j].x() - delaunay_vertex[i].x());
         _y = (delaunay_vertex[j].y() - delaunay_vertex[i].y());
         _rr = (_x*_x)+(_y*_y);
         if(_rr < _duplicate_tolerance){
           if(_delete_duplicate){
             duplicate++;
           }else if(_average_duplicate){
             duplicate++;
             _duplicate_value++;
             delaunay_vertex[i][2] += delaunay_vertex[j][2];
           }
           delaunay_vertex[j] = delaunay_vertex[_number_vertices];
           delaunay_vertex.null(_number_vertices);
           delaunay_vertex.delete_last();
           j--;
           _number_vertices--;
        }
      }
      if(_average_duplicate && (_duplicate_value>1))
      delaunay_vertex[i][2] = delaunay_vertex[i][2]/(real)_duplicate_value;
    }
    // sort the samples in x coordinate
    for(i=1;i<_number_vertices;i++){
      for(j=i+1;j<=_number_vertices;j++){
        if(delaunay_vertex[j].x() < delaunay_vertex[i].x()){
          t3p        = delaunay_vertex[i];
          delaunay_vertex[i] = delaunay_vertex[j];
          delaunay_vertex[j] = t3p;
        }
      }
    }
  }
#ifdef _SHOW_MESSAGE_
  for(i=1;i<delaunay_vertex.size();i++){
    printf("(%f,%f)=%f\n",delaunay_vertex[i].x(),delaunay_vertex[i].y(),delaunay_vertex[i].z());
  }
  printf("vertex size = %i - duplicates = %i\n",delaunay_vertex.size(),duplicate);
#endif
#ifdef _DEBUG
  printf("Performing a Delaunay triangulation\n");
#endif
  // perform a Delaunay triangulation
  delaunay_xz_();
}

void fl_contour::initialize_levels_(){
  _delta_x    = (xmax - xmin) / xcells;
  _delta_y    = (ymax - ymin) / ycells;
  dq    = sqrt(_delta_x*_delta_x+_delta_y*_delta_y);
  sin1  = _delta_x/dq;
  sin2  = _delta_y/dq;
#ifdef _SHOW_MESSAGE_
  printf("initialize_levels_\n");
#endif
}

// set levels in "_lvls"
void fl_contour::set_z_levels_(const gm_uint _lv){
  _lvls = _lv;
#ifdef _SHOW_MESSAGE_
  printf("set_lvls\n");
#endif
}

// set graphic mode
void fl_contour::set_view_mode_(const int mode){
  _view_mode = mode;
#ifdef _SHOW_MESSAGE_
  printf("set graph mode\n");
#endif
}

// set display mode
void fl_contour::set_method_mode_(const int mode){
  _view_method = mode;
#ifdef _SHOW_MESSAGE_
  printf("set calc mode\n");
#endif
}

void fl_contour::eval_contour_map(bool _op){
  gm_uint i;
  gm_uint _lvl, tp;
  float _lv;
  update_palette_();
  if(_op)line_3d_contour.clear();
  line_2d_contour.clear();
  for(_lvl=0; _lvl<=_lvls; _lvl++){ // Contour _lvl
    _lv = _lvl*dz;
    for(i=0; i<square_grid.size(); i++){
      if(eval_cell_curvature_(square_grid[i])){
        set_2d_triangle_(square_grid[i],1);
        tp=triangle_type_(_lv);
        if(tp){
          line_interpolation_(_lv, tp);
          if(_op) add_3d_line_(_line,_lv+zmin);
          else add_2d_line_(_line,_lv+zmin);
        }
        set_2d_triangle_(square_grid[i],2);
        tp=triangle_type_(_lv);
        if(tp){
          line_interpolation_(_lv, tp);
          if(_op) add_3d_line_(_line,_lv+zmin);
          else add_2d_line_(_line,_lv+zmin);
        }
      }else{
        set_2d_triangle_(square_grid[i],3);
        tp=triangle_type_(_lv);
        if(tp){
          line_interpolation_(_lv, tp);
          if(_op) add_3d_line_(_line,_lv+zmin);
          else add_2d_line_(_line,_lv+zmin);
        }
        set_2d_triangle_(square_grid[i],4);
        tp=triangle_type_(_lv);
        if(tp){
          line_interpolation_(_lv, tp);
          if(_op) add_3d_line_(_line,_lv+zmin);
          else add_2d_line_(_line,_lv+zmin);
        }
      }
    }
    for(i=0; i<_triangle_mesh.size(); i++){
      _Trg = _triangle_mesh[i];
      if(_type_triangle_list[i]==1){
        //set_2d_triangle_(square_grid[i],1);
        tp=triangle_type_(_lv);
        if(tp){
          line_interpolation_(_lv, tp);
          if(_op) add_3d_line_(_line,_lv+zmin);
          else add_2d_line_(_line,_lv+zmin);
        }
      }else if(_type_triangle_list[i]==2){
        //set_2d_triangle_(square_grid[i],2);
        tp=triangle_type_(_lv);
        if(tp){
          line_interpolation_(_lv, tp);
          if(_op) add_3d_line_(_line,_lv+zmin);
          else add_2d_line_(_line,_lv+zmin);
        }
      }else if(_type_triangle_list[i]==3){
        //set_2d_triangle_(square_grid[i],3);
        tp=triangle_type_(_lv);
        if(tp){
          line_interpolation_(_lv, tp);
          if(_op) add_3d_line_(_line,_lv+zmin);
          else add_2d_line_(_line,_lv+zmin);
        }
      }else if(_type_triangle_list[i]==4){
        //set_2d_triangle_(square_grid[i],4);
        tp=triangle_type_(_lv);
        if(tp){
          line_interpolation_(_lv, tp);
          if(_op) add_3d_line_(_line,_lv+zmin);
          else add_2d_line_(_line,_lv+zmin);
        }
      }
    }
  }
#ifdef _SHOW_MESSAGE_
  printf("eval_contour_map()\n");
#endif
}

void fl_contour::eval_color_map(){
  float _lv;
  gm_uint i;
  triangle_color_mesh.clear();
  update_palette_();
  for(gm_uint _lvl=0; _lvl<=_lvls; _lvl++){
    // Contour _lvl
    _lv = _lvl*dz;
    for(i=0; i<square_grid.size(); i++){
      if(eval_cell_curvature_(square_grid[i])){
        set_2d_triangle_(square_grid[i],1);
        add_3d_triangle_();
        set_2d_triangle_(square_grid[i],2);
        add_3d_triangle_();
      }else{
        set_2d_triangle_(square_grid[i],3);
        add_3d_triangle_();
        set_2d_triangle_(square_grid[i],4);
        add_3d_triangle_();
      }
    }
    for(i=0; i<_triangle_mesh.size(); i++){
      _Trg = _triangle_mesh[i];
      if(_type_triangle_list[i]==1){
        add_3d_triangle_();
      }else if(_type_triangle_list[i]==2){
        add_3d_triangle_();
      }else if(_type_triangle_list[i]==3){
        add_3d_triangle_();
      }else if(_type_triangle_list[i]==4){
        add_3d_triangle_();
      }
    }
  }
#ifdef _SHOW_MESSAGE_
  printf("eval_color_map\n");
#endif
}

void fl_contour::save_xyz_submesh(char* _f){
  vector<gm_real> _vr;
  TVector<gm_real> _vc(3);
  TMatrix<gm_real> _mcut(0,3);
  for(unsigned int l=0; l<square_grid.size(); l++){
    for(unsigned int l2=0; l2<4; l2++){
      //cout<<"x="<<square_grid[l][l2].x();
      //cout<<"y="<<square_grid[l][l2].y();
      //cout<<"z="<<square_grid[l][l2].z()<<endl;
      //_vr = square_grid[l][l2].get_vec();
      //_vc = vcast(_vr);
      _vc[0]=square_grid[l][l2].x();
      _vc[1]=square_grid[l][l2].y();
      _vc[2]=square_grid[l][l2].z();
      if(_vc[0] >= _nor_lim[0][0] && _vc[0] <= _nor_lim[0][1] && _vc[1] >= _nor_lim[1][0] && _vc[1] <= _nor_lim[1][1])
        _mcut.add_row(_vc);
      //cout<<"v="<<_vc;
        //else
      //cout<<_nor_lim;
        //cout<<square_grid[l][l2];
    }
  }
  //cout<<"mesh="<<_mcut;
  gm_uint _num_rows = _mcut.rows();
  gm_uint _num_cols = _mcut.cols();
  gm_real _dup_tol = 0.001;
  for(gm_uint _i=0;_i<_num_rows-1;_i++){
    for(gm_uint _j=_i+1;_j<_num_rows;_j++){
      if(fabs(_mcut[_j][0] - _mcut[_i][0])<_dup_tol){
        if(fabs(_mcut[_j][1] - _mcut[_i][1])<_dup_tol){
          _mcut[_j] = _mcut[_num_rows-1];
          _num_rows--;
          _j--;
        }
      }
    }
  }
  _mcut.redim(_num_rows,_num_cols);
  _mcut = mscalexy(_mcut,_m_grid_limits,_ss);
  _mcut.write_file(_f);
}

bool fl_contour::eval_cell_curvature_(Geom_Edges<gm_real,4,3>& grid_cell){
    return ( fabs( fabs(grid_cell.z(0))+fabs(grid_cell.z(2)) ) > fabs( fabs(grid_cell.z(1))+fabs(grid_cell.z(3))) );
}

// overall inverse distance interpolation
void fl_contour::inverse_distance_interp(void){
  gm_uint i, j, k;
  gm_real hfi, sh, hi, _dx, _dy, _i_hi;
  square_grid.clear();
  _triangle_mesh.clear();
  _type_triangle_list.clear();
  Geom_Vertex<gm_real,2> pxy;
  TMatrix<gm_real> f(xcells+1,ycells+1);
  for (j=0;j<=ycells;j++){
    pxy.y() = ymin + (gm_real) j * _delta_y;
    for (i=0;i<=xcells;i++){
      pxy.x() = xmin + (gm_real) i * _delta_x;
      sh = 0;
      hfi = 0;
      for(k=1;k<=_number_vertices;k++){
        _dx = delaunay_vertex[k].x()-pxy.x();
        _dy = delaunay_vertex[k].y()-pxy.y();
        hi = sqrt((_dx*_dx)+(_dy*_dy));
        if( hi > GM_EPSILON ){
         _i_hi = 1.0/hi;
         sh  += _i_hi;
         hfi += _i_hi* (delaunay_vertex[k].z()-0.5);
        }
      }
      f[i][j]= hfi/sh;
    }
  }
  zmin = f.get_min();
  zmax = f.get_max();
  //initialize_levels_();
#ifdef _DEBUG
  printf("Start\n");
  printf("zMin = %f\n",zmin);
  printf("zMax = %f\n",zmax);
#endif
  for (j=0;j<ycells;j++){
    //printf("j%i\n",j);
    _vsquare[0].y() = ymin + j * _delta_y;
    _vsquare[1].y() = _vsquare[0].y();
    _vsquare[2].y() = _vsquare[0].y() + _delta_y;
    _vsquare[3].y() = _vsquare[2].y();
    for (i=0;i<xcells;i++){
      //printf("i%i\n",i);
      _vsquare[0].x() = xmin + i * _delta_x;
      _vsquare[1].x() = _vsquare[0].x() + _delta_x;
      _vsquare[2].x() = _vsquare[1].x();
      _vsquare[3].x() = _vsquare[0].x();

      _vsquare[0].z() = f[i][j];
      _vsquare[1].z() = f[i+1][j];
      _vsquare[2].z() = f[i+1][j+1];
      _vsquare[3].z() = f[i][j+1];
      square_grid.add(_vsquare);
    }
  }
#ifdef _SHOW_MESSAGE_
  printf("Number of vertex: %i  \n",square_grid.size());
  printf("InvDistance\n");
#endif
}

// nearest weighted interpolation
void fl_contour::nearest_neighbor_interp(void){
    gm_uint i, j, k, _tp, corner;
    real sh, _vhi, hfi, _vx, _vy;
    square_grid.clear();
    _triangle_mesh.clear();
    _type_triangle_list.clear();
#ifdef _DEBUG
    printf("Performing nearest Interpolation\n");
#endif
    zmin = _m_input_data.get_col(2).get_min()-0.5;
    zmax = _m_input_data.get_col(2).get_max()-0.5;
    //dz = (zmax - zmin) / _lvls;
#ifdef _DEBUG
    printf("zMin = %f\n",zmin);
    printf("zMax = %f\n",zmax);
    printf("_lvls = %i\n",_lvls);
    //printf("dz = %f\n",dz);
#endif
    for(j=0;j<=ycells;j++){
        _vsquare[0].y() = (ymin+j*_delta_y);
        _vsquare[1].y() = (_vsquare[0].y());
        _vsquare[2].y() = (_vsquare[0].y()+_delta_y);
        _vsquare[3].y() = (_vsquare[2].y());
        for(i=0;i<=xcells;i++){
          _vsquare[0].x() = (xmin+i*_delta_x);
          _vsquare[1].x() = (_vsquare[0].x()+_delta_x);
          _vsquare[2].x() = (_vsquare[1].x());
      _vsquare[3].x() = (_vsquare[0].x());
          for(corner=0;corner<4;corner++){
                _found[corner] = false;
    for(k=0;k<_number_triangles;k++){
        _p1 = delaunay_vertex[_vtriangle[k].x()];
        _p2 = delaunay_vertex[_vtriangle[k].y()];
        _p3 = delaunay_vertex[_vtriangle[k].z()];
        sh=0;
        hfi=0;
        if(is_point_in_triangle_(_vsquare[corner].x(), _vsquare[corner].y())){
      _vx=_vsquare[corner].x()-_p1.x();
      _vy=_vsquare[corner].y()-_p1.y();
      _vhi = sqrt(_vx*_vx+_vy*_vy);
      if( _vhi > GM_EPSILON ){
          sh += 1.0/_vhi;
          hfi+= _p1.z()/_vhi;
      }
      _vx=_vsquare[corner].x()-_p2.x();
      _vy=_vsquare[corner].y()-_p2.y();
      _vhi = sqrt(_vx*_vx+_vy*_vy);
      if( _vhi > GM_EPSILON ){
          sh += 1.0/_vhi;
          hfi+= _p2.z()/_vhi;
      }
      _vx=_vsquare[corner].x()-_p3.x();
      _vy=_vsquare[corner].y()-_p3.y();
      _vhi = sqrt(_vx*_vx+_vy*_vy);
      if( _vhi > GM_EPSILON ){
          sh += 1.0/_vhi;
          hfi+= _p3.z()/_vhi;
      }
      _vsquare[corner].z()=(hfi/sh)-0.5;
      zmin = mini(zmin,_vsquare[corner].z());
      zmax = maxi(zmax,_vsquare[corner].z());
      _found[corner]=true;
      break;
        }
    }
      }
      if(_found[0] && _found[1] && _found[2] && _found[3]){
    square_grid.add(_vsquare);
      }else if( (_tp=is_triange_in_triangle_()) ){
    set_2d_triangle_(_vsquare, _tp);
    add_3d_triangle_mesh_();
    _type_triangle_list.push_back(_tp);
      }
  }
    }
#ifdef _SHOW_MESSAGE_
    printf("Number of vertex: %i  \n",square_grid.size());
    printf("Lineal()\n");
#endif
}

// linear interpolation
void fl_contour::lineal_interpolation(void){
  gm_uint i, j, k, _tp, corner;
  square_grid.clear();
  _triangle_mesh.clear();
  _type_triangle_list.clear();
  zmax = 0;
  zmin = 0;
#ifdef _DEBUG
  printf("Lineal Interpolation\n");
#endif
  zmin = _m_input_data.get_col(2).get_min()-0.5;
  zmax = _m_input_data.get_col(2).get_max()-0.5;
#ifdef _DEBUG
  printf("zMin = %f\n",zmin);
  printf("zMax = %f\n",zmax);
  printf("_lvls = %i\n",_lvls);
  //printf("dz = %f\n",dz);
#endif
  for(j=0;j<=ycells;j++){
    _vsquare[0].y() = (ymin+j*_delta_y);
    _vsquare[1].y() = (_vsquare[0].y());
    _vsquare[2].y() = (_vsquare[0].y()+_delta_y);
    _vsquare[3].y() = (_vsquare[2].y());
    for(i=0;i<=xcells;i++){
      _vsquare[0].x() = (xmin+i*_delta_x);
      _vsquare[1].x() = (_vsquare[0].x()+_delta_x);
      _vsquare[2].x() = (_vsquare[1].x());
      _vsquare[3].x() = (_vsquare[0].x());
      for(corner=0;corner<4;corner++){
        _found[corner] = false;
        for(k=0;k<_number_triangles;k++){
          _p1 = delaunay_vertex[_vtriangle[k].x()];
          _p2 = delaunay_vertex[_vtriangle[k].y()];
          _p3 = delaunay_vertex[_vtriangle[k].z()];
          if(is_point_in_triangle_(_vsquare[corner].x(), _vsquare[corner].y())){
            _vsquare[corner].z()=eval_plane_point_(_vsquare[corner].x(),_vsquare[corner].y())-0.5;
            zmin = mini(zmin,_vsquare[corner].z());
            zmax = maxi(zmax,_vsquare[corner].z());
            _found[corner]=true;
            break;
          }
        }
      }
      if(_found[0] && _found[1] && _found[2] && _found[3])
        square_grid.add(_vsquare);
      else if( (_tp=is_triange_in_triangle_()) ){
        set_2d_triangle_(_vsquare, _tp);
        add_3d_triangle_mesh_();
        _type_triangle_list.push_back(_tp);
      }
    }
  }
#ifdef _SHOW_MESSAGE_
  printf("Number of vertex: %i  \n",square_grid.size());
  printf("Lineal()\n");
#endif
}

//
void fl_contour::eval_big_triangle_(){
  gm_real dx,dy,dmax,xmid,ymid;
  dx = xmax - xmin;
  dy = ymax - ymin;
  dmax = (dx > dy) ? dx : dy;
  xmid = 0.5 * (xmax + xmin);
  ymid = 0.5 * (ymax + ymin);
  //
  delaunay_vertex[_number_vertices + 1].x() = xmid - 2.0 * dmax;
  delaunay_vertex[_number_vertices + 1].y() = ymid - dmax;
  delaunay_vertex[_number_vertices + 1].z() = 0.0;
  delaunay_vertex[_number_vertices + 2].x() = xmid;
  delaunay_vertex[_number_vertices + 2].y() = ymid + 2.0 * dmax;
  delaunay_vertex[_number_vertices + 2].z() = 0.0;
  delaunay_vertex[_number_vertices + 3].x() = xmid + 2.0 * dmax;
  delaunay_vertex[_number_vertices + 3].y() = ymid - dmax;
  delaunay_vertex[_number_vertices + 3].z() = 0.0;

#ifdef _SHOW_MESSAGE_
  printf("eval_big_triangle_\n");
  printf("(%f,%f,%f)\n",delaunay_vertex[_number_vertices + 1].x(), delaunay_vertex[_number_vertices + 1].y(), delaunay_vertex[_number_vertices + 1].z());
  printf("(%f,%f,%f)\n",delaunay_vertex[_number_vertices + 2].x(), delaunay_vertex[_number_vertices + 2].y(), delaunay_vertex[_number_vertices + 2].z());
  printf("(%f,%f,%f)\n",delaunay_vertex[_number_vertices + 3].x(), delaunay_vertex[_number_vertices + 3].y(), delaunay_vertex[_number_vertices + 3].z());
#endif
}

// incremental Delaunay list of triangles in clockwise order.
void fl_contour::delaunay_xz_(void){
  vector<gm_uint> _valid;
  gm_uint nvertex, x, y, z;
  gm_uint  i,j,k;
  bool inside;
  //
  Geom_Triangle<gm_real> tt;
  Geom_Vertex<gm_real,3> _p1, _p2, _p3, _point;
  Geom_Vertex<gm_short,2> sp1, sp2, sp3, pnull;
  Geom_Edges<gm_short ,0 ,2> vertex;
  Geom_Edges<gm_real,3,2> p;
  //
  vertex.resize(EMAX+1000);
  pnull.point2d(0,0);
  _valid.resize(trimax);
  eval_big_triangle_();
  _vtriangle[0].x() = _number_vertices + 1;
  _vtriangle[0].y() = _number_vertices + 2;
  _vtriangle[0].z() = _number_vertices + 3;
  _valid[0] = false;
  _number_triangles = 1;
  // insert each point at due time into the triangle mesh
  for(i=1;i<=_number_vertices;i++){
      _point.pointnd(delaunay_vertex[i]);
      nvertex = 0;
      j = 0;
      do{
    if (!_valid[j]){
  x = _vtriangle[j].x();
  y = _vtriangle[j].y();
  z = _vtriangle[j].z();
  // discrete triangle coordinates into discrete triangle array
  sp1.point2d(x,y); sp2.point2d(y,z); sp3.point2d(z,x);
  // real triangle coordinates into gm_real coordinates array
  _p1=delaunay_vertex[x]; _p2=delaunay_vertex[y]; _p3=delaunay_vertex[z];
  tt.data(_p1, _p2, _p3);
  // lies _point inside tt circle.
  inside = eval_circumcircle_( tt, _point);
  if(_cc.x() + _radio < _point.x())
    _valid[j] = true;
  //if(_cc.y() + _radio < _point.y())
    //_valid[j] = true;
  if (inside>0){
      // check if we haven't exceeded the vertex list size
      if (nvertex+3 > EMAX) exit(0);
      vertex.add3(nvertex,sp1,sp2,sp3);
      //
      nvertex += 3;
      _vtriangle[j] = _vtriangle[_number_triangles-1];
      _valid[j] = _valid[_number_triangles-1];
      j--;
      _number_triangles--;
  }
    }
      j++;
  }while(j < _number_triangles);

  for(j=0;j<nvertex-1;j++){
      for(k=j+1;k<nvertex;k++){
    if(vertex.comp_xy(j,k) && vertex.comp_xy(k,j)){
        vertex.null2d(j);
        vertex.null2d(k);
    }
    if(vertex.comp_xx(j,k) && vertex.comp_yy(k,j)){
        vertex.null2d(j);
        vertex.null2d(k);
    }

      }
  }
  for(j=0;j<nvertex;j++){
      if(vertex.nonull(j)){
    if(_number_triangles > trimax) exit(0);
    _number_triangles++;
    _vtriangle[_number_triangles-1].point3d(vertex[j].x(),vertex[j].y(),i);
    _valid[_number_triangles-1] = false;
      }
  }
    }
    i = 0;
    do{
  if((_vtriangle[i].x() > _number_vertices) || (_vtriangle[i].y() > _number_vertices) || (_vtriangle[i].z() > _number_vertices)){
      _vtriangle[i] = _vtriangle[_number_triangles-1];
      i--;
      _number_triangles--;
  }
  i++;
    }while (i < _number_triangles);
#ifdef _SHOW_MESSAGE_
  printf("Vertex = %i, triangles = %i\n",delaunay_vertex.size(),_number_triangles);
  printf("Delaunay... bye bye\n");
#endif
}

bool fl_contour::eval_circumcircle_(Geom_Triangle<gm_real>& o, Geom_Vertex<gm_real,3>& p){
    gm_real dx, dy, rsqr, drsqr;
    // look if these points are not collinear
    if (o.collinear()){
      printf("colinear points\n");
      return false;
    }
    eval_centroide_(o);
    rsqr = eval_square_radius_(o);
    _radio = sqrt(rsqr);
    dx = p.x() - _cc.x();
    dy = p.y() - _cc.y();
    drsqr = dx*dx + dy*dy;
    return (drsqr <= rsqr) ? true : false;
}

void fl_contour::eval_centroide_(Geom_Triangle<gm_real>& _o){
    gm_real  m1,m2,mx1,mx2,my1,my2;
    if (fabs(_o.y(1)-_o.y(0)) <= GM_EPSILON){
        m2    =  -(_o.x(2)-_o.x(1)) / (_o.y(2)-_o.y(1));
        mx2   = 0.5 * (_o.x(1) + _o.x(2));
        my2   = 0.5 * (_o.y(1) + _o.y(2));
        _cc[0] = 0.5 * (_o.x(1) + _o.x(0));
        _cc[1] = m2 * (_cc[0] - mx2) + my2;
    }else if (fabs(_o.y(2)-_o.y(1)) <= GM_EPSILON){
        m1    =  -(_o.x(1)-_o.x(0)) / (_o.y(1)-_o.y(0));
        mx1   = 0.5 * (_o.x(0) + _o.x(1));
        my1   = 0.5 * (_o.y(0) + _o.y(1));
        _cc[0] = 0.5 * (_o.x(2) + _o.x(1));
        _cc[1] = m1 * (_cc[0] - mx1) + my1;
    }else{
        m1    =  -(_o.x(1)-_o.x(0)) / (_o.y(1)-_o.y(0));
        m2    =  -(_o.x(2)-_o.x(1)) / (_o.y(2)-_o.y(1));
        mx1   = 0.5 * (_o.x(0) + _o.x(1));
        mx2   = 0.5 * (_o.x(1) + _o.x(2));
        my1   = 0.5 * (_o.y(0) + _o.y(1));
        my2   = 0.5 * (_o.y(1) + _o.y(2));
        _cc[0] = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
        _cc[1] = m1 * (_cc[0] - mx1) + my1;
    }
    // Sat Nov 26 11:03:31 MST 2011: bug caught by Long To
    //_cc.z() = 0;
}

gm_real fl_contour::eval_square_radius_(Geom_Triangle<gm_real>& _o){
    gm_real dx, dy;
    dx = _o.x(1) - _cc.x();
    dy = _o.y(1) - _cc.y();
    return (dx*dx+dy*dy);
}

// return the type of intersects of a triangle with a plane.
gm_uint fl_contour::triangle_type_(gm_real lb){
    for(gm_uint j=0; j<3; j++)
  if((_Trg[j].z()-zmin) > lb)
      l[j]=2;
  else if((_Trg[j].z()-zmin) < lb)
      l[j]=0;
  else
      l[j]=1;
    return typet[l[0]][l[1]][l[2]];
}

// returns true if a point lies inside a triangle
gm_bool fl_contour::is_point_in_triangle_(gm_real xp,gm_real yp){
    int side1,side2,side3;
    side1 = which_side_(xp,yp,_p1.x(),_p1.y(),_p2.x(),_p2.y());
    side2 = which_side_(xp,yp,_p2.x(),_p2.y(),_p3.x(),_p3.y());
    side3 = which_side_(xp,yp,_p3.x(),_p3.y(),_p1.x(),_p1.y());
    if (side1 == 0 && side2 == 0)
        return(true);
    if (side1 == 0 && side3 == 0)
        return(true);
    if (side2 == 0 && side3 == 0)
        return(true);
    if (side1 == 0 && (side2 == side3))
        return(true);
    if (side2 == 0 && (side1 == side3))
        return(true);
    if (side3 == 0 && (side1 == side2))
        return(true);
    if ((side1 == side2) && (side2 == side3))
        return(true);
    return(false);
}

// returns true if an small traingle lies inside a big triangle
gm_uint fl_contour::is_triange_in_triangle_(void){
    gm_uint _t;
    if(_found[0] && _found[1] && _found[2])      _t = 1;
    else if(_found[2] && _found[3] && _found[0]) _t = 2;
    else if(_found[1] && _found[2] && _found[3]) _t = 3;
    else if(_found[0] && _found[1] && _found[3]) _t = 4;
    else _t = 0;
    return _t;
}

gm_real fl_contour::eval_plane_point_(gm_real xp,gm_real yp){
    gm_real a,b,c,d;
    a = _p1.y() * (_p2.z() - _p3.z()) + _p2.y() * (_p3.z() - _p1.z()) + _p3.y() * (_p1.z() - _p2.z());
    b = _p1.z() * (_p2.x() - _p3.x()) + _p2.z() * (_p3.x() - _p1.x()) + _p3.z() * (_p1.x() - _p2.x());
    c = _p1.x() * (_p2.y() - _p3.y()) + _p2.x() * (_p3.y() - _p1.y()) + _p3.x() * (_p1.y() - _p2.y());
    d = - _p1.x() * (_p2.y() * _p3.z() - _p3.y() * _p2.z()) - _p2.x() * (_p3.y() * _p1.z() - _p1.y() * _p3.z()) - _p3.x() * (_p1.y() * _p2.z() - _p2.y() * _p1.z());

    if (fabs(c) > GM_EPSILON)
        return(- (a * xp + b * yp + d) / c);
    else
        return(0.0);
}

// which side of a line a point lies
int fl_contour::which_side_(gm_real x,gm_real y,gm_real x1,gm_real y1,gm_real x2,gm_real y2){
    gm_real dist;
    dist = (y - y1) * (x2 - x1) - (y2 - y1) * (x - x1);
    if (dist > 0)
        return(-1);
    else if (dist < 0)
        return(1);
    else
  return(0);
}

// interpolation Subrutines
void fl_contour::set_2d_triangle_(Geom_Edges<gm_real, 4, 3>& p, int _tp){
    _Trg[0]=p[sidet[_tp-1][0]];
    _Trg[1]=p[sidet[_tp-1][1]];
    _Trg[2]=p[sidet[_tp-1][2]];
}

void fl_contour::line_interpolation_(gm_real lv, gm_uint c){
    if (c!=0) {
  switch (c) {
      // case 1 - line between vertices 1 and 2
      case 1:
    _line[0].x()=_Trg[0].x();
    _line[0].y()=_Trg[0].y();
    _line[1].x()=_Trg[1].x();
    _line[1].y()=_Trg[1].y();
    break;
      // case 2 - line between vertices 2 and 3
      case 2:
    _line[0].x()=_Trg[1].x();
    _line[0].y()=_Trg[1].y();
    _line[1].x()=_Trg[2].x();
    _line[1].y()=_Trg[2].y();
    break;
      // case 3 - line between vertices 3 and 1
      case 3:
    _line[0].x()=_Trg[2].x();
    _line[0].y()=_Trg[2].y();
    _line[1].x()=_Trg[0].z();
    _line[1].y()=_Trg[0].y();
    break;
      // case 4 - line between vertex 1 and side 2-3
      case 4:
    _line[0].x()=_Trg[0].x();
    _line[0].y()=_Trg[0].y();
    _line[1]=linear_interpolation_(0,1,lv,zmin);
    break;
      // case 5 - line between vertex 2 and side 3-1
      case 5:
    _line[0].x()=_Trg[1].x();
    _line[0].y()=_Trg[1].y();
    _line[1]=linear_interpolation_(0,2,lv,zmin);
    break;
      // case 6 - line between vertex 3 and side 1-2
      case 6:
    _line[0].x()=_Trg[0].z();
    _line[0].y()=_Trg[0].y();
    _line[1]=linear_interpolation_(0,2,lv,zmin);
    break;
      // case 7 - line between sides 1-2 and 2-3
      case 7:
    _line[0]=linear_interpolation_(0,2,lv,zmin);
    _line[1]=linear_interpolation_(0,1,lv,zmin);
     break;
       // case 8 - line between sides 2-3 and 3-1
      case 8:
    _line[0]=linear_interpolation_(1,2,lv,zmin);
    _line[1]=linear_interpolation_(0,2,lv,zmin);
    break;
      // case 9 - line between sides 3-1 and 1-2
      case 9:
    _line[0]=linear_interpolation_(0,1,lv,zmin);
    _line[1]=linear_interpolation_(2,1,lv,zmin);
    break;
      default:
    break;
  }
    }
}

Geom_Vertex<gm_real,2> fl_contour::linear_interpolation_(gm_int _i, gm_int _j, gm_real _lv, gm_real _zmn){
  Geom_Vertex<gm_real,2> _vtx;
  register gm_real _z = (_lv-(_Trg[_i].z()-zmin))/(_Trg[_j].z()-_Trg[_i].z());
  _vtx.x()=_Trg[_i].x()+(_Trg[_j].x()-_Trg[_i].x())*_z;
  _vtx.y()=_Trg[_i].y()+(_Trg[_j].y()-_Trg[_i].y())*_z;
  return _vtx;
}

//
void fl_contour::add_2d_line_(Geom_Edges<gm_real, 2, 2>& _l, gm_real _lv){
  line_2d_contour.push_back(_lv);
}

void fl_contour::add_3d_line_(Geom_Edges<gm_real, 2, 2>& _l, gm_real _lv){
  Geom_Edges<gm_real,2,3>  line_3d;
  line_3d[0].x() = _line[0].x();
  line_3d[0].y() = _line[0].y();
  line_3d[0].z() = (_lv);
  line_3d[1].x() = _line[1].x();
  line_3d[1].y() = _line[1].y();
  line_3d[1].z() = (_lv);
  line_2d_contour.push_back(_lv);
  line_3d_contour.add(line_3d);
}

void fl_contour::add_3d_triangle_(void){
  triangle_color_mesh.push_back(_Trg);
}

void fl_contour::add_3d_triangle_mesh_(void){
  _triangle_mesh.push_back(_Trg);
}

void fl_contour::set_palette(gm_uint u){
  _cpalette.set_color(u);
}

void fl_contour::update_palette_(void){
  dz = (zmax - zmin) / _lvls;
  _cpalette.initialize(zmin,zmax,_lvls);
  _cpalette.update_palette_();
}

// submesh cutting utils
void fl_contour::set_submesh_limits(int _i, int _j, real _r){ 
    _cut_lim[_i][_j]=_r;
    //cout<<"Cut mesh ="<<_cut_lim;  }; 
}

void fl_contour::normalize_submesh_limits(void){
  _nor_lim = _cut_lim;
  _nor_lim.transpose();
  _nor_lim = mnormxy(_nor_lim,_m_grid_limits,_ss);
  _nor_lim.transpose();
  _nor_lim.redim(2,2);
  //cout<<"Normalized Limits ="<<_nor_lim;
}

gm_real fl_contour::get_submesh_limits(int _i, int _j){
  return _cut_lim[_i][_j];
};

TMatrix<gm_real> fl_contour::get_normalize_submesh_limits(void){
  return _nor_lim;
}

/// RUBISH


