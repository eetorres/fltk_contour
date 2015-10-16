// Latest modification 21/03/2008
// ========================================================================
// version: 0.1
// name:    fl_palette.cxx
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


#include <fl_palette.h>

fl_palette::fl_palette(void){
    _ipalette = 0;
}

void fl_palette::set(double w){
    xdelta1 = w/w0;
    xdelta2 = xdelta1+xdelta1;
    x0=xdelta1;
    x1=xdelta1+x0;
    x2=xdelta1+x1;
    x3=xdelta1+x2;
    x4=xdelta1+x3;
    x5=xdelta1+x4;
    ydelta1 = 255.0;
}

void fl_palette::set_color(unsigned int u){ 
    _ipalette = u;
    switch(_ipalette){
	case 0: 
	    w0 = 1.0;//Linear
	    break;
	default:
	    w0 = 6.0;//HSV
	    break;
    }
};

void fl_palette::update_palette_(void){
    double lv;
    _color_palette.clear();
#ifdef _SHOW_MESSAGE_
    printf("_lvls = %i\n",_lvls);
    printf("dz = %f\n",dz);
#endif
    dz=(zmax-zmin)/_lvls;
    set(zmax-zmin);
    for(unsigned int _p=0; _p<(_lvls+1); _p++){
	lv = (_p*dz);
	rgb = palette_selection_(lv);
	rgb.r/=255.0;
	rgb.g/=255.0;
        rgb.b/=255.0;
	//printf("r=%f, g=%f, b=%f\n",rgb.r,rgb.g,rgb.b);
	_color_palette.push_back(rgb);
    }
}

gm_rgb fl_palette::palette_selection_(double _x){
    switch(_ipalette){
	case 0: 
	    return linear_palette_(_x);
	    break;
	case 1: 
	    return terrain_palette_(_x);
	    break;
	case 2:  
	    return earth_palette_(_x);
	    break;
	case 3: 
	    return hsv_palette_(_x);
	    break;
	default: 
	    return rgb_palette_(_x);
	    break;
    }
}

gm_rgb fl_palette::get_color(double val){
    unsigned int c = (unsigned int)(val/dz);
    return _color_palette[c];
}

double fl_palette::color_interpolation_(double dx, double dy, double _cl){
    return (dy/dx)*_cl;
}

//// Palettes

// HSV palette
gm_rgb fl_palette::hsv_palette_(double val){
    double rcf1=100.0;
    double bcf2=200.0;
    double gcf3=255.0;
    ydelta1 = 255.0;
    if(val >= 0 && val <= x0){ //1
	rgb.r=rcf1;
	rgb.g=0.0;
	rgb.b=color_interpolation_(xdelta1, bcf2, val);
    }else if(val > x0 && val <= x1){ //2
	rgb.r=rcf1-color_interpolation_(xdelta1, rcf1, val-x0);
	rgb.g=0.0;
	rgb.b=bcf2;
    }else if(val > x1 && val <= x2){ // 3
	rgb.r=0.0;
	rgb.g=color_interpolation_(xdelta1, gcf3, val-x1);
	rgb.b=bcf2;
    }else if(val > x2 && val <= x3){ // 4
	rgb.r=0.0;
	rgb.g=gcf3;
	rgb.b=bcf2-color_interpolation_(xdelta1, bcf2, val-x2);
    }else if(val > x3 && val <= x4){ // 5
	rgb.r=color_interpolation_(xdelta1, rcf1, val-x3);
	rgb.g=gcf3;
	rgb.b=0.0;
    }else if(val > x4 && val <= x5){ // 6
	rgb.r=rcf1;
	rgb.g=gcf3-color_interpolation_(xdelta1, gcf3, val-x4);
	rgb.b=0.0;
    }
    return rgb;
}

// Earth palette
gm_rgb fl_palette::earth_palette_(double val){
    cf1.r=250.0; cf1.g=250.0;  cf1.b=255.0;
    ydelta1 = 255.0;
    if(val >= 0 && val <= x0){ // 1
	rgb.r=0.0;
	rgb.g=0.0;
	rgb.b=color_interpolation_(xdelta1, cf1.b, val);
    }else if(val >= x0 && val <= x1){ // 2
	rgb.r=0.0;
	rgb.g=color_interpolation_(xdelta1, cf1.g, val-x0);
	rgb.b=0.0;
    }else if(val > x1 && val <= x2){ //3
	rgb.r=0.0;
	rgb.g=cf1.g;
	rgb.b=0.0;
    }else if(val > x2 && val <= x3){ //4
	rgb.r=color_interpolation_(xdelta1, cf1.r, val-x2);
	rgb.g=cf1.g;
	rgb.b=0.0;
    }else if(val > x3 && val <= x4){ //5
	rgb.r=0.5*cf1.r+0.25*color_interpolation_(xdelta1, cf1.r, val-x3);
	rgb.g=0.5*cf1.r+0.25*color_interpolation_(xdelta1, cf1.g, val-x3);
	rgb.b=0.5*cf1.r+0.25*color_interpolation_(xdelta1, cf1.b, val-x3);
    }else if(val > x4){ //6
	rgb.r=0.75*cf1.r+0.25*color_interpolation_(xdelta1, cf1.r, val-x4);
	rgb.g=0.75*cf1.g+0.25*color_interpolation_(xdelta1, cf1.g, val-x4);
	rgb.b=0.75*cf1.b+0.25*color_interpolation_(xdelta1, cf1.b, val-x4);
    }
    return rgb;
}

// Terrain palette
gm_rgb fl_palette::terrain_palette_(double val){
    cf1.r=250.0; cf1.g=250.0;  cf1.b=255.0;
    ydelta1 = 255.0;
    if(val >= 0 && val <= x0){ // 1
	rgb.r=color_interpolation_(xdelta1, 184.0, val);
	rgb.g=color_interpolation_(xdelta1, 134.0, val);
	rgb.b=color_interpolation_(xdelta1, 11.2, val);
    }else if(val >= x0 && val <= x1){ // 2
	rgb.r=0.0;
	rgb.g=0.6*color_interpolation_(xdelta1, cf1.g, val-x0);
	rgb.b=0.0;
    }else if(val > x1 && val <= x2){ //3
	rgb.r=0.0;
	rgb.g=0.6*cf1.g;
	rgb.b=0.0;
    }else if(val > x2 && val <= x3){ //4
	rgb.r=34.0+0.3*color_interpolation_(xdelta1, cf1.r, val-x2);
	rgb.g=0.6*cf1.g;
	rgb.b=34.0+0.3*color_interpolation_(xdelta1, cf1.b, val-x2);;
    }else if(val > x3 && val <= x4){ //5
	rgb.r=0.5*cf1.r+0.25*color_interpolation_(xdelta1, cf1.r, val-x3);
	rgb.g=0.5*cf1.r+0.25*color_interpolation_(xdelta1, cf1.g, val-x3);
	rgb.b=0.5*cf1.r+0.25*color_interpolation_(xdelta1, cf1.b, val-x3);
    }else if(val > x4){ //6
	rgb.r=0.75*cf1.r+0.25*color_interpolation_(xdelta1, cf1.r, val-x4);
	rgb.g=0.75*cf1.g+0.25*color_interpolation_(xdelta1, cf1.g, val-x4);
	rgb.b=0.75*cf1.b+0.25*color_interpolation_(xdelta1, cf1.b, val-x4);
    }
    return rgb;
}

// RGB palette
gm_rgb fl_palette::rgb_palette_(double val){
    cf1.r=250.0; cf1.g=250.0;  cf1.b=255.0;
    ydelta1 = 255.0;
    if(val >= 0 && val <= x0){ // 1
	rgb.r=0.0;
	rgb.g=0.0;
	rgb.b=color_interpolation_(xdelta1, cf1.b, val);
    }else if(val >= x0 && val <= x1){ // 2
	rgb.r=0.0;
	rgb.g=color_interpolation_(xdelta1, cf1.g, val-x0);
	rgb.b=cf1.b;
    }else if(val > x1 && val <= x2){ //3
	rgb.r=0.0;
	rgb.g=cf1.g;
	rgb.b=cf1.b-color_interpolation_(xdelta1, cf1.b, val-x1);
    }else if(val > x2 && val <= x3){ //4
	rgb.r=color_interpolation_(xdelta1, cf1.r, val-x2);
	rgb.g=cf1.g;
	rgb.b=0.0;
    }else if(val > x3 && val <= x4){ //5
	rgb.r=cf1.r;
	rgb.g=cf1.g-color_interpolation_(xdelta1, cf1.g, val-x3);
	rgb.b=0.0;
    }else if(val > x4){ //6
	rgb.r=cf1.r-color_interpolation_(xdelta1, cf1.r, val-x4);
	rgb.g=0.0;
	rgb.b=0.0;
    }
    return rgb;
}

// Linear palette
gm_rgb fl_palette::linear_palette_(double val){
    double rydelta1, gydelta1, bydelta1;
	
    cf1.r=20.0;  cf1.g=20.0;  cf1.b=20.0;
    cf2.r=255.0; cf2.g=255.0; cf2.b=255.0;
    ydelta1 = 255.0;
	
    rydelta1 = cf2.r-cf1.r;
    gydelta1 = cf2.g-cf1.g;
    bydelta1 = cf2.b-cf1.b;
    if(val >= 0 && val <= x0){ // 1
	rgb.r=cf1.r+color_interpolation_(xdelta1, rydelta1, val);
	rgb.g=cf1.g+color_interpolation_(xdelta1, gydelta1, val);
	rgb.b=cf1.b+color_interpolation_(xdelta1, bydelta1, val);
    }
    return rgb;
}



/// RUBISH


