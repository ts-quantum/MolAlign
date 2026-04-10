#include "colors.inc"
#include "textures.inc"
#include "shapes3.inc"
#include "stones1.inc"

global_settings{ assumed_gamma 1.0 }
#default{ finish{ ambient 0.1 diffuse 0.9 }}

#declare phi=0.5;
#declare h=5;
#declare r=15;
#declare Camera_1 = camera { /*ultra_wide_angle*/ angle 55  
							location < cos(phi)*r,h ,sin(phi)*r>
							right x*image_width/image_height
							look_at < 1.5,2.0,1.0 > }

camera{Camera_1}

  
light_source { <1500, 2500, 0> color rgb <1,1,1> }
  
plane{<0,1,0>, 1 hollow  // 
      
        texture{ pigment {color rgb<0.1,0.3,0.75>*0.7}
                 #if (version = 3.7 )  finish {emission 1 diffuse 0}
                 #else                 finish { ambient 1 diffuse 0}
                 #end 
               } // end texture 1

        texture{ pigment{ bozo turbulence 0.75
                          octaves 6  omega 0.7 lambda 2 
                          color_map {
                          [0.0  color rgb <0.95, 0.95, 0.95> ]
                          [0.05  color rgb <1, 1, 1>*1.25 ]
                          [0.15 color rgb <0.85, 0.85, 0.85> ]
                          [0.55 color rgbt <1, 1, 1, 1>*1 ]
                          [1.0 color rgbt <1, 1, 1, 1>*1 ]
                          } // end color_map 
                         translate< 3, 0,-1>
                         scale <0.3, 0.4, 0.2>*3
                        } // end pigment
                 #if (version = 3.7 )  finish {emission 1 diffuse 0}
                 #else                 finish { ambient 1 diffuse 0}
                 #end 
               } // end texture 2
       scale 10000 
     } //-------------------------------------------------------------
 
// ground fog at the horizon -----------------------------------------
fog{ fog_type   2
     distance   1000
     color      rgb<1,1,1>*0.9
     fog_offset 0.1
     fog_alt    20
     turbulence 1.8
   } //---------------------------------------------------------------
 
// ground ------------------------------------------------------------
  
plane{ <0,1,0>, 0 
    texture{ Silver_Metal 
              } // end of texture
     } // end of plane

// coordinate system
#declare axis = union {
cylinder {
	<0,0,0>,
	<1,0,0>,
	0.01
	scale <10,1,1>
	translate <-5,0,0>
}
cone {
	<0,0,0>, 0.05
	<1,0,0>, 0
	translate<5,0,0>
}
}

#declare coordinatesystem = union {
object {axis rotate x*0 } 
 text   { ttf "arial.ttf",  "x",  0.05,  0 
          rotate<0,-90,0> scale 0.75 translate <5,0.2,0.2>}

object {axis rotate z*90 }
 text   { ttf "arial.ttf",  "z",  0.05,  0   // y axis in povray coordinates
          rotate y*-90 scale 0.75 translate <0,5,0.2>}

object {axis rotate y*-90 scale z*1 }
 text   { ttf "arial.ttf",  "y",  0.05,  0   // z axis in povray coordinates 
          rotate y*-90 scale 0.75 translate <-0.2,0.2,5>}
}


//object {coordinatesystem scale 0.7 rotate <0,0,0> translate <0,2,0> 
//    pigment {color NavyBlue} } //draw coordinate system
    
#include "irc.inc"
#declare K=clock;

object { molecule[K]  scale 0.75 rotate x*60 rotate y*-90 translate <5,3,2>}

