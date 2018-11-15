#version 400
layout(location = 0) in vec3 vertex_position;
layout(location = 1) in vec3 vertex_color;
smooth out vec3 thePosition; // Interpolate position among fragments
uniform mat4 ORTHO;
uniform float iGlobalTime; 
out vec3 color;

vec2 hash( vec2 p )   
{   
   p=vec2(dot(p, vec2(127.1, 311.7)), dot(p, vec2(269.5, 183.3)));   
     
   return fract(sin(p)*18.5453);  
}  

vec2 voronoi( in vec2 x )  
{  
    vec2 n = floor( x );              // cell(n)  
    vec2 f = fract( x );              // 當前像素在cell space的坐標    
   vec3 m = vec3( 8. );               // 影響每個cell的大小，影響背景?色  
   // 訪問相鄰的9個cell  
    for( int j=-1; j<=1; j++ )  
    {  
       for( int i=-1; i<=1; i++ )  
       {  
           vec2  g = vec2( float(i), float(j) );  // 鄰近的 cell id offset  
           // n+g 鄰近的 cell(n+g) 的?機像素坐標 o (cell space)  
           vec2  o = hash( n + g );   // 影響cell的?色               
           // 
           vec2  r = g - f + (0.5+0.5*sin(iGlobalTime+6.2831*o));  
           //vec2  r = g - f + o;     // cell(n+g)[o] - cell(n)[f]   
           // 
           float d = dot( r, r );  
           // 保存更小的d  
           if( d<m.x )  
           {  
              m = vec3( d, o );  
           }  
       }  
    }  
    return vec2( sqrt(m.x), m.y+m.z );  
}  

void main () {

		color = vertex_color;
		gl_Position = ORTHO*vec4 (vertex_position, 1.0);
		thePosition = vertex_position; 

		//vec2 iResolution = vec2(400.0, 400.0);
		//vec2 p = thePosition.xy/max(iResolution.x,iResolution.y);  
	    // computer voronoi patterm  
		//vec2 c = voronoi( (14.0+6.0*sin(0.2*iGlobalTime))*p );   
		 //float  n=(c.y+c.x)/2;
		//float h = clamp(n/1.0,0.0,1.0);
		//thePosition.z+=h;
		
}
