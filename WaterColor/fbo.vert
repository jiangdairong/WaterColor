#version 400
layout(location = 0) in vec3 vertex_position;
layout(location = 1) in vec3 vertex_position1;
out vec2 UV;

        
void main () {
	gl_Position = vec4 (vertex_position, 1.0);
	
	UV = ( vertex_position.xy+vec2(1,1))/2.0;
	
	
}
