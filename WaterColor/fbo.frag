#version 400
in vec2 UV;
out vec3 color;

uniform sampler2D renderedTexture;
uniform sampler2D renderedTexture1;
uniform sampler2D  renderedTexture2;
uniform sampler2D renderedTexture3;



void main(){
   color = texture( renderedTexture, UV).xyz ;
   color= texture(renderedTexture1 , UV).xyz ;

}
