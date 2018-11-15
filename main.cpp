#define GLEW_STATIC
#define _USE_MATH_DEFINES 
#include <GL/glew.h> // include GLEW and new version of GL on Windows
#include <GLFW/glfw3.h> // GLFW helper library
#include <GL/glui.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <cstdio>
#include <cmath>
#include <Windows.h>
#include <gl/glut.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
//#include "tbb/tbb.h"

#include "paper.h"
#include "Watercolor.h"
#include "PerlinNoise.h"

using namespace std;
#define HM_SIZE_X 200 // Dimensions of our heightmap
#define HM_SIZE_Y 200

glm::vec3 k;
glm::vec3 s;
glm::vec3 R[HM_SIZE_X*HM_SIZE_Y];
glm::vec3 T[HM_SIZE_X*HM_SIZE_Y];
glm::vec3 R0;
glm::vec3 T0;
glm::vec3 R1;
glm::vec3 T1;
glm::vec3 map_k[HM_SIZE_Y*HM_SIZE_X][2];
glm::vec3 map_s[HM_SIZE_Y*HM_SIZE_X][2];
// Table of pigments  
// from Computer-Generated Watercolor. Cassidy et al. 
// K is absorption. S is scattering 
// a 
#define K_QuinacridoneRose glm::vec3(0.22, 1.47, 0.57) 
#define S_QuinacridoneRose  glm::vec3(0.05, 0.003, 0.03) 
// b 
#define K_IndianRed  glm::vec3(0.46, 1.07, 1.50) 
#define S_IndianRed  glm::vec3(1.28, 0.38, 0.21) 
// c 
#define K_CadmiumYellow  glm::vec3(0.10, 0.36, 3.45) 
#define S_CadmiumYellow  glm::vec3(0.97, 0.65, 0.007) 
// d 
#define K_HookersGreen  glm::vec3(1.62, 0.61, 1.64) 
#define S_HookersGreen  glm::vec3(0.01, 0.012, 0.003) 
// e 
#define K_CeruleanBlue  glm::vec3(1.52, 0.32, 0.25) 
#define S_CeruleanBlue  glm::vec3(0.06, 0.26, 0.40) 
// f 
#define K_BurntUmber  glm::vec3(0.74, 1.54, 2.10) 
#define S_BurntUmber  glm::vec3(0.09, 0.09, 0.004) 
// g 
#define K_CadmiumRed  glm::vec3(0.14, 1.08, 1.68) 
#define S_CadmiumRed  glm::vec3(0.77, 0.015, 0.018) 
// h 
#define K_BrilliantOrange  glm::vec3(0.13, 0.81, 3.45) 
#define S_BrilliantOrange  glm::vec3(0.009, 0.007, 0.01) 
// i 
#define K_HansaYellow  glm::vec3(0.06, 0.21, 1.78) 
#define S_HansaYellow glm::vec3(0.50, 0.88, 0.009) 
// j 
#define K_PhthaloGreen glm::vec3(1.55, 0.47, 0.63) 
#define S_PhthaloGreen glm::vec3(0.01, 0.05, 0.035) 
// k 
#define K_FrenchUltramarine glm::vec3(0.86, 0.86, 0.06) 
#define S_FrenchUltramarine glm::vec3(0.005, 0.005, 0.09) 
// l 
#define K_InterferenceLilac glm::vec3(0.08, 0.11, 0.07) 
#define S_InterferenceLilac glm::vec3(1.25, 0.42, 1.43) 

unsigned int uiVBOHeightmapData; // Here are stored heightmap data (vertices)
unsigned int uiVBOIndices; // And here indices for rendering heightmap

unsigned int uiVAOHeightmap; // One VAO for heightmap
unsigned int colors_vbo = 0;	//宣告頂點顏色vbo

//fbo
GLuint FramebufferName = 0;
GLuint renderedTexture[4];
GLuint depthrenderbuffer;




glm::vec3 vHeightmapData[HM_SIZE_X*HM_SIZE_Y];
float colors[HM_SIZE_X*HM_SIZE_Y*3];
float paint[3];
float mix[3];
/*	fluid capacity*/
float cmax=0.2f;
float cmin=0.0f;
/***************/
float r_bound[HM_SIZE_X*HM_SIZE_Y];
float l_bound[HM_SIZE_X*HM_SIZE_Y];
float u_bound[HM_SIZE_X*HM_SIZE_Y];
float d_bound[HM_SIZE_X*HM_SIZE_Y];
float viscosity = 0.1;
float viscous_drag = 0.01;
float wp[HM_SIZE_X*HM_SIZE_Y]; //water pressure
float newUs[HM_SIZE_X][HM_SIZE_Y][2];
float pigment[HM_SIZE_X*HM_SIZE_Y][12];
float Primepigment[HM_SIZE_X*HM_SIZE_Y][12];
float M_wet[HM_SIZE_X][HM_SIZE_Y];
float M_flow[HM_SIZE_X*HM_SIZE_Y];
float Mask[25];
float epsilon = 0.1;
#define tolerance  0.0001
float layer[HM_SIZE_Y*HM_SIZE_X][2];
float fHeights[HM_SIZE_X*HM_SIZE_Y];
int iIndices[HM_SIZE_Y*(HM_SIZE_X - 1) * 2 + HM_SIZE_X - 1];
int order = 0;
bool flag;

/*	fluid capacity*/
float c[HM_SIZE_X][HM_SIZE_Y];
float weight[11][11] = { 0 };
float blurred[HM_SIZE_X][HM_SIZE_Y];
/***************/
//Paper *g_Paper;
//WaterColor *g_WaterColor;

////




//int load_height_normal()
//{
//	g_Surface->write_bmp_height_("Sand_Height_Map.bmp", g_Surface->m_GridNum, g_Surface->m_GridNum);
//
//	glEnable(GL_TEXTURE_2D);
//	textureheight = auxDIBImageLoadA("Sand_Height_Map.bmp");
//	glBindTexture(GL_TEXTURE_2D, 4);
//	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
//	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, textureheight->sizeX, textureheight->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, textureheight->data);
//
//	g_Surface->write_bmp_height_("Water_In_Sand_Map.bmp", g_Surface->m_GridNum, g_Surface->m_GridNum);
//	glEnable(GL_TEXTURE_2D);
//	textureheight = auxDIBImageLoadA("Water_In_Sand_Map.bmp");
//	glBindTexture(GL_TEXTURE_2D, 10);
//	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
//	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, textureheight->sizeX, textureheight->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, textureheight->data);
//
//	g_Surface->write_bmp_height_("Height_Map.bmp", g_Surface->m_GridNum, g_Surface->m_GridNum);
//
//	glEnable(GL_TEXTURE_2D);
//	textureheight = auxDIBImageLoadA("Height_Map.bmp");
//	glBindTexture(GL_TEXTURE_2D, 9);
//	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
//	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, textureheight->sizeX, textureheight->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, textureheight->data);
//
//	return 1;
//}
void loadFile(const char* filename, std::string &string)
{
	std::ifstream fp(filename);
	if (!fp.is_open()){
		std::cout << "Open <" << filename << "> error." << std::endl;
		return;
	}

	char temp[300];
	while (!fp.eof()){
		fp.getline(temp, 300);
		string += temp;
		string += '\n';
	}

	fp.close();
}
bool OutOfBound(int x,int y)
{
	if (x >= 0 && x < HM_SIZE_X &&y >= 0 && y < HM_SIZE_Y)
		return false;
	else 
		return true;
}

void ColorChanged(GLFWwindow* window, int key, int scancode, int action, int mode)
{
	if (key == GLFW_KEY_A && action == GLFW_PRESS)
	{
		order = 0;
		k = K_QuinacridoneRose;
		s = S_QuinacridoneRose;
	}
	if (key == GLFW_KEY_B && action == GLFW_PRESS)
	{
		order = 1;
		k = K_IndianRed;
		s = S_IndianRed;
	}
	if (key == GLFW_KEY_C && action == GLFW_PRESS)
	{
		order = 2;
		k = K_CadmiumYellow;
		s = S_CadmiumYellow;
	}
	if (key == GLFW_KEY_D && action == GLFW_PRESS)
	{
		order = 3;
		k = K_HookersGreen;
		s = S_HookersGreen;
	}
	if (key == GLFW_KEY_E && action == GLFW_PRESS)
	{
		order = 4;
		k = K_CeruleanBlue;
		s = S_CeruleanBlue;
	}
	if (key == GLFW_KEY_F && action == GLFW_PRESS)
	{
		order = 5;
		k = K_BurntUmber;
		s = S_BurntUmber;
	}
	if (key == GLFW_KEY_G && action == GLFW_PRESS)
	{
		order = 6;
		k = K_CadmiumRed;
		s = S_CadmiumRed;
	}
	if (key == GLFW_KEY_H && action == GLFW_PRESS)
	{
		order = 7;
		k = K_BrilliantOrange;
		s = S_BrilliantOrange;
	}
	if (key == GLFW_KEY_I && action == GLFW_PRESS)
	{
		order = 8;
		k = K_HansaYellow;
		s = S_HansaYellow;
	}
	if (key == GLFW_KEY_J && action == GLFW_PRESS)
	{
		order = 9;
		k = K_PhthaloGreen;
		s = S_PhthaloGreen;
	}
	if (key == GLFW_KEY_K && action == GLFW_PRESS)
	{
		order = 10;
		k = K_FrenchUltramarine;
		s = S_FrenchUltramarine;
	}
	if (key == GLFW_KEY_L && action == GLFW_PRESS)
	{
		order = 11;
		k = K_InterferenceLilac;
		s = S_InterferenceLilac;
	}
}


float rand_FloatRange(float a, float b) {
	return ((b - a)*((float)rand() / RAND_MAX)) + a;
}
float CalVelocity_U(int i,int j)
{
	//float u = (r_bound[j*HM_SIZE_X + i] - l_bound[j*HM_SIZE_X + i]) / 2;
	float u = (r_bound[j*HM_SIZE_X + i] + l_bound[j*HM_SIZE_X + i]) / 2;
	return u;
}
float CalVelocity_V(int i,int j)
{
	float v = (u_bound[j*HM_SIZE_X + i] - d_bound[j*HM_SIZE_X + i]) / 2;
	return v;
}
float thirdStepUV(float A, float B, int i, int j) {
	
	float t;
	try {
		t = r_bound[j*HM_SIZE_X + i] + 1 * (A - viscosity*B + wp[j*HM_SIZE_X + i] - wp[j*HM_SIZE_X + (i + 1)] - viscous_drag*r_bound[j*HM_SIZE_X + i]);
	}
	catch (...) {
		t = r_bound[j*HM_SIZE_X + i] + 1 * (A - viscosity*B + wp[j*HM_SIZE_X + i] - wp[j*HM_SIZE_X + i] - viscous_drag*r_bound[j*HM_SIZE_X + i]);
	}
	return t;
}
float sixthStepUV(float A, float B, int i, int j) {
	// r_bound = vi(j+.5) + timechange(A - mewB + pij - p(i+1)j - kv(i+.5)j;
	float s;
	try {
		s = u_bound[j*HM_SIZE_X + i] + 1 * (A - viscosity*B + wp[j*HM_SIZE_X + i] - wp[(j + 1)*HM_SIZE_X + i] - viscous_drag* u_bound[j*HM_SIZE_X + i]);
	}
	catch (...) {
		s = u_bound[j*HM_SIZE_X + i] + 1 * (A - viscosity*B + wp[j*HM_SIZE_X + i] - wp[j*HM_SIZE_X + i] - viscous_drag* u_bound[j*HM_SIZE_X + i]);
		
	}
	return s;
}
void update(){
	for (int i = 0; i < HM_SIZE_X; i++)
	{
		for (int j = 0; j < HM_SIZE_Y; j++)
		{
			r_bound[j*HM_SIZE_X + i] = newUs[j][i][0];
			if (i + 1 < HM_SIZE_X)
				l_bound[j*HM_SIZE_X + (i + 1)] = 1 * newUs[j][i][0];
			//l_bound[j*HM_SIZE_X + (i + 1)] = -1 * newUs[j][i][0];
			else
			{
				;
			}
			u_bound[j*HM_SIZE_X + i] = newUs[j][i][1];
			if (j + 1 < HM_SIZE_Y)
				d_bound[(j + 1)*HM_SIZE_X + i] = 1 * newUs[j][i][1];
				//d_bound[(j + 1)*HM_SIZE_X + i] = -1 * newUs[j][i][1];
			else
			{
				;
			}
		}
	}
}
void initialize()
{
	for (int i = 0; i < HM_SIZE_X; i++)
	{
		for (int j = 0; j < HM_SIZE_Y; j++)
		{
			r_bound[j*HM_SIZE_X+i] = 0;
			l_bound[j*HM_SIZE_X + i] = 0;
			u_bound[j*HM_SIZE_X + i] = 0;
			d_bound[j*HM_SIZE_X + i] = 0;
			M_wet[j][i] = 0.0f;
			M_flow[j*HM_SIZE_X + i] = 0.0f;
			layer[j*HM_SIZE_X+i][0] = 0.0f;
			layer[j*HM_SIZE_X+i][1] = 0.0f;
			wp[j*HM_SIZE_X + i] = 0.0f;
			for (int n = 0; n < 12; n++)
			{
				pigment[j*HM_SIZE_X + i][n] = 0;
				Primepigment[j*HM_SIZE_X + i][n] = 0;
			}
		}
	}
	//Gaussian Blur mask
	//int s = 3;
	//for (int i = 0; i <5; i++)
	//{
	//	for (int j = 0; j < 5; j++)
	//	{
	//		Mask[j * 5 + i] = exp(-(pow((i - 2), 2) + pow((j - 2), 2)) / (2 * s*s)) / (2 * M_PI*s*s);
	//	}
	//}
	//float sum = 0;
	//for (int i = 0; i < 25; i++)
	//{
	//	sum += Mask[i];
	//}
	//for (int i = 0; i < 25; i++)
	//{
	//	Mask[i] = Mask[i] / sum;
	//}

	/******GAUSSIAN***********/
	/*
	int mradius = 5;

	
	float sigma = 1.5;
	
	for (int i = -5; i <= mradius; i++)
	{
		for (int j = -5; j <= mradius; j++)
		{
			weight[i + 5][j + 5] = (1 / (2 * PI*pow(sigma, 2)))*(exp(-(i*i + j*j) / (2 * pow(sigma, 2))));
			cout << weight[i + 5][j + 5] << "   ";
		}
		cout << endl;
	}
	*/

	/***************************/




}
void KM(int index, glm::vec3 k, glm::vec3 s, float x){
	glm::vec3 a = (k + s) / s;
	glm::vec3 b = sqrt(a*a - glm::vec3(1.0));
	glm::vec3 bSx = b*s*glm::vec3(x);
	glm::vec3 sinh_bSx = sinh(bSx);
	glm::vec3 c = a*sinh_bSx + b*cosh(bSx);
	R[index] = sinh_bSx / c;
	T[index] = b / c;

}
void CompositeLayers(int index)
{
	glm::vec3 tmp = glm::vec3(1.0) / (glm::vec3(1.0) - R0*R1);
	R[index] = R0 + T0 * T0 * R1 * tmp;
	T[index] = T0 * T1 * tmp;

}

void UpdateHeightMap(int index, float add,glm::vec3 temp_k,glm::vec3 temp_s)
{
	//increase  water 
	
	float mix[3];
	//pigment
	pigment[index][order] += add;
	//cout << pigment[index][order] << endl;
	//water height
	vHeightmapData[index].z += add;
	if (vHeightmapData[index].z - fHeights[index] < 0) vHeightmapData[index].z = fHeights[index];
	wp[index] = (vHeightmapData[index].z - fHeights[index])*0.0098;
	
	//wp[index] += add*0.0098;
	M_wet[(int)(index / HM_SIZE_X)][index%HM_SIZE_X] = 1;
	glBindBuffer(GL_ARRAY_BUFFER, uiVBOHeightmapData);
	glBufferSubData(GL_ARRAY_BUFFER, sizeof(glm::vec3)*(index), sizeof(glm::vec3), &vHeightmapData[index]);
	if (add>0.0f){
		glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
		if (map_k[index][0] != temp_k && layer[index][0] > 0.0f)	//different color with layer[0] && layer[0] has pigment =>add layer[1]
		{
			//	layer[index][1] += add;
			layer[index][1] = pigment[index][order];
			KM(index, map_k[index][0], map_s[index][0], layer[index][0]);
			R0 = R[index];
			T0 = T[index];
			KM(index, temp_k, temp_s, layer[index][1]);
			R1 = R[index];
			T1 = T[index];
			CompositeLayers(index);
		}
		else{	//same color with layer[0] and layer[0]
			//layer[index][0] += add;
			layer[index][0] = pigment[index][order];
			map_k[index][0] = temp_k;
			map_s[index][0] = temp_s;
			KM(index, temp_k, temp_s, layer[index][0]);

		}
		mix[0] = T[index].x + R[index].x;
		mix[1] = T[index].y + R[index].y;
		mix[2] = T[index].z + R[index].z;
		glBufferSubData(GL_ARRAY_BUFFER, sizeof(GLfloat)*(index)* 3, sizeof(GLfloat)* 3, mix);
	}
}
void decreaseHeight(int x, int y, float decrease)
{
	//if (vHeightmapData[y*HM_SIZE_X + x].z - fHeights[y*HM_SIZE_X + x] > 0){
		pigment[y*HM_SIZE_X + x][order] -= decrease;
		vHeightmapData[y *  HM_SIZE_X + x].z -= decrease;
		wp[y*HM_SIZE_X + x] = (vHeightmapData[y*HM_SIZE_X + x].z - fHeights[y*HM_SIZE_X + x])*0.0098;
		glBindBuffer(GL_ARRAY_BUFFER, uiVBOHeightmapData);
		glBufferSubData(GL_ARRAY_BUFFER, sizeof(glm::vec3)*(y *  HM_SIZE_X + x), sizeof(glm::vec3), &vHeightmapData[y *  HM_SIZE_X + x]);
	//}
}
void updateVelocities()
{
	float maxVelocity=0.0f;
	for (int i = 0; i < HM_SIZE_X; i++)
	{
		for (int j = 0; j<HM_SIZE_Y; j++)
		{
			int index = j*HM_SIZE_X + i;		
			//if (i + 1< HM_SIZE_X)
			//	r_bound[j*HM_SIZE_X + i] = r_bound[j*HM_SIZE_X + i] + fHeights[j*HM_SIZE_X + i] - fHeights[j*HM_SIZE_X + i + 1];
			//if (i - 1 >= 0)
			//	l_bound[j*HM_SIZE_X + i] = l_bound[j*HM_SIZE_X + i] + fHeights[j*HM_SIZE_X + i] - fHeights[j*HM_SIZE_X + i - 1];
			//if (j + 1< HM_SIZE_Y)
			//	u_bound[j*HM_SIZE_X + i] = u_bound[j*HM_SIZE_X + i] + fHeights[j*HM_SIZE_X + i] - fHeights[(j + 1)*HM_SIZE_X + i];
			//if (j - 1 >= 0)
			//	d_bound[j*HM_SIZE_X + i] = d_bound[j*HM_SIZE_X + i] + fHeights[j*HM_SIZE_X + i] - fHeights[(j - 1)*HM_SIZE_X + i];
			if (i + 1< HM_SIZE_X)
				r_bound[index] = fHeights[index] - fHeights[index + 1];	 //u(i+.5)j 
			if (i - 1 >= 0)
				l_bound[index] = -(fHeights[index] - fHeights[index - 1]);	//u(i - .5)j
			if (j + 1< HM_SIZE_Y)
				u_bound[index] =  fHeights[index] - fHeights[(j + 1)*HM_SIZE_X + i];		//vi(j+.5)
			if (j - 1 >= 0)
				d_bound[index] = -(fHeights[index] - fHeights[(j - 1)*HM_SIZE_X + i]);	//vi(j - .5)

		}
	}
	float A;
	float B;

	for (int i = 0; i < HM_SIZE_X; i++)
	{
		for (int j = 0; j < HM_SIZE_Y; j++)
		{
			if (i + 1 < HM_SIZE_X)
				A = pow(CalVelocity_U(i, j), 2) - pow(CalVelocity_U(i + 1, j), 2) + r_bound[j*HM_SIZE_X + i] * (d_bound[j*HM_SIZE_X + i]) - r_bound[j*HM_SIZE_X + i] * u_bound[j*HM_SIZE_X + i];
				//A = pow(CalVelocity_U(i, j), 2) - pow(CalVelocity_U(i + 1, j), 2) + r_bound[j*HM_SIZE_X + i] * (-d_bound[j*HM_SIZE_X + i]) - r_bound[j*HM_SIZE_X + i] * u_bound[j*HM_SIZE_X + i];
			else
				A = pow(CalVelocity_U(i, j), 2) - pow(CalVelocity_U(i, j), 2) + r_bound[j*HM_SIZE_X + i] * (d_bound[j*HM_SIZE_X + i]) - r_bound[j*HM_SIZE_X + i] * u_bound[j*HM_SIZE_X + i];
				//A = pow(CalVelocity_U(i, j), 2) - pow(CalVelocity_U(i, j), 2) + r_bound[j*HM_SIZE_X + i] * (-d_bound[j*HM_SIZE_X + i]) - r_bound[j*HM_SIZE_X + i] * u_bound[j*HM_SIZE_X + i];
			if (i + 1 < HM_SIZE_X && j + 1 < HM_SIZE_Y &&  j - 1 >= 0 && i - 1 >= 0)		//B = u(i+1.5)j + u(i-.5)j + u(i+.5)(j+1) + u(i+.5)(j-1) - 4u(i+.5)j
				B = r_bound[j*HM_SIZE_X + (i + 1)] + (l_bound[j*HM_SIZE_X + i]) + r_bound[(j + 1)*HM_SIZE_X + i] + r_bound[(j - 1)*HM_SIZE_X + i] - 4 * r_bound[j*HM_SIZE_X + i];
				//B = r_bound[j*HM_SIZE_X + (i+1)] + (-l_bound[j*HM_SIZE_X + i]) + r_bound[(j + 1)*HM_SIZE_X + i] + r_bound[(j - 1)*HM_SIZE_X + i] - 4 * r_bound[j*HM_SIZE_X + i];
			else
				B = r_bound[j*HM_SIZE_X + i] + l_bound[j*HM_SIZE_X + i] + r_bound[j*HM_SIZE_X + i] + r_bound[j*HM_SIZE_X + i] - 4 * r_bound[j*HM_SIZE_X + i];
				//B = r_bound[j*HM_SIZE_X + i] + (-l_bound[j*HM_SIZE_X + i]) + r_bound[j*HM_SIZE_X + i] + r_bound[j*HM_SIZE_X + i] - 4 * r_bound[j*HM_SIZE_X + i];
			
			float t = thirdStepUV(A, B, i, j);

			if (j + 1<HM_SIZE_Y)
				A = pow(CalVelocity_V(i, j), 2) - pow(CalVelocity_V(i, j + 1), 2) + l_bound[j*HM_SIZE_X + i] * u_bound[j*HM_SIZE_X + i] - r_bound[j*HM_SIZE_X + i] * u_bound[j*HM_SIZE_X + i];
				//A = pow(CalVelocity_V(i, j), 2) - pow(CalVelocity_V(i, j + 1), 2) + (-l_bound[j*HM_SIZE_X + i]) * u_bound[j*HM_SIZE_X + i] - r_bound[j*HM_SIZE_X + i] * u_bound[j*HM_SIZE_X + i];
			else
				A = pow(CalVelocity_V(i, j), 2) - pow(CalVelocity_V(i, j), 2) + l_bound[j*HM_SIZE_X + i] * u_bound[j*HM_SIZE_X + i] - r_bound[j*HM_SIZE_X + i] * u_bound[j*HM_SIZE_X + i];
				//A = pow(CalVelocity_V(i, j), 2) - pow(CalVelocity_V(i, j), 2) + (-l_bound[j*HM_SIZE_X + i]) * u_bound[j*HM_SIZE_X + i] - r_bound[j*HM_SIZE_X + i] * u_bound[j*HM_SIZE_X + i];
			if (i + 1 < HM_SIZE_X && j + 1<HM_SIZE_Y && i - 1 >= 0 && j - 1 >= 0)
				B = u_bound[j*HM_SIZE_X + (i + 1)] + d_bound[j*HM_SIZE_X + i] + u_bound[(j + 1)*HM_SIZE_X + i] + u_bound[(j - 1)*HM_SIZE_X + i] - 4 * u_bound[j*HM_SIZE_X + i];
				//B = u_bound[j*HM_SIZE_X + (i + 1)] + (-d_bound[j*HM_SIZE_X + i]) + u_bound[(j + 1)*HM_SIZE_X + i] + u_bound[(j - 1)*HM_SIZE_X + i] - 4 * u_bound[j*HM_SIZE_X + i];
			else
				B = u_bound[j*HM_SIZE_X + i] + u_bound[j*HM_SIZE_X + i] + u_bound[j*HM_SIZE_X + i] + u_bound[j*HM_SIZE_X + i] - 4 * u_bound[j*HM_SIZE_X + i];
				//B = u_bound[j*HM_SIZE_X + i] + u_bound[j*HM_SIZE_X + i] + u_bound[j*HM_SIZE_X + i] + u_bound[j*HM_SIZE_X + i] - 4 * u_bound[j*HM_SIZE_X + i];
			
			float s = sixthStepUV(A, B, i, j);
			newUs[j][i][0] = t;
			newUs[j][i][1] = s;

		}
	}

	//Update all new boundaries
	for (int i = 0; i < HM_SIZE_X; i++)
	{
		for (int j = 0; j < HM_SIZE_Y; j++)
		{
			r_bound[j*HM_SIZE_X + i] = newUs[j][i][0];
			if (i + 1 < HM_SIZE_X)
				l_bound[j*HM_SIZE_X + (i + 1)] = 1 * newUs[j][i][0];
				//l_bound[j*HM_SIZE_X + (i + 1)] = -1 * newUs[j][i][0];
			else
			{
				;
			}
			u_bound[j*HM_SIZE_X + i] = newUs[j][i][1];
			if (j + 1 < HM_SIZE_Y)
				d_bound[(j + 1)*HM_SIZE_X + i] = 1 * newUs[j][i][1];
				//d_bound[(j+1)*HM_SIZE_X + i] = -1 * newUs[j][i][1];
			else
			{;}
		}
	}
	//not in wet area velocity to zero
	for (int i = 0; i < HM_SIZE_X; i++)
	{
		for (int j = 0; j < HM_SIZE_Y; j++)
		{
			//if (wp[j*HM_SIZE_X + i] == 0)
			if (M_wet[j][i]==0)
			{
				r_bound[j*HM_SIZE_X + i] = 0;
				l_bound[j*HM_SIZE_X + i] = 0;
				u_bound[j*HM_SIZE_X + i] = 0;
				d_bound[j*HM_SIZE_X + i] = 0;
			}

		}
	}

}
void relaxDivergence()
{
	float deltamax = 0;
	int index;
	float delta=0;
	for (int t = 0; t < 50;t++){
	
		// (u',v') = (u,v)
		for (int i = 0; i < HM_SIZE_X; i++)
		{
			for (int j = 0; j < HM_SIZE_Y; j++)

			{
				newUs[j][i][0] = (r_bound[j*HM_SIZE_X + i] + l_bound[j*HM_SIZE_X + i]) / 2;
				newUs[j][i][1] = (u_bound[j*HM_SIZE_X + i] + d_bound[j*HM_SIZE_X + i]) / 2;
				//newUs[j][i][0] = (r_bound[j*HM_SIZE_X + i] - l_bound[j*HM_SIZE_X + i]) / 2;
				//newUs[j][i][1] = (u_bound[j*HM_SIZE_X + i] - d_bound[j*HM_SIZE_X + i]) / 2;
			}
		}
		// end of  (u',v') = (u,v)

		deltamax = 0;
		for (int i = 0; i < HM_SIZE_X; i++)
		{
			for (int j = 0; j < HM_SIZE_Y; j++)
			{
				index = j*HM_SIZE_X + i;
				if (wp[index] > 0.00f){
					//delta = epsilon*(r_bound[j*HM_SIZE_X + i] - (-l_bound[j*HM_SIZE_X + i]) + u_bound[j*HM_SIZE_X + i] - (-d_bound[j*HM_SIZE_X + i]));
					delta = epsilon*(r_bound[j*HM_SIZE_X + i] - l_bound[j*HM_SIZE_X + i] + u_bound[j*HM_SIZE_X + i] - d_bound[j*HM_SIZE_X + i]);
					wp[index] += delta;
					newUs[j][i][0] += 1 * delta;
					try{
						newUs[j][i - 1][0] = newUs[j][i - 1][0] - 1 * delta;
					}
					catch (...) {}
					newUs[j][i][1] += 1 * delta;
					try{
						newUs[j - 1][i][1] = newUs[j - 1][i][1] -1 * delta;
					}
					catch (...) {}
					deltamax = max((float)abs(1 * delta), deltamax);
				}
			}
		}//end for

	//update
		update();
	//update end

		if (deltamax <=1*tolerance ) {
			break;
		}
		//cout << "time out" << endl;
	}
	//end of 50 times
}
void edgeDarkening()
{
	for (int i = 0; i < HM_SIZE_X; i++)
	{
		for (int j = 0; j < HM_SIZE_Y ; j++)
		{
			float sum = 0;

			for (int m = i - 5; m < i + 5; i++)
			{
				for (int n = j - 5; n < j + 5; j++)
				{
					try{
						sum += M_wet[m][n] * weight[m + 5][n + 5];
					}
					catch (...){ sum += M_wet[i][j] * weight[m + 5][n + 5]; }
				}
			}
			blurred[i][j] = sum;
		}
	}
}

void flowOutward()
{   //calculate M'=M_flow
	//for (int x = 0; x < HM_SIZE_X; x++)
	//{
	//	for (int y = 0; y < HM_SIZE_Y; y++)
	//	{
	//		float temp=0;
	//		for (int i = -2; i < 3; i++)
	//		{
	//			for (int j = -2; j < 3; j++)
	//			{
	//				if (y+j>=0 && y+j<HM_SIZE_Y &&x+i>=0 && x+i<HM_SIZE_X)
	//					temp += M_wet[(y + j)][ (x + i)] * Mask[(i + 2)*(j + 2)];
	//				else{
	//					temp += M_wet[y][x] * Mask[(i + 2)*(j + 2)];
	//				}
	//			}
	//		}
	//		M_flow[y*HM_SIZE_X + x] = temp;
	//	}
	//}
	////end of result
	for (int i = 0; i < HM_SIZE_X; i++)
	{
		for (int j = 0; j < HM_SIZE_Y; j++)
		{
			float sum = 0;

			for (int m = i - 5; m < i + 5; i++)
			{
				for (int n = j - 5; n < j + 5; j++)
				{
					try{
						sum += M_wet[m][n] * weight[m + 5][n + 5];
					}
					catch (...){ sum += M_wet[i][j] * weight[m + 5][n + 5]; }
				}
			}
			blurred[i][j] = sum;
		}
	}

	//calculate new p
	for (int i = 0; i < HM_SIZE_X; i++)
	{
		for (int j = 0; j < HM_SIZE_Y; j++)
		{
			if (wp[j*HM_SIZE_X + i] - 0.1*(1 - M_flow[j*HM_SIZE_X + i])*M_wet[j][i]>0.01f)  {		//decrease water
				wp[j*HM_SIZE_X + i] = wp[j*HM_SIZE_X + i] - 0.1*(1 - M_flow[j*HM_SIZE_X + i])*M_wet[j][i];
				//decreaseHeight(i, j, 0.1*(1 - M_flow[j*HM_SIZE_X + i])*M_wet[j][i] / 9.80);
			}
			//if (wp[j*HM_SIZE_X + i] - 0.1*(1 - blurred[j][i])*M_wet[j][i]>0.01f)  {		//decrease water
			//	wp[j*HM_SIZE_X + i] = wp[j*HM_SIZE_X + i] - 0.1*(1 - blurred[j][i])*M_wet[j][i];
			//	decreaseHeight(i, j, 0.1*(1 - blurred[j][i])*M_wet[j][i] / 9.80);
			//}
		}
	}
	//end of new pressure




}

void movePigments() {
	int index;
	for (int n = 0; n < 1; n++)
	{
		for (int i = 0; i < HM_SIZE_X; i++)
		{
			for (int j = 0; j < HM_SIZE_Y; j++)
			{
				index = j*HM_SIZE_X + i;
				Primepigment[index][n] = pigment[index][n];
			}
		}
	}
	
	//cal delta time
	float deltatime = 0;
	for (int i = 0; i < HM_SIZE_X; i++)
	{
		for (int j = 0; j < HM_SIZE_Y; j++)
		{
			if (r_bound[j*HM_SIZE_X + i]>deltatime)
				deltatime = r_bound[j*HM_SIZE_X + i];
			if (l_bound[j*HM_SIZE_X + i]>deltatime)
				deltatime = l_bound[j*HM_SIZE_X + i];
			if (u_bound[j*HM_SIZE_X + i]>deltatime)
				deltatime = u_bound[j*HM_SIZE_X + i];
			if (d_bound[j*HM_SIZE_X + i]>deltatime)
				deltatime = d_bound[j*HM_SIZE_X + i];
		}
	}


	for (int n = 0; n < 1; n++)
	{
		for (int i = 0; i < HM_SIZE_X; i++)
		{
			for (int j = 0; j < HM_SIZE_Y; j++)
			{
				index = j*HM_SIZE_X + i;
				try{
					Primepigment[index + 1][n] +=  0.001*max(0, r_bound[j*HM_SIZE_X + i] * pigment[index][n]);
				}
				catch (...){}
				try{
					Primepigment[index - 1][n] +=  0.001*max(0, l_bound[j*HM_SIZE_X + i] * pigment[index][n]);
				}
				catch (...){}
				try{
					Primepigment[(j + 1)*HM_SIZE_X + i][n] += 0.001*max(0, u_bound[j*HM_SIZE_X + i] * pigment[index][n]);
				}
				catch (...){}
				try{
					Primepigment[(j - 1)*HM_SIZE_X + i][n] +=0.001*max(0, d_bound[j*HM_SIZE_X + i] * pigment[index][n]);
				}
				catch (...){}
				Primepigment[index][n] = Primepigment[index][n] - 0.001*max(0, r_bound[j*HM_SIZE_X + i] * pigment[index][n]) + 0.001*max(0, l_bound[j*HM_SIZE_X + i] * pigment[index][n]) + 0.001* max(0, u_bound[j*HM_SIZE_X + i] * pigment[index][n]) + 0.001*max(0, d_bound[j*HM_SIZE_X + i] * pigment[index][n]);
			}
		}
	}
	//0.001->0.01
	for (int n = 0; n < 1; n++)
	{
		for (int i = 0; i < HM_SIZE_X; i++)
		{
			for (int j = 0; j < HM_SIZE_Y; j++)
			{
			    index = j*HM_SIZE_X + i;
				pigment[index][n] = Primepigment[index][n];		
			}
		}
	}
}

void click(GLFWwindow *window,float xpos_scale, float ypos_scale)
{
	flag = true;
	int x_int, y_int;
	x_int = (int)(xpos_scale);
	y_int = (int)(ypos_scale);
	int x = x_int;
	int y = y_int;
	int index;
	int i = 0; int j = 0;
	int  radius=20;
	//cout << x << "   " << y << endl;
	for (int i = -radius; i <= radius; i++)
	{
		for (int j = -radius; j <= radius; j++)
		{	
			x = x_int + i;
			y = y_int + j;
			
			if (!OutOfBound(x,y) && (i*i + j*j) <= radius*radius)
			{
				index =(y*HM_SIZE_X +x);
				//cout << index<<endl;
				UpdateHeightMap(index,0.15f,k,s);
				M_wet[y_int + j][x_int + i] = 1;
			}
		}
	}
}
void updateCanvas()
{
	updateVelocities();
	relaxDivergence();
	//flowOutward();
	//movePigments();

	
}
void Fbo(){
	

	// create a framebuffer object, you need to delete them when program exits
	glGenFramebuffersEXT(1, &FramebufferName);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, FramebufferName);//应用绑定FBO
	/*multitexture*/
	for (int i = 0; i < 4; i++){
		//// The texture we're going to render to
		glGenTextures(1, &renderedTexture[i]);
		// "Bind" the newly created texture : all future texture functions will modify this texture
		glBindTexture(GL_TEXTURE_2D, renderedTexture[i]);
		// Give an empty image to OpenGL ( the last "0" )
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 800, 800, 0, GL_RGB, GL_UNSIGNED_BYTE, 0);
		// Poor filtering. Needed !
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		// Set "renderedTexture0" as our colour attachement #0
		glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0 + i, renderedTexture[i], 0);
	}


	



	/************************************/
	// The depth buffer
	/*glGenRenderbuffers(1, &depthrenderbuffer);
	glBindRenderbuffer(GL_RENDERBUFFER, depthrenderbuffer);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, HM_SIZE_X,  HM_SIZE_Y);
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthrenderbuffer);*/
	


	/*****/

	

	
	// Set the list of draw buffers.
	GLenum DrawBuffers[4] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1,GL_COLOR_ATTACHMENT2, GL_COLOR_ATTACHMENT3 };
	glDrawBuffers(4, DrawBuffers); // "1" is the size of DrawBuffers
	/****/
	//for (int i = 0; i < 1;i++)
	//	glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0+i, textureId[i], i);
	//static const GLenum fboBuffs[] = { GL_COLOR_ATTACHMENT0/*, GL_COLOR_ATTACHMENT1, */};
	//glDrawBuffers(1, fboBuffs);
	




	GLenum status;
	status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	switch (status)
	{
	case GL_FRAMEBUFFER_COMPLETE_EXT:
		cout << "good ";
		break;
	default:
		;
	}
}


int main(){


	if (!glfwInit()) {  //檢查GLFW是否正常被初始化
		fprintf(stderr, "ERROR: could not start GLFW3\n");
		return 1;
	}
	
	
	float fSizeX = HM_SIZE_X, fSizeY = HM_SIZE_Y;
	int wSizeX = 800, wSizeY = 800;
	double xpos, ypos;
	float scale;
	scale =(float)( fSizeX / wSizeX);

	/*
	g_Paper = new Paper(fSizeX, fSizeY);
	g_WaterColor = new WaterColor(fSizeX, fSizeY,4);
	*/

	GLFWwindow* window = glfwCreateWindow(wSizeX, wSizeY, "PAPER", NULL, NULL);
	GLFWcursor *cursor = glfwCreateStandardCursor(GLFW_ARROW_CURSOR);
	glfwSetCursor(window, cursor);
	glfwSetCursorPos(window, 0, 0);
	glfwGetFramebufferSize(window, &wSizeX, &wSizeY);
	

	if (!window) {  //檢查視窗是否正確建立
		fprintf(stderr, "ERROR: could not open window with GLFW3\n");
		glfwTerminate();
		return 1;
	}
	glfwMakeContextCurrent(window);  //將我們的視窗視為openGL當前的context



	glewExperimental = GL_TRUE;
	glewInit();		//初始化GLEW

	const GLubyte* renderer = glGetString(GL_RENDERER);
	const GLubyte* version = glGetString(GL_VERSION);
	printf("Renderer: %s\n", renderer);
	printf("OpenGL version supported %s\n", version);

	/*paper generation*/
	float **map = createMap();
	for (int i = 0; i < HM_SIZE_X; i++)
	{
		for (int j = 0; j < HM_SIZE_Y; j++)
		{		 
			fHeights[i*HM_SIZE_X + j] = fabs(map[i][j]);
		}
	}
	/*********/
	/*	fluid capacity*/
	for (int i = 0; i < HM_SIZE_X; i++)
	{
		for (int j = 0; j < HM_SIZE_Y; j++)
		{
			c[i][j] = map[i][j] * (cmax - cmin) + cmin;
		}
	}
	/**************/
	//edgeDarkening();


	GLint maxbuffers;		//value is  8
	glGetIntegerv(GL_MAX_DRAW_BUFFERS, &maxbuffers);
	

	glEnable(GL_DEPTH_TEST);	//深度測試，讓近的東西蓋掉遠的
	glEnable(GL_TEXTURE_2D);
	glGenVertexArrays(1, &uiVAOHeightmap); // Create one VAO
	glGenBuffers(1, &uiVBOHeightmapData); // One VBO for data


	glBindVertexArray(uiVAOHeightmap);
	glBindBuffer(GL_ARRAY_BUFFER, uiVBOHeightmapData);

	/*for (int i = 0; i < HM_SIZE_X*HM_SIZE_Y; i++)
	{
		fHeights[i] = rand_FloatRange(0, 0.005);
	}*/
	for (int k = 0; k < HM_SIZE_X*HM_SIZE_Y * 3; k++)
	{
		colors[k] = 0.9999999f;
	}

	for (int i = 0; i < HM_SIZE_X*HM_SIZE_Y; i++)
	{
		float column = float(i%HM_SIZE_X), row = float(i / HM_SIZE_X);
		vHeightmapData[i] = glm::vec3(
			/*-fSizeX / 2 +*/ fSizeX*column / float(HM_SIZE_X - 1), // X Coordinate
			/*-fSizeZ / 2 +*/ fSizeY*row / float(HM_SIZE_Y - 1),	// Z Coordinate
			fHeights[i]									// Y Coordinate (it's height)			
			);
	}

	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3)*HM_SIZE_X*HM_SIZE_Y, vHeightmapData, GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(0);





	//add
	//color
	glGenBuffers(1, &colors_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);	//綁定vbo
	glBufferData(GL_ARRAY_BUFFER, sizeof(colors), colors, GL_STATIC_DRAW);
	////將陣列colors的值傳給vbo
	glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(1);

	glGenBuffers(1, &uiVBOIndices);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, uiVBOIndices);

	Fbo();

	
	////////
	unsigned int vIndex = -1;
	int ox = 0, oz = HM_SIZE_X;
	int width = HM_SIZE_X;


	for (int i = 0; i < HM_SIZE_X - 1; i++){
		for (int j = 0; j < HM_SIZE_Y; j++)
		{
			iIndices[++vIndex] = ox + i*width + j;
			iIndices[++vIndex] = oz + i*width + j;
		}
		if (i != HM_SIZE_Y - 2)
			iIndices[++vIndex] = HM_SIZE_X*HM_SIZE_Y;
	}

	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(iIndices), iIndices, GL_STATIC_DRAW);
	glEnable(GL_PRIMITIVE_RESTART);
	glPrimitiveRestartIndex(HM_SIZE_X*HM_SIZE_Y);

	//////fbo//////

	static const GLfloat g_quad_vertex_buffer_data[] = {
		-(wSizeX / 2), -(wSizeY / 2), 0.00f,
		(wSizeX / 2), -(wSizeY / 2), 0.00f,
		-(wSizeX / 2), (wSizeY / 2), 0.00f,
		-(wSizeX / 2), (wSizeY / 2), 0.00f,
		(wSizeX / 2), -(wSizeY / 2), 0.00f,
		(wSizeX / 2), (wSizeY / 2), 0.00f,
	};
	
	//static const GLfloat g_quad_vertex_buffer_data[] = {
	//	-(HM_SIZE_X / 2), -(HM_SIZE_Y / 2), 0.00f,
	//	(HM_SIZE_X / 2), -(HM_SIZE_Y / 2), 0.00f,
	//	-(HM_SIZE_X / 2), (HM_SIZE_Y / 2), 0.00f,
	//	-(HM_SIZE_X / 2), (HM_SIZE_Y / 2), 0.00f,
	//	(HM_SIZE_X / 2), -(HM_SIZE_Y / 2), 0.00f,
	//	(HM_SIZE_X / 2), (HM_SIZE_Y / 2), 0.00f,
	//};
	//static const GLfloat g_quad_vertex_buffer_data[] = {
	//	-200, -200, 0.00f,
	//	0, -200, 0.00f,
	//	-200, 0, 0.00f,
	//	-200,0, 0.00f,
	//	0, -200, 0.00f,
	//	0, 0, 0.00f,
	//};
	static const GLfloat g_quad_vertex_buffer_data1[] = {
		0, 0, 0.00f,
		200, 0, 0.00f,
		0, 200, 0.00f,
		0, 200, 0.00f,
		200, 0, 0.00f,
		200, 200, 0.00f,
	};
	GLuint quad_vertexbuffer0;
	glGenBuffers(1, &quad_vertexbuffer0);
	glBindBuffer(GL_ARRAY_BUFFER, quad_vertexbuffer0);
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_quad_vertex_buffer_data), g_quad_vertex_buffer_data, GL_STATIC_DRAW);

	GLuint quad_vertexbuffer1;
	glGenBuffers(1, &quad_vertexbuffer1);
	glBindBuffer(GL_ARRAY_BUFFER, quad_vertexbuffer1);
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_quad_vertex_buffer_data1), g_quad_vertex_buffer_data1, GL_STATIC_DRAW);

	

	//shader
	GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);		//宣告VertexShader的ID
	GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);	//宣告FragmentShader的ID

	std::string source;
	loadFile("shader.vert", source);
	const char *vertex_shader = source.c_str();			//把std::string結構轉換成const char*

	glShaderSource(VertexShaderID, 1, &vertex_shader, NULL);		//指明VertexShader的原始碼來源
	glCompileShader(VertexShaderID);								//編譯VertexShader
	source = "";
	loadFile("shader.frag", source);
	const char *fragment_shader = source.c_str();

	glShaderSource(FragmentShaderID, 1, &fragment_shader, NULL);	//指明FragmentShader的原始碼來源
	glCompileShader(FragmentShaderID);								//編譯FragmentShader

	unsigned int shader_programme = glCreateProgram();		//建立shader program
	glAttachShader(shader_programme, VertexShaderID);		//把shader attach到shader program上
	glAttachShader(shader_programme, FragmentShaderID);
	glLinkProgram(shader_programme);						//綁定shader program

	////fbo shader/////
	GLuint FboVertexShaderID = glCreateShader(GL_VERTEX_SHADER);		//宣告VertexShader的ID
	GLuint FboFragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);	//宣告FragmentShader的ID
	source = "";
	loadFile("fbo.vert", source);
	const char *fbovertex_shader = source.c_str();			//把std::string結構轉換成const char*

	glShaderSource(FboVertexShaderID, 1, &fbovertex_shader, NULL);		//指明VertexShader的原始碼來源
	glCompileShader(FboVertexShaderID);								//編譯VertexShader
	source = "";
	loadFile("fbo.frag", source);
	const char *fbofragment_shader = source.c_str();

	glShaderSource(FboFragmentShaderID, 1, &fbofragment_shader, NULL);	//指明FragmentShader的原始碼來源
	glCompileShader(FboFragmentShaderID);								//編譯FragmentShader

	unsigned int fboshader_programme = glCreateProgram();		//建立shader program
	glAttachShader(fboshader_programme, FboVertexShaderID);		//把shader attach到shader program上
	glAttachShader(fboshader_programme, FboFragmentShaderID);
	glLinkProgram(fboshader_programme);						//綁定shader program



	//GLFW
	GLint matrix_location = glGetUniformLocation(shader_programme, "ORTHO");
	

	initialize();
	


	k = K_QuinacridoneRose;
	s = S_QuinacridoneRose;
	order = 0;
	while (!glfwWindowShouldClose(window)) {
		glBindFramebuffer(GL_FRAMEBUFFER, FramebufferName);
	
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glViewport(0, 0, wSizeX, wSizeY);
		
		glUseProgram(shader_programme);							//使用shader program
		//my add
		glm::mat4 ortho = glm::mat4(1);
		ortho *= glm::ortho<float>(0.0f, HM_SIZE_X, HM_SIZE_Y, 0, -450.0f,450.0f);
		glUniformMatrix4fv(matrix_location, 1, GL_FALSE, &ortho[0][0]);//傳送uniform值給shader

		glBindBuffer(GL_ARRAY_BUFFER, uiVBOHeightmapData);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(0);

		glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);	//綁定vbo
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
		glEnableVertexAttribArray(1);

		glBindVertexArray(uiVAOHeightmap);		//綁定vao
		glDrawElements(GL_TRIANGLE_STRIP, HM_SIZE_Y*(HM_SIZE_X - 1) * 2 + HM_SIZE_X - 1, GL_UNSIGNED_INT, 0);		//繪圖


		//handle pigment
		if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
		{
			glfwGetCursorPos(window, &xpos, &ypos);
			if (xpos*scale < HM_SIZE_X && xpos*scale >= 0.0f && ypos*scale < HM_SIZE_Y&& ypos*scale >= 0.0f)
			{				
				click(window, xpos*scale, ypos*scale);
			}
		}
		glfwSetKeyCallback(window, ColorChanged);
		updateCanvas();
		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);
		glUseProgram(0);																//取消綁定shader program

		// Render to the screen
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		// Render on the whole framebuffer, complete from the lower left corner to the upper right

		glViewport(0, 0, wSizeX, wSizeY);
		// Clear the screen
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glUseProgram(fboshader_programme);

		// Bind our texture in Texture Unit 0
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, renderedTexture[0]);
		GLint texture_location = glGetUniformLocation(fboshader_programme, "renderedTexture");
		glUniform1i(texture_location, 0);

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, renderedTexture[1]);
		texture_location = glGetUniformLocation(fboshader_programme, "renderedTexture1");
		glUniform1i(texture_location, 0);
		
		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, renderedTexture[2]);
		texture_location = glGetUniformLocation(fboshader_programme, "renderedTexture2");
		glUniform1i(texture_location, 1);
		
		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, renderedTexture[3]);
		texture_location = glGetUniformLocation(fboshader_programme, "renderedTexture3");
		glUniform1i(texture_location, 2);
		


	


		
		// 1rst attribute buffer : vertices
		glEnableVertexAttribArray(0);

		glBindBuffer(GL_ARRAY_BUFFER, quad_vertexbuffer0);
		glVertexAttribPointer(
			0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
			3,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
			);

		// Draw the triangles !
		glDrawArrays(GL_TRIANGLES, 0, 6); // 2*3 indices starting at 0 -> 2 triangles
		glDisableVertexAttribArray(0);

		glEnableVertexAttribArray(1);

		glBindBuffer(GL_ARRAY_BUFFER, quad_vertexbuffer1);
		glVertexAttribPointer(
			0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
			3,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
			);

		// Draw the triangles !
		glDrawArrays(GL_TRIANGLES, 0, 6); // 2*3 indices starting at 0 -> 2 triangles
		glDisableVertexAttribArray(1);
		
		glfwSwapBuffers(window);										//刷新畫面
		glfwPollEvents();
	
	}
	glDeleteShader(VertexShaderID);
	glDeleteShader(FragmentShaderID);
	glDeleteProgram(shader_programme);
	glDeleteFramebuffers(1, &FramebufferName);
	//glDeleteTextures(1, &renderedTexture);
	for (int i = 0; i < 4;i++)
		glDeleteTextures(1, &renderedTexture[i]);
	glDeleteBuffers(1, &quad_vertexbuffer0);
	//glDeleteBuffers(1, &quad_vertexbuffer1);
	glfwTerminate();
	
	

	return 0;
} 