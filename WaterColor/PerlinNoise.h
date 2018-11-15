/*****************************************FILENAME DEFINITION***/
const char oname[16] = "test.bmp";

/*******************************************STRUCT DEFINITION***/
struct color{
	//v[0]=red, v[1]=green, v[2]=blue
	unsigned char v[3];

	color(unsigned char r, unsigned char g, unsigned char b){
		v[0] = r;
		v[1] = g;
		v[2] = b;
	}
};

/********************************************SIZE DEFINITIONS***/
const unsigned hgrid =  200,//x dimension of the grid
								vgrid =  200;//y dimension of the grid

/*****************************************FUNCTION PROTOTYPES***/
float **createMap();
float random(float r);
color lerp(color c1, color c2, float value);//LERP = Linear intERPolation

void fillMap(float **map, float &min, float &max);// this is the algorithm part (the interesting part)
void printMap(float map[][vgrid], float min, float max);//bitmap part
void printPage(time_t beginning, time_t end);//webpage part

//the remaining functions are all used by the Perlin Noise function, not the template
int   myfloor(float value);
float dotproduct(float grad[], float x, float y);
float lerp(float left, float right, float amount);
float fade(float x);