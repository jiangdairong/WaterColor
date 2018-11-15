#include "paper.h"
using namespace std;
Paper::Paper(float a_fWidth, float a_fHeigh)
{
	HeightMap = (float**)malloc(GridHeight*sizeof(float));
}
