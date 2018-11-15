#include "Watercolor.h"
using namespace std;

WaterColor::WaterColor(int i, int j, int numTex)
{
	DiffusionActive = false;
	mDx = i;//length of x-axis 1080 720 540
	mDy = j;//length of y-axis 720  480 360
	TexNum = numTex;

	mHeightMap = new HeightMap(mDx, mDy, numTex);
	LastHeight = new float*[mDx];
	TempHeight = new float*[mDx];
	ActiveState = new int*[mDx];
	ActiveWaterState = new int*[mDx];
	IsHeightChg = new bool*[mDx];
	for (int i = 0; i<mDx; ++i)
	{
		LastHeight[i] = new float[mDy];
		TempHeight[i] = new float[mDy];
		ActiveState[i] = new int[mDy];
		ActiveWaterState[i] = new int[mDy];
		IsHeightChg[i] = new bool[mDy];
	}

	//readScript("Script.txt");
	//readHeightMap("output.png");
	//readHeightMap("hand0.bmp");
	//readHeightMap("sand_qq.bmp");
	//TestHeight();

	//$$$$$$$$$$$$$$$$$$$$$$$
	// 參數設定
	//$$$$$$$$$$$$$$$$$$$$$$$
	//GravityAcc = -16.0f; // 重力加速度 G
	//TimeStep = 0.03f; // 時間 T
	//BounceFactor = 0.1f; // 反彈後的速度scale factor
	ParticleIndexer = 0;
	ParticleRadius = 0.25f;
	AddHeight = 8.0f;

	CMcount = 0;
	CMMAX = 30;
	IsParticleExist = false;

	//IsLimitHeight = true;
	IsLimitHeight = false;
	printf("IsLimitHeight : %d\n", IsLimitHeight);

	ThresholdDiff = tan(34.0f*PI / 180.0f);// 休止角34度
	printf("ThresholdDiff : %f\n", ThresholdDiff);

	Sharpness = 0.6f;
	shakeFrames = 0;
	maxheight = 0.0f;
	waterMaxHeight = 0.0f;
	clearActiveState(3);
	clearActiveWaterState(0);

	//printf("SetUpdateRegion");
	//setUpdateRegion();
	printf("...ok\n");
}
WaterColor::~WaterColor()
{
	for (int i = 0; i<mDx; ++i)
	{
		delete[] LastHeight[i];
		delete[] TempHeight[i];
		delete[] ActiveState[i];
		delete[] IsHeightChg[i];
	}
	delete mHeightMap;
}

void WaterColor::clearActiveState(int flag)
{
	//printf("ClearActiveState(%d)\n", flag);
	for (int i = 0; i<mDx; ++i)
	for (int j = 0; j<mDy; ++j)
	{
		ActiveState[i][j] = flag;
	}
}
void WaterColor::clearActiveWaterState(int flag)
{
	//printf("ClearActiveWaterState(%d)\n", flag);
	for (int i = 0; i<mDx; ++i)
	for (int j = 0; j<mDy; ++j)
	{
		ActiveWaterState[i][j] = flag;
	}
}

void WaterColor::InitializeHeightMap(float** heightmap, int height, int width)
{
	for (int i = 0; i<height; i++)
	for (int j = 0; j<width; j++)
		heightmap[i][j] = 0.0f;
}




void WaterColor::readHeightMap(const char* filename)
{
	int width, height;
	unsigned char* HeightMap_rgb_data;
	unsigned char* HeightMap_gray_data;
	FILE* load;

	// bmp data
	width = (int)mDx;
	height = (int)mDy;

	// allocate buffer
	HeightMap_rgb_data = (unsigned char*)malloc(width * height * 3);
	HeightMap_gray_data = (unsigned char*)malloc(width * height);
	// open and read texture data
	load = fopen(filename, "rb");
	if (load == NULL)
	{
		printf("LOAD HEIGHT MAP %s ERROR\n", filename);
	}
	else
	{
		printf("LOAD HEIGHT MAP %s OK\n", filename);
	}
	fseek(load, 54, SEEK_SET);//跳過前面header
	fread(HeightMap_rgb_data, width*height * 3, 1, load);
	fclose(load);
	float max = 0.0f;
	float min = 1024.0f;// 255*3 最多也才700多
	int j = 0;
	for (int i = 0; i<width*height * 3; i += 3, j++)
	{
		float sum_rgb = HeightMap_rgb_data[i + 2] + HeightMap_rgb_data[i + 1] + HeightMap_rgb_data[i];
		HeightMap_gray_data[j] = sum_rgb / 3;
		if (HeightMap_gray_data[j] > max) max = HeightMap_gray_data[j];
		if (HeightMap_gray_data[j] < min) min = HeightMap_gray_data[j];
	}

	for (int i = 0; i<width*height; ++i)
	{
		int row = (i) / mDx;
		int col = (i) % mDy;
		float gray = HeightMap_gray_data[i];

		if (gray <= min + 16)
		{
			mHeightMap->lattices[row][col].height = 0.0f;
			continue;
		}
		// log10(1) = 0, log10(10) = 1
		// 想要高度介於0 - maxh
		// 則logN須為0-1的值
		// 則N需為1-10的值
		// 所以乘(1-灰階值)*9+1
		mHeightMap->lattices[row][col].height = mHeightMap->MaxHeight*log10((1.0f - gray / 255)*9.0f + 1.0f);
		Texel texel;
		texel.u = row / (float)mDx;
		texel.v = col / (float)mDy;
		mHeightMap->lattices[row][col].addLayer(mHeightMap->lattices[row][col].height, texel, 2);
	}
	free(HeightMap_rgb_data);
	free(HeightMap_gray_data);
}