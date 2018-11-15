#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>

#include "MyStructures.h"

class WaterColor
{
public:
	// Resolution
	int mDx;
	int mDy;

	HeightMap* mHeightMap;
	float** LastHeight; //  recording for undo
	float** TempHeight; // for collapse
	int** ActiveState; // 紀錄哪些lattice的狀態, 
	// 0:不用更新 1:要檢查collapse 2:要做粒子處理 3:強迫重繪
	int** ActiveWaterState;
	bool** IsHeightChg; // 記錄這回合是否有改變高度, 若只用ActiveState, 
	// 則一回合內, 有些格子會因為它附近的鄰居將高度分給他, 
	// 造成他也將高度再分給別人, 
	// Collapse
	float ThresholdDiff; // if height difference > threshold, then sand go from high to low

	// SandParticle System

	// 為了控管記憶體, 一開始可以create最大值的particle, 開成一個static array 
	// 每個particle記錄是否在使用中, 
	// 要增加particle時, 須要先檢查剩下的非使用中的particle夠不夠用

	float GravityAcc; // 重力加速度 G
	float TimeStep; // 時間 T
	float BounceFactor; // 反彈後的速度scale factor
	int ParticleIndexer; // 追蹤粒子的位置用
	float AddHeight; // 新增時沙的初始高度
	float ParticleRadius; // 沙顆粒的預設半徑
	bool IsParticleExist; // true則run AdvanceParticle, false則不用
	int ParticleCount; // AdvanceParticle時統計目前共有多少顆粒子

	// Read Height Map From Script
	float InitHeight;
	bool IsNoise;

	// Check Max Height 
	bool IsLimitHeight;
	int CMcount;
	int CMMAX;
	float maxheight;
	float waterMaxHeight;
	// Local Updating
	//myMath::myPoint mUpdateLeftTop;
	//myMath::myPoint mUpdateRightBottom;
	//myMath::myPoint mUpdateWaterLeftTop;
	//myMath::myPoint mUpdateWaterRightBottom;
	bool DiffusionActive;

	// Rendering 
	float Sharpness;
	int TexNum;
	int shakeFrames;
	int shakeX;
	int shakeY;
	float EyePos[3]; // for Shadow (目前沒在用)
	float LightPos[3]; // for Shadow (目前沒在用)
	//AUX_RGBImageRec* mTextureheight;
	//AUX_RGBImageRec* mWaterTexture;
	void clearActiveState(int flag);
	void clearActiveWaterState(int flag);
	void setUpdateRegion();
	void InitializeHeightMap(float** heightmap, int height, int width);
	void readHeightMap(const char* filename);
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	// Constructor and Destructor
	WaterColor(int i, int j, int numTex);
	~WaterColor();
};