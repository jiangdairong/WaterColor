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
	int** ActiveState; // ��������lattice�����A, 
	// 0:���Χ�s 1:�n�ˬdcollapse 2:�n���ɤl�B�z 3:�j����ø
	int** ActiveWaterState;
	bool** IsHeightChg; // �O���o�^�X�O�_�����ܰ���, �Y�u��ActiveState, 
	// �h�@�^�X��, ���Ǯ�l�|�]�������񪺾F�~�N���פ����L, 
	// �y���L�]�N���צA�����O�H, 
	// Collapse
	float ThresholdDiff; // if height difference > threshold, then sand go from high to low

	// SandParticle System

	// ���F���ްO����, �@�}�l�i�Hcreate�̤j�Ȫ�particle, �}���@��static array 
	// �C��particle�O���O�_�b�ϥΤ�, 
	// �n�W�[particle��, ���n���ˬd�ѤU���D�ϥΤ���particle��������

	float GravityAcc; // ���O�[�t�� G
	float TimeStep; // �ɶ� T
	float BounceFactor; // �ϼu�᪺�t��scale factor
	int ParticleIndexer; // �l�ܲɤl����m��
	float AddHeight; // �s�W�ɨF����l����
	float ParticleRadius; // �F���ɪ��w�]�b�|
	bool IsParticleExist; // true�hrun AdvanceParticle, false�h����
	int ParticleCount; // AdvanceParticle�ɲέp�ثe�@���h�����ɤl

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
	float EyePos[3]; // for Shadow (�ثe�S�b��)
	float LightPos[3]; // for Shadow (�ثe�S�b��)
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