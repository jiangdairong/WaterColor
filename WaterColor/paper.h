#include <vector>
#include <queue>
using namespace std;
class Paper{
public:
	/************** Water Droplet ****************/
	//vector<Droplet> DropletCntr;
	float ** HeightMap;
	//int ** IDMap;
	//bool **isStaticMap;

	Paper(float a_fWidth, float a_fHeight);
	//void Initialize_WaterDropletSimulation();
	//void CreateDroplets();
	//void newDroplet(Particle p, float parVel);
	////**** 5.2.1 Apply external force - Gravity ****//
	//void AddGravity();
	////**** Update position of each particle every frame ****//
	//void IntegratePos();
	////*** Assign surface's material, rough / smooth surface  ***//
	//// if it is smooth, constant Cij of each grid is similar to each other // 
	//void AssignSurfaceMaterial();
	//void AddDroplets(int type);
	//void AddBigDroplets();

	////save normal & height map
	//int write_bmp_normal_(const char *filename, int width, int height);
	//int write_bmp_height_(const char *filename, int width, int height);
	////*** Pass height map of current frame to rendering pipeline to do rendering. ***//
	//int write_bmp_height(AUX_RGBImageRec *ndata, int width, int height, bool reset);

	//Droplets
	////**** 5.2.2 Choose Direction of movement. ****//
	//void ChooseDirection(int curDropIndex);
	////*** Merge two droplets ***//
	//void MergeDroplets(int DropletIndex_I, int DropletIndex_M);
	//int  Merge(int ID1, int ID2);
	////*** 5.2.6 Merging detection. ***//
	//void MergeDetect();
	////1 With no water ahead decelerates and stops 
	////2 when mass is less than dynamic critical mass
	//void DecelerateAndStop(int curDropIndex);
	//int  CheckPlaneBound(CVector3D CurPos, int Direction, int MoveCount);

	//Residual water
	//void CreateResidualParticle(int a_Index);

	//Height Map
	//*** 5.2.4 Compute height map. ***//
	void ComputeHeightMap();
	//*** 5.2.5 Update height map. ***//
	void SmoothHeightMap();
	void ReduceHeightMap();
	//*** droplet can move more than 2 grids, interpolate intermediate steps ***//
	//void DrawDrop(int i);
	//void DrawHeightMap();
	//int   CheckDirType(CVector3D vel);
	//bool  OutOfBound(int x, int y);
	////*** 5.2.3 Assign shape of droplets. ***//
	//void  SelectShapeFromDatabase(Droplet &curDroplet, float ParentVel);

	//float computeMass();
	//CVector3D ComputeNormal(int i, int j);

	//vector<int> RecordDropletIndex;
	//vector<int> IDCanUse;

	//Constant  m_Constant;
	//CVector3D m_DirArray[5];

	int m_GridNum;
	int m_ParticleNum;
	float TimerAcc;
	float maxheight;

	int saveCount;
	void WaterDropletSimulation();
	int GridWidth, GridHeight;

	void WaterDropletCreateMesh(int cnt);

	/************** Water Droplet ****************/

	/************** Water Art ****************/
	//time_t timer, start,  end;


	//brush mBrush;

	//CVector3D mBrushCenter;

	//CVector3D ScreenPos;
	//vector<WaterArea> WaterAreaList;
	//Grid** GridIDMap;
	//GridPtr** GridPtrIDMap;

	//int ActiveParticles;

	//vector<WaterAreaID> WaterAreaIDList;


	//int ** StaticMap;
	//int ** WaterAreaIDMap, **ParticleIDMap; // look up water area IDMap first, then look up particle IDMap
	//int ** WaterAreaIDContainerMap;
	//int ** BrushMap;

	//Follower* mFollower;
	//void WaterArtSimulation(bool);
	//void ResetWaterArtSimulation();
	//void Initialize_WaterArtSimulation();
	//void AddWaterArea();

	//void AddForce(CVector3D velocity, float ratio);
	////void ApplyFriction();
	//void AdherePlane();
	//void ApplyViscosityImpulse();
	//void Intergrate();
	//void SavePos();
	//void SaveVel();
	//void AdjustSpring();
	//void DoubleDensityRelaxation();
	//void ApplySpringDisplacements();

	//void ApplyCohesiveForce();
	//void HandleCollision();

	//void WaterAreaMergeDetect();
	//int  MergeWaterArea(int, int);

	//void SmoothHeight();
	//void ReduceHeight();
	//void ComputeHeight(bool isClear);
	//void ComputeHeightVolume();
	//void ComputeSphereHeight(int WaterAreaID, int ParticleID);

	////math
	//void ComputeGapHeight(Particle pi, Particle pj);
	//vector<CVector3D> ComputeTangentPointsSameRadius(CVector3D P1, CVector3D P2, float r);
	//vector<CVector3D> ComputeTangentPointsOfTwoSpheres(CVector3D P1, CVector3D P2, float r0, float r1, CVector3D& intersectionPoint);
	//float ComputeTriangleArea(CVector3D A, CVector3D B, CVector3D C);
	//bool isInsideQuad(vector<CVector3D> quad, CVector3D P);
	//bool isInsideRegion(CVector3D P, vector<CVector3D> tangentPointsList, CVector3D intersectionPoint, CVector3D P1, CVector3D P2);
	//bool existOuterTangents(vector<CVector3D> tangentPointsList);
	//void ComputeSearchBox(vector<CVector3D> PointsList, float* min, float* max);
	//CVector3D ComputeIntersectionOfPointToLine(CVector3D A, CVector3D B, CVector3D P);
	//float InterpolateHeight(Particle pi, Particle pj, CVector3D P);

	////math, for same radius of particles
	//void ComputeGapHeight_SameRadius(Particle pi, Particle pj);
	//vector<CVector3D> ComputeTangentPointsOfTwoSpheres_SameRadius(CVector3D P1, CVector3D P2, float r, CVector3D& intersectionPoint);
	//float InterpolateHeight_SameRadius(CVector3D A, CVector3D B, CVector3D P, float r);


	////visual funciton for debugging
	//void VisualFunction();
	//void DrawWaterAreaIDMap();
	//void DrawBrushCenter();
	//void DrawParticles();
	//void ShowSpring();
	//void DrawParticleIDMap();
	//void DrawWaterAreaContainerIDMap();
	//int ComputeParticleNumber();
	//void WaterArtCreateMesh(int cnt, bool isPressRender);
	//void printParticlePos();
	//int WriteHeightMapBmp(const char *filename, int width, int height);
	//CVector3D ComputeNormalY(int i, int j);

	////Brush
	//void CollisionDetectWithBrush();
	//void AdhereToBrush();
	//void AdjustHydrogenBonds(bool showfpsFlag);
	//void DrawHydrogenBonds();
	//void DrawHydrogenBondsToBrush();
	//void PullToBrush();
	//bool isRightButtonDown;
	//void ReleaseBrush();
	//void DeleteAdhesionHydrogenBonds();
	//void ApplyWaterPlaneFriction();

	////GUI
	//void SaveCurrentSurface(int cnt, string fileName);
	//void LoadSurface(int cnt, string fileName);
	//int saveCount_surface;
	//bool startExport, startConstructSurface;
	//void ExportFrames();
	//void ConstructSurface();
	//void ImportSurface(int cnt);
	//int frameCnt, constructCnt;
	//void WritePovRayScript(int cnt);
	//bool importHasBeenStopped;
	//bool isParticleSameRadius;
	//bool isStopExport;


	//vector<WaterArea> WaterAreaList_backupList;
	//int backupStep;
	//void BackupOneStep();
	//void Undo();
	//void Redo();
	//void Render(); // pov ray
	//AUX_RGBImageRec* mTextureheight;
	//void clearRenderedImage();
	//// for redo undo
	//int currentStep;

	////spatial
	//bool isUsingSpatial;
	//vector<Bucket> BucketListSpring;
	//float cellSize;
	//int SpatialWidth, SpatialHeight;
	//void UpdateBuckets();
	//void UpdateSpringMapSpatial();
};