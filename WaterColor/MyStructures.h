#include <vector>
#include "MyMath.h"

struct Texel{
	float u;
	float v;
	Texel()
	{
		u = 0.0f;
		v = 0.0f;
	}
	Texel& operator=(Texel const& copy)
	{
		u = copy.u;
		v = copy.v;
		return *this;
	}
};
struct SandLayer{
	int texIndex;
	float height;
	Texel texel;
	SandLayer()
	{
		height = 0.0f;
	}
};
struct SandParticle{
	float velocity[3];
	float position[3];
	float radius;
	//int index; // particle created index, for debugging
	SandParticle(){}
	~SandParticle(){}
	Texel texel;
	int texIndex;
};
struct Lattice{
	Lattice()
	{
		height = 0.0f;
		texIndex = 0;
		ParticleSize = 0;
		ParticleHeight = 0.0f;
		normal[0] = normal[1] = normal[2] = 0.0f;
	}
	~Lattice()
	{
		particle.clear();
	}
	float height;
	float normal[3];

	vector<SandParticle> particle;
	int ParticleSize;
	float ParticleHeight;

	vector<SandLayer> vLayers; // vector
	SandLayer* aLayers; // array

	Texel texel;
	int texIndex;

	int highestParticleIndex()
	{
		float maxH = particle[0].position[1];
		int maxI = 0;
		for (int i = 1; i<(int)particle.size(); ++i)
		{
			float h = particle[i].position[1];
			if (h > maxH)
			{
				maxI = i;
				maxH = h;
			}
		}
		return maxI;
	}
	void addLayerByTexIndex(float _height, Texel _texel, int _texIndex)
	{
		aLayers[_texIndex].texel = _texel;
		aLayers[_texIndex].height += _height;
		texIndex = _texIndex;
		texel = _texel;
	}
	void addLayer(float _height, Texel _texel, int _texIndex)
	{
		// case 0: 這個lattice上沒有任何沙子
		// case 1: 這個lattice上有沙子, 且texIndex一樣
		// case 2: 這個lattice上有沙子, 但texIndex不一樣
		if (vLayers.size() == 0)
		{
			SandLayer layer;
			layer.height = _height;
			layer.texel = _texel;
			layer.texIndex = _texIndex;
			vLayers.push_back(layer);
			return;
		}
		int topLayer = vLayers.size() - 1;
		int texIndexTop = vLayers[topLayer].texIndex;
		if (texIndexTop == _texIndex)
		{
			vLayers[topLayer].height += _height;
			vLayers[topLayer].texel = _texel;
		}
		else
		{
			SandLayer layer;
			layer.height = _height;
			layer.texel = _texel;
			layer.texIndex = _texIndex;
			vLayers.push_back(layer);
		}
	}
};



class HeightMap{
public:
	int mDx;
	int mDy;
	float MaxHeight;
	Lattice** lattices;
	HeightMap(int x, int y, int numTex)
	{
		MaxHeight = 2.5f;
		mDx = x;
		mDy = y;
		lattices = new Lattice*[mDx];
		for (int i = 0; i<mDx; ++i)
		{
			lattices[i] = new Lattice[mDy];
			for (int j = 0; j<mDy; ++j)
			{
				lattices[i][j].ParticleSize = 0;
				lattices[i][j].aLayers = new SandLayer[numTex];
				for (int k = 0; k<numTex; ++k)
				{
					lattices[i][j].aLayers[k].height = 0.0f;
				}
				//lattices[i][j].normal[0] = 0.0f;
				//lattices[i][j].normal[1] = 1.0f;
				//lattices[i][j].normal[2] = 0.0f;
			}
		}

	}
	~HeightMap(){}
	void setHeightAssign(float h, bool IsNoise)
	{
		for (int i = 0; i<mDx; ++i)
		{
			for (int j = 0; j<mDy; ++j)
			{
				if (!IsNoise)
				{
					lattices[i][j].height = h;
				}
				else
				{
					lattices[i][j].height = h *((rand() % 10) / 10.0f);
				}
			}
		}
	}
	void computeNormal(int row, int col)
	{
		float v0[3] = { -2.0f, lattices[row - 1][col].height - lattices[row + 1][col].height, 0.0f };
		float v1[3] = { 0.0f, lattices[row][col + 1].height - lattices[row][col - 1].height, 2.0f };
		myMath::Cross3X3(v0, v1, lattices[row][col].normal);

		// noise 
		// cX + cY +cZ = 0.1
		//float changeX = rand()%11/2000.0f; // cX = [0.0, 0.05]
		//float changeY = rand()%11/2000.0f; // cY = [0.0, 0.05]
		//float changeZ = 0.1f - changeX - changeY;
		//lattices[row][col].normal[0] += changeX;
		//lattices[row][col].normal[1] += changeY;
		//lattices[row][col].normal[2] += changeZ;
	}
	void computeTotalLatticeNormal()
	{
		for (int i = 1; i<mDx - 1; ++i)
		for (int j = 1; j<mDy - 1; ++j)
		{
			computeNormal(i, j);
		}
	}
};