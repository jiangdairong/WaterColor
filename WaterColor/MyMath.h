#include <math.h>
#include <vector>
#define PI 3.14159265f
using namespace std;
namespace myMath{
	struct myPoint{
		int x;
		int y;
		int index;
		int thickness;
		myPoint(int a, int b)
		{
			x = a;
			y = b;
			thickness = 1;
		}
		myPoint()
		{
			x = 0;
			y = 0;
			thickness = 1;
		}
		myPoint operator= (const myPoint & p)
		{
			x = p.x;
			y = p.y;
			thickness = p.thickness;
			return *this;
		}
		myPoint operator- (const myPoint & p)
		{
			x = x - p.x;
			y = y - p.y;
			return *this;
		}
		myPoint operator+ (const myPoint & p)
		{
			x = x + p.x;
			y = y + p.y;
			return *this;
		}
		myPoint operator* (const int & s)
		{
			x = x*s;
			y = y*s;
			return *this;
		}
		bool operator< (const myPoint & p) const
		{
			if (x < p.x && y < p.y)
				return true;
			else
				return false;
		}
		float DistanceTo(const myPoint& p)
		{
			float mDx = (float)x - p.x;
			float mDy = (float)y - p.y;
			return sqrt(mDx*mDx + mDy*mDy);
		}
	};
	struct FPoint{
		float x;
		float y;
		FPoint(int a, int b)
		{
			x = float(a);
			y = float(b);
		}
		FPoint()
		{
			x = 0.0f;
			y = 0.0f;
		}
		void operator= (const FPoint & p)
		{
			x = p.x;
			y = p.y;
		}
		void normalize()
		{
			float l = (float)sqrt(x*x + y*y);
			x = x / l;
			y = y / l;
		}
		float length()
		{
			return sqrt(x*x + y*y);
		}
		void scale(float f)
		{
			x = x*f;
			y = y*f;
		}
	};
	struct Line{
		myPoint p0;
		myPoint p1;
		Line(int x0, int y0, int x1, int y1)
		{
			p0.x = x0;
			p0.y = y0;
			p1.x = x1;
			p1.y = y1;
		}
		Line(myPoint pa, myPoint pb)
		{
			p0.x = pa.x;
			p0.y = pa.y;
			p1.x = pb.x;
			p1.y = pb.y;
		}
		float length()
		{
			float dx = p1.x - p0.x;
			float dy = p1.y - p0.y;
			return sqrt(dx*dx + dy*dy);
		}
	};
	struct Vector2{
		// 因為除法會造成誤差
		// 所以用加法跟乘法來作一些向量的計算
		// 是比較好的做法
		float x, y;
		Vector2()
		{
			x = 0;
			y = 0;
		}
		Vector2(int x0, int y0, int x1, int y1)
		{
			x = x1 - x0;
			y = y1 - y0;
		}
		void operator= (const Vector2 & v)
		{
			x = v.x;
			y = v.y;
		}
		void operator= (const myPoint & p)
		{
			x = p.x;
			y = p.y;
		}
		void normalize()
		{
			float l = length();
			if (l == 0.0f)
			{
				return;
			}
			else
			{
				x = x / l;
				y = y / l;
			}
		}
		float norm()
		{
			return x*x + y*y;
		}
		float length()
		{
			return sqrt(x*x + y*y);
		}
		float radFromXaxis()
		{
			normalize();
			return atan2(x, y);
		}
		float degFromXaxis()
		{
			float rad = radFromXaxis();
			return rad*(180 / 3.141592653589793238);
		}
		static float dot(const Vector2& v0, const Vector2& v1)
		{
			return v0.x * v1.x + v0.y * v1.y;
		}
		static float cross(const Vector2& v0, const Vector2& v1)
		{
			// 外積運算，回傳純量（除去方向）
			// 二維的話就像是求det([v0 v1])
			return v0.x * v1.y - v0.y * v1.x;
		}
		static float distance(const myPoint& p, const Line& l)
		{
			// 面積除以底得高
			Vector2 v0, v1;
			v0.x = p.x - l.p0.x;
			v0.y = p.y - l.p0.y;
			v1.x = l.p1.x - l.p0.x;
			v1.y = l.p1.y - l.p0.y;

			return abs(cross(v0, v1)) / v1.length();
		}
		static myPoint project(const myPoint& pOut, const Line& l)
		{
			myPoint pProj;
			Vector2 v0, v1, v2;
			v0.x = pOut.x - l.p0.x;
			v0.y = pOut.y - l.p0.y;
			v1.x = l.p1.x - l.p0.x;
			v1.y = l.p1.y - l.p0.y;
			v2.x = v1.x * dot(v0, v1) / v1.norm();
			v2.y = v1.y * dot(v0, v1) / v1.norm();
			pProj.x = l.p0.x + v2.x;
			pProj.y = l.p0.y + v2.y;
			return pProj;
		}
	};
	template<class T>
	inline T Square(T x)
	{
		return x*x;
	}
	template<class T>
	T DistanceSquare(const T *a, const T *b)
	{
		T t = a[0] - b[0];
		T u = a[1] - b[1];
		T v = a[2] - b[2];
		return (t*t + u*u + v*v);
	}
	template<class T>
	T DistanceSquare2d(const T *a, const T *b)
	{
		T t = a[0] - b[0];
		T u = a[1] - b[1];
		return (t*t + u*u);
	}
	template<class T>
	T TriangleArea(const T *a, const T *b)
	{
		T t = a[0] * b[1] - a[1] * b[0];
		T u = a[1] * b[2] - a[2] * b[1];
		T v = a[2] * b[0] - a[0] * b[2];
		return 0.5f*sqrt(t*t + u*u + v*v);
	}
	template<class T>
	inline T TriangleArea2d(const T *a, const T *b)
	{
		return 0.5f*abs(a[0] * b[1] - a[1] * b[0]);
	}
	template<class T>
	inline T VectorLength(const T *a)
	{
		return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	}
	template<class T>
	inline T VectorLength2d(const T *a)
	{
		return sqrt(a[0] * a[0] + a[1] * a[1]);
	}
	template<class T>
	inline T Dot(const T *a, const T *b)
	{
		return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	}
	template<class T>
	inline T Dot2(const T *a, const T *b)
	{
		return a[0] * b[0] + a[1] * b[1];
	}
	template<class T>
	T Sign(T n)
	{
		if (n < 0)
			return -1;
		else
			return 1;
	}
	inline void Normalize(float *v)
	{
		float l = VectorLength(v);
		v[0] /= l;
		v[1] /= l;
		v[2] /= l;
	}
	inline void Normalize2d(float *v)
	{
		float l = VectorLength2d(v);
		v[0] /= l;
		v[1] /= l;
	}
	inline void Cross3X3(const float *a, const float*b, float *out)
	{
		out[0] = a[1] * b[2] - a[2] * b[1];
		out[1] = -(a[0] * b[2] - a[2] * b[0]);
		out[2] = a[0] * b[1] - a[1] * b[0];
		Normalize(out);
	}
	inline float RandomSign()
	{
		int r = rand() % 2;
		if (r == 0)
			return 1.0f;
		else
			return -1.0f;
	}
	inline void ComputeCenter(const vector<myPoint>* vecP, myPoint* rCenter)
	{
		int pSize = vecP->size();

		rCenter->x = 0;
		rCenter->y = 0;
		for (int i = 0; i<pSize; ++i)
		{
			rCenter->x += vecP->at(i).x;
			rCenter->y += vecP->at(i).y;
		}
		rCenter->x /= pSize;
		rCenter->y /= pSize;
	}
	inline bool IsVectorSameDir(const float* a, const float* b, const float* c)
	{
		//兩兩做內積, 若兩向量夾角 < 90度則為正號, 也就是同向
		//又我們得到的三個法向量, 扣掉小誤差, 只有兩種可能, 朝上或朝下, 所以可以用dot的正負號判斷
		//利用漸進式的判斷, 可讓速度變快
		float d = Dot(a, b);
		if (d > 0)
		{
			float e = Dot(b, c);
			if (e > 0)
			{
				float f = Dot(c, a);
				if (f > 0)
					return true;
			}
		}
		return false;
	}
	inline bool IsInTriangle(const float* va, const float* vb, const float* vc, const float* vp)
	{
		/**
		*      a
		*      \
		*     /|\
		*    / p \   /
		*   / / \ \ /
		b-------c
		*/
		float ap[3];
		float bp[3];
		float cp[3];

		float ab[3];
		float bc[3];
		float ca[3];

		float n0[3];
		float n1[3];
		float n2[3];

		ap[0] = vp[0] - va[0];
		ap[1] = vp[1] - va[1];
		ap[2] = vp[2] - va[2];
		bp[0] = vp[0] - vb[0];
		bp[1] = vp[1] - vb[1];
		bp[2] = vp[2] - vb[2];
		cp[0] = vp[0] - vc[0];
		cp[1] = vp[1] - vc[1];
		cp[2] = vp[2] - vc[2];

		ab[0] = vb[0] - va[0];
		ab[1] = vb[1] - va[1];
		ab[2] = vb[2] - va[2];
		bc[0] = vc[0] - vb[0];
		bc[1] = vc[1] - vb[1];
		bc[2] = vc[2] - vb[2];
		ca[0] = va[0] - vc[0];
		ca[1] = va[1] - vc[1];
		ca[2] = va[2] - vc[2];

		Cross3X3(ap, ab, n0);
		Cross3X3(bp, bc, n1);
		Cross3X3(cp, ca, n2);

		if (IsVectorSameDir(n0, n1, n2))
			return true;

		return false;
	}
	inline bool IsVec3Equal(const float* va, const float* vb)
	{
		if (va[0] == vb[0] && va[1] == vb[1] && va[2] == vb[2])
		{
			return true;
		}
		else
			return false;
	}
	template<class T>
	inline T pow2(T n)
	{
		return n*n;
	}
};