#pragma once

#include "Core.h"
#include "Sampling.h"

class Ray
{
public:
	Vec3 o;
	Vec3 dir;
	Vec3 invDir;
	Ray()
	{
	}
	Ray(Vec3 _o, Vec3 _d)
	{
		init(_o, _d);
	}
	void init(Vec3 _o, Vec3 _d)
	{
		o = _o;
		dir = _d;
		invDir = Vec3(1.0f / dir.x, 1.0f / dir.y, 1.0f / dir.z);
	}
	Vec3 at(const float t) const
	{
		return (o + (dir * t));
	}
};

class Plane
{
	
	Ray ray;

public:
	Vec3 n;
	float d;
	void init(Vec3& _n, float _d)
	{
		n = _n;
		d = _d;
	}
	// Add code here
	bool rayIntersect(Ray& r, float& t)
	{

		float denom = Dot(n, r.dir); 

		if (fabs(denom) < 1e-6f) {
			return false;
		}

		t = -(Dot(n, r.o) + d) / denom; 

		return (t >= 0); 
	}

};

#define EPSILON 0.001f

class Triangle
{
public:
	Vertex vertices[3];
	Vec3 e1; // Edge 1
	Vec3 e2; // Edge 2
	Vec3 n; // Geometric Normal
	float area; // Triangle area
	float d; // For ray triangle if needed
	unsigned int materialIndex;
	void init(Vertex v0, Vertex v1, Vertex v2, unsigned int _materialIndex)
	{
		materialIndex = _materialIndex;
		vertices[0] = v0;
		vertices[1] = v1;
		vertices[2] = v2;
		e1 = vertices[2].p - vertices[1].p;
		e2 = vertices[0].p - vertices[2].p;
		n = e1.cross(e2).normalize();
		area = e1.cross(e2).length() * 0.5f;
		d = Dot(n, vertices[0].p);
	}
	Vec3 centre() const
	{
		return (vertices[0].p + vertices[1].p + vertices[2].p) / 3.0f;
	}
	// Add code here
	bool rayIntersect(const Ray& r, float& t, float& u, float& v) const
	{
		
		//float denom = Dot(n, r.dir);
		//if (denom == 0) { return false; }
		//t = (d - Dot(n, r.o)) / denom;
		//if (t < 0) { return false; }
		//Vec3 p = r.at(t);
		//float invArea = 1.0f / Dot(e1.cross(e2), n);
		//u = Dot(e1.cross(p - vertices[1].p), n) * invArea;
		//if (u < 0 || u > 1.0f) { return false; }
		//v = Dot(e2.cross(p - vertices[2].p), n) * invArea;
		//if (v < 0 || (u + v) > 1.0f) { return false; }
		//return true;

		const Vec3& p0 = vertices[0].p;
		const Vec3& p1 = vertices[1].p;
		const Vec3& p2 = vertices[2].p;

		Vec3 edge1 = p1 - p0;
		Vec3 edge2 = p2 - p0;
		Vec3 h = r.dir.cross(edge2); 
		float a = Dot(edge1, h);
		if (fabs(a) < EPSILON) return false;

		float f = 1.0f / a;
		Vec3 s = r.o - p0;
		u = f * Dot(s, h);
		if (u < 0.0f || u > 1.0f) return false;

		Vec3 q = s.cross(edge1);
		v = f * Dot(r.dir, q);
		if (v < 0.0f || u + v > 1.0f) return false;

		t = f * Dot(edge2, q);
		if (t > EPSILON)
			return true;

		return false;

	}
	void interpolateAttributes(const float alpha, const float beta, const float gamma, Vec3& interpolatedNormal, float& interpolatedU, float& interpolatedV) const
	{
		interpolatedNormal = vertices[0].normal * alpha + vertices[1].normal * beta + vertices[2].normal * gamma;
		interpolatedNormal = interpolatedNormal.normalize();
		interpolatedU = vertices[0].u * alpha + vertices[1].u * beta + vertices[2].u * gamma;
		interpolatedV = vertices[0].v * alpha + vertices[1].v * beta + vertices[2].v * gamma;
	}

	// Add code here
	//Vec3 sample(Sampler* sampler, float& pdf)
	//{
	//	float r1 = sample->next();
	//	//return Vec3(0, 0, 0);
	//}

	Vec3 sample(Sampler* sampler, float& pdf)
	{
		float r1 = sampler->next();
		float r2 = sampler->next();

		float sqrtR1 = sqrt(r1);
		float b0 = 1.0f - sqrtR1;
		float b1 = r2 * sqrtR1;
		float b2 = 1.0f - b0 - b1;

		pdf = 1.0f / area;

		return vertices[0].p * b0 + vertices[1].p * b1 + vertices[2].p * b2;
	}

	Vec3 gNormal()
	{
		return (n * (Dot(vertices[0].normal, n) > 0 ? 1.0f : -1.0f));
	}
};


class AABB
{
public:
	Vec3 max;
	Vec3 min;
	AABB()
	{
		reset();
	}
	void reset()
	{
		max = Vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
		min = Vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	}
	void extend(const Vec3 p)
	{
		max = Max(max, p);
		min = Min(min, p);
	}


	bool rayAABB(const Ray& r, float& t) {
		float tMin, tMax;
		bool hit = rayIntersect(r, tMin, tMax);
		t = tMin;
		return hit;
	}

	bool rayAABB(const Ray& r) {
		float tMin, tMax;
		return rayIntersect(r, tMin, tMax);
	}


	// Add code here

	bool rayIntersect(const Ray& ray, float& tMin, float& tMax) const {
		tMin = 0.0f;
		tMax = FLT_MAX;

		for (int i = 0; i < 3; ++i) {
			float invD = 1.0f / ray.dir[i];
			float t0 = (min[i] - ray.o[i]) * invD;
			float t1 = (max[i] - ray.o[i]) * invD;

			if (invD < 0.0f)
				std::swap(t0, t1);

			tMin = std::max(tMin, t0);
			tMax = std::min(tMax, t1);

			if (tMax < tMin)
				return false;
		}
		return true;
	}

	bool rayAABB(const Ray& r, float& tMin, float& tMax)
	{
		return rayIntersect(r, tMin, tMax);
	}


	void extend(const AABB& box)
	{
		extend(box.min);
		extend(box.max);
	}



	float area()
	{
		Vec3 size = max - min;
		return ((size.x * size.y) + (size.y * size.z) + (size.x * size.z)) * 2.0f;
	}

};

class Sphere
{
public:
	Vec3 centre;
	float radius;
	void init(Vec3& _centre, float _radius)
	{
		centre = _centre;
		radius = _radius;
	}
	// Add code here
	bool rayIntersect(Ray& r, float& t)
	{
		//return false;
		Vec3 L = centre - r.o;
		float tca = Dot(L, r.dir);
		if (tca < 0) return false;

		float d2 = Dot(L, L) - tca * tca;
		float radius2 = radius * radius;
		if (d2 > radius2) return false;

		float thc = sqrt(radius2 - d2);
		float t0 = tca - thc;
		float t1 = tca + thc;

		if (t0 > EPSILON) {
			t = t0;
			return true;
		}
		else if (t1 > EPSILON) {
			t = t1;
			return true;
		}

		return false;
	}
};

struct IntersectionData
{
	unsigned int ID;
	float t;
	float alpha;
	float beta;
	float gamma;
};

#define MAXNODE_TRIANGLES 8
#define TRAVERSE_COST 1.0f
#define TRIANGLE_COST 2.0f
#define BUILD_BINS 32


class BVHNode
{
public:
	AABB bounds;
	BVHNode* l;
	BVHNode* r;

	int startIndex;
	int endIndex;

	BVHNode()
		: l(nullptr)
		, r(nullptr)
		, startIndex(0)
		, endIndex(0)
	{
	}

	static AABB triangleBounds(const Triangle& tri)
	{
		AABB box;
		box.reset();
		box.extend(tri.vertices[0].p);
		box.extend(tri.vertices[1].p);
		box.extend(tri.vertices[2].p);
		return box;
	}

	void buildRecursive(std::vector<Triangle>& triangles, int start, int end)
	{
		bounds.reset();
		for (int i = start; i < end; i++) {
			bounds.extend(triangles[i].vertices[0].p);
			bounds.extend(triangles[i].vertices[1].p);
			bounds.extend(triangles[i].vertices[2].p);
		}

		int numTriangles = end - start;
		if (numTriangles <= MAXNODE_TRIANGLES)
		{
			this->startIndex = start;
			this->endIndex = end;
			return;
		}


		Vec3 size = bounds.max - bounds.min;
		int axis = 0;
		if (size.y > size.x && size.y > size.z)
			axis = 1;
		else if (size.z > size.x && size.z > size.y)
			axis = 2;

		std::sort(triangles.begin() + start, triangles.begin() + end,
			[axis](const Triangle& a, const Triangle& b) {
				float ca = (axis == 0) ? a.centre().x : (axis == 1) ? a.centre().y : a.centre().z;
				float cb = (axis == 0) ? b.centre().x : (axis == 1) ? b.centre().y : b.centre().z;
				return ca < cb;
			});

		int n = numTriangles;
		std::vector<AABB> leftBounds(n);
		std::vector<AABB> rightBounds(n);
		leftBounds[0] = triangleBounds(triangles[start]);
		for (int i = 1; i < n; i++) {
			leftBounds[i] = leftBounds[i - 1];
			leftBounds[i].extend(triangleBounds(triangles[start + i]));
		}
		rightBounds[n - 1] = triangleBounds(triangles[end - 1]);
		for (int i = n - 2; i >= 0; i--) {
			rightBounds[i] = rightBounds[i + 1];
			rightBounds[i].extend(triangleBounds(triangles[start + i]));
		}

		float totalArea = bounds.area();
		float bestCost = FLT_MAX;
		int bestSplit = -1;

		for (int i = 1; i < n; i++) {
			float leftArea = leftBounds[i - 1].area();
			float rightArea = rightBounds[i].area();
			int leftCount = i;
			int rightCount = n - i;
			float cost = 1.0f + (leftArea * leftCount + rightArea * rightCount) / totalArea;
			if (cost < bestCost) {
				bestCost = cost;
				bestSplit = i;
			}
		}

		if (bestCost >= static_cast<float>(numTriangles)) {
			this->startIndex = start;
			this->endIndex = end;
			return;
		}

		int mid = start + bestSplit;
		l = new BVHNode();
		r = new BVHNode();
		l->buildRecursive(triangles, start, mid);
		r->buildRecursive(triangles, mid, end);
	}

	void build(std::vector<Triangle>& inputTriangles, std::vector<Triangle>& triangles)
	{
		triangles = inputTriangles;
		buildRecursive(triangles, 0, static_cast<int>(triangles.size()));
	}

	void traverse(const Ray& ray, const std::vector<Triangle>& triangles, IntersectionData& intersection)
	{
		float tBox;
		if (!bounds.rayAABB(ray, tBox)) {
			return;
		}

		if (!l && !r)
		{
			for (int i = startIndex; i < endIndex; i++)
			{
				float t, u, v;
				if (triangles[i].rayIntersect(ray, t, u, v) && t > 1e-4f && t < intersection.t)
				{
					intersection.t = t;
					intersection.alpha = u;
					intersection.beta = v;
					intersection.gamma = 1 - u - v;
					intersection.ID = i;
				}
			}
			return;
		}

		if (l) l->traverse(ray, triangles, intersection);
		if (r) r->traverse(ray, triangles, intersection);
	}

	IntersectionData traverse(const Ray& ray, const std::vector<Triangle>& triangles)
	{
		IntersectionData intersection;
		intersection.t = FLT_MAX;
		traverse(ray, triangles, intersection);
		return intersection;
	}

	bool traverseVisible(const Ray& ray, const std::vector<Triangle>& triangles, float maxT)
	{
		float tBox;
		float tMin, tMax;
		if (!bounds.rayIntersect(ray, tMin, tMax)) {
			return true;
		}

		if (tMin > maxT) {
			return true;
		}

		if (!l && !r)
		{
			for (int i = startIndex; i < endIndex; i++)
			{
				float t, u, v;
				if (triangles[i].rayIntersect(ray, t, u, v) && t > 1e-4f && t < maxT)
				{
					return false;
				}
			}
			return true;
		}

		bool leftVis = l ? l->traverseVisible(ray, triangles, maxT) : true;
		bool rightVis = r ? r->traverseVisible(ray, triangles, maxT) : true;
		return leftVis && rightVis;
	}

};
