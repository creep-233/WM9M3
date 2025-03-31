//#pragma once
//
//#include "Core.h"
//#include "Sampling.h"
//#include "Geometry.h"
//#include "Imaging.h"
//#include "Materials.h"
//#include "Lights.h"
//
//class Camera
//{
//public:
//	Matrix projectionMatrix;
//	Matrix inverseProjectionMatrix;
//	Matrix camera;
//	Matrix cameraToView;
//	float width = 0;
//	float height = 0;
//	Vec3 origin;
//	Vec3 viewDirection;
//	float Afilm;
//	void init(Matrix ProjectionMatrix, int screenwidth, int screenheight)
//	{
//		projectionMatrix = ProjectionMatrix;
//		inverseProjectionMatrix = ProjectionMatrix.invert();
//		width = (float)screenwidth;
//		height = (float)screenheight;
//		float Wlens = (2.0f / ProjectionMatrix.a[1][1]);
//		float aspect = ProjectionMatrix.a[0][0] / ProjectionMatrix.a[1][1];
//		float Hlens = Wlens * aspect;
//		Afilm = Wlens * Hlens;
//	}
//	void updateView(Matrix V)
//	{
//		camera = V;
//		cameraToView = V.invert();
//		origin = camera.mulPoint(Vec3(0, 0, 0));
//		viewDirection = inverseProjectionMatrix.mulPointAndPerspectiveDivide(Vec3(0, 0, 1));
//		viewDirection = camera.mulVec(viewDirection);
//		viewDirection = viewDirection.normalize();
//	}
//	// Add code here
//	Ray generateRay(float x, float y)
//	{
//		
//		float xprime = x / width;
//		float yprime = 1.0f - (y / height);
//		xprime = (xprime * 2.0f) - 1.0f;
//		yprime = (yprime * 2.0f) - 1.0f;
//
//		Vec3 dir(xprime, yprime, 1.0f);
//
//		dir = inverseProjectionMatrix.mulPoint(dir);
//		dir = camera.mulVec(dir);
//		dir = dir.normalize();
//
//		return Ray(origin, dir);
//
//
//		//Vec3 dir(0, 0, 1);
//		//return Ray(origin, dir);
//	}
//
//	bool projectOntoCamera(const Vec3& p, float& x, float& y)
//	{
//		Vec3 pview = cameraToView.mulPoint(p);
//		Vec3 pproj = projectionMatrix.mulPointAndPerspectiveDivide(pview);
//		x = (pproj.x + 1.0f) * 0.5f;
//		y = (pproj.y + 1.0f) * 0.5f;
//		if (x < 0 || x > 1.0f || y < 0 || y > 1.0f)
//		{
//			return false;
//		}
//		x = x * width;
//		y = 1.0f - y;
//		y = y * height;
//		return true;
//	}
//
//};
//
//class Scene
//{
//public:
//	std::vector<Triangle> triangles;
//	std::vector<BSDF*> materials;
//	std::vector<Light*> lights;
//	Light* background = NULL;
//	BVHNode* bvh = NULL;
//	Camera camera;
//	AABB bounds;
//	void build()
//	{
//		// Add BVH building code here
//		
//		// Do not touch the code below this line!
//		// Build light list
//		for (int i = 0; i < triangles.size(); i++)
//		{
//			if (materials[triangles[i].materialIndex]->isLight())
//			{
//				AreaLight* light = new AreaLight();
//				light->triangle = &triangles[i];
//				light->emission = materials[triangles[i].materialIndex]->emission;
//				lights.push_back(light);
//			}
//		}
//	}
//	IntersectionData traverse(const Ray& ray)
//	{
//		IntersectionData intersection;
//		intersection.t = FLT_MAX;
//		for (int i = 0; i < triangles.size(); i++)
//		{
//			float t;
//			float u;
//			float v;
//			if (triangles[i].rayIntersect(ray, t, u, v))
//			{
//				if (t < intersection.t)
//				{
//					intersection.t = t;
//					intersection.ID = i;
//					intersection.alpha = u;
//					intersection.beta = v;
//					intersection.gamma = 1.0f - (u + v);
//				}
//			}
//		}
//		return intersection;
//	}
//	Light* sampleLight(Sampler* sampler, float& pmf)
//	{
//		//float r1 = sampler->next();
//		//pmf = 1.0f / (float)Lights.size();
//		//return lights[std::min(int)(r1 * Lights.size()), (int)(Lights.size() - 1))];
//
//			// 1. 生成随机数 r1
//		float r1 = sampler->next();
//
//		// 2. 确保光源数量不为 0
//		int numLights = lights.size();
//		if (numLights == 0) return nullptr;
//
//		// 3. 计算 PMF
//		pmf = 1.0f / numLights;
//
//		// 4. 计算光源索引，并确保索引不越界
//		int lightIndex = std::min((int)(r1 * numLights), numLights - 1);
//
//		// 5. 返回选中的光源
//		return lights[lightIndex];
//
//	}
//	// Do not modify any code below this line
//	void init(std::vector<Triangle> meshTriangles, std::vector<BSDF*> meshMaterials, Light* _background)
//	{
//		for (int i = 0; i < meshTriangles.size(); i++)
//		{
//			triangles.push_back(meshTriangles[i]);
//			bounds.extend(meshTriangles[i].vertices[0].p);
//			bounds.extend(meshTriangles[i].vertices[1].p);
//			bounds.extend(meshTriangles[i].vertices[2].p);
//		}
//		for (int i = 0; i < meshMaterials.size(); i++)
//		{
//			materials.push_back(meshMaterials[i]);
//		}
//		background = _background;
//		if (background->totalIntegratedPower() > 0)
//		{
//			lights.push_back(background);
//		}
//	}
//	bool visible(const Vec3& p1, const Vec3& p2)
//	{
//		Ray ray;
//		Vec3 dir = p2 - p1;
//		float maxT = dir.length() - (2.0f * EPSILON);
//		dir = dir.normalize();
//		ray.init(p1 + (dir * EPSILON), dir);
//		
//		
//		//return bvh->traverseVisible(ray, triangles, maxT);
//
//		for (int i = 0; i < triangles.size(); i++)
//		{
//			float t;
//			float u;
//			float v;
//			if (triangles[i].rayIntersect(ray, t, u, v))
//			{
//				if (t < maxT)
//				{
//					return false;
//				}
//			}
//		}
//		return true;
//
//	}
//	Colour emit(Triangle* light, ShadingData shadingData, Vec3 wi)
//	{
//		return materials[light->materialIndex]->emit(shadingData, wi);
//	}
//	ShadingData calculateShadingData(IntersectionData intersection, Ray& ray)
//	{
//		ShadingData shadingData = {};
//		if (intersection.t < FLT_MAX)
//		{
//			shadingData.x = ray.at(intersection.t);
//			shadingData.gNormal = triangles[intersection.ID].gNormal();
//			triangles[intersection.ID].interpolateAttributes(intersection.alpha, intersection.beta, intersection.gamma, shadingData.sNormal, shadingData.tu, shadingData.tv);
//			shadingData.bsdf = materials[triangles[intersection.ID].materialIndex];
//			shadingData.wo = -ray.dir;
//			if (shadingData.bsdf->isTwoSided())
//			{
//				if (Dot(shadingData.wo, shadingData.sNormal) < 0)
//				{
//					shadingData.sNormal = -shadingData.sNormal;
//				}
//				if (Dot(shadingData.wo, shadingData.gNormal) < 0)
//				{
//					shadingData.gNormal = -shadingData.gNormal;
//				}
//			}
//			shadingData.frame.fromVector(shadingData.sNormal);
//			shadingData.t = intersection.t;
//		} else
//		{
//			shadingData.wo = -ray.dir;
//			shadingData.t = intersection.t;
//		}
//		return shadingData;
//	}
//};






#pragma once

#include "Core.h"
#include "Sampling.h"
#include "Geometry.h"
#include "Imaging.h"
#include "Materials.h"
#include "Lights.h"

class Camera
{
public:
	Matrix projectionMatrix;
	Matrix inverseProjectionMatrix;
	Matrix camera;
	Matrix cameraToView;
	float width = 0;
	float height = 0;
	Vec3 origin;
	void init(Matrix ProjectionMatrix, int screenwidth, int screenheight)
	{
		projectionMatrix = ProjectionMatrix;
		inverseProjectionMatrix = ProjectionMatrix.invert();
		width = (float)screenwidth;
		height = (float)screenheight;
	}
	void updateView(Matrix V)
	{
		camera = V;
		cameraToView = V.invert();
		origin = camera.mulPoint(Vec3(0, 0, 0));
	}
	// Add code here
	//Ray generateRay(float x, float y)
	//{
	//	Vec3 dir(0, 0, 1);
	//	return Ray(origin, dir);
	//}
	Ray generateRay(float x, float y)
	{
		
		float xprime = x / width;
		float yprime = 1.0f - (y / height);
		xprime = (xprime * 2.0f) - 1.0f;
		yprime = (yprime * 2.0f) - 1.0f;

		Vec3 dir(xprime, yprime, 1.0f);

		dir = inverseProjectionMatrix.mulPoint(dir);
		dir = camera.mulVec(dir);
		dir = dir.normalize();

		return Ray(origin, dir);


		//Vec3 dir(0, 0, 1);
		//return Ray(origin, dir);
	}
	bool projectOntoCamera(const Vec3& p, float& x, float& y)
	{
		Vec3 pview = cameraToView.mulPoint(p);
		Vec3 pproj = projectionMatrix.mulPointAndPerspectiveDivide(pview);
		x = (pproj.x + 1.0f) * 0.5f;
		y = (pproj.y + 1.0f) * 0.5f;
		if (x < 0 || x > 1.0f || y < 0 || y > 1.0f)
		{
			return false;
		}
		x = x * width;
		y = 1.0f - y;
		y = y * height;
		return true;
	}
};

class Scene
{
public:
	std::vector<Triangle> triangles;
	std::vector<Triangle> originalTriangles; // 你的原始 mesh
	std::vector<Triangle> bvhTriangles;      // 要传给 BVH 并由其内部使用的副本

	std::vector<BSDF*> materials;
	std::vector<Light*> lights;
	Light* background = NULL;
	BVHNode* bvh = NULL;
	Camera camera;
	AABB bounds;
	void build()
	{
		// Add BVH building code here


		if (!bvh)
			bvh = new BVHNode();

		// 调用新的 build 函数（注意，这里要声明一个 trianglesCopy 用来排序）
		bvh->build(triangles, triangles); // triangles 作为原始数据和副本都传进去

		// Do not touch the code below this line!
		// Build light list
		for (int i = 0; i < triangles.size(); i++)
		{
			if (materials[triangles[i].materialIndex]->isLight())
			{
				AreaLight* light = new AreaLight();
				light->triangle = &triangles[i];
				light->emission = materials[triangles[i].materialIndex]->emission;
				lights.push_back(light);
			}
		}
	}
	//IntersectionData traverse(const Ray& ray)
	//{
	//	IntersectionData intersection;
	//	intersection.t = FLT_MAX;
	//	for (int i = 0; i < triangles.size(); i++)
	//	{
	//		float t;
	//		float u;
	//		float v;
	//		if (triangles[i].rayIntersect(ray, t, u, v))
	//		{
	//			if (t < intersection.t)
	//			{
	//				intersection.t = t;
	//				intersection.ID = i;
	//				intersection.alpha = u;
	//				intersection.beta = v;
	//				intersection.gamma = 1.0f - (u + v);
	//			}
	//		}
	//	}
	//	return intersection;
	//}
	// 
	
	IntersectionData traverse(const Ray& ray)
	{
		if (bvh)
			return bvh->traverse(ray, triangles); 

		// fallback
		IntersectionData intersection;
		intersection.t = FLT_MAX;
		for (int i = 0; i < triangles.size(); i++)
		{
			float t, u, v;
			if (triangles[i].rayIntersect(ray, t, u, v) && t < intersection.t)
			{
				intersection.t = t;
				intersection.ID = i;
				intersection.alpha = u;
				intersection.beta = v;
				intersection.gamma = 1.0f - u - v;
			}
		}
		return intersection;
	}


	//Light* sampleLight(Sampler* sampler, float& pmf)
	//{
	//	return NULL;
	//}
	// 
	 Light* sampleLight(Sampler* sampler, float& pmf)
	{
		//float r1 = sampler->next();
		//pmf = 1.0f / (float)Lights.size();
		//return lights[std::min(int)(r1 * Lights.size()), (int)(Lights.size() - 1))];

			// 1. 生成随机数 r1
		float r1 = sampler->next();

		// 2. 确保光源数量不为 0
		int numLights = lights.size();
		if (numLights == 0) return nullptr;

		// 3. 计算 PMF
		pmf = 1.0f / numLights;

		// 4. 计算光源索引，并确保索引不越界
		int lightIndex = std::min((int)(r1 * numLights), numLights - 1);

		// 5. 返回选中的光源
		return lights[lightIndex];

	}
	// Do not modify any code below this line
	void init(std::vector<Triangle> meshTriangles, std::vector<BSDF*> meshMaterials, Light* _background)
	{
		for (int i = 0; i < meshTriangles.size(); i++)
		{
			triangles.push_back(meshTriangles[i]);
			bounds.extend(meshTriangles[i].vertices[0].p);
			bounds.extend(meshTriangles[i].vertices[1].p);
			bounds.extend(meshTriangles[i].vertices[2].p);
		}
		for (int i = 0; i < meshMaterials.size(); i++)
		{
			materials.push_back(meshMaterials[i]);
		}
		background = _background;
		if (background->totalIntegratedPower() > 0)
		{
			lights.push_back(background);
		}
	}
	bool visible(const Vec3& p1, const Vec3& p2)
	{
		Ray ray;
		Vec3 dir = p2 - p1;
		float maxT = dir.length() - (2.0f * EPSILON);
		dir = dir.normalize();
		ray.init(p1 + (dir * EPSILON), dir);
		return bvh->traverseVisible(ray, triangles, maxT);
		//for (int i = 0; i < triangles.size(); i++)
		//{
		//	float t;
		//	float u;
		//	float v;
		//	if (triangles[i].rayIntersect(ray, t, u, v))
		//	{
		//		if (t < maxT)
		//		{
		//			return false;
		//		}
		//	}
		//}
		//return true;
	}
	Colour emit(Triangle* light, ShadingData shadingData, Vec3 wi)
	{
		return materials[light->materialIndex]->emit(shadingData, wi);
	}
	ShadingData calculateShadingData(IntersectionData intersection, Ray& ray)
	{
		ShadingData shadingData = {};
		if (intersection.t < FLT_MAX)
		{
			shadingData.x = ray.at(intersection.t);
			shadingData.gNormal = triangles[intersection.ID].gNormal();
			triangles[intersection.ID].interpolateAttributes(intersection.alpha, intersection.beta, intersection.gamma, shadingData.sNormal, shadingData.tu, shadingData.tv);
			shadingData.bsdf = materials[triangles[intersection.ID].materialIndex];
			shadingData.wo = -ray.dir;
			if (shadingData.bsdf->isTwoSided())
			{
				if (Dot(shadingData.wo, shadingData.sNormal) < 0)
				{
					shadingData.sNormal = -shadingData.sNormal;
				}
				if (Dot(shadingData.wo, shadingData.gNormal) < 0)
				{
					shadingData.gNormal = -shadingData.gNormal;
				}
			}
			shadingData.frame.fromVector(shadingData.sNormal);
			shadingData.t = intersection.t;
		}
		else
		{
			shadingData.wo = -ray.dir;
			shadingData.t = intersection.t;
		}
		return shadingData;
	}

	Light* getLightFromHit(const IntersectionData& isect)
	{
		if (isect.t >= FLT_MAX) return nullptr;
		for (auto* light : lights)
		{
			if (!light->isArea()) continue;

			AreaLight* area = dynamic_cast<AreaLight*>(light);
			if (area && area->triangle == &triangles[isect.ID])
			{
				return light;
			}
		}
		return nullptr;
	}
};