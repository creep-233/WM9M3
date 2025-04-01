#pragma once

#include "Core.h"
#include "Sampling.h"
#include "Geometry.h"
#include "Imaging.h"
#include "Materials.h"
#include "Lights.h"

const int MAX_DEPTH = 5;


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

	float Afilm;
	Vec3 viewDirection;

	void init(Matrix ProjectionMatrix, int screenwidth, int screenheight)
	{
		projectionMatrix = ProjectionMatrix;
		inverseProjectionMatrix = ProjectionMatrix.invert();
		width = (float)screenwidth;
		height = (float)screenheight;

		Afilm = width * height;
	}
	void updateView(Matrix V)
	{
		camera = V;
		cameraToView = V.invert();
		origin = camera.mulPoint(Vec3(0, 0, 0));

		viewDirection = camera.mulVec(Vec3(0, 0, -1)).normalize();
	}
	// Add code here

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
	std::vector<Triangle> originalTriangles; 
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


		bvh->build(triangles, triangles); 

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

		float r1 = sampler->next();
		int numLights = lights.size();
		if (numLights == 0) return nullptr;
		pmf = 1.0f / numLights;
		int lightIndex = std::min((int)(r1 * numLights), numLights - 1);
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


	struct VPL {
		ShadingData shadingData;
		Colour Le;
	};

	std::vector<VPL> vpls;

	void traceVPLs(Sampler* sampler, int N_VPLs)
	{
		vpls.clear();

		for (int i = 0; i < N_VPLs; ++i)
		{
			float pmf;
			Light* light = sampleLight(sampler, pmf);
			if (!light || pmf <= 0.0f) continue;

			float pdfPos;
			Vec3 p = light->samplePositionFromLight(sampler, pdfPos);

			float pdfDir;
			Vec3 wi = light->sampleDirectionFromLight(sampler, pdfDir);

	
			ShadingData fakeSD;
			fakeSD.x = p;

			Vec3 n = light->normal(fakeSD, wi);
			if (Dot(wi, n) <= 0.0f) continue;

			Colour Le = light->evaluate(-wi);
			if (Le.Lum() == 0.0f) continue;

			float weight = pmf * pdfPos * pdfDir;
			if (weight <= EPSILON) continue;

			Colour initLe = Le * Dot(wi, n);
			initLe = initLe / (weight * float(N_VPLs) + EPSILON);

			Ray ray(p + wi * EPSILON, wi);
			Colour pathThroughput = initLe;

			traceVPLPath(ray, pathThroughput, sampler, 0);
		}
	}




	void traceVPLPath(Ray r, Colour pathThroughput, Sampler* sampler, int depth)
	{
		if (depth > MAX_DEPTH) return;

		IntersectionData isect = traverse(r);
		if (isect.t == FLT_MAX) return;

		ShadingData sd = calculateShadingData(isect, r);
		if (sd.bsdf == nullptr || sd.bsdf->isPureSpecular()) return; 

		VPL vpl;
		vpl.shadingData = sd;
		vpl.Le = pathThroughput;
		vpls.push_back(vpl);

		// Russian Roulette
		float pRR = std::min(0.9f, pathThroughput.Lum());
		if (sampler->next() > pRR) return;
		pathThroughput = pathThroughput / pRR;

		// Sample next direction
		float pdf;
		Colour f;
		Vec3 wi = sd.bsdf->sample(sd, sampler, f, pdf);
		if (pdf <= 0.0f || f.Lum() <= 0.0f) return;

		float cosTheta = std::max(Dot(wi, sd.sNormal), 0.0f);
		pathThroughput = pathThroughput * f * cosTheta / pdf;

		r.init(sd.x + wi * EPSILON, wi);
		traceVPLPath(r, pathThroughput, sampler, depth + 1);
	}
};