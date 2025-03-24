#pragma once

#include "Core.h"
#include "Sampling.h"
#include "Geometry.h"
#include "Imaging.h"
#include "Materials.h"
#include "Lights.h"
#include "Scene.h"
#include "GamesEngineeringBase.h"
#include <thread>
#include <functional>

const int MAX_DEPTH = 5;

class RayTracer
{
public:
	Scene* scene;
	GamesEngineeringBase::Window* canvas;
	Film* film;
	MTRandom *samplers;
	std::thread **threads;
	int numProcs;
	void init(Scene* _scene, GamesEngineeringBase::Window* _canvas)
	{
		scene = _scene;
		canvas = _canvas;
		film = new Film();
		film->init((unsigned int)scene->camera.width, (unsigned int)scene->camera.height, /*new BoxFilter()*/ new GaussianFilter());
		SYSTEM_INFO sysInfo;
		GetSystemInfo(&sysInfo);
		numProcs = sysInfo.dwNumberOfProcessors;
		threads = new std::thread*[numProcs];
		samplers = new MTRandom[numProcs];
		clear();
	}
	void clear()
	{
		film->clear();
	}



	Colour computeDirect(ShadingData shadingData, Sampler* sampler)
	{
		if (shadingData.bsdf->isPureSpecular() == true)
		{
			return Colour(0.0f, 0.0f, 0.0f);
		}
		// Sample a light
		float pmf;
		Light* light = scene->sampleLight(sampler, pmf);
		// Sample a point on the light
		float pdf;
		Colour emitted;
		Vec3 p = light->sample(shadingData, sampler, emitted, pdf);
		if (light->isArea())
		{
			// Calculate GTerm
			Vec3 wi = p - shadingData.x;
			float l = wi.lengthSq();
			wi = wi.normalize();
			float GTerm = (max(Dot(wi, shadingData.sNormal), 0.0f) * max(-Dot(wi, light->normal(shadingData, wi)), 0.0f)) / l;
			if (GTerm > 0)
			{
				// Trace
				if (scene->visible(shadingData.x, p))
				{
					// Shade
					return shadingData.bsdf->evaluate(shadingData, wi) * emitted * GTerm / (pmf * pdf);
				}
			}
		}
		else
		{
			// Calculate GTerm
			Vec3 wi = p;
			float GTerm = max(Dot(wi, shadingData.sNormal), 0.0f);
			if (GTerm > 0)
			{
				// Trace
				if (scene->visible(shadingData.x, shadingData.x + (p * 10000.0f)))
				{
					// Shade
					return shadingData.bsdf->evaluate(shadingData, wi) * emitted * GTerm / (pmf * pdf);
				}
			}
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}
	Colour pathTrace(Ray& r, Colour& pathThroughput, int depth, Sampler* sampler, bool canHitLight = true)
	{
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
			{
				if (canHitLight == true)
				{
					return pathThroughput * shadingData.bsdf->emit(shadingData, shadingData.wo);
				}
				else
				{
					return Colour(0.0f, 0.0f, 0.0f);
				}
			}
			Colour direct = pathThroughput * computeDirect(shadingData, sampler);
			if (depth > MAX_DEPTH)
			{
				return direct;
			}
			float russianRouletteProbability = min(pathThroughput.Lum(), 0.9f);
			if (sampler->next() < russianRouletteProbability)
			{
				pathThroughput = pathThroughput / russianRouletteProbability;
			}
			else
			{
				return direct;
			}
			Colour bsdf;
			float pdf;
			Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
			pdf = SamplingDistributions::cosineHemispherePDF(wi);
			wi = shadingData.frame.toWorld(wi);
			bsdf = shadingData.bsdf->evaluate(shadingData, wi);
			pathThroughput = pathThroughput * bsdf * fabsf(Dot(wi, shadingData.sNormal)) / pdf;
			r.init(shadingData.x + (wi * EPSILON), wi);
			return (direct + pathTrace(r, pathThroughput, depth + 1, sampler, shadingData.bsdf->isPureSpecular()));
		}
		return scene->background->evaluate(shadingData, r.dir);
	}
	Colour direct(Ray& r, Sampler* sampler)
	{
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
			{
				return shadingData.bsdf->emit(shadingData, shadingData.wo);
			}
			return computeDirect(shadingData, sampler);
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}

	Colour albedo(Ray& r)
	{
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
			{
				return shadingData.bsdf->emit(shadingData, shadingData.wo);
			}
			return shadingData.bsdf->evaluate(shadingData, Vec3(0, 1, 0));
		}
		return scene->background->evaluate(shadingData, r.dir);
	}
	Colour viewNormals(Ray& r)
	{
		IntersectionData intersection = scene->traverse(r);
		if (intersection.t < FLT_MAX)
		{
			ShadingData shadingData = scene->calculateShadingData(intersection, r);
			return Colour(fabsf(shadingData.sNormal.x), fabsf(shadingData.sNormal.y), fabsf(shadingData.sNormal.z));
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}
	void render()
	{
		film->incrementSPP();
		for (unsigned int y = 0; y < film->height; y++)
		{
			for (unsigned int x = 0; x < film->width; x++)
			{
				float px = x + 0.5f;
				float py = y + 0.5f;
				Ray ray = scene->camera.generateRay(px, py);
				//Colour col = viewNormals(ray);
				//Colour col = albedo(ray);
				Colour col = direct(ray, samplers);

				//int threadId = 0; // �������ʱֻ�õ��߳�
				//Colour col = direct(ray, &samplers[threadId]);

				//Colour throughput(1.0f, 1.0f, 1.0f);
				//Colour col = pathTrace(ray, throughput, 0, &samplers[threadId]);

				film->splat(px, py, col);
				unsigned char r = (unsigned char)(col.r * 255);
				unsigned char g = (unsigned char)(col.g * 255);
				unsigned char b = (unsigned char)(col.b * 255);

				film->tonemap(x, y, r, g, b);

				canvas->draw(x, y, r, g, b);
			}
		}
	}
	int getSPP()
	{
		return film->SPP;
	}
	void saveHDR(std::string filename)
	{
		film->save(filename);
	}
	void savePNG(std::string filename)
	{
		stbi_write_png(filename.c_str(), canvas->getWidth(), canvas->getHeight(), 3, canvas->getBackBuffer(), canvas->getWidth() * 3);
	}
};

//#pragma once
//
//#include "Core.h"
//#include "Sampling.h"
//#include "Geometry.h"
//#include "Imaging.h"
//#include "Materials.h"
//#include "Lights.h"
//#include "Scene.h"
//#include "GamesEngineeringBase.h"
//#include <thread>
//#include <functional>
//
//const int MAX_DEPTH = 5;
//
//class RayTracer
//{
//public:
//	Scene* scene;
//	GamesEngineeringBase::Window* canvas;
//	Film* film;
//	MTRandom* samplers;
//	std::thread** threads;
//	int numProcs;
//	void init(Scene* _scene, GamesEngineeringBase::Window* _canvas)
//	{
//		scene = _scene;
//		canvas = _canvas;
//		film = new Film();
//		film->init((unsigned int)scene->camera.width, (unsigned int)scene->camera.height, new GaussianFilter());
//		SYSTEM_INFO sysInfo;
//		GetSystemInfo(&sysInfo);
//		numProcs = sysInfo.dwNumberOfProcessors;
//		threads = new std::thread * [numProcs];
//		samplers = new MTRandom[numProcs];
//		clear();
//	}
//	void clear()
//	{
//		film->clear();
//	}
//
//	Colour computeDirect(ShadingData shadingData, Sampler* sampler)
//	{
//		if (shadingData.bsdf->isPureSpecular())
//			return Colour(0.0f, 0.0f, 0.0f);
//		float pmf;
//		Light* light = scene->sampleLight(sampler, pmf);
//		float pdf;
//		Colour emitted;
//		Vec3 p = light->sample(shadingData, sampler, emitted, pdf);
//		if (light->isArea())
//		{
//			Vec3 wi = p - shadingData.x;
//			float l = wi.lengthSq();
//			wi = wi.normalize();
//			float GTerm = (max(Dot(wi, shadingData.sNormal), 0.0f) * max(-Dot(wi, light->normal(shadingData, wi)), 0.0f)) / l;
//			if (GTerm > 0 && scene->visible(shadingData.x, p))
//				return shadingData.bsdf->evaluate(shadingData, wi) * emitted * GTerm / (pmf * pdf);
//		}
//		else
//		{
//			Vec3 wi = p;
//			float GTerm = max(Dot(wi, shadingData.sNormal), 0.0f);
//			if (GTerm > 0 && scene->visible(shadingData.x, shadingData.x + (p * 10000.0f)))
//				return shadingData.bsdf->evaluate(shadingData, wi) * emitted * GTerm / (pmf * pdf);
//		}
//		return Colour(0.0f, 0.0f, 0.0f);
//	}
//
//	Colour pathTrace(Ray& r, Colour& pathThroughput, int depth, Sampler* sampler, bool canHitLight = true)
//	{
//		IntersectionData intersection = scene->traverse(r);
//		ShadingData shadingData = scene->calculateShadingData(intersection, r);
//		if (shadingData.t < FLT_MAX)
//		{
//			if (shadingData.bsdf->isLight())
//				return canHitLight ? pathThroughput * shadingData.bsdf->emit(shadingData, shadingData.wo) : Colour(0.0f);
//			Colour direct = pathThroughput * computeDirect(shadingData, sampler);
//			if (depth > MAX_DEPTH) return direct;
//			float rr = min(pathThroughput.Lum(), 0.9f);
//			if (sampler->next() >= rr) return direct;
//			pathThroughput = pathThroughput / rr;
//			Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
//			float pdf = SamplingDistributions::cosineHemispherePDF(wi);
//			wi = shadingData.frame.toWorld(wi);
//			Colour bsdf = shadingData.bsdf->evaluate(shadingData, wi);
//			pathThroughput *= bsdf * fabsf(Dot(wi, shadingData.sNormal)) / pdf;
//			r.init(shadingData.x + (wi * EPSILON), wi);
//			return direct + pathTrace(r, pathThroughput, depth + 1, sampler, shadingData.bsdf->isPureSpecular());
//		}
//		return scene->background->evaluate(shadingData, r.dir);
//	}
//
//	void render()
//	{
//		film->incrementSPP();
//		auto renderTile = [&](int startY, int endY, int threadId) {
//			for (int y = startY; y < endY; ++y) {
//				for (unsigned int x = 0; x < film->width; ++x) {
//					float px = x + 0.5f;
//					float py = y + 0.5f;
//					Ray ray = scene->camera.generateRay(px, py);
//					Colour throughput(1.0f);
//					Colour col = pathTrace(ray, throughput, 0, &samplers[threadId]);
//					film->splat(px, py, col);
//				}
//			}
//			};
//		int slice = film->height / numProcs;
//		for (int i = 0; i < numProcs; ++i) {
//			int startY = i * slice;
//			int endY = (i == numProcs - 1) ? film->height : startY + slice;
//			threads[i] = new std::thread(renderTile, startY, endY, i);
//		}
//		for (int i = 0; i < numProcs; ++i) {
//			threads[i]->join();
//			delete threads[i];
//		}
//	}
//
//	int getSPP() { return film->SPP; }
//	void saveHDR(std::string filename) { film->save(filename); }
//	void savePNG(std::string filename) {
//		stbi_write_png(filename.c_str(), canvas->getWidth(), canvas->getHeight(), 3, canvas->getBackBuffer(), canvas->getWidth() * 3);
//	}
//};
//
