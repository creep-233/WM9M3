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
		if (shadingData.bsdf->isPureSpecular())
			return Colour(0.0f,0.0f,0.0f);

		Colour result = Colour(0.0f, 0.0f, 0.0f);


		// Light sampling path

		float pmf;
		Light* light = scene->sampleLight(sampler, pmf);
		if (light && pmf > 0.0f)
		{
			float lightPdf;
			Colour Le;

			Vec3 wi;
			float G = 1.0f;

			if (light->isArea())
			{
				Vec3 p = light->sample(shadingData, sampler, Le, lightPdf);
				wi = (p - shadingData.x);
				float l2 = wi.lengthSq();
				wi = wi.normalize();
				G = (max(Dot(wi, shadingData.sNormal), 0.0f) * max(-Dot(wi, light->normal(shadingData, wi)), 0.0f)) / (l2 + EPSILON);

				if (G > 0.0f && scene->visible(shadingData.x, p))
				{
					float bsdfPdf = shadingData.bsdf->PDF(shadingData, wi);
					float weight = (lightPdf * pmf) / (lightPdf * pmf + bsdfPdf + EPSILON);
					Colour f = shadingData.bsdf->evaluate(shadingData, wi);
					result = result + f * Le * G * weight / (lightPdf * pmf);
				}
			}
			else
			{
				wi = light->sample(shadingData, sampler, Le, lightPdf); // environment
				G = max(Dot(wi, shadingData.sNormal), 0.0f);
				if (G > 0.0f && scene->visible(shadingData.x, shadingData.x + wi * 10000.0f))
				{
					float bsdfPdf = shadingData.bsdf->PDF(shadingData, wi);
					float weight = (lightPdf * pmf) / (lightPdf * pmf + bsdfPdf + EPSILON);
					Colour f = shadingData.bsdf->evaluate(shadingData, wi);
					result = result + f * Le * G * weight / (lightPdf * pmf);
				}
			}
		}


		// BSDF sampling path

		float bsdfPdf;
		Colour f;
		Vec3 wi_bsdf = shadingData.bsdf->sample(shadingData, sampler, f, bsdfPdf);

		if (bsdfPdf > 0.0f)
		{
			Ray shadowRay(shadingData.x + shadingData.sNormal * EPSILON, wi_bsdf);
			IntersectionData isect = scene->traverse(shadowRay);
			if (isect.t < FLT_MAX)
			{
				ShadingData lightHit = scene->calculateShadingData(isect, shadowRay);
				if (lightHit.bsdf && lightHit.bsdf->isLight())
				{
					Colour Le = lightHit.bsdf->emit(lightHit, -wi_bsdf);
					float lightPdf = 0.0f;
					Light* hitLight = scene->getLightFromHit(isect);
					if (hitLight)
						lightPdf = hitLight->PDF(shadingData, wi_bsdf);

					float weight = bsdfPdf / (bsdfPdf + lightPdf * pmf + EPSILON);
					float cosTheta = max(Dot(wi_bsdf, shadingData.sNormal), 0.0f);
					Colour contrib = f * Le * cosTheta * weight / bsdfPdf;
					result = result + contrib;
				}
			}
		}

		return result;
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
			//Colour bsdf;
			//float pdf;
			//Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
			//pdf = SamplingDistributions::cosineHemispherePDF(wi);
			//wi = shadingData.frame.toWorld(wi);
			//bsdf = shadingData.bsdf->evaluate(shadingData, wi);
			//pathThroughput = pathThroughput * bsdf * fabsf(Dot(wi, shadingData.sNormal)) / pdf;
			//r.init(shadingData.x + (wi * EPSILON), wi);
			//return (direct + pathTrace(r, pathThroughput, depth + 1, sampler, shadingData.bsdf->isPureSpecular()));

			float pdf;
			Colour bsdf;
			Vec3 wi;

			// 对于 specular 材质，必须调用 sample
			if (shadingData.bsdf->isPureSpecular())
			{
				wi = shadingData.bsdf->sample(shadingData, sampler, bsdf, pdf);
				if (pdf <= 0.0f || bsdf.Lum() <= 0.0f)
					return direct;

				pathThroughput = pathThroughput * bsdf;
				r.init(shadingData.x + wi * EPSILON, wi);
				return pathTrace(r, pathThroughput, depth + 1, sampler, true);
			}
			else
			{
				// 继续原有的 cosineSampleHemisphere 路径
				wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
				pdf = SamplingDistributions::cosineHemispherePDF(wi);
				wi = shadingData.frame.toWorld(wi);
				bsdf = shadingData.bsdf->evaluate(shadingData, wi);
				if (pdf <= 0.0f || bsdf.Lum() <= 0.0f)
					return direct;

				pathThroughput = pathThroughput * bsdf * fabsf(Dot(wi, shadingData.sNormal)) / pdf;
				r.init(shadingData.x + (wi * EPSILON), wi);
				return direct + pathTrace(r, pathThroughput, depth + 1, sampler, false);
			}

		}
		return scene->background->evaluate(r.dir);
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
		return scene->background->evaluate(r.dir);
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


	template<typename T>
	inline T clamp(T x, T a, T b)
	{
		return x < a ? a : (x > b ? b : x);
	}

	void render()
	{
		film->incrementSPP();

	//	auto renderBlock = [&](int threadId, int yStart, int yEnd) {
	//		for (int y = yStart; y < yEnd; ++y)
	//		{
	//			for (unsigned int x = 0; x < film->width; ++x)
	//			{
	//				float px = x + 0.5f;
	//				float py = y + 0.5f;
	//				Ray ray = scene->camera.generateRay(px, py);

	//				//Colour col = direct(ray, &samplers[threadId]);
	//				Colour pathThroughput(1.0f, 1.0f, 1.0f);
	//				Colour col = pathTrace(ray, pathThroughput, 0, &samplers[threadId]);
	//				film->splat(px, py, col);

	//				unsigned char r = (unsigned char)(col.r * 255);
	//				unsigned char g = (unsigned char)(col.g * 255);
	//				unsigned char b = (unsigned char)(col.b * 255);

	//				film->tonemap(x, y, r, g, b);
	//				canvas->draw(x, y, r, g, b);
	//			}
	//		}
	//		};

		auto renderBlock = [&](int threadId, int yStart, int yEnd) {
			for (int y = yStart; y < yEnd; ++y)
			{
				for (unsigned int x = 0; x < film->width; ++x)
				{
					float px = x + 0.5f;
					float py = y + 0.5f;

					Colour accum = Colour(0.0f, 0.0f, 0.0f);
					float sumLum = 0.0f;
					float sumLum2 = 0.0f;

					//Adaptive Sampling
					//This will affect the running speed
					const int maxSamples = 8;
					const int minSamples = 2;
					const float varianceThreshold = 0.001f;

					int s = 0;
					for (; s < maxSamples; ++s)
					{
						Ray ray = scene->camera.generateRay(px, py);

						//RayTracing
						//Colour sample = direct(ray, &samplers[threadId]);
						//PathTracing
						Colour pathThroughput(1.0f, 1.0f, 1.0f);
						Colour sample = pathTrace(ray, pathThroughput, 0, &samplers[threadId]);

						accum = accum + sample;
						float lum = sample.Lum();
						sumLum += lum;
						sumLum2 += lum * lum;

						if (s >= minSamples)
						{
							float mean = sumLum / (s + 1);
							float var = (sumLum2 / (s + 1)) - (mean * mean);
							if (var < varianceThreshold)
								break;  // stop early
						}
					}

					Colour final = accum / float(s + 1);
					film->splat(px, py, final);
					if (final.Lum() < 1e-4f) {
						// 可以标记黑色点
						final = Colour(1.0f, 0.0f, 1.0f); // 粉色警告色
					}


					unsigned char r = (unsigned char)(clamp(final.r, 0.0f, 1.0f) * 255);
					unsigned char g = (unsigned char)(clamp(final.g, 0.0f, 1.0f) * 255);
					unsigned char b = (unsigned char)(clamp(final.b, 0.0f, 1.0f) * 255);

					film->tonemap(x, y, r, g, b);
					canvas->draw(x, y, r, g, b);
				}
			}
			};


		int blockSize = film->height / numProcs;
		std::vector<std::thread> threadPool;

		for (int i = 0; i < numProcs; ++i)
		{
			int yStart = i * blockSize;
			int yEnd = (i == numProcs - 1) ? film->height : yStart + blockSize;
			threadPool.emplace_back(renderBlock, i, yStart, yEnd);
		}

		for (auto& t : threadPool)
		{
			t.join();
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


