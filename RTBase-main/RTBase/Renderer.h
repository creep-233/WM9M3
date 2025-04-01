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
#include <OpenImageDenoise/oidn.hpp>


class RayTracer
{
public:
	Scene* scene;
	GamesEngineeringBase::Window* canvas;
	Film* film;
	MTRandom *samplers;
	std::thread **threads;
	int numProcs;

	// OIDN
	float* colorBuffer;
	float* normalBuffer;
	float* albedoBuffer;
	float* outputBuffer;



	void init(Scene* _scene, GamesEngineeringBase::Window* _canvas)
	{

		scene = _scene;
		canvas = _canvas;

		int width = scene->camera.width;
		int height = scene->camera.height;
		colorBuffer = new float[width * height * 3];
		normalBuffer = new float[width * height * 3];
		albedoBuffer = new float[width * height * 3];
		outputBuffer = new float[width * height * 3];


		film = new Film();
		film->init((unsigned int)scene->camera.width, (unsigned int)scene->camera.height, /*new BoxFilter()*/ new GaussianFilter());
		SYSTEM_INFO sysInfo;
		GetSystemInfo(&sysInfo);
		numProcs = sysInfo.dwNumberOfProcessors;
		threads = new std::thread*[numProcs];
		samplers = new MTRandom[numProcs];
		clear();
		scene->traceVPLs(&samplers[0], 1000);
	}
	void clear()
	{
		film->clear();
	}


	void denoiseAndSave(std::string filename)
	{
		int width = scene->camera.width;
		int height = scene->camera.height;

	
		oidn::DeviceRef device = oidn::newDevice();
		device.setErrorFunction([](void*, oidn::Error code, const char* message) {
			std::cerr << "OIDN error: " << message << std::endl;
			});
		device.commit();


		oidn::BufferRef colorBuf = device.newBuffer(width * height * 3 * sizeof(float));
		oidn::BufferRef albedoBuf = device.newBuffer(width * height * 3 * sizeof(float));
		oidn::BufferRef normalBuf = device.newBuffer(width * height * 3 * sizeof(float));
		oidn::BufferRef outputBuf = device.newBuffer(width * height * 3 * sizeof(float));


		std::memcpy(colorBuf.getData(), colorBuffer, width * height * 3 * sizeof(float));
		std::memcpy(albedoBuf.getData(), albedoBuffer, width * height * 3 * sizeof(float));
		std::memcpy(normalBuf.getData(), normalBuffer, width * height * 3 * sizeof(float));

	
		oidn::FilterRef filter = device.newFilter("RT");
		filter.setImage("color", colorBuf, oidn::Format::Float3, width, height);
		filter.setImage("albedo", albedoBuf, oidn::Format::Float3, width, height);
		filter.setImage("normal", normalBuf, oidn::Format::Float3, width, height);
		filter.setImage("output", outputBuf, oidn::Format::Float3, width, height);
		filter.set("hdr", true);
		filter.commit();


		filter.execute();


		std::memcpy(outputBuffer, outputBuf.getData(), width * height * 3 * sizeof(float));

	
		for (int y = 0; y < height; ++y)
		{
			for (int x = 0; x < width; ++x)
			{
				int idx = (y * width + x) * 3;
				float r = outputBuffer[idx + 0];
				float g = outputBuffer[idx + 1];
				float b = outputBuffer[idx + 2];


				unsigned char cr = (unsigned char)(min(r, 1.0f) * 255.0f);
				unsigned char cg = (unsigned char)(min(g, 1.0f) * 255.0f);
				unsigned char cb = (unsigned char)(min(b, 1.0f) * 255.0f);
				canvas->draw(x, y, cr, cg, cb);
			}
		}

		savePNG(filename);
	}

	void connectToCamera(Vec3 p, Vec3 n, Colour col)
	{
		float px, py;
		if (!scene->camera.projectOntoCamera(p, px, py))
			return;

		Vec3 camDir = scene->camera.viewDirection.normalize();
		Vec3 dirToCam = (scene->camera.origin - p).normalize();

		float cosTheta = Dot(camDir, dirToCam);
		if (cosTheta <= 0.0f) return;

		float geometry = 1.0f / (scene->camera.Afilm * powf(cosTheta, 4.0f));
		Colour importance = col * geometry;

		// splat onto film
		film->splat(px, py, importance);
	}

	void lightTrace(Sampler* sampler)
	{
		float pmf;
		Light* light = scene->sampleLight(sampler, pmf);
		if (!light) return;

		float pdfPos;
		Vec3 p = light->samplePositionFromLight(sampler, pdfPos);
		float pdfDir;
		Vec3 wi = light->sampleDirectionFromLight(sampler, pdfDir);
		if (pdfPos <= 0.0f || pdfDir <= 0.0f) return;

		Colour Le = light->evaluate(-wi) / pdfPos;
		Colour throughput = Le;


		//connectToCamera(p, light->normal(p, wi), throughput);
		ShadingData shadingData;
		shadingData.x = p; 
		connectToCamera(p, light->normal(shadingData, wi), throughput);

	
		Ray r(p + wi * EPSILON, wi);
		lightTracePath(r, throughput, Le, sampler);
	}

	void lightTracePath(Ray& r, Colour pathThroughput, Colour Le, Sampler* sampler)
	{
		for (int depth = 0; depth < MAX_DEPTH; ++depth)
		{
			IntersectionData isect = scene->traverse(r);
			if (isect.t >= FLT_MAX)
				break;

			ShadingData shadingData = scene->calculateShadingData(isect, r);
			Vec3 nextWi;
			Colour f;
			float pdf;

			if (!shadingData.bsdf)
				break;

	
			connectToCamera(shadingData.x, shadingData.sNormal, pathThroughput * shadingData.bsdf->emit(shadingData, shadingData.wo));

			if (shadingData.bsdf->isPureSpecular())
			{
				nextWi = shadingData.bsdf->sample(shadingData, sampler, f, pdf);
			}
			else
			{
				nextWi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
				pdf = SamplingDistributions::cosineHemispherePDF(nextWi);
				nextWi = shadingData.frame.toWorld(nextWi);
				f = shadingData.bsdf->evaluate(shadingData, nextWi);
			}

			if (pdf <= 0.0f || f.Lum() <= 0.0f)
				break;

			float cosTheta = max(0.0f, Dot(nextWi, shadingData.sNormal));
			pathThroughput = pathThroughput * f * cosTheta / pdf;

			// Russian Roulette
			float rrProb = min(0.9f, pathThroughput.Lum());
			if (sampler->next() > rrProb)
				break;
			pathThroughput = pathThroughput / rrProb;

			r = Ray(shadingData.x + nextWi * EPSILON, nextWi);
		}
	}




	//Colour computeDirect(ShadingData shadingData, Sampler* sampler)
	//{
	//	if (shadingData.bsdf->isPureSpecular())
	//		return Colour(0.0f,0.0f,0.0f);

	//	Colour result = Colour(0.0f, 0.0f, 0.0f);


	//	// Light sampling path

	//	float pmf;
	//	Light* light = scene->sampleLight(sampler, pmf);
	//	if (light && pmf > 0.0f)
	//	{
	//		float lightPdf;
	//		Colour Le;

	//		Vec3 wi;
	//		float G = 1.0f;

	//		if (light->isArea())
	//		{
	//			Vec3 p = light->sample(shadingData, sampler, Le, lightPdf);
	//			wi = (p - shadingData.x);
	//			float l2 = wi.lengthSq();
	//			wi = wi.normalize();
	//			G = (max(Dot(wi, shadingData.sNormal), 0.0f) * max(-Dot(wi, light->normal(shadingData, wi)), 0.0f)) / (l2 + EPSILON);

	//			if (G > 0.0f && scene->visible(shadingData.x, p))
	//			{
	//				float bsdfPdf = shadingData.bsdf->PDF(shadingData, wi);
	//				float weight = (lightPdf * pmf) / (lightPdf * pmf + bsdfPdf + EPSILON);
	//				Colour f = shadingData.bsdf->evaluate(shadingData, wi);
	//				result = result + f * Le * G * weight / (lightPdf * pmf);
	//			}
	//		}
	//		else
	//		{
	//			wi = light->sample(shadingData, sampler, Le, lightPdf); // environment
	//			G = max(Dot(wi, shadingData.sNormal), 0.0f);
	//			if (G > 0.0f && scene->visible(shadingData.x, shadingData.x + wi * 10000.0f))
	//			{
	//				float bsdfPdf = shadingData.bsdf->PDF(shadingData, wi);
	//				float weight = (lightPdf * pmf) / (lightPdf * pmf + bsdfPdf + EPSILON);
	//				Colour f = shadingData.bsdf->evaluate(shadingData, wi);
	//				result = result + f * Le * G * weight / (lightPdf * pmf);
	//			}
	//		}
	//	}


	//	// BSDF sampling path

	//	float bsdfPdf;
	//	Colour f;
	//	Vec3 wi_bsdf = shadingData.bsdf->sample(shadingData, sampler, f, bsdfPdf);

	//	if (bsdfPdf > 0.0f)
	//	{
	//		Ray shadowRay(shadingData.x + shadingData.sNormal * EPSILON, wi_bsdf);
	//		IntersectionData isect = scene->traverse(shadowRay);
	//		if (isect.t < FLT_MAX)
	//		{
	//			ShadingData lightHit = scene->calculateShadingData(isect, shadowRay);
	//			if (lightHit.bsdf && lightHit.bsdf->isLight())
	//			{
	//				Colour Le = lightHit.bsdf->emit(lightHit, -wi_bsdf);
	//				float lightPdf = 0.0f;
	//				Light* hitLight = scene->getLightFromHit(isect);
	//				if (hitLight)
	//					lightPdf = hitLight->PDF(shadingData, wi_bsdf);

	//				float weight = bsdfPdf / (bsdfPdf + lightPdf * pmf + EPSILON);
	//				float cosTheta = max(Dot(wi_bsdf, shadingData.sNormal), 0.0f);
	//				Colour contrib = f * Le * cosTheta * weight / bsdfPdf;
	//				result = result + contrib;
	//			}
	//		}
	//	}

	//	return result;
	//}

	Colour computeDirect(ShadingData shadingData, Sampler* sampler)
	{
		if (shadingData.bsdf->isPureSpecular())
			return Colour(0.0f, 0.0f, 0.0f);

		Colour result = Colour(0.0f, 0.0f, 0.0f);

		// === Light sampling path ===
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
					Colour f = shadingData.bsdf->evaluate(shadingData, wi);

		
					if (lightPdf * pmf > 1e-4f && f.Lum() > 1e-6f)
					{
						float weight = (lightPdf * pmf) / (lightPdf * pmf + bsdfPdf + EPSILON);
						Colour contrib = f * Le * G * weight / (lightPdf * pmf);

			
						if (contrib.Lum() < 20.0f) 
							result = result + contrib;
					}
				}
			}
			else
			{
				wi = light->sample(shadingData, sampler, Le, lightPdf); // env map
				G = max(Dot(wi, shadingData.sNormal), 0.0f);

				if (G > 0.0f && scene->visible(shadingData.x, shadingData.x + wi * 10000.0f))
				{
					float bsdfPdf = shadingData.bsdf->PDF(shadingData, wi);
					Colour f = shadingData.bsdf->evaluate(shadingData, wi);

				
					if (lightPdf * pmf > 1e-4f && f.Lum() > 1e-6f)
					{
						float weight = (lightPdf * pmf) / (lightPdf * pmf + bsdfPdf + EPSILON);
						Colour contrib = f * Le * G * weight / (lightPdf * pmf);

				
						if (contrib.Lum() < 20.0f)
							result = result + contrib;
					}
				}
			}
		}

		// === BSDF sampling path ===
		float bsdfPdf;
		Colour f;
		Vec3 wi_bsdf = shadingData.bsdf->sample(shadingData, sampler, f, bsdfPdf);

		if (bsdfPdf > 1e-4f && f.Lum() > 1e-6f)
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

					if (Le.Lum() > 1e-6f && lightPdf > 1e-4f)
					{
						float weight = bsdfPdf / (bsdfPdf + lightPdf * pmf + EPSILON);
						float cosTheta = max(Dot(wi_bsdf, shadingData.sNormal), 0.0f);
						Colour contrib = f * Le * cosTheta * weight / bsdfPdf;

				
						if (contrib.Lum() < 20.0f)
							result = result + contrib;
					}
				}
			}
		}

		return result;
	}



	Colour computeDirectFromVPLs(ShadingData shadingData)
	{
		if (shadingData.bsdf->isPureSpecular())
			return Colour(0.0f, 0.0f, 0.0f);

		Colour result = Colour(0.0f, 0.0f, 0.0f);

		for (const Scene::VPL& vpl : scene->vpls)
		{
			Vec3 wi = vpl.shadingData.x - shadingData.x;
			float d2 = wi.lengthSq();
			wi = wi.normalize();

			float cosCamera = max(Dot(wi, shadingData.sNormal), 0.0f);
			float cosVPL = max(Dot(-wi, vpl.shadingData.sNormal), 0.0f);
			float G = cosCamera * cosVPL / (d2 + EPSILON);
			if (G <= 0.0f) continue;

			if (!scene->visible(shadingData.x, vpl.shadingData.x))
				continue;

			Colour f1 = shadingData.bsdf->evaluate(shadingData, wi);
			Colour f2 = vpl.shadingData.bsdf->evaluate(vpl.shadingData, -wi);

			result = result + vpl.Le * f1 * f2 * G;
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
			//Colour direct = pathThroughput * computeDirectFromVPLs(shadingData);

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

				wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
				pdf = SamplingDistributions::cosineHemispherePDF(wi);
				wi = shadingData.frame.toWorld(wi);
				bsdf = shadingData.bsdf->evaluate(shadingData, wi);
				if (pdf <= 0.0f || bsdf.Lum() <= 0.0f)
					return direct;

				//pathThroughput = pathThroughput * bsdf * fabsf(Dot(wi, shadingData.sNormal)) / pdf;
				float cosTheta = fabsf(Dot(wi, shadingData.sNormal));

				// 护栏：下限
				if (pdf > 1e-4f && bsdf.Lum() > 1e-6f && cosTheta > 1e-4f)
				{
					pathThroughput = pathThroughput * bsdf * cosTheta / pdf;

					// 护栏：上限 clamp
					if (pathThroughput.Lum() < 100.0f)
					{
						r.init(shadingData.x + wi * EPSILON, wi);
						return direct + pathTrace(r, pathThroughput, depth + 1, sampler, false);
					}
					else
					{
						// discard explosive path
						return direct;
					}
				}
				else
				{
					return direct;
				}

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
			//return computeDirectFromVPLs(shadingData);

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

		//scene->traceVPLs(&samplers[0], 1000);

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

					int idx = (y * film->width + x) * 3;


					colorBuffer[idx + 0] = final.r;
					colorBuffer[idx + 1] = final.g;
					colorBuffer[idx + 2] = final.b;


					Ray ray = scene->camera.generateRay(px, py);
					Colour normal = viewNormals(ray);
					normalBuffer[idx + 0] = normal.r;
					normalBuffer[idx + 1] = normal.g;
					normalBuffer[idx + 2] = normal.b;


					Colour albedoVal = albedo(ray);
					albedoBuffer[idx + 0] = albedoVal.r;
					albedoBuffer[idx + 1] = albedoVal.g;
					albedoBuffer[idx + 2] = albedoVal.b;


					film->splat(px, py, final);
					if (final.Lum() < 1e-4f) {

						final = Colour(1.0f, 0.0f, 1.0f); 
					}


					unsigned char r = (unsigned char)(clamp(final.r, 0.0f, 1.0f) * 255);
					unsigned char g = (unsigned char)(clamp(final.g, 0.0f, 1.0f) * 255);
					unsigned char b = (unsigned char)(clamp(final.b, 0.0f, 1.0f) * 255);

					final.r = powf(clamp(final.r, 0.0f, 1.0f), 1.0f / 2.2f);
					final.g = powf(clamp(final.g, 0.0f, 1.0f), 1.0f / 2.2f);
					final.b = powf(clamp(final.b, 0.0f, 1.0f), 1.0f / 2.2f);


					film->tonemap(x, y, r, g, b);
					canvas->draw(x, y, r, g, b);
				}
			}
			};

		//LightTracing
		//auto renderBlock = [&](int threadId, int yStart, int yEnd) {
		//	const int samplesPerThread = 100000; 
		//	for (int i = 0; i < samplesPerThread; ++i)
		//	{
		//		lightTrace(&samplers[threadId]);
		//	}
		//	};



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


