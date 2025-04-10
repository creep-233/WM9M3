#pragma once

#include "Core.h"
#include "Geometry.h"
#include "Materials.h"
#include "Sampling.h"

#pragma warning( disable : 4244)

class SceneBounds
{
public:
	Vec3 sceneCentre;
	float sceneRadius;
};

class Light
{
public:
	virtual Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& emittedColour, float& pdf) = 0;

	virtual Colour evaluate(const Vec3& wi) = 0;
	virtual float PDF(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual bool isArea() = 0;
	virtual Vec3 normal(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual float totalIntegratedPower() = 0;
	virtual Vec3 samplePositionFromLight(Sampler* sampler, float& pdf) = 0;
	virtual Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf) = 0;
};

class AreaLight : public Light
{
public:
	Triangle* triangle = NULL;
	Colour emission;
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& emittedColour, float& pdf)
	{
		emittedColour = emission;
		return triangle->sample(sampler, pdf);
	}
	Colour evaluate(const Vec3& wi)
	{
		if (Dot(wi, triangle->gNormal()) < 0)
		{
			return emission;
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		return 1.0f / triangle->area;
	}
	bool isArea()
	{
		return true;
	}
	Vec3 normal(const ShadingData& shadingData, const Vec3& wi)
	{
		return triangle->gNormal();
	}
	float totalIntegratedPower()
	{
		return (triangle->area * emission.Lum());
	}
	Vec3 samplePositionFromLight(Sampler* sampler, float& pdf)
	{
		return triangle->sample(sampler, pdf);
	}
	Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf)
	{
		// Add code to sample a direction from the light
		Vec3 wi = Vec3(0, 0, 1);
		pdf = 1.0f;
		Frame frame;
		frame.fromVector(triangle->gNormal());
		return frame.toWorld(wi);
	}

};

class BackgroundColour : public Light
{
public:
	Colour emission;
	BackgroundColour(Colour _emission)
	{
		emission = _emission;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 wi = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
		pdf = SamplingDistributions::uniformSpherePDF(wi);
		reflectedColour = emission;
		return wi;
	}
	Colour evaluate(const Vec3& wi)
	{
		return emission;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		return SamplingDistributions::uniformSpherePDF(wi);
	}
	bool isArea()
	{
		return false;
	}
	Vec3 normal(const ShadingData& shadingData, const Vec3& wi)
	{
		return -wi;
	}
	float totalIntegratedPower()
	{
		return emission.Lum() * 4.0f * M_PI;
	}
	Vec3 samplePositionFromLight(Sampler* sampler, float& pdf)
	{
		Vec3 p = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
		p = p * use<SceneBounds>().sceneRadius;
		p = p + use<SceneBounds>().sceneCentre;
		pdf = 4 * M_PI * use<SceneBounds>().sceneRadius * use<SceneBounds>().sceneRadius;
		return p;
	}
	Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf)
	{
		Vec3 wi = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
		pdf = SamplingDistributions::uniformSpherePDF(wi);
		return wi;
	}
};

class EnvironmentMap : public Light
{
public:
	Texture* env;
	std::vector<float> marginalDist; 
	std::vector<std::vector<float>> conditionalDists; 

	EnvironmentMap(Texture* _env)
	{
		env = _env;

		marginalDist.resize(env->height);
		conditionalDists.resize(env->height);

		for (int y = 0; y < env->height; ++y) {
			std::vector<float>& row = conditionalDists[y];
			row.resize(env->width);
			float rowSum = 0.0f;

			float theta = M_PI * (y + 0.5f) / env->height;
			float sinTheta = sinf(theta);

			for (int x = 0; x < env->width; ++x) {
				Colour c = env->texels[y * env->width + x];
				float lum = c.Lum();
				float weight = lum * sinTheta;
				row[x] = weight;
				rowSum += weight;
			}
			marginalDist[y] = rowSum;

			// Normalize conditional distribution
			if (rowSum > 0) {
				for (float& val : row) val /= rowSum;
			}
		}

		// Normalize marginal distribution
		float marginalSum = 0.0f;
		for (float val : marginalDist) marginalSum += val;
		if (marginalSum > 0) {
			for (float& val : marginalDist) val /= marginalSum;
		}
	}
	//Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	//{
	//	// Assignment: Update this code to importance sampling lighting based on luminance of each pixel
	//	//Vec3 wi = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
	//	//pdf = SamplingDistributions::uniformSpherePDF(wi);
	//	//reflectedColour = evaluate(wi);
	//	//return wi;

	//	float rv = sampler->next();
	//	float ru = sampler->next();

	//	// Sample row (v) from marginal distribution
	//	int y = 0;
	//	float acc = 0.0f;
	//	for (; y < env->height - 1; ++y) {
	//		acc += marginalDist[y];
	//		if (rv <= acc) break;
	//	}
	//	rv = (float(y) + sampler->next()) / env->height;

	//	// Sample column (u) from conditional distribution
	//	const std::vector<float>& row = conditionalDists[y];
	//	int x = 0;
	//	acc = 0.0f;
	//	for (; x < env->width - 1; ++x) {
	//		acc += row[x];
	//		if (ru <= acc) break;
	//	}
	//	ru = (float(x) + sampler->next()) / env->width;

	//	reflectedColour = env->sample(ru, rv);

	//	// Convert to world direction
	//	float theta = rv * M_PI;
	//	float phi = ru * 2.0f * M_PI;
	//	float sinTheta = sinf(theta);
	//	Vec3 wi = Vec3(sinTheta * cosf(phi), cosf(theta), sinTheta * sinf(phi));

	//	// PDF: p(u, v) / (sinTheta * 2π²)
	//	float marginalPDF = marginalDist[y];
	//	float conditionalPDF = conditionalDists[y][x];
	//	pdf = marginalPDF * conditionalPDF * env->width * env->height / (2.0f * M_PI * M_PI * sinTheta + EPSILON);

	//	return wi;

	//	std::cout << "Sampled direction: " << wi.x << " " << wi.y << " " << wi.z << ", PDF: " << pdf << "\n";


	//}

	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Importance sampling environment map using tabulated distribution
		float rv = sampler->next();
		float ru = sampler->next();

		//Sample marginal (v - theta direction)
		int y = 0;
		float acc = 0.0f;
		for (; y < env->height - 1; ++y) {
			acc += marginalDist[y];
			if (rv <= acc) break;
		}
		rv = (float(y) + sampler->next()) / env->height;


		const std::vector<float>& row = conditionalDists[y];
		int x = 0;
		acc = 0.0f;
		for (; x < env->width - 1; ++x) {
			acc += row[x];
			if (ru <= acc) break;
		}
		ru = (float(x) + sampler->next()) / env->width;


		reflectedColour = env->sample(ru, rv);

		// Convert (u,v) to direction vector 
		float theta = rv * M_PI;
		float phi = ru * 2.0f * M_PI;
		float sinTheta = sinf(theta);
		Vec3 wi = Vec3(
			sinTheta * cosf(phi),
			cosf(theta),
			sinTheta * sinf(phi)
		);

		// Compute PDF in solid angle
		float marginalPDF = marginalDist[y];
		float conditionalPDF = conditionalDists[y][x];
		pdf = marginalPDF * conditionalPDF * env->width * env->height / (2.0f * M_PI * M_PI * sinTheta + EPSILON);

		return wi;
	}

	Colour evaluate(const Vec3& wi)
	{
		float u = atan2f(wi.z, wi.x);
		u = (u < 0.0f) ? u + (2.0f * M_PI) : u;
		u = u / (2.0f * M_PI);
		float v = acosf(wi.y) / M_PI;
		return env->sample(u, v);
	}
	//float PDF(const ShadingData& shadingData, const Vec3& wi)
	//{
	//	// Assignment: Update this code to return the correct PDF of luminance weighted importance sampling
	//	//return SamplingDistributions::uniformSpherePDF(wi);

	//	float theta = acosf(wi.y);
	//	float phi = atan2f(wi.z, wi.x);
	//	if (phi < 0) phi += 2.0f * M_PI;
	//	float u = phi / (2.0f * M_PI);
	//	float v = theta / M_PI;

	//	int x = std::min(int(u * env->width), env->width - 1);
	//	int y = std::min(int(v * env->height), env->height - 1);

	//	float sinTheta = sinf(theta);
	//	float marginalPDF = marginalDist[y];
	//	float conditionalPDF = conditionalDists[y][x];
	//	return marginalPDF * conditionalPDF * env->width * env->height / (2.0f * M_PI * M_PI * sinTheta + EPSILON);

	//}

	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Convert direction to spherical coordinates (theta, phi)
		float theta = acosf(wi.y);
		float phi = atan2f(wi.z, wi.x);
		if (phi < 0.0f) phi += 2.0f * M_PI;

		// Convert to (u, v) in [0,1]
		float u = phi / (2.0f * M_PI);
		float v = theta / M_PI;

		// Discrete bin index
		int x = std::min(int(u * env->width), env->width - 1);
		int y = std::min(int(v * env->height), env->height - 1);

		float sinTheta = sinf(theta);
		if (sinTheta < 1e-6f) return 0.0f;  // Prevent division by zero

		// Lookup tabulated PDF
		float marginalPDF = marginalDist[y];
		float conditionalPDF = conditionalDists[y][x];

		// Convert to solid angle density
		return marginalPDF * conditionalPDF * env->width * env->height / (2.0f * M_PI * M_PI * sinTheta + EPSILON);
	}


	bool isArea()
	{
		return false;
	}
	Vec3 normal(const ShadingData& shadingData, const Vec3& wi)
	{
		return -wi;
	}
	float totalIntegratedPower()
	{
		float total = 0;
		for (int i = 0; i < env->height; i++)
		{
			float st = sinf(((float)i / (float)env->height) * M_PI);
			for (int n = 0; n < env->width; n++)
			{
				total += (env->texels[(i * env->width) + n].Lum() * st);
			}
		}
		total = total / (float)(env->width * env->height);
		return total * 4.0f * M_PI;
	}
	Vec3 samplePositionFromLight(Sampler* sampler, float& pdf)
	{
		// Samples a point on the bounding sphere of the scene. Feel free to improve this.
		Vec3 p = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
		p = p * use<SceneBounds>().sceneRadius;
		p = p + use<SceneBounds>().sceneCentre;
		pdf = 1.0f / (4 * M_PI * SQ(use<SceneBounds>().sceneRadius));
		return p;
	}
	//Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf)
	//{
	//	// Replace this tabulated sampling of environment maps
	//	Vec3 wi = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
	//	pdf = SamplingDistributions::uniformSpherePDF(wi);
	//	return wi;
	//}

	Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf)
	{
		float rv = sampler->next();
		float ru = sampler->next();

		//Sample marginal (v - theta direction)
		int y = 0;
		float acc = 0.0f;
		for (; y < env->height - 1; ++y) {
			acc += marginalDist[y];
			if (rv <= acc) break;
		}
		rv = (float(y) + sampler->next()) / env->height;

		//Sample conditional (u - phi direction)
		const std::vector<float>& row = conditionalDists[y];
		int x = 0;
		acc = 0.0f;
		for (; x < env->width - 1; ++x) {
			acc += row[x];
			if (ru <= acc) break;
		}
		ru = (float(x) + sampler->next()) / env->width;

		// Convert to spherical direction
		float theta = rv * M_PI;
		float phi = ru * 2.0f * M_PI;
		float sinTheta = sinf(theta);
		Vec3 wi = Vec3(
			sinTheta * cosf(phi),
			cosf(theta),
			sinTheta * sinf(phi)
		);

		// Compute PDF in solid angle
		float marginalPDF = marginalDist[y];
		float conditionalPDF = conditionalDists[y][x];
		pdf = marginalPDF * conditionalPDF * env->width * env->height / (2.0f * M_PI * M_PI * sinTheta + EPSILON);

		return wi;
	}


};