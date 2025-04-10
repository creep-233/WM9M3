#pragma once

#include "Core.h"
#include <random>
#include <algorithm>

class Sampler
{
public:
	virtual float next() = 0;
};

class MTRandom : public Sampler
{
public:
	std::mt19937 generator;
	std::uniform_real_distribution<float> dist;
	MTRandom(unsigned int seed = 1) : dist(0.0f, 1.0f)
	{
		generator.seed(seed);
	}
	float next()
	{
		return dist(generator);
	}
};

// Note all of these distributions assume z-up coordinate system
class SamplingDistributions
{
public:
	static Vec3 uniformSampleHemisphere(float r1, float r2) {

		float theta = acos(r1);          
		float phi = 2.0f * M_PI * r2;     

		float x = sin(theta) * cos(phi);  
		float y = sin(theta) * sin(phi);
		float z = cos(theta);            

		return Vec3(x, y, z);

	}
	
	static float uniformHemispherePDF(const Vec3 wi) {
		return 1.0f / (2.0f * M_PI);
	}

	static Vec3 cosineSampleHemisphere(float r1, float r2)
	{
		// Add code here
		//return Vec3(0, 0, 1);
		float r = sqrtf(r1);
		float theta = 2.0f * M_PI * r2;

		float x = r * cosf(theta);
		float y = r * sinf(theta);
		float z = sqrtf(std::max(0.0f, 1.0f - x * x - y * y));  

		return Vec3(x, y, z);
	}
	static float cosineHemispherePDF(const Vec3 wi)
	{
		// Add code here
		//return 1.0f;
		if (wi.z <= 0.0f) return 0.0f;  
		return wi.z / M_PI;
	}
	static Vec3 uniformSampleSphere(float r1, float r2)
	{
		float theta = acos(1 - 2 * r1); 
		float phi = 2.0f * M_PI * r2;   

		float x = sin(theta) * cos(phi);  
		float y = sin(theta) * sin(phi);
		float z = cos(theta);

		return Vec3(x, y, z);
	}
	static float uniformSpherePDF(const Vec3& wi)
	{
		return 1.0f / (4.0f * M_PI);
	}

	static Vec3 GGXSample(float alpha, float u1, float u2)
	{
		// Trowbridge-Reitz GGX importance sampling of the half-vector h
		float phi = 2.0f * M_PI * u1;

		// Transform uniform random variable u2 into the distribution of h.z
		float cosTheta = sqrtf((1.0f - u2) / (1.0f + (alpha * alpha - 1.0f) * u2));
		float sinTheta = sqrtf(1.0f - cosTheta * cosTheta);

		return Vec3(cosf(phi) * sinTheta, sinf(phi) * sinTheta, cosTheta);
	}



};

