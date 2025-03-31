#pragma once

#include "Core.h"
#include "Imaging.h"
#include "Sampling.h"

#pragma warning( disable : 4244)

class BSDF;

class ShadingData
{
public:
	Vec3 x;
	Vec3 wo;
	Vec3 sNormal;
	Vec3 gNormal;
	float tu;
	float tv;
	Frame frame;
	BSDF* bsdf;
	float t;
	ShadingData() {}
	ShadingData(Vec3 _x, Vec3 n)
	{
		x = _x;
		gNormal = n;
		sNormal = n;
		bsdf = NULL;
	}
};

class ShadingHelper
{
public:


	template<typename T>
	static T clamp(T x, T a, T b)
	{
		return x < a ? a : (x > b ? b : x);
	}

	static float fresnelDielectric(float cosTheta, float iorInt, float iorExt)
	{
		// Calculate refracted direction
		// Calculate Fresnel (|| and T)
		// return average^2
		// Add code here
		//return 1.0f;



		cosTheta = ShadingHelper::clamp(cosTheta, -1.0f, 1.0f);
		bool entering = cosTheta > 0.0f;
		float etaI = entering ? iorExt : iorInt;
		float etaT = entering ? iorInt : iorExt;

		float sinThetaT2 = (etaI / etaT) * sqrtf(std::max(0.0f, 1 - cosTheta * cosTheta));
		if (sinThetaT2 >= 1.0f) return 1.0f; // Total internal reflection

		float cosThetaT = sqrtf(std::max(0.0f, 1 - sinThetaT2 * sinThetaT2));

		float rParallel = ((etaT * cosTheta) - (etaI * cosThetaT)) /
			((etaT * cosTheta) + (etaI * cosThetaT));
		float rPerpendicular = ((etaI * cosTheta) - (etaT * cosThetaT)) /
			((etaI * cosTheta) + (etaT * cosThetaT));
		return 0.5f * (rParallel * rParallel + rPerpendicular * rPerpendicular);
	}
	static Colour fresnelConductor(float cosTheta, Colour ior, Colour k)
	{
		// Add code here
		//return Colour(1.0f, 1.0f, 1.0f);

		cosTheta = ShadingHelper::clamp(cosTheta, 0.0f, 1.0f);
		Colour cos2 = Colour(cosTheta * cosTheta, cosTheta * cosTheta, cosTheta * cosTheta);
		Colour sin2 = Colour(1.0f, 1.0f, 1.0f) - cos2;

		Colour ior2 = ior * ior;
		Colour k2 = k * k;

		Colour t0 = ior2 - k2 - sin2;
		Colour a2plusb2 = sqrt(t0 * t0 + ior2 * k2 * 4.0f);
		Colour t1 = a2plusb2 + cos2;
		Colour t2 = Colour(2.0f * cosTheta, 2.0f * cosTheta, 2.0f * cosTheta) * sqrt(a2plusb2);

		Colour Rs = (t1 - t2) / (t1 + t2);
		Colour Rp = Rs;

		return (Rs + Rp) * 0.5f;

	}
	static float lambdaGGX(Vec3 wi, float alpha)
	{
		// Add code here
		//return 1.0f;
		float cosTheta = fabsf(wi.z);
		if (cosTheta < 1e-4f) return 0.0f;
		float tanTheta = sqrtf(1.0f - cosTheta * cosTheta) / cosTheta;
		float a = 1.0f / (alpha * tanTheta);
		if (a >= 1.6f) return 0.0f;
		return (1.0f - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181f * a * a);
	}
	static float Gggx(Vec3 wi, Vec3 wo, float alpha)
	{
		// Add code here
		//return 1.0f;
		return 1.0f / (1.0f + lambdaGGX(wi, alpha) + lambdaGGX(wo, alpha));
	}
	static float Dggx(Vec3 h, float alpha)
	{
		// Add code here
		//return 1.0f;
		float cosThetaH = h.z;
		float cosThetaH2 = cosThetaH * cosThetaH;
		float alpha2 = alpha * alpha;
		float denom = (cosThetaH2 * (alpha2 - 1.0f) + 1.0f);
		return alpha2 / (M_PI * denom * denom);
	}
};

class BSDF
{
public:
	Colour emission;
	virtual Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf) = 0;
	virtual Colour evaluate(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual float PDF(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual bool isPureSpecular() = 0;
	virtual bool isTwoSided() = 0;
	bool isLight()
	{
		return emission.Lum() > 0 ? true : false;
	}
	void addLight(Colour _emission)
	{
		emission = _emission;
	}
	Colour emit(const ShadingData& shadingData, const Vec3& wi)
	{
		return emission;
	}
	virtual float mask(const ShadingData& shadingData) = 0;
};


class DiffuseBSDF : public BSDF
{
public:
	Texture* albedo;
	DiffuseBSDF() = default;
	DiffuseBSDF(Texture* _albedo)
	{
		albedo = _albedo;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Add correct sampling code here
		//Vec3 wi = Vec3(0, 1, 0);
		//pdf = 1.0f;
		//reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		//wi = shadingData.frame.toWorld(wi);
		//return wi;

		Vec3 wiLocal = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		Vec3 wiWorld = shadingData.frame.toWorld(wiLocal);
		pdf = wiLocal.z / M_PI;
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;

		return wiWorld;

	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Add correct PDF code here
		//return 1.0f;
		Vec3 wiLocal = shadingData.frame.toLocal(wi);

		if (wiLocal.z <= 0.0f)
			return 0.0f;

		return wiLocal.z / M_PI;
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class MirrorBSDF : public BSDF
{
public:
	Texture* albedo;
	MirrorBSDF() = default;
	MirrorBSDF(Texture* _albedo)
	{
		albedo = _albedo;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Replace this with Mirror sampling code
		//Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		//pdf = wi.z / M_PI;
		//reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		//wi = shadingData.frame.toWorld(wi);
		//return wi;

		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo); // 出射方向
		Vec3 wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);   // 理想镜面反射（z 轴法线）

		// 如果反射方向在表面下方（理论上不会发生，但可以加保险）
		if (wiLocal.z <= 0.0f) {
			pdf = 0.0f;
			reflectedColour = Colour(0.0f, 0.0f, 0.0f);
			return Vec3(0.0f, 0.0f, 0.0f);
		}

		pdf = 1.0f; // delta 分布，不是连续分布
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv);
		return shadingData.frame.toWorld(wiLocal);

	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Mirror evaluation code
		//return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		return Colour(0.0f, 0.0f, 0.0f);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Mirror PDF
		//Vec3 wiLocal = shadingData.frame.toLocal(wi);
		//return SamplingDistributions::cosineHemispherePDF(wiLocal);
		return 0.0f;
	}
	bool isPureSpecular()
	{
		return true;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};


class ConductorBSDF : public BSDF
{
public:
	Texture* albedo;
	Colour eta;
	Colour k;
	float alpha;
	ConductorBSDF() = default;
	ConductorBSDF(Texture* _albedo, Colour _eta, Colour _k, float roughness)
	{
		albedo = _albedo;
		eta = _eta;
		k = _k;
		alpha = 1.62142f * sqrtf(roughness);
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Replace this with Conductor sampling code
		//Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		//pdf = wi.z / M_PI;
		//reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		//wi = shadingData.frame.toWorld(wi);
		//return wi;

		{
			// Convert wo to local space
			Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
			if (wo.z <= 0.0f)
			{
				pdf = 0.0f;
				reflectedColour = Colour(0.0f, 0.0f, 0.0f);
				return Vec3(0.0f, 0.0f, 0.0f);
			}

			// Sample microfacet normal h using GGX distribution
			Vec3 h = SamplingDistributions::GGXSample(alpha, sampler->next(), sampler->next());
			if (Dot(wo, h) <= 0.0f)
			{
				pdf = 0.0f;
				reflectedColour = Colour(0.0f, 0.0f, 0.0f);
				return Vec3(0.0f, 0.0f, 0.0f);
			}

			// Reflect wo about h to get wi
			Vec3 wi = Reflect(wo, h);
			if (wi.z <= 0.0f)
			{
				pdf = 0.0f;
				reflectedColour = Colour(0.0f, 0.0f, 0.0f);
				return Vec3(0.0f, 0.0f, 0.0f);
			}

			// Compute components
			Colour F = ShadingHelper::fresnelConductor(Dot(wi, h), eta, k);
			float D = ShadingHelper::Dggx(h, alpha);
			float G = ShadingHelper::Gggx(wi, wo, alpha);
			float denom = 4.0f * wo.z * wi.z;

			// Final reflected BSDF
			Colour f = albedo->sample(shadingData.tu, shadingData.tv) * F * D * G / denom;

			// PDF: GGX sampling PDF
			float Jh = 1.0f / (4.0f * Dot(wo, h)); // Jacobian
			float D_pdf = D * h.z / (4.0f * Dot(wo, h)); // half-vector sampling pdf
			pdf = D_pdf;

			reflectedColour = f;

			return shadingData.frame.toWorld(wi);
		}
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Conductor evaluation code
		//return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		if (wi.z <= 0.0f || woLocal.z <= 0.0f) return Colour(0.0f, 0.0f, 0.0f);

		Vec3 h = (wi + woLocal).normalize();
		if (Dot(wi, h) <= 0.0f) return Colour(0.0f, 0.0f, 0.0f);

		Colour F = ShadingHelper::fresnelConductor(Dot(wi, h), eta, k);
		float D = ShadingHelper::Dggx(h, alpha);
		float G = ShadingHelper::Gggx(wi, woLocal, alpha);
		float denom = 4.0f * wi.z * woLocal.z;

		Colour f = albedo->sample(shadingData.tu, shadingData.tv) * F * D * G / denom;
		return f;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Conductor PDF
		//Vec3 wiLocal = shadingData.frame.toLocal(wi);
		//return SamplingDistributions::cosineHemispherePDF(wiLocal);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		if (wi.z <= 0.0f || woLocal.z <= 0.0f) return 0.0f;

		Vec3 h = (wi + woLocal).normalize();
		float D = ShadingHelper::Dggx(h, alpha);
		float pdf = D * h.z / (4.0f * Dot(woLocal, h));
		return pdf;
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class GlassBSDF : public BSDF
{
public:
	Texture* albedo;
	float intIOR;
	float extIOR;
	GlassBSDF() = default;
	GlassBSDF(Texture* _albedo, float _intIOR, float _extIOR)
	{
		albedo = _albedo;
		intIOR = _intIOR;
		extIOR = _extIOR;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		if (woLocal.z == 0.0f) {
			pdf = 0.0f;
			reflectedColour = Colour(0.0f, 0.0f, 0.0f);
			return Vec3(0.0f,0.0f,0.0f);
		}

		float etaI = extIOR;
		float etaT = intIOR;
		bool entering = woLocal.z > 0.0f;
		if (!entering) std::swap(etaI, etaT);

		float eta = etaI / etaT;
		float cosThetaI = fabsf(woLocal.z);
		float sin2ThetaT = eta * eta * (1.0f - cosThetaI * cosThetaI);

		// Total internal reflection
		if (sin2ThetaT >= 1.0f) {
			Vec3 wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv);
			pdf = 1.0f;
			return shadingData.frame.toWorld(wiLocal);
		}

		float cosThetaT = sqrtf(std::max(0.0f, 1.0f - sin2ThetaT));
		float Fr = ShadingHelper::fresnelDielectric(cosThetaI, intIOR, extIOR);

		// Russian roulette: decide reflection or refraction
		if (sampler->next() < Fr) {
			// Reflect
			Vec3 wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * Fr;
			pdf = Fr;
			return shadingData.frame.toWorld(wiLocal);
		}
		else {
			// Refract
			Vec3 wiLocal;
			float sign = entering ? 1.0f : -1.0f;
			wiLocal.x = -woLocal.x * eta;
			wiLocal.y = -woLocal.y * eta;
			wiLocal.z = sign * cosThetaT;

			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * (1.0f - Fr);
			pdf = 1.0f - Fr;
			return shadingData.frame.toWorld(wiLocal);
		}
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Glass evaluation code
		//return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		return Colour(0.0f,0.0f,0.0f);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with GlassPDF
		//Vec3 wiLocal = shadingData.frame.toLocal(wi);
		//return SamplingDistributions::cosineHemispherePDF(wiLocal);
		return 0.0f;
	}
	bool isPureSpecular()
	{
		return true;
	}
	bool isTwoSided()
	{
		return false;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class DielectricBSDF : public BSDF
{
public:
	Texture* albedo;
	float intIOR;
	float extIOR;
	float alpha;
	DielectricBSDF() = default;
	DielectricBSDF(Texture* _albedo, float _intIOR, float _extIOR, float roughness)
	{
		albedo = _albedo;
		intIOR = _intIOR;
		extIOR = _extIOR;
		alpha = 1.62142f * sqrtf(roughness);
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		//// Replace this with Dielectric sampling code
		//Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		//pdf = wi.z / M_PI;
		//reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		//wi = shadingData.frame.toWorld(wi);
		//return wi;
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);  // 出射方向
		bool entering = wo.z > 0.0f;

		float etaI = entering ? extIOR : intIOR;
		float etaT = entering ? intIOR : extIOR;
		float eta = etaI / etaT;

		// 法线方向
		Vec3 n = Vec3(0.0f, 0.0f, 1.0f);
		if (!entering) n = -n;

		// 计算 cosθ
		float cosTheta = Dot(wo, n);
		cosTheta = ShadingHelper::clamp(cosTheta, -1.0f, 1.0f);

		// Fresnel 反射率
		float Fr = ShadingHelper::fresnelDielectric(cosTheta, intIOR, extIOR);

		if (sampler->next() < Fr) {
			// 反射分支
			Vec3 wi = Reflect(wo, n);
			pdf = Fr;
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * Fr / std::abs(wi.z);
			return shadingData.frame.toWorld(wi);
		}
		else {
			// 折射分支
			float sin2ThetaI = std::max(0.0f, 1.0f - cosTheta * cosTheta);
			float sin2ThetaT = eta * eta * sin2ThetaI;

			if (sin2ThetaT >= 1.0f) {
				// 全反射，作为保底
				Vec3 wi = Reflect(wo, n);
				pdf = 1.0f;
				reflectedColour = albedo->sample(shadingData.tu, shadingData.tv);
				return shadingData.frame.toWorld(wi);
			}

			float cosThetaT = sqrtf(std::max(0.0f, 1.0f - sin2ThetaT));
			Vec3 wi = (wo * -eta + n * (eta * cosTheta - cosThetaT));
			pdf = 1.0f - Fr;
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * (1.0f - Fr) / std::abs(wi.z);
			return shadingData.frame.toWorld(wi);
		}
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Dielectric evaluation code
		//return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		return Colour(0.0f, 0.0f, 0.0f);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Dielectric PDF
		//Vec3 wiLocal = shadingData.frame.toLocal(wi);
		//return SamplingDistributions::cosineHemispherePDF(wiLocal);
		return 0.0f;
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return false;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class OrenNayarBSDF : public BSDF
{
public:
	Texture* albedo;
	float sigma;
	OrenNayarBSDF() = default;
	OrenNayarBSDF(Texture* _albedo, float _sigma)
	{
		albedo = _albedo;
		sigma = _sigma;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Replace this with OrenNayar sampling code
		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		pdf = wi.z / M_PI;
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		wi = shadingData.frame.toWorld(wi);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with OrenNayar evaluation code
		//return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);

		if (wi.z <= 0.0f || woLocal.z <= 0.0f)
			return Colour(0.0f, 0.0f, 0.0f);

		// Convert sigma (in degrees) to radians and precompute sigma²
		float sigmaRad = sigma * (M_PI / 180.0f);
		float sigma2 = sigmaRad * sigmaRad;

		// Oren-Nayar coefficients
		float A = 1.0f - (sigma2 / (2.0f * (sigma2 + 0.33f)));
		float B = 0.45f * sigma2 / (sigma2 + 0.09f);

		float sinThetaI = sqrtf(std::max(0.0f, 1.0f - wi.z * wi.z));
		float sinThetaO = sqrtf(std::max(0.0f, 1.0f - woLocal.z * woLocal.z));

		float maxCos = 0.0f;
		if (sinThetaI > 1e-4f && sinThetaO > 1e-4f) {
			float sinPhiI = wi.y / sinThetaI;
			float cosPhiI = wi.x / sinThetaI;
			float sinPhiO = woLocal.y / sinThetaO;
			float cosPhiO = woLocal.x / sinThetaO;
			float dCos = cosPhiI * cosPhiO + sinPhiI * sinPhiO;
			maxCos = std::max(0.0f, dCos);
		}

		float sinAlpha = std::min(sinThetaI, sinThetaO);
		float tanBeta = std::max(sinThetaI, sinThetaO) / std::max(1e-4f, std::min(wi.z, woLocal.z));

		float orenNayar = A + B * maxCos * sinAlpha * tanBeta;

		Colour reflectance = albedo->sample(shadingData.tu, shadingData.tv);

		return reflectance * orenNayar / M_PI;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with OrenNayar PDF
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		return SamplingDistributions::cosineHemispherePDF(wiLocal);
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class PlasticBSDF : public BSDF
{
public:
	Texture* albedo;
	float intIOR;
	float extIOR;
	float alpha;
	PlasticBSDF() = default;
	PlasticBSDF(Texture* _albedo, float _intIOR, float _extIOR, float roughness)
	{
		albedo = _albedo;
		intIOR = _intIOR;
		extIOR = _extIOR;
		alpha = 1.62142f * sqrtf(roughness);
	}
	float alphaToPhongExponent()
	{
		return (2.0f / SQ(std::max(alpha, 0.001f))) - 2.0f;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Replace this with Plastic sampling code
		//Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		//pdf = wi.z / M_PI;
		//reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		//wi = shadingData.frame.toWorld(wi);
		//return wi;

		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		if (wo.z <= 0.0f) {
			pdf = 0.0f;
			reflectedColour = Colour(0.0f, 0.0f, 0.0f);
			return Vec3(0.0f, 0.0f, 0.0f);
		}

		// Fresnel 计算
		float F = ShadingHelper::fresnelDielectric(wo.z, intIOR, extIOR);

		// 随机采样决定是镜面还是漫反射
		if (sampler->next() < F) {
			// 镜面反射
			Vec3 wi = Vec3(-wo.x, -wo.y, wo.z);  // perfect mirror
			if (wi.z <= 0.0f) {
				pdf = 0.0f;
				reflectedColour = Colour(0.0f, 0.0f, 0.0f);
				return Vec3(0.0f, 0.0f, 0.0f);
			}
			pdf = F;
			//reflectedColour = Colour(F) * albedo->sample(shadingData.tu, shadingData.tv);
			reflectedColour = Colour(F, F, F) * albedo->sample(shadingData.tu, shadingData.tv);
			return shadingData.frame.toWorld(wi);
		}
		else {
			// 漫反射
			Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
			pdf = (1.0f - F) * wi.z / M_PI;
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * (1.0f - F) / M_PI;
			return shadingData.frame.toWorld(wi);
		}

	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Plastic evaluation code
		//return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);

		if (wi.z <= 0.0f || woLocal.z <= 0.0f)
			return Colour(0.0f, 0.0f, 0.0f);

		// 检查是否为镜面反射方向
		Vec3 mirror = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
		if (fabsf(wi.x - mirror.x) < 1e-4f &&
			fabsf(wi.y - mirror.y) < 1e-4f &&
			fabsf(wi.z - mirror.z) < 1e-4f)
		{
			float F = ShadingHelper::fresnelDielectric(woLocal.z, intIOR, extIOR);
			return albedo->sample(shadingData.tu, shadingData.tv) * F;
		}

		// Diffuse lobe
		float F = ShadingHelper::fresnelDielectric(woLocal.z, intIOR, extIOR);
		Colour diffuse = albedo->sample(shadingData.tu, shadingData.tv) * (1.0f - F) / M_PI;
		return diffuse;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Plastic PDF
		//Vec3 wiLocal = shadingData.frame.toLocal(wi);
		//return SamplingDistributions::cosineHemispherePDF(wiLocal);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);

		if (wi.z <= 0.0f || woLocal.z <= 0.0f)
			return 0.0f;

		float F = ShadingHelper::fresnelDielectric(woLocal.z, intIOR, extIOR);

		// 如果是镜面方向，返回 delta pdf
		Vec3 mirror = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
		if (fabsf(wi.x - mirror.x) < 1e-4f &&
			fabsf(wi.y - mirror.y) < 1e-4f &&
			fabsf(wi.z - mirror.z) < 1e-4f)
		{
			return F;
		}

		// 否则是 diffuse pdf
		return (1.0f - F) * wi.z / M_PI;
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class LayeredBSDF : public BSDF
{
public:
	BSDF* base;
	Colour sigmaa;
	float thickness;
	float intIOR;
	float extIOR;
	LayeredBSDF() = default;
	LayeredBSDF(BSDF* _base, Colour _sigmaa, float _thickness, float _intIOR, float _extIOR)
	{
		base = _base;
		sigmaa = _sigmaa;
		thickness = _thickness;
		intIOR = _intIOR;
		extIOR = _extIOR;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Add code to include layered sampling
		//return base->sample(shadingData, sampler, reflectedColour, pdf);
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		if (wo.z <= 0.0f)
		{
			pdf = 0.0f;
			reflectedColour = Colour(0.0f, 0.0f, 0.0f);
			return Vec3(0.0f, 0.0f, 0.0f);
		}

		// Fresnel to decide reflect or transmit
		float F = ShadingHelper::fresnelDielectric(wo.z, intIOR, extIOR);
		float r = sampler->next();

		if (r < F)
		{
			// Perfect specular reflection at air-coating interface
			Vec3 wi = Vec3(-wo.x, -wo.y, wo.z);
			reflectedColour = Colour(F, F, F);
			pdf = F;
			return shadingData.frame.toWorld(wi);
		}
		else
		{
			// Transmission through coating and interaction with base
			ShadingData baseData = shadingData;
			baseData.wo = shadingData.wo;

			Vec3 wi = base->sample(baseData, sampler, reflectedColour, pdf);
			Vec3 wiLocal = shadingData.frame.toLocal(wi);

			if (wiLocal.z <= 0.0f || pdf <= 0.0f)
			{
				pdf = 0.0f;
				reflectedColour = Colour(0.0f, 0.0f, 0.0f);
				return Vec3(0.0f, 0.0f, 0.0f);
			}

			// Compute attenuation from coating (2 passes through)
			float cosThetaI = std::abs(wo.z);
			float cosThetaO = std::abs(wiLocal.z);
			//Colour attenuation = exp(-sigmaa * thickness * (1.0f / cosThetaI + 1.0f / cosThetaO));
			Colour absorption = sigmaa * thickness * (1.0f / cosThetaI + 1.0f / cosThetaO);
			Colour attenuation = Colour(
				exp(-absorption.r),
				exp(-absorption.g),
				exp(-absorption.b)
			);


			//reflectedColour *= attenuation * (1.0f - F);
			reflectedColour = reflectedColour * attenuation * (1.0f - F);
			pdf *= (1.0f - F);

			return wi;
		}
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Add code for evaluation of layer
		//return base->evaluate(shadingData, wi);
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);

		if (wo.z <= 0.0f || wiLocal.z <= 0.0f)
			return Colour(0.0f, 0.0f, 0.0f);

		// Check for perfect mirror reflection
		Vec3 mirror = Vec3(-wo.x, -wo.y, wo.z);
		if ((wiLocal - mirror).length() < 1e-4f)
		{
			float F = ShadingHelper::fresnelDielectric(wo.z, intIOR, extIOR);
			return Colour(F, F, F);
		}

		// Otherwise, evaluate base BSDF and apply attenuation
		float F = ShadingHelper::fresnelDielectric(wo.z, intIOR, extIOR);
		Colour baseEval = base->evaluate(shadingData, wi);

		float cosThetaI = std::abs(wo.z);
		float cosThetaO = std::abs(wiLocal.z);
		//Colour attenuation = exp(-sigmaa * thickness * (1.0f / cosThetaI + 1.0f / cosThetaO));
		Colour absorption = sigmaa * thickness * (1.0f / cosThetaI + 1.0f / cosThetaO);
		Colour attenuation = Colour(
			exp(-absorption.r),
			exp(-absorption.g),
			exp(-absorption.b)
		);


		return baseEval * attenuation * (1.0f - F);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Add code to include PDF for sampling layered BSDF
		//return base->PDF(shadingData, wi);
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);

		if (wo.z <= 0.0f || wiLocal.z <= 0.0f)
			return 0.0f;

		// Check for specular path
		Vec3 mirror = Vec3(-wo.x, -wo.y, wo.z);
		if ((wiLocal - mirror).length() < 1e-4f)
		{
			float F = ShadingHelper::fresnelDielectric(wo.z, intIOR, extIOR);
			return F;
		}

		float F = ShadingHelper::fresnelDielectric(wo.z, intIOR, extIOR);
		return base->PDF(shadingData, wi) * (1.0f - F);
	}
	bool isPureSpecular()
	{
		return base->isPureSpecular();
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return base->mask(shadingData);
	}
};
