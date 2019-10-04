/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core/core.h"

TR_NAMESPACE_BEGIN

/**
 * Modified Phong reflectance model
 */
struct PhongBSDF : BSDF {

    std::unique_ptr<Texture < v3f>> specularReflectance;
    std::unique_ptr<Texture < v3f>> diffuseReflectance;
    std::unique_ptr<Texture < float>> exponent;
    float specularSamplingWeight;
    float scale;

    PhongBSDF(const WorldData& scene, const Config& config, const size_t& matID) : BSDF(scene, config, matID) {
        const tinyobj::material_t& mat = scene.materials[matID];

        if (mat.specular_texname.empty())
            specularReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.specular)));
        else
            specularReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.specular_texname));

        if (mat.diffuse_texname.empty())
            diffuseReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.diffuse)));
        else
            diffuseReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.diffuse_texname));

        exponent = std::unique_ptr<Texture<float>>(new ConstantTexture1f(mat.shininess));

        //get scale value to ensure energy conservation
        v3f maxValue = specularReflectance->getMax() + diffuseReflectance->getMax();
        float actualMax = max(max(maxValue.x, maxValue.y), maxValue.z);
        scale = actualMax > 1.0f ? 0.99f * (1.0f / actualMax) : 1.0f;

        float dAvg = getLuminance(diffuseReflectance->getAverage() * scale);
        float sAvg = getLuminance(specularReflectance->getAverage() * scale);
        specularSamplingWeight = sAvg / (dAvg + sAvg);

        components.push_back(EGlossyReflection);
        components.push_back(EDiffuseReflection);

        combinedType = 0;
        for (unsigned int component : components)
            combinedType |= component;
    }

    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }

    v3f eval(const SurfaceInteraction& i) const override {
        v3f val(0.f);

        // TODO(A2): Implement this
		// 1. Test for front-facing incoming and outgoing rays;
		// 2. Evaluate the Phong BRDF in Equation (1);
		// 3. Return this value multiplied by the cosine foreshortening factor.

		if (Frame::cosTheta(i.wi) >= 0 && Frame::cosTheta(i.wo) >= 0) {
			// Perform step 2. 
			// Multiply p_d and p_s by PhongBSDF::scale to follow the energy law. 
			v3f p_d = diffuseReflectance->eval(worldData, i) * PhongBSDF::scale;
			v3f p_s = specularReflectance->eval(worldData, i) * PhongBSDF::scale;

			// Compute the Phong exponent. Do similar evaluations as diffuse cases. 
			float n = exponent->eval(worldData, i);

			// Compute the perfect specular direction w_r by using PhongBSDF::reflect();
			v3f w_r = PhongBSDF::reflect(i.wi);

			// Compute the cosine of alpha. 
			// alpha is the angle between w_r and lighting direction. Cosine becomes the dot product of the two vectors. 
			// a * b = |a||b|cos(theta). -> cos(theta) = a*b / (|a||b|).
			float cosine = glm::dot(i.wo, w_r) / (glm::length(i.wo) * glm::length(w_r));

			// TODO: cosine = 0 when cosine < 0?
			// Overwrite cosine by cosine with exponential. 
			cosine = pow(cosine, n);

			if (cosine >= 0) {
				val = (p_d / M_PI + (p_s * (n + 2) * cosine) / (2 * M_PI)) * Frame::cosTheta(glm::normalize(i.wi));
			}
			else {
				val = p_d / M_PI * Frame::cosTheta(glm::normalize(i.wi));
			}
		}
		else {
			// Return black. 
			val = v3f(0.f, 0.f, 0.f);
		}

        return val;
    }

    float pdf(const SurfaceInteraction& i) const override {
        float pdf = 0.f;

        // TODO(A3): Implement this

        return pdf;
    }

    v3f sample(SurfaceInteraction& i, Sampler& sampler, float* pdf) const override {
        v3f val(0.f);

        // TODO(A3): Implement this

        return val;
    }

    std::string toString() const override { return "Phong"; }
};

TR_NAMESPACE_END