/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once
#include <random>
#include "bsdfs/phong.h"

TR_NAMESPACE_BEGIN

/**
 * Reflection occlusion integrator
 */
struct ROIntegrator : Integrator {

    float m_exponent;

    explicit ROIntegrator(const Scene& scene) : Integrator(scene) {
        m_exponent = scene.config.integratorSettings.ro.exponent;
    }

    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }


    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f Li(0.f);

        // TODO(A3): Implement this
		// High level similar to AO. 
		// No need to limit the max shadow ray length. 
		float val = 0.f;

		SurfaceInteraction i;
		bool hit = scene.bvh->intersect(ray, i);
		if (hit) {
			v2f sample = sampler.next2D();

			// Create the view vector mirror reflection. 
			// this wo_r should be the center of the cosine sampling lobe. 
			v3f wo_r = glm::normalize(i.frameNs.toWorld(reflect(i.wo)));
			// Rotate the sampled incident ray to match the correct lobe center. 
			// Use Frame structure to perform the rotation. 
			Frame lobeCenterFrame = Frame(wo_r);

			// Sample the incident rays around the center of the lobe. 
			// Warp the incident ray
			// Calculate the pdf here as this ray will be modified later. 
			i.wi = Warp::squareToPhongLobe(sample, m_exponent);
			float pdf = Warp::squareToPhongLobePdf(i.wi, m_exponent);

			// Calculate cos(alpha). 
			// alpha is the angle b/w mirror reflection of the *view vector about the normal* with the incident direction. 
			// cosine alpha is always greater than 0. 
			float cosineAlpha = lobeCenterFrame.cosTheta(i.wi);
			if (cosineAlpha < 0) {
				cosineAlpha = 0;
			}

			// Transform the sampled incident ray to the world coordinate. 
			i.wi = glm::normalize(lobeCenterFrame.toWorld(i.wi));

			// Construct shadow ray
			Ray shadowRay(i.p, i.wi);
			bool visible = !scene.bvh->intersect(shadowRay, i);

			if (visible) {
				val = pow(cosineAlpha, m_exponent) * (m_exponent + 2) / (2 * M_PI * pdf);
			}
		}

		Li = v3f(val);

        return Li;
    }
};

TR_NAMESPACE_END