/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <core/platform.h>
#include <core/integrator.h>

TR_NAMESPACE_BEGIN

/**
 * Surface normal integrator.
 */
struct NormalIntegrator : Integrator {
    explicit NormalIntegrator(const Scene& scene) : Integrator(scene) { }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        // HINT: Use the scene.bvh->intersect method. It's definition is in src/accel.h
        // TODO(A1): Implement this
		
		SurfaceInteraction interaction;
		bool hit = scene.bvh->intersect(ray, interaction);
		v3f color;

		if (!hit) {
			// No hit, color is black.
			color = v3f(0.f, 0.f, 0.f);
		}
		else {
			// Retrieve the normal at the hit point. 
			color = glm::abs(interaction.frameNs.n);
		}
		return color;
    }
};

TR_NAMESPACE_END