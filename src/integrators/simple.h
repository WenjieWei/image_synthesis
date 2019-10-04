/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Simple direct illumination integrator.
 */
struct SimpleIntegrator : Integrator {
    explicit SimpleIntegrator(const Scene& scene) : Integrator(scene) { }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f Li(0.f);

        // TODO(A2): Implement this
		// Perform ray tracing and check if an intersection occured. 
		SurfaceInteraction act;
		bool hit = scene.bvh->intersect(ray, act);

		if (hit) {
			// Retrieve the light position and intensity. 
			
			v3f pointLight_position = scene.getFirstLightPosition();
			v3f pointLight_intensity = scene.getFirstLightIntensity();

			// Compute the distance between the intersection point and the light source. 
			// World space to local coordinates. *wi should be a direction.
			// Retrieve intersected surface material. 
			float r = glm::length(pointLight_position - act.p);
			act.wi = glm::normalize(act.frameNs.toLocal(pointLight_position - act.p));
			v3f material = getBSDF(act)->eval(act);

			// 1.4: Construct shadow ray. 
			Ray shadowRay(act.p, glm::normalize(pointLight_position - act.p));
			shadowRay.max_t = glm::length(pointLight_position - act.p);
			bool visible = !scene.bvh->intersect(shadowRay, act);

			// Calculate BRDF.
			if (visible) {
				Li = (pointLight_intensity / (pow(r, 2))) * material;
			}
		} 
		else {
			Li = v3f(0.f, 0.f, 0.f);
		}

        return glm::abs(Li);
    }
};

TR_NAMESPACE_END