/*
	This file is part of TinyRender, an educative rendering system.

	Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
	Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Path tracer integrator
 */
	struct PathTracerIntegrator : Integrator {
	explicit PathTracerIntegrator(const Scene& scene) : Integrator(scene) {
		m_isExplicit = scene.config.integratorSettings.pt.isExplicit;
		m_maxDepth = scene.config.integratorSettings.pt.maxDepth;
		m_rrDepth = scene.config.integratorSettings.pt.rrDepth;
		m_rrProb = scene.config.integratorSettings.pt.rrProb;
	}

	v3f renderImplicit(const Ray& ray, Sampler& sampler, SurfaceInteraction& hit) const {
		v3f Li(0.f);

		// TODO(A5): Implement this
		std::vector<v3f> bsdfVal;	// Vector to hold the values of every point's bsdf. 

		for (int depth = 0; depth <= m_maxDepth; depth++) {
			// If the current interaction point is the light source
			// Verify that it is the front face of the light.
			// If yes, multiply with the contributions in the vector to compute traced contributions
			// Otherwise, bounce normally if max depth is not arrived. 
			bool frontFacing = hit.wo.z > 0;	// This assumes that the light is always pointing downwards. TODO: change this?
			if (getEmission(hit) != v3f(0.f) && frontFacing) {
				v3f val(1.f);

				for (int i = 0; i < bsdfVal.size(); i++) {
					val *= bsdfVal[i];
				}
				Li = val * getEmission(hit);

				break;
			}
			else {
				// The current interaction point is not a light source
				// Evaluate the BSDF at this point and push into the vector to store the bsdf value. 
				float pdf;
				v3f contribution = getBSDF(hit)->sample(hit, sampler, &pdf);
				v3f wiW = glm::normalize(hit.frameNs.toWorld(hit.wi));

				if (m_maxDepth > 0) {
					Ray shadowRay(hit.p, wiW);

					// Make sure that the shadow ray hits something and
					// Most importantly: update the interaction point
					if (scene.bvh->intersect(shadowRay, hit)) {
						bsdfVal.push_back(contribution);
					}
					else {
						// The reflected ray didn't hit anything in the scene.
						// Return black. 
						return v3f(0.f);
					}
				}
			}
		}

		return clampBelow(Li, 0.f);
	}

	/**
	 *	This function performs direct illumination shading for the hit point using emitter area sampling.
	 *	The implementation is very similar to the polygonal integrator.
	 *	High level: sample the emitter by sampleEmitterByPosition, and use (pe - hit.p) as the incident ray
	 *	Compute geometry jacobian terms, bsdf, and the emitter contribution at the shading point.
	 *
	 *	Arguments:
	 *	- Ray& ray: the current ray that's being traced;
	 *	- Sampler& sampler: rand number sampler.
	 *	- SurfaceInteraction& hit: the interaction point.
	 *	Returns:
	 *	- The direct illumination value at the shading point.
	 */
	v3f directShading(const Ray& ray, Sampler& sampler, SurfaceInteraction& hit) const {
		v3f Ldir(0.f);

		float emPdf, pdf_A;
		v3f pe, ne;
		size_t id = selectEmitter(sampler.next(), emPdf);
		const Emitter& em = getEmitterByID(id);

		sampleEmitterPosition(sampler, em, ne, pe, pdf_A);
		v3f wiW = pe - hit.p;
		hit.wi = hit.frameNs.toLocal(wiW);

		if (glm::dot(glm::normalize(ne), glm::normalize(-wiW)) >= 0) {
			// Compute Jacobian of the geometry term
			float cosTheta_o = glm::abs(glm::dot(ne, glm::normalize(-wiW)));
			float Jacobian = cosTheta_o / glm::distance2(hit.p, pe);

			SurfaceInteraction shadowInteraction;
			Ray shadowRay(hit.p, wiW);

			if (scene.bvh->intersect(shadowRay, shadowInteraction)) {
				float pdf = emPdf * pdf_A;
				if (pdf != 0) {
					Ldir = getBSDF(hit)->eval(hit) * Jacobian * getEmission(shadowInteraction) / pdf;
				}
			}
		}

		return Ldir;
	}

	v3f shade(const Ray& ray, Sampler& sampler, SurfaceInteraction& hit, int depth) const {
		v3f Li(0.f), Ldir(0.f), Lind(0.f);
		float rr_factor = 1;

		// Check the condition to stop the recursion
		if (m_maxDepth == -1) {
			// RR recursion termination. 
			if (depth >= m_rrDepth) {
				if (sampler.next() >= m_rrProb) return v3f(0.f);
				else rr_factor = m_rrProb;
			}
		}
		else {
			// Regular recursion termination, check if depth has reached the max depth
			if (depth >= m_maxDepth) return v3f(0.f);

			if (m_rrDepth > 0 && depth >= m_rrDepth) {
				if (sampler.next() >= m_rrProb) return v3f(0.f);
				else rr_factor = m_rrProb;
			}
		}

		if (scene.bvh->intersect(ray, hit)) {
			bool frontFacing = hit.wo.z > 0;	// This assumes that the light is always pointing downwards. TODO: change this?
			if (getEmission(hit) != v3f(0.f) && frontFacing)
				return getEmission(hit);
			else {
				// Compute direct illumination
				Ldir = directShading(ray, sampler, hit);

				// Compute indirect illumination
				// Generate the shadow ray by sampling the upper hemisphere, BSDF::sample. 
				// Avoid double counting by checking if the sampled ray hits an emitter. 
				while (1) {
					float pdf;
					v3f val = getBSDF(hit)->sample(hit, sampler, &pdf);

					v3f wiW = hit.frameNs.toWorld(hit.wi);
					Ray shadowRay(hit.p, glm::normalize(wiW));
					SurfaceInteraction shadowInteraction;

					if (scene.bvh->intersect(shadowRay, shadowInteraction)) {
						if (getEmission(shadowInteraction) != v3f(0.f)) {
							continue;
						}
						else {
							Lind += val * shade(shadowRay, sampler, shadowInteraction, depth + 1);
							break;
						}
					}
					else {
						break;
					}
				}
			}
		}

		Li = Ldir + Lind / rr_factor;

		return Li;
	}

	v3f renderExplicit(const Ray& ray, Sampler& sampler, SurfaceInteraction& hit) const {
		v3f Li(0.f);

		// TODO(A5): Implement this
		if (m_maxDepth == 0) {
			if (scene.bvh->intersect(ray, hit))
				Li = getEmission(hit);
		}
		else {
			int depth = 0;
			Li = clampBelow(shade(ray, sampler, hit, depth), 0.f);
		}

		return Li;
	}


	v3f render(const Ray& ray, Sampler& sampler) const override {
		Ray r = ray;
		SurfaceInteraction hit;

		if (scene.bvh->intersect(r, hit)) {
			if (m_isExplicit)
				return this->renderExplicit(ray, sampler, hit);
			else
				return this->renderImplicit(ray, sampler, hit);
		}
		return v3f(0.0);
	}

	int m_maxDepth;     // Maximum number of bounces
	int m_rrDepth;      // When to start Russian roulette
	float m_rrProb;     // Russian roulette probability
	bool m_isExplicit;  // Implicit or explicit
};

TR_NAMESPACE_END
