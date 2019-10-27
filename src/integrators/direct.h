/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Direct illumination integrator with MIS
 */
struct DirectIntegrator : Integrator {
    explicit DirectIntegrator(const Scene& scene) : Integrator(scene) {
        m_emitterSamples = scene.config.integratorSettings.di.emitterSamples;
        m_bsdfSamples = scene.config.integratorSettings.di.bsdfSamples;
        m_samplingStrategy = scene.config.integratorSettings.di.samplingStrategy;
    }

    static inline float balanceHeuristic(float nf, float fPdf, float ng, float gPdf) {
        float f = nf * fPdf, g = ng * gPdf;
        return f / (f + g);
    }

    void sampleSphereByCosineHemisphere(const p2f& sample,
                                        const v3f& n,
                                        const p3f& pShading,
                                        const v3f& emitterCenter,
                                        float emitterRadius,
                                        v3f& wiW,
                                        float& pdf) const {
        // TODO(A3): Implement this
		wiW = Warp::squareToCosineHemisphere(sample);
		pdf = Warp::squareToCosineHemispherePdf(wiW);
    }

    void sampleSphereByArea(const p2f& sample,
                            const p3f& pShading,
                            const v3f& emitterCenter,
                            float emitterRadius,
                            v3f& pos,
                            v3f& ne,
                            v3f& wiW,
                            float& pdf) const {
        // TODO(A3): Implement this
		// samples on the emitter is a vector with radius of the emitter radius.
		v3f lightSample = Warp::squareToUniformSphere(sample) * emitterRadius;
		// shift the samples according to the emitter center. 
		pos = lightSample + emitterCenter;
		// wiW = x -> y
		wiW = glm::normalize(pos - pShading);
		ne = glm::normalize(lightSample);

		pdf = INV_FOURPI / pow(emitterRadius, 2);
    }

    void sampleSphereBySolidAngle(const p2f& sample,
                                  const p3f& pShading,
                                  const v3f& emitterCenter,
                                  float emitterRadius,
                                  v3f& wiW,
                                  float& pdf) const {
        // TODO(A3): Implement this
    }

    v3f renderArea(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);
		// TODO(A3): Implement this

		SurfaceInteraction i;
		bool hit = scene.bvh->intersect(ray, i);

		if (hit) {
			// Render the light emitter source
			if (getEmission(i) != v3f(0.f)) {
				// intersection point is on emitter. 
				// render the emitter.
				size_t emId = getEmitterIDByShapeID(i.shapeID);
				Emitter em = getEmitterByID(emId);
				Lr = em.getRadiance();
			}
			else {
				for (int j = 0; j < m_bsdfSamples; j++) {
					// Sampling function requires the following variables:
					// sample, pShading, emitterCenter, emitterRadius, pos, ne, wiW, pdf 
					// Retrieve emitter parameters
					float emPdf;
					size_t id = selectEmitter(sampler.next(), emPdf);
					const Emitter& em = getEmitterByID(id);
					v3f emCenter = scene.getShapeCenter(em.shapeID);
					float emRadius = scene.getShapeRadius(em.shapeID);
					float pdf;

					v3f pos(0.f);
					v3f ne(0.f);

					sampleSphereByArea(sampler.next2D(), i.p, emCenter, emRadius, pos, ne, i.wi, pdf);

					// Check direct illumination visibility
					SurfaceInteraction shadowInteraction;

					//if(scene.bvh->intersect())
				}
			}
		}
		
		return Lr;
    }

    v3f renderCosineHemisphere(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);		
		
        // TODO(A3): Implement this
		SurfaceInteraction i;
		bool hit = scene.bvh->intersect(ray, i);
		if (hit) {
			// If there is an interaction
			// Check if the intersection point is an emitter point. 
			// getEmission(i) should be nonzero if it is on an emitter. 
			if (getEmission(i) != v3f(0.f)) {
				// intersection point is on emitter. 
				// render the emitter.
				
				size_t emId = getEmitterIDByShapeID(i.shapeID);
				Emitter em = getEmitterByID(emId);
				Lr = em.getRadiance();
			}
			else {
				// Perform MC Estimation
				for (int j = 0; j < m_bsdfSamples; j++) {
					// Sampling function requires the following variables:
					// sample, n, pShading, emitterCenter, emitterRadius, wiW, pdf&. 
					// Retrieve emitter parameters
					float emPdf;
					size_t id = selectEmitter(sampler.next(), emPdf);
					const Emitter& em = getEmitterByID(id);
					v3f emCenter = scene.getShapeCenter(em.shapeID);
					float emRadius = scene.getShapeRadius(em.shapeID);
					float pdf;

					// Call the sampling function to set wiW and pdf. 
					sampleSphereByCosineHemisphere(sampler.next2D(), i.frameNs.n,
						i.p, emCenter, emRadius, i.wi, pdf);

					Ray shadowRay(i.p, i.frameNs.toWorld(glm::normalize(i.wi)));
					SurfaceInteraction shadowInteraction;

					if (scene.bvh->intersect(shadowRay, shadowInteraction)) {
						Lr += getBSDF(i)->eval(i) * getEmission(shadowInteraction) / pdf;
					}
				}

				Lr = Lr / m_bsdfSamples;
			}
		}

        return Lr;
    }

    v3f renderBSDF(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

        // TODO(A3): Implement this
		SurfaceInteraction i;
		bool hit = scene.bvh->intersect(ray, i);
		if (hit) {
			// If there is an interaction
			// Check if the intersection point is an emitter point. 
			// getEmission(i) should be nonzero if it is on an emitter. 
			if (getEmission(i) != v3f(0.f)) {
				// intersection point is on emitter. 
				// render the emitter.

				size_t emId = getEmitterIDByShapeID(i.shapeID);
				Emitter em = getEmitterByID(emId);
				Lr = em.getRadiance();
			}
			else {
				for (int j = 0; j < m_bsdfSamples; j++) {
					// test the diffuse BSDF. 
					float pdf;

					v3f val = getBSDF(i)->sample(i, sampler, &pdf);

					Ray shadowRay(i.p, glm::normalize(i.frameNs.toWorld(i.wi)));
					SurfaceInteraction shadowInteraction;

					// If the shadow ray hits an emitter, integrate the BSDF. 
					if (scene.bvh->intersect(shadowRay, shadowInteraction)) {
						Lr += val * getEmission(shadowInteraction);
					}
				}

				Lr /= m_bsdfSamples;
			}
		}
        return Lr;
    }

    v3f renderSolidAngle(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

        // TODO(A3): Implement this

        return Lr;
    }

    v3f renderMIS(const Ray& ray, Sampler& sampler) const {

        v3f Lr(0.f);

        // TODO(A4): Implement this

        return Lr;
    }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        if (m_samplingStrategy == ESamplingStrategy::EMIS)
            return this->renderMIS(ray, sampler);
        else if (m_samplingStrategy == ESamplingStrategy::EArea)
            return this->renderArea(ray, sampler);
        else if (m_samplingStrategy == ESamplingStrategy::ESolidAngle)
            return this->renderSolidAngle(ray, sampler);
        else if (m_samplingStrategy == ESamplingStrategy::ECosineHemisphere)
            return this->renderCosineHemisphere(ray, sampler);
        else
            return this->renderBSDF(ray, sampler);
    }

    size_t m_emitterSamples;     // Number of emitter samples
    size_t m_bsdfSamples;        // Number of BSDF samples
    ESamplingStrategy m_samplingStrategy;   // Sampling strategy to use
};

TR_NAMESPACE_END