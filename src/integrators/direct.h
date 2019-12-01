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
		// arguments:
		//	- sample: 2D canonical sample point
		//	- pShading: i.p
		//	- pos: position of the sample in world coordinates
		//	- ne: normal of emitter point
		ne = Warp::squareToUniformSphere(sample);
		pos = emitterCenter + emitterRadius * ne;
		wiW = glm::normalize(pos - pShading);
		pdf = INV_FOURPI / (pow(emitterRadius, 2));
    }

    void sampleSphereBySolidAngle(const p2f& sample,
                                  const p3f& pShading,
                                  const v3f& emitterCenter,
                                  float emitterRadius,
                                  v3f& wiW,
                                  float& pdf) const {
        // TODO(A3): Implement this
		float dist = glm::distance(pShading, emitterCenter);
		float theta = acos(dist / sqrt(pow(dist, 2) + pow(emitterRadius, 2)));

		// Construct a frame about x->y
		Frame coneFrame = Frame(glm::normalize(emitterCenter - pShading));

		wiW = coneFrame.toWorld(Warp::squareToUniformCone(sample, cos(theta)));
		pdf = Warp::squareToUniformConePdf(cos(theta));
	}

    v3f renderArea(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

        // TODO(A3): Implement this
		SurfaceInteraction i;
		bool hit = scene.bvh->intersect(ray, i);
		if (hit) {
			// Check if intersection is an emitter.
			if (getEmission(i) != v3f(0.f)) {
				// Render the emitter if nonzero
				size_t emId = getEmitterIDByShapeID(i.shapeID);
				Emitter em = getEmitterByID(emId);
				Lr = em.getRadiance();
			}
			else {
				// The interaction point is not on the emitter. 
				// Do BSDF shading according to the by area method. 
				for (int j = 0; j < m_emitterSamples; j++) {
					// sampleArea function signature:
					// sample, pShading, emitterCenter, emitterRadius, pos, ne, wiW, pdf
					float emPdf, pdf_A;
					v3f pos, ne;

					size_t id = selectEmitter(sampler.next(), emPdf);
					const Emitter& em = getEmitterByID(id);
					v3f emCenter = scene.getShapeCenter(em.shapeID);
					float emRadius = scene.getShapeRadius(em.shapeID);

					sampleSphereByArea(sampler.next2D(), i.p, emCenter, emRadius, pos, ne, i.wi, pdf_A);

					// After the sampling, all vectors should be in world coordinates. 
					// wiW is already normalized. 
					// Trace the shadow ray to perform direct illumination MC. 
					SurfaceInteraction shadowInteraction;
					Ray shadowRay(i.p, i.wi);
					if (scene.bvh->intersect(shadowRay, shadowInteraction)) {
						// If (x->y) dot ne is negative
						// That means the sampled point is behind the sphere, should be discarded. 
						if (glm::dot((emCenter - i.p), ne) >= 0) {
							// do the integration
							float pdf_Omega = pdf_A * glm::distance2(pos, i.p) / glm::dot(i.wi, ne);

							// Convert i.wi to local coord to do the bsdf evaluation. 
							i.wi = i.frameNs.toLocal(i.wi);
							Lr += getBSDF(i)->eval(i) * getEmission(shadowInteraction) / (pdf_Omega * emPdf);
						}
					}
				}

				Lr /= m_emitterSamples;
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
		SurfaceInteraction i;
		bool hit = scene.bvh->intersect(ray, i);
		if (hit) {
			if (getEmission(i) != v3f(0.f)) {
				// Interaction point is an emitter. 
				size_t emId = getEmitterIDByShapeID(i.shapeID);
				Emitter em = getEmitterByID(emId);
				Lr = em.getRadiance();
			}
			else {
				for (int j = 0; j < m_emitterSamples; j++) {
					// Sampling function requires the following variables:
					// sample, pShading, emitterCenter, emitterRadius, wiW, pdf
					float emPdf;
					size_t id = selectEmitter(sampler.next(), emPdf);
					const Emitter& em = getEmitterByID(id);
					v3f emCenter = scene.getShapeCenter(em.shapeID);
					float emRadius = scene.getShapeRadius(em.shapeID);
					float pdf;

					sampleSphereBySolidAngle(sampler.next2D(), i.p, emCenter, emRadius, i.wi, pdf);

					Ray shadowRay(i.p, i.wi);
					SurfaceInteraction shadowInteraction;
					if (scene.bvh->intersect(shadowRay, shadowInteraction)) {
						i.wi = i.frameNs.toLocal(i.wi);
						Lr += getBSDF(i)->eval(i) * getEmission(shadowInteraction) / (pdf * emPdf);
					}
				}

				Lr /= m_emitterSamples;
			}
		}

        return Lr;
    }

    v3f renderMIS(const Ray& ray, Sampler& sampler) const {

        v3f Lr(0.f);

        // TODO(A4): Implement this
		SurfaceInteraction i;
		bool hit = scene.bvh->intersect(ray, i);

		if (hit) {
			if (getEmission(i) != v3f(0.f)) {
				// Render the emitter if the emission value is not zero. 
				size_t emId = getEmitterIDByShapeID(i.shapeID);
				Emitter em = getEmitterByID(emId);
				Lr = em.getRadiance();
			}
			else {
				// Intersection point is not on the emitter.
				// Perform the MC Integration for solid angle first. 

				v3f val_emitter(0.f), val_bsdf(0.f);
				float pdf_bsdf = 0.f, pdf_emitter = 0.f;

				for (int j = 0; j < m_emitterSamples; j++) {
					// Sampling function requires the following variables:
					// sample, pShading, emitterCenter, emitterRadius, wiW, pdf
					float emPdf;
					size_t id = selectEmitter(sampler.next(), emPdf);
					const Emitter& em = getEmitterByID(id);
					v3f emCenter = scene.getShapeCenter(em.shapeID);
					float emRadius = scene.getShapeRadius(em.shapeID);

					sampleSphereBySolidAngle(sampler.next2D(), i.p, emCenter, emRadius, i.wi, pdf_emitter);

					Ray shadowRay(i.p, i.wi);
					SurfaceInteraction shadowInteraction;
					if (scene.bvh->intersect(shadowRay, shadowInteraction) && getEmission(shadowInteraction) != v3f(0.f)) {
						i.wi = i.frameNs.toLocal(i.wi);
						pdf_bsdf = getBSDF(i)->pdf(i);
						v3f val = getBSDF(i)->eval(i);

						float ws = balanceHeuristic(m_emitterSamples, (pdf_emitter * emPdf), m_bsdfSamples, pdf_bsdf);
						if ((pdf_emitter * emPdf) != 0) {
							val_emitter += getEmission(shadowInteraction) * val * ws / (pdf_emitter * emPdf);
						}
					}
				}

				// Perform MC Integration for BSDF. 
				for (int k = 0; k < m_bsdfSamples; k++) {
					v3f val = getBSDF(i)->sample(i, sampler, &pdf_bsdf);
					i.wi = i.frameNs.toWorld(i.wi);

					Ray shadowRay(i.p, i.wi);
					SurfaceInteraction shadowInteraction;

					if (scene.bvh->intersect(shadowRay, shadowInteraction) && getEmission(shadowInteraction) != v3f(0.f)) {
						// Retrieve the parameters of the emitter that the shadow ray intersects. 
						float emPdf = 1.0f / scene.emitters.size();
						Emitter em = getEmitterByID(getEmitterIDByShapeID(shadowInteraction.shapeID));
						v3f emCenter = scene.getShapeCenter(em.shapeID);
						float emRadius = scene.getShapeRadius(em.shapeID);

						float dist = glm::distance(i.p, emCenter);
						float thetaMax = acos(dist / sqrt(pow(dist, 2) + pow(emRadius, 2)));
						pdf_emitter = Warp::squareToUniformConePdf(cos(thetaMax));

						// Call the balance heuristic function to determine the weight of both samplings. 
						float ws = balanceHeuristic(m_bsdfSamples, pdf_bsdf, m_emitterSamples, (pdf_emitter * emPdf));
						val_bsdf += val * getEmission(shadowInteraction) * ws;
					}
				}

				// Double check to avoid the divide by zero exception. 
				if (m_emitterSamples == 0) {
					Lr = val_bsdf / m_bsdfSamples;
				}
				else if (m_bsdfSamples == 0) {
					Lr = val_emitter / m_emitterSamples;
				}
				else {
					Lr = val_emitter / m_emitterSamples + val_bsdf / m_bsdfSamples;
				}
			}
		}

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