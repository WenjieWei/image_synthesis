/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <tiny_obj_loader.h>
#define RAY_EPS_CV 1e-5 // Use when setting min and max dist for ray in control variates code
TR_NAMESPACE_BEGIN

/**
 * Direct illumination integrator for polygonal light sources
 * Follows Arvo '94.
 */
struct PolygonalIntegrator : Integrator {

	float m_alpha;             // Control variates "strength"
	size_t m_visSamples;       // # of samples to estimate h - alpha*g
	bool m_traceShadows;       // Trace shadows or not
	EPolygonalMethod m_method; // Method to use (Arvo, or control variates)

	std::vector<std::vector<v3f>> m_triangles; // Data structure to store triangles

    explicit PolygonalIntegrator(const Scene& scene) : Integrator(scene) {
        m_alpha = scene.config.integratorSettings.poly.alpha;
        m_visSamples = scene.config.integratorSettings.poly.visSamples;
        m_traceShadows = scene.config.integratorSettings.poly.traceShadows;
        m_method = scene.config.integratorSettings.poly.method;

		/**
		 * 1) Get # of triangles on emitter
		 * 2) Store vertices in m_triangles
		 */
		// TODO(A4): Implement this
		size_t shapeID = scene.getFirstLight();
		auto shape = scene.worldData.shapes[shapeID];
		int numTriangles = shape.mesh.indices.size() / 3;

		for (int i = 0; i < numTriangles; i++) {
			std::vector<v3f> tri(3);	// Create empty triangle

			v3f vertexPosition;
			for (int j = 0; j < 3; j++) {
				vertexPosition = scene.getObjectVertexPosition(shapeID, j + i * 3);
				tri[j] = vertexPosition;
			}
			m_triangles.push_back(tri);
		}
    }

    /// Reflect
    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }

    /**
     * === PHONG BONUS ONLY ===
     * Compute the following integral:
     *    T(a, b, n, x) = \int_0^x [a \cos(\theta) + b \sin(\theta)]ˆn d\theta
     * Uses a recurrent relation (see Snyder's note, 1996)
     *
     * Series function:
     *    T_sum(a, b, n, x) = \sum_{i=0}ˆ{(n-1)/2} T(a, b, 2i+1, x)
     * assuming n is _odd_
     */
    float cosineSinePowerIntegralSum(float a, float b, int exp, float theta) const {
        if (exp % 2 == 0) exp += 1; // Make exponent odd
        float Tsum = 0.f;

		// Implementing this function may be useful if you attempt the bonus

        // TODO(A4): Implement this

        return Tsum;
    }

    /**
     * Compute edge (v1--v2) contribution
	 * The exp term is only needed if you attempt the bonus, otherwise, you can ignore it
     */
    float getEdgeContrib(const v3f& v1, const v3f& v2, const SurfaceInteraction& i, int exp = 0) const {
        float contrib = 0.f;

        // TODO(A4): Implement this
		float Theta;
		v3f Gamma(0.f);

		v3f r1 = v1 - i.p;
		v3f r2 = v2 - i.p;

		Theta = acos(glm::dot(glm::normalize(r1), glm::normalize(r2)));
		Gamma = glm::normalize(glm::cross(r2, r1));
		contrib = Theta * (glm::dot(Gamma, i.frameNs.n));

        return contrib;
    }
	   

    /// Direct illumination using Arvo '94 analytic solution for polygonal lights
    v3f renderAnalytic(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

        // TODO(A4): Implement this
		Emitter em = scene.emitters[0];

		SurfaceInteraction i;
		bool hit = scene.bvh->intersect(ray, i);
		if (hit) {
			if (getEmission(i) != v3f(0.f)) {
				size_t emId = getEmitterIDByShapeID(i.shapeID);
				Emitter em = getEmitterByID(emId);
				Lr = em.getRadiance();
			}
			else {
				v3f Ei(0.f);
				Emitter em = scene.emitters[0];
				float contribution = 0.f;	// Sum of each edge's contribution. 

				for (int j = 0; j < m_triangles.size(); j++) {
					// loop through triangle vector and compute Equation 3.
					std::vector<v3f> tri = m_triangles[j];

					// Loop through the triangle vertices and compute the contribution of each edge. 
					float Theta = 0.f;
					v3f Gamma(0.f);
					for (int k = 0; k < tri.size(); k++) {
						v3f v1(0.f), v2(0.f);
						if (k != tri.size() - 1) {
							v1 = tri[k];
							v2 = tri[k + 1];
						}
						else {
							v1 = tri[k];
							v2 = tri[0];
						}

						contribution += getEdgeContrib(v1, v2, i);
					}
				}

				Ei = INV_TWOPI * (em.getPower() / em.area) * contribution;
				v3f wiW = i.wi;         // A backup of i.wi
				i.wi = v3f(0.f, 0.f, 1.f);

				v3f rho_d = getBSDF(i)->eval(i);
				Lr = clampBelow(rho_d * Ei, 0);
			}
		}

        return Lr;
    }

    /**
     * Stand-alone estimator for h - alpha*g (with primary ray)
     * Trace a primary ray, check for emitter hit, and then call `estimateVisDiff()`
     * Used by polygonal render pass
     */
    v3f estimateVisDiffRealTime(const Ray& ray, Sampler& sampler, const Emitter& em) {
        v3f D(0.f);

        SurfaceInteraction hit;
        if (!scene.bvh->intersect(ray, hit)) return D;

        const BSDF* bsdf = getBSDF(hit);
        if (bsdf->isEmissive()) return D;

        hit.wi = v3f(0, 0, 1); // Trick to get 1/pi * albedo without cosine term
        D = estimateVisDiff(sampler, hit, em);

        return D;
    }

    /// Stand-alone estimator for h - alpha*g (without primary ray)
	/// Use RAY_EPS_CV when setting min and max dist for shadow ray
    v3f estimateVisDiff(Sampler& sampler, SurfaceInteraction& i, const Emitter& em) const {
        v3f sum(0.f);

        // TODO(A4): Implement this
		v3f pe, ne;
		float pdf;

		sampleEmitterPosition(sampler, em, ne, pe, pdf);

		i.wi = pe - i.p;
		if (glm::dot(glm::normalize(ne), glm::normalize(-i.wi)) >= 0) {
			v3f emission(0.f), g(0.f), h(0.f);

			float cosTheta_o = glm::abs(glm::dot(ne, glm::normalize(-i.wi)));
			float Jacobian = cosTheta_o / glm::distance2(i.p, pe);

			v3f wiW = i.wi;
			SurfaceInteraction shadowInteraction;
			Ray shadowRay(i.p, wiW);
			shadowRay.min_t = RAY_EPS_CV;

			i.wi = i.frameNs.toLocal(i.wi);
			// evaluate h(x)
			if (scene.bvh->intersect(shadowRay, shadowInteraction)) {
				h = getBSDF(i)->eval(i) * Jacobian * getEmission(shadowInteraction);
			}

			// evaluate g(x)
			g = getBSDF(i)->eval(i) * Jacobian * scene.getFirstLightIntensity();

			if (pdf != 0) {
				sum = (h - m_alpha * g) / pdf;
			}
		}

        return sum;
    }

    /// Control variates using Arvo '94 for direct illumination; ray trace shadows
	
    v3f renderControlVariates(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

        // TODO(A4): Implement this
		SurfaceInteraction i;
		bool hit = scene.bvh->intersect(ray, i);
		if (hit) {
			if (getEmission(i) != v3f(0.f)) {
				size_t emId = getEmitterIDByShapeID(i.shapeID);
				Emitter em = getEmitterByID(emId);
				Lr = em.getRadiance();
			}
			else {
				// Analytic part
				v3f Ei(0.f), G(0.f), H(0.f);
				Emitter em = scene.emitters[0];

				// Perform Monte Carlo for h-ag. 
				for (int j = 0; j < m_visSamples; j++) {
					// Do the control variance. 
					H += estimateVisDiff(sampler, i, em);
				}
				H = H / m_visSamples;
				G = renderAnalytic(ray, sampler);

				Lr = m_alpha * G + H;
			}
		}

		return clampBelow(Lr, 0);
    }

    /// Direct illumination using surface area sampling
    v3f renderArea(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

        // TODO(A4): Implement this
		SurfaceInteraction i;
		bool hit = scene.bvh->intersect(ray, i);
		if (hit) {
			if (getEmission(i) != v3f(0.f)) {
				size_t emId = getEmitterIDByShapeID(i.shapeID);
				Emitter em = getEmitterByID(emId);
				Lr = em.getRadiance();
			}
			else {
				Emitter em = scene.emitters[0];
				for (int j = 0; j < m_visSamples; j++) {
					v3f pe, ne; // Point on emitter and Surface normal at emitter point. 
					float pdf;

					sampleEmitterPosition(sampler, em, ne, pe, pdf);

					i.wi = pe - i.p;
					if (glm::dot(glm::normalize(ne), glm::normalize(-i.wi)) >= 0) {
						// Use illumination equation (2) in a3
						v3f emission(0.f);

						if (m_traceShadows) {
							SurfaceInteraction shadowInteraction;
							Ray shadowRay(i.p, i.wi);
							if (scene.bvh->intersect(shadowRay, shadowInteraction)) {
								emission = getEmission(shadowInteraction);
							}
						}
						else {
							emission = v3f(1.f);
						}

						// Front facing. 
						// Calculatte the geometry terms: 04-DirectIllumination-I&II, page 39.
						float cosTheta_o = glm::abs(glm::dot(ne, glm::normalize(-i.wi)));
						float Jacobian = cosTheta_o / glm::distance2(i.p, pe);

						i.wi = i.frameNs.toLocal(i.wi);
						if (pdf != 0) {
							if (m_traceShadows) {
								Lr = getBSDF(i)->eval(i) * Jacobian * emission / pdf;
							}
							else {
								Lr = Lr = getBSDF(i)->eval(i) * Jacobian * scene.getFirstLightIntensity() / pdf;
							}
						}
					}
				}
			}
		}

        return Lr;
    }

    /// Branch to corresponding method
    v3f render(const Ray& ray, Sampler& sampler) const override {
        switch (m_method) {
            case EPolygonalMethod::ESurfaceArea:
                return PolygonalIntegrator::renderArea(ray, sampler);
                break;
            case EPolygonalMethod::EControlVariates:
                return PolygonalIntegrator::renderControlVariates(ray, sampler);
                break;
            default:
                return PolygonalIntegrator::renderAnalytic(ray, sampler);
                break;
        }
    }

};

TR_NAMESPACE_END