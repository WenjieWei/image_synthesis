/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Ambient occlusion integrator
 */
struct AOIntegrator : Integrator {

	// Use this in your switch statement to select the sampling type 
	ESamplingType m_samplingStrategy;

    explicit AOIntegrator(const Scene& scene) : Integrator(scene) { 
		m_samplingStrategy = scene.config.integratorSettings.ao.sampling_type;
	}

    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f Li(0.f);
		
		/*
		Use the m_sampling_type variable to set wi and the corresponding pdf 
		appropriately for sphere, hemisphere, or cosine sampling.

		You can use a switch statement or an if/else block.

		The m_sampling_type variable is an enum. The different values of the enum 
		can be accessed through:
		ESamplingType::ESpherical
		ESamplingType::EHemispherical
		ESamplingType::ECosineHemispherical
		*/
		
        // TODO(A3): Implement this
		// 1. Intersect the eye rays with the scene geometry
		SurfaceInteraction i;
		bool hit = scene.bvh->intersect(ray, i);
		// If there is an intersection
		if (hit) {
			// Compute the MC AO estimate at the shading point. 
			v2f sample = sampler.next2D();
			// Determine which warping function to use. 
			ESamplingType m_samping_type;

			// TODO: Figure out how to dynamically use sampling techniques. 
			//i.wi = Warp::squareToUniformSphere(sample);
			i.wi = Warp::squareToUniformHemisphere(sample);
			//i.wi = Warp::squareToCosineHemisphere(sample);
			
			// Construct the shadow ray
			// First retrieve the bounding sphere and its radius.

			// Use a half of this value as the max traverse distance of the shadow ray.
			BSphere BSphere = scene.aabb.getBSphere();
			float maxShadowRayLength = BSphere.radius / 2;
			Ray shadowRay(i.p, glm::normalize(i.frameNs.toWorld(i.wi)));
			shadowRay.max_t = maxShadowRayLength;
			
			bool visible = !scene.bvh->intersect(shadowRay, i);
			// Calculate the ambient occlusion
			// refer to page ~84 on direct illumination I & II
			if (visible) {
				// albedo rho is set to 1 for ambient occlusion. 
				float albedo = 1.f;

				if (Frame::cosTheta(i.wi) >= 0) {
					//Li = v3f(albedo) * Frame::cosTheta(i.wi) / (M_PI * pdf);
					Li = v3f(albedo) * Frame::cosTheta(i.wi) / M_PI;
				}
				else {
					Li = v3f(0.f);
				}
			}
		}
		//float pdf = Warp::squareToUniformSpherePdf();
		float pdf = Warp::squareToUniformHemispherePdf(Li);
		//float pdf = Warp::squareToCosineHemispherePdf(i.wi);

        return Li / pdf;
    }
};

TR_NAMESPACE_END