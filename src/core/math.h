/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

inline float safeSqrt(float v){
    return std::sqrt(std::max(float(0), v));
}

/**
 * Computes barycentric coordinates.
 */
template<class T>
inline T barycentric(const T& a, const T& b, const T& c, const float u, const float v) {
    return a * (1 - u - v) + b * u + c * v;
}

/**
 * Restricts a value to a given interval.
 */
template<class T>
inline T clamp(T v, T min, T max) {
    return std::min(std::max(v, min), max);
}

/**
 * Checks if vector is zero.
 */
inline bool isZero(const v3f v) {
    return glm::dot(v, v) < Epsilon;
}

/**
 * Generates coordinate system.
 */
inline void coordinateSystem(const v3f& a, v3f& b, v3f& c) {
    if (std::abs(a.x) > std::abs(a.y)) {
        float invLen = 1.f / std::sqrt(a.x * a.x + a.z * a.z);
        c = v3f(a.z * invLen, 0.f, -a.x * invLen);
    } else {
        float invLen = 1.f / std::sqrt(a.y * a.y + a.z * a.z);
        c = v3f(0.f, a.z * invLen, -a.y * invLen);
    }
    b = glm::cross(c, a);
}

/**
 * Converts RGB value to luminance.
 */
inline float getLuminance(const v3f& rgb) {
    return glm::dot(rgb, v3f(0.212671f, 0.715160f, 0.072169f));
}

/**
 * Pseudo-random sampler (Mersenne Twister 19937) structure.
 */
struct Sampler {
    std::mt19937 g;
    std::uniform_real_distribution<float> d;
    explicit Sampler(int seed) {
        g = std::mt19937(seed);
        d = std::uniform_real_distribution<float>(0.f, 1.f);
    }
    float next() { return d(g); }
    p2f next2D() { return {d(g), d(g)}; }
    void setSeed(int seed) {
        g.seed(seed);
        d.reset();
    }
};

/**
 * 1D discrete distribution.
 */
struct Distribution1D {
    std::vector<float> cdf{0};
    bool isNormalized = false;

    inline void add(float pdfVal) {
        cdf.push_back(cdf.back() + pdfVal);
    }

    size_t size() {
        return cdf.size() - 1;
    }

    float normalize() {
        float sum = cdf.back();
        for (float& v : cdf) {
            v /= sum;
        }
        isNormalized = true;
        return sum;
    }

    inline float pdf(size_t i) const {
        assert(isNormalized);
        return cdf[i + 1] - cdf[i];
    }

    int sample(float sample) const {
        assert(isNormalized);
        const auto it = std::upper_bound(cdf.begin(), cdf.end(), sample);
        return clamp(int(distance(cdf.begin(), it)) - 1, 0, int(cdf.size()) - 2);
    }
};


/**
 * Warping functions.
 */
namespace Warp {


inline v3f squareToUniformSphere(const p2f& sample) {
    v3f v(0.f);
    // TODO(A3): Implement this
	// p42 of 03-Monte-Carlo-I
	// zeta_1 is the first rand, x coordinates of the sample
	// zeta_2 is the second rand, y coordinates of the sample.
	v.z = 2 * sample.x - 1;

	float r = sqrt(1 - pow(v.z, 2));
	float phi = 2 * M_PI * sample.y;

	v.x = r * cos(phi);
	v.y = r * sin(phi);

    return v;
}

inline float squareToUniformSpherePdf() {
    float pdf = 0.f;
    // TODO(A3): Implement this
	// pdf * surface area = 1.
	// area = 4 * pi * r^2. 
	// uniform sphere, so r = 1
	float r = 1.f;
	pdf = 1 / (4 * M_PI * pow(r, 2));

    return pdf;
}

inline v3f squareToUniformHemisphere(const p2f& sample) {
    v3f v(0.f);
    // TODO(A3): Implement this
	// Hemisphere is just a half of the sphere. 
	// hemisphere is aligned with the positive z axis, so simply take the inverse if z is negative. 
	v.z = 2 * sample.x - 1;
	if (v.z < 0) {
		v.z = -v.z;
	}

	float r = sqrt(1 - pow(v.z, 2));
	float phi = 2 * M_PI * sample.y;

	v.x = r * cos(phi);
	v.y = r * sin(phi);

    return v;
}

inline float squareToUniformHemispherePdf(const v3f& v) {
    float pdf = 0.f;
    // TODO(A3): Implement this
	// pdf * surface area = 1.
	// surface area = 2 * pi * r^2. 
	// uniform sphere, so r = 1
	float r = 1.f;
	pdf = 1 / (2 * M_PI * pow(r, 2));

    return pdf;
}

inline v2f squareToUniformDiskConcentric(const p2f& sample) {
    v2f v(0.f);
    // TODO(A3): Implement this
	// Using Nusselt Analog. 
	// PBRT 13.6.2
	float r, phi;

	r = sample.x;
	phi = 2 * M_PI * sample.y;

	v.x = cos(phi) * sqrt(r);
	v.y = sin(phi) * sqrt(r);

    return v;
}

inline v3f squareToCosineHemisphere(const p2f& sample) {
    v3f v(0.f);
    // TODO(A3): Implement this
	// Warp to the uniform disk
	v2f p = squareToUniformDiskConcentric(sample);
	v.x = p.x;
	v.y = p.y;
	v.z = sqrt(1 - pow(p.x, 2) - pow(p.y, 2));

    return v;
}

inline float squareToCosineHemispherePdf(const v3f& v) {
    float pdf = 0.f;
    // TODO(A3): Implement this
	// PBRT: Returns a weight of cos(theta) / pi. 
	// cos(theta) is the value of v.z
	pdf = v.z / M_PI;

    return pdf;
}

inline v3f squareToPhongLobe(const p2f& sample, float exponent) {
    v3f v(0.f);
    // TODO(A3): Implement this
	// Formula taken from: 
	// http://www.cim.mcgill.ca/~derek/ecse689_a3.html?fbclid=IwAR1Wf5g0tY7ElS8COAVg7-QDorbvBCg9PIy3ftH3P3l_SqsiLNyPvVxP_k0
	float phi, theta;

	phi = 2 * M_PI * sample.y;
	theta = acos(pow(sample.x, 1 / (exponent + 1)));

	v.x = sin(theta) * cos(phi);
	v.y = sin(theta) * cos(phi);
	v.z = cos(theta);

    return v;
}

inline float squareToPhongLobePdf(const v3f& v, float exponent) {
    float pdf = 0.f;
    // TODO(A3): Implement this
	float theta = acos(v.z / sqrt(pow(v.x, 2) + pow(v.y, 2) + pow(v.z, 2)));
	pdf = (exponent + 2) * pow(cos(theta), exponent) / (2 * M_PI);

    return pdf;
}

inline v2f squareToUniformTriangle(const p2f& sample) {
    v2f v(0.f);
    float u = std::sqrt(1.f - sample.x);
    v = {1 - u, u * sample.y};
    return v;
}

inline v3f squareToUniformCone(const p2f& sample, float cosThetaMax) {
    v3f v(0.f);
    // TODO(A3): Implement this
    return v;
}

inline float squareToUniformConePdf(float cosThetaMax) {
    float pdf = 0.f;
    // TODO(A3): Implement this
    return pdf;
}

}

TR_NAMESPACE_END