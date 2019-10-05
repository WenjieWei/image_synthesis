/**
 *  Fragment shader for the phong reflection model. 
 *
 *  Shading logic combines the functionality of PhongBRDF::eval() and
 *  SimpleIntegrator. 
 *  
 *  Use the glUniform variables created in renderpasses/simple.h
 *  Use the logic from bsdf/phong.h
 *  All computations are still done in world space. 
 */

#version 330 core

// Define shader uniforms. 
uniform vec3 camPos;			// camera position in world space. 
uniform vec3 lightPos;          // light position in world space. 
uniform vec3 lightIntensity;
uniform vec3 rho_d;             // diffuse albedo.
uniform vec3 rho_s;             // specular albedo. 
uniform float exponent;          // phong exponent. 

in vec3 vNormal;
in vec3 vPosition;
out vec3 color;

// Define the value of pi. 
#define M_PI       3.14159265358979323846f

void main() {
    // Compute the distance from the light source to the intersection point. 
    // Calculate the intensity decay first. 
    float dist = distance(lightPos, vPosition);
    vec3 energy = lightIntensity / pow(dist, 2);

    // Compute the incident and view vectors.
    vec3 wo = normalize(camPos - vPosition);
    vec3 wi = normalize(lightPos - vPosition);

    // Calculate the cosine of the incident angle. 
    // Cosine is calculated by: cos(theta) = a * b / (|a||b|).
    float cosTheta = dot(wi, vNormal);

    // Compute the perfect specular reflection direction wr. 
    // wr = 2n(n * v) - v.
    vec3 wr = 2 * vNormal * dot(vNormal, wi) - wi;

    // Compute the cosine forshortening factor of cosine alpha. 
    // alpha is the angle between wr and lighting direction.    
    // Overwrite this term with cos^n(alpha)
    float cosAlpha = dot(wo, wr) / (length(wo) * length(wr));
    cosAlpha = pow(cosAlpha, exponent);

    if (cosAlpha > 0) {
        color = (rho_d / M_PI + (rho_s * (exponent + 2) * cosAlpha) / (2 * M_PI)) * cosTheta;
        color = color * energy;
    } else {
        color = (rho_d / M_PI * cosTheta);
        color = color * energy;
    }
}
