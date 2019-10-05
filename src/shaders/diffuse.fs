/**
 *  Fragment shader for the diffuse reflection model. 
 *
 *  Shading logic combines the functionality of DiffuseBRDF::eval() and
 *  SimpleIntegrator. 
 *  
 *  Use the glUniform variables created in renderpasses/simple.h
 *  Use the logic from bsdf/diffuse.h
 *  All computations are still done in world space. 
 */

#version 330 core

// Define shader uniforms
uniform vec3 camPos;			// camera position in world space
uniform vec3 lightPos;          // light position in world space. 
uniform vec3 lightIntensity;
uniform vec3 albedo;

in vec3 vNormal;
in vec3 vPosition;
out vec3 color;

// Define the value of pi. 
#define M_PI       3.14159265358979323846f

void main() {
    // val = (albedo->eval(worldData, i) / M_PI) * Frame::cosTheta(glm::normalize(i.wi));
    // Compute the distance from the light source to the intersection point. 
    // Calculate the intensity decay: I / r^2. 
    float dist = distance(lightPos, vPosition);
    vec3 energy = lightIntensity / pow(dist, 2);
    
    // Compute the incident ray vector to the intersection point. 
    // Intersection point = vertex position?
    vec3 wi = normalize(lightPos - vPosition);

    // Calculate the cosine of the incident angle. 
    // Cosine is calculated by: cos(theta) = a * b / (|a||b|).
    float cosTheta = dot(wi, vNormal);

    // Calculate the color. 
    color = energy * albedo * cosTheta / M_PI;
}
