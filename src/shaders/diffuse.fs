/**
 *  Fragment shader for the diffuse material. 
 *
 *  Shading logic combines the functionality of DiffuseBRDF::eval() and
 *  SimpleIntegrator. 
 *  
 *  Use the glUniform variables created in renderpasses/simple.h.
 */

#version 330 core

// Define shader uniforms
uniform vec3 cameraPos;         // camera position in world space
uniform vec3 lightPos;          // light position in world space. 
uniform vec3 lightIntensity;
uniform vec3 albedo;

in vec3 vNormal;
in vec3 vPosition;
out vec3 color;

// Define the value of pi. 
#define M_PI       3.14159265358979323846f

void main() {
    
}
