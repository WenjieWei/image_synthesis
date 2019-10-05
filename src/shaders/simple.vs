/**
 *	The new vertex shader file. 
 *	
 *	Start definitions are copied from normal.vs from A1 as-is. 
 */

// All shading computations to be computed in world space. 
#version 330 core

uniform mat4 model;         // Model matrix.
uniform mat4 view;          // View matrix.
uniform mat4 projection;    // Projection matrix. 
uniform mat4 normalMat;     // Normal matrix. 

layout(location = 0) in vec3 position;  // Vertex positions in world space
layout(location = 1) in vec3 normal;    // Vertex normals in world space. 

out vec3 vPosition;
out vec3 vNormal;

void main() {
	// High level: 
	// 1. Take Triangles,
	// 2. Figure out which pixels they cover. 
	// 3. Perform computations over the pixels. 

    gl_Position = projection * view * model * vec4(position, 1.f);
    // vNormal abs(normal);

    // Calculate the out variables. 
    vPosition = vec3(model * vec4(position, 1.f));
    vNormal = vec3(normalMat * model * vec4(normal, 0.f));
}
