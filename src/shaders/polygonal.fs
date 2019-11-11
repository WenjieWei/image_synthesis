/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#version 330 core


#define PI 3.14159265359
#define MAX_NUM_EMITTER_TRIANGLES 40 // Max number of emitter triangles allowed (tuned for A4)
uniform float emitterVertices[MAX_NUM_EMITTER_TRIANGLES*3*3]; // Need to specify a size at compile time (max size is 512 floats)

uniform int nbTriangles; // Use nbTriangles to index the buffer correctly
uniform vec3 lightIrradiance;
uniform vec3 albedo;
uniform vec2 windowSize; // [width, height]

uniform sampler2D cvTerm; // Creates a 2D texture we can sample from (range x,y = [0,1])

in vec3 vNormal;
in vec3 vPos;

out vec3 color;

// Compute edge (v1--v2) contribution
float getEdgeContrib(vec3 v1, vec3 v2, vec3 pos) {
	// Adapt your getEdgeContrib code from the offline part
	float value = 0.f;
	// TODO(A4): Implement this
	float Theta;
	vec3 Gamma;

	vec3 r1 = v1 - pos;
	vec3 r2 = v2 - pos;

	Theta = acos(dot(normalize(r1), normalize(r2)));
	Gamma = normalize(cross(r2, r1));
	value = Theta * dot(Gamma, vNormal);

	return value;
}


void main()
{	
	// 1) Extract vertices of triangles from `emitterVertices` buffer using `nbTriangles`
	// 2) Calculate G term
	// 3) Subtract modification term for G after extracting it from texture (use built-in `texture()` function)
	//	    e.g. `vec3 delta = texture(cvTerm, coords).xyz;`

	color = vec3(0);
	vec3 G = vec3(0);

    // TODO(A4): Implement this
	float contribution = 0;
	vec3 Ei = vec3(0);
	vec3 v1 = vec3(0); 
	vec3 v2 = vec3(0);

	// emitterVertices store the complete list of coordinate information. (x, y, z values in float)
	// first three variables: one vertex
	// first nine variables: one triangle. 
	for (int i = 0; i < nbTriangles; i++) {
		for (int j = 0; j < 3; j++) {
			if (j != 2){
				v1.x = emitterVertices[i * 9 + j * 3];
				v1.y = emitterVertices[i * 9 + j * 3 + 1];
				v1.z = emitterVertices[i * 9 + j * 3 + 2];

				v2.x = emitterVertices[i * 9 + j * 3 + 3];
				v2.y = emitterVertices[i * 9 + j * 3 + 4];
				v2.z = emitterVertices[i * 9 + j * 3 + 5];

			} else {
				// for example:
				// the v1 is the last vertex of the first triangle. 
				// x, y, z are 6, 7, 8
				// x, y, z for v2 are 0, 1, 2. 
				v1.x = emitterVertices[i * 9 + j * 3];
				v1.y = emitterVertices[i * 9 + j * 3 + 1];
				v1.z = emitterVertices[i * 9 + j * 3 + 2];

				v2.x = emitterVertices[i * 9 + j * 3 - 6];
				v2.y = emitterVertices[i * 9 + j * 3 - 5];
				v2.z = emitterVertices[i * 9 + j * 3 - 4];
			}

			contribution += getEdgeContrib(v1, v2, vPos);
		}
	}

	G = albedo * lightIrradiance * contribution / (2 * PI * PI);

	vec2 textureCoord = vec2(0);
	textureCoord.x = gl_FragCoord.x / windowSize.x;
	textureCoord.y = 1 - gl_FragCoord.y / windowSize.y;

	vec3 alpha = -texture(cvTerm, textureCoord).xyz;

	color = G + alpha;
}

