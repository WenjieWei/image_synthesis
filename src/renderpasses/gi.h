/*
	This file is part of TinyRender, an educative rendering system.

	Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
	Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <core/core.h>
#include "core/renderpass.h"
#include "tiny_obj_loader.h"
#include "integrators/path.h"

TR_NAMESPACE_BEGIN

/**
 * Global Illumination baking renderpass.
 */
	struct GIPass : RenderPass {
	GLuint shader{ 0 };

	GLuint modelMatUniform{ 0 };
	GLuint viewMatUniform{ 0 };
	GLuint projectionMatUniform{ 0 };

	int m_samplePerVertex;

	std::unique_ptr<PathTracerIntegrator> m_ptIntegrator;

	explicit GIPass(const Scene& scene) : RenderPass(scene) {
		m_ptIntegrator = std::unique_ptr<PathTracerIntegrator>(new PathTracerIntegrator(scene));
		m_ptIntegrator->m_maxDepth = scene.config.integratorSettings.gi.maxDepth;
		m_ptIntegrator->m_rrProb = scene.config.integratorSettings.gi.rrProb;
		m_ptIntegrator->m_rrDepth = scene.config.integratorSettings.gi.rrDepth;
		m_samplePerVertex = scene.config.integratorSettings.gi.samplesByVertex;
	}

	virtual void buildVBO(size_t objectIdx) override {
		GLObject& obj = objects[objectIdx];

		// TODO(A5): Implement this
		// Set the shading point directly to the scene vertex location.
		obj.nVerts = scene.getObjectNbVertices(objectIdx);
		// Resize vertices vector to contain all attributes of all vertices
		//obj.vertices.resize(obj.nVerts * N_ATTR_PER_VERT);

		// Iterate over all vertices.
		for (int vertexIdx = 0; vertexIdx < obj.nVerts; vertexIdx++) {
			SurfaceInteraction i = SurfaceInteraction();

			// Define all the parameters of i. 
			// First get the position and normal of the current vertex.
			v3f vPos = scene.getObjectVertexPosition(objectIdx, vertexIdx);
			v3f vNorm = scene.getObjectVertexNormal(objectIdx, vertexIdx);
			Frame nFrame = Frame(glm::normalize(vNorm));

			// List of attributes of surface interaction:
			// v3f p, wo, wi
			// float t, u, v
			// size_t shapeID, primID
			// Frame frameNs, frameNg
			// int matID;
			// unsigned int sampledComponent, sampledType
			// Create geometry and shading frames. 
			i.frameNg = nFrame;
			i.frameNs = nFrame;
			// Create arbitrary wo wi and shift the vertex position. 
			i.wo = v3f(0.f, 0.f, 1.f);
			i.p = vPos + vNorm * Epsilon;
			// Assign primitive, shape, and material ID
			i.primID = scene.getPrimitiveID(vertexIdx);
			i.shapeID = objectIdx;
			i.matID = scene.getMaterialID(objectIdx, i.primID);

			// Start tracing the ray from the vertex we want to shade
			Ray ray(i.p, vNorm);
			Sampler sampler = Sampler(260685967);

			// Compute RGB values
			v3f rgb(0.f);
			for (int j = 0; j < m_samplePerVertex; j++) {
				rgb += m_ptIntegrator->renderExplicit(ray, sampler, i);
			}
			rgb /= m_samplePerVertex;

			// Assign the position and RGB data to objVertices. 
			obj.vertices.push_back(vPos.x);
			obj.vertices.push_back(vPos.y);
			obj.vertices.push_back(vPos.z);
			obj.vertices.push_back(rgb.x);
			obj.vertices.push_back(rgb.y);
			obj.vertices.push_back(rgb.z);
		}

		// VBO
		glGenVertexArrays(1, &obj.vao);
		glBindVertexArray(obj.vao);

		glGenBuffers(1, &obj.vbo);
		glBindBuffer(GL_ARRAY_BUFFER, obj.vbo);
		glBufferData(GL_ARRAY_BUFFER,
			sizeof(GLfloat) * obj.nVerts * N_ATTR_PER_VERT,
			(GLvoid*)(&obj.vertices[0]),
			GL_STATIC_DRAW);
	}

	bool init(const Config& config) override {
		RenderPass::init(config);

		// Create shader
		GLuint vs = compileShader("gi.vs", GL_VERTEX_SHADER);
		GLuint fs = compileShader("gi.fs", GL_FRAGMENT_SHADER);
		shader = compileProgram(vs, fs);
		glDeleteShader(vs);
		glDeleteShader(fs);

		// Create uniforms
		modelMatUniform = GLuint(glGetUniformLocation(shader, "model"));
		viewMatUniform = GLuint(glGetUniformLocation(shader, "view"));
		projectionMatUniform = GLuint(glGetUniformLocation(shader, "projection"));

		// Create vertex buffers
		objects.resize(scene.worldData.shapes.size());
		for (size_t i = 0; i < objects.size(); i++) {
			buildVBO(i);
			buildVAO(i);
		}

		return true;
	}

	void cleanUp() override {
		// Delete vertex buffers
		for (size_t i = 0; i < objects.size(); i++) {
			glDeleteBuffers(1, &objects[i].vbo);
			glDeleteVertexArrays(1, &objects[i].vao);
		}

		RenderPass::cleanUp();
	}

	void render() override {
		glBindFramebuffer(GL_FRAMEBUFFER, postprocess_fboScreen);
		glClearColor(0.f, 0.f, 0.f, 1.f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glEnable(GL_DEPTH_TEST);

		// TODO(A5): Implement this
		// shader is already defined in init().
		glUseProgram(shader);

		// Update camera
		glm::mat4 model, view, projection;
		camera.Update();
		camera.GetMatricies(projection, view, model);

		// Pass uniforms
		glUniformMatrix4fv(modelMatUniform, 1, GL_FALSE, &(modelMat[0][0]));
		glUniformMatrix4fv(viewMatUniform, 1, GL_FALSE, &(view[0][0]));
		glUniformMatrix4fv(projectionMatUniform, 1, GL_FALSE, &(projection[0][0]));

		// Draw
		for (auto& object : objects) {
			// Bind vertex array of current object/
			// Draw its triangles.
			// Bind vertex array to 0.
			glBindVertexArray(object.vao);
			glDrawArrays(GL_TRIANGLES, 0, object.nVerts);
			glBindVertexArray(0);
		}

		RenderPass::render();
	}

};

TR_NAMESPACE_END
