#include "Buffer.h"

using namespace std;

Buffer::Buffer(const std::vector<size_t>& attrSizes, const GLfloat* vertices, size_t vtxSize, const GLuint* indices, size_t idxSize) {
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, vtxSize, (const void*)vertices, GL_STATIC_DRAW);

	size_t stride = 0, numFloat = 0;
	vector<size_t> offsets(attrSizes.size());
	for (size_t i = 0; i < attrSizes.size(); ++i) {
		offsets[i] = stride;
		size_t x = attrSizes[i] * sizeof(GLfloat);
		stride += x;
		numFloat += attrSizes[i];
	}
	numVert = vtxSize / sizeof(GLfloat) / numFloat;

	for (size_t i = 0; i < attrSizes.size(); ++i) {
		glVertexAttribPointer(i, attrSizes[i], GL_FLOAT, GL_FALSE, stride, (const void*)offsets[i]);
		glEnableVertexAttribArray(i);
	}

	if (indices != nullptr) {
		useElement = true;
		glGenBuffers(1, &ebo);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, idxSize, (const void*)indices, GL_STATIC_DRAW);
		// Do NOT unbind the EBO while a VAO is active as the bound element buffer object IS stored in the VAO; keep the EBO bound.
		//glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		numIdx = idxSize / sizeof(GLuint);
	} else {
		useElement = false;
		numIdx = 0;
	}

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}

void Buffer::bind() {
	glBindVertexArray(vao);
}

void Buffer::unbind() {
	glBindVertexArray(0);
}

void Buffer::draw(GLenum mode) {
	bind();
	if (useElement)
		glDrawElements(mode, numIdx, GL_UNSIGNED_INT, 0);
	else
		glDrawArrays(mode, 0, numVert);
}

Buffer::~Buffer() {
	glDeleteVertexArrays(1, &vao);
	glDeleteBuffers(1, &vbo);
	if (useElement)
		glDeleteBuffers(1, &ebo);
}
