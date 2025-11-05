#pragma once

#include <gl/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "../BoundingBox/BoundingBox.h"

struct MeshData {
    glm::vec3 position;
    glm::vec2 texCoords;
    glm::vec3 normal;

    MeshData(glm::vec3 pos, glm::vec2 tex, glm::vec3 norm)
        : position(pos), texCoords(tex), normal(norm) {}
    MeshData(){}
};



class Mesh {
public:
    Mesh(std::vector<MeshData> vertices, std::vector<unsigned short> indices);
    Mesh(std::string filePath, int meshNum);
    ~Mesh();
    
	std::vector<MeshData> m_vertices;
    std::vector<unsigned short> m_indices;
	BoundingBox innerBB = BoundingBox(glm::vec3(0), glm::vec3(0));
	BoundingBox outerBB = BoundingBox(glm::vec3(0), glm::vec3(0));

private:
    void MakeMesh(std::vector<MeshData> vertices, std::vector<unsigned short> indices);

    GLuint m_vertextBuffer;
    GLuint m_indexBuffer;
public:
    void DrawMesh();
};