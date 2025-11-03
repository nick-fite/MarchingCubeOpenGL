#pragma once
#include "../Shaders/ShadersProgram/ShadersProgram.h"
#include "../Texture/Texture.h"
#include <glm/gtc/matrix_transform.hpp>
#include <vector>

class Material {
private:
    ShaderProgram* m_shaderProgram = nullptr;

    std::vector<GLuint> m_textureUniforms;
    std::vector<Texture*> m_textures;

    std::vector<GLuint> m_matrixUniforms;
    std::vector<glm::mat4> m_matrices;
public:
    Material(ShaderProgram* shaderProgram);
    ~Material();
    void SetTexture(char* name, Texture* texture);
    void SetMatrix(char* name, glm::mat4 matrix);

    void Bind();
    void UnBind();
};