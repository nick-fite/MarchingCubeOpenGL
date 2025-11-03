#pragma once
#include "../Shader.h"
#include <iostream>

class ShaderProgram
{
private:
    Shader* m_vertShader = nullptr;
    Shader* m_fragShader = nullptr;

    GLuint m_shaderProgram;

    bool m_built = false;

    unsigned int m_refCount = 0;
public:
    ShaderProgram();
    ~ShaderProgram();
    GLuint GetGLShaderProgram();
    void AttachShader(Shader* shader);
    void Bind();
    void Unbind();
    void IncRefCount();
    void DecRefCount();
};