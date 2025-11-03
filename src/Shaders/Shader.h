#pragma once

#include <gl/glew.h>
#include <glfw/glfw3.h>
#include <string>
#include <iostream>
#include <fstream>
class Shader {
private:
    GLuint m_shader;
    GLenum m_type;

    unsigned int m_refCount = 0;

public:
    Shader(std::string filepath, GLenum shaderType);
    ~Shader();

    GLuint GetGLShader();
    GLenum GetGLShaderType();

    bool InitFromFile(std::string filePath, GLenum shaderType);
    bool InitFromString(std::string shaderCode, GLenum shaderType);

    void IncRefCount();
    void DecRefCount();
};