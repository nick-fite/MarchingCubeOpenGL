#pragma once
#include <gl/glew.h>
#include <GLFW/glfw3.h>
#include "../Shaders/Shader.h"
#include <FreeImage.h>

class Texture {
private:
    GLuint m_texture;
    unsigned int m_refCount = 0;
public:
    Texture(char* filePath);
    ~Texture();
    void IncRefCount();
    void DecRefCount();
    GLuint GetGLTexture();
};