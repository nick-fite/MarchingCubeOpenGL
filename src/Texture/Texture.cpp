#include "Texture.h"
Texture::Texture(char* filePath)
{
    FIBITMAP* bitmap = FreeImage_Load(FreeImage_GetFileType(filePath), filePath);
    FIBITMAP* bitmap32 = FreeImage_ConvertTo32Bits(bitmap);

    glGenTextures(1, &m_texture);

    glBindTexture(GL_TEXTURE_2D, m_texture);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, FreeImage_GetWidth(bitmap32), FreeImage_GetHeight(bitmap32), 0, GL_BGRA, GL_UNSIGNED_BYTE, static_cast<void*>(FreeImage_GetBits(bitmap32)));

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glBindTexture(GL_TEXTURE_2D, 0);

    FreeImage_Unload(bitmap);
    FreeImage_Unload(bitmap32);

}

Texture::~Texture()
{
    glDeleteTextures(1, &m_texture);
}

void Texture::IncRefCount() 
{
    m_refCount++;

}

void Texture::DecRefCount()
{
    m_refCount--;
    if(m_refCount == 0)
    {
        delete this;
    }
}

GLuint Texture::GetGLTexture()
{
    return m_texture;

}