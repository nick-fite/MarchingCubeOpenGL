
#include "ShadersProgram.h"
ShaderProgram::ShaderProgram()
{
    m_shaderProgram = glCreateProgram();
}
ShaderProgram::~ShaderProgram()
{
    glDeleteProgram(m_shaderProgram);

    if(m_vertShader != nullptr)
        m_vertShader->DecRefCount();
    if(m_fragShader != nullptr)
        m_fragShader->DecRefCount();
}
GLuint ShaderProgram::GetGLShaderProgram()
{
    return m_shaderProgram;
}

void ShaderProgram::AttachShader(Shader* shader)
{
    Shader** currentShader;
    switch(shader->GetGLShaderType())
    {
        case GL_VERTEX_SHADER:
            currentShader = &m_vertShader;
            break;
        
        case GL_FRAGMENT_SHADER:
            currentShader = &m_fragShader;
            break;
        default: return;
    }
    shader->IncRefCount();
    
    if(*currentShader != nullptr)
        (*currentShader)->DecRefCount();
    
    *currentShader = shader;

    if(shader->GetGLShader() != 0)
    {
        glAttachShader(m_shaderProgram, shader->GetGLShader());
        m_built = false;
    }
    else
    {
        std::cout << "failed to attachm not initialized" << std::endl;
    }
}
void ShaderProgram::Bind()
{
    if(!m_built)
    {
        glLinkProgram(m_shaderProgram);
        m_built = true;
    }
    glUseProgram(m_shaderProgram);
}

void ShaderProgram::Unbind()
{
    glUseProgram(0);
}

void ShaderProgram::IncRefCount()
{
    m_refCount++;
}

void ShaderProgram::DecRefCount()
{
    
    m_refCount--;
    if(m_refCount == 0)
        delete this;
}