#include "Material.h"
Material::Material(ShaderProgram* shaderProgram)
{
    shaderProgram->IncRefCount();
    m_shaderProgram = shaderProgram;
}

Material::~Material()
{
    if(m_shaderProgram != nullptr)
        m_shaderProgram->DecRefCount();

    for(size_t i = 0, textureLen = m_textures.size(); i < textureLen; i++)
    {
        m_textures[i]->DecRefCount();
    }

}

void Material::SetTexture(char* name, Texture* texture)
{
    m_shaderProgram->Bind();

    GLint uniform = glGetUniformLocation(m_shaderProgram->GetGLShaderProgram(), name);

    if(uniform == -1)
    {
        //std::cout << "Uniform: " << name << " not found in shader program." << std::endl;
        return;
    }

    texture->IncRefCount();

    for(size_t i = 0, uniformsLen = m_textureUniforms.size(); i < uniformsLen; i++)
    {
        if(m_textureUniforms[i] == uniform)
        {
            m_textures[i]->DecRefCount();
            m_textures[i] = texture;
            return;
        }
    }
    m_textureUniforms.push_back(uniform);
    m_textures.push_back(texture);
}

void Material::SetMatrix(char* name, glm::mat4 matrix)
{
    m_shaderProgram->Bind();

    GLuint uniform = glGetUniformLocation(m_shaderProgram->GetGLShaderProgram(), name);

    if(uniform == -1)
    {
        //std::cout << "Uniform: " << name << " not found in shader program." << std::endl;
        return;
    }

    for(size_t i = 0, matrixUniformsLen = m_matrixUniforms.size(); i < matrixUniformsLen; i++)
    {
        if(m_matrixUniforms[i] == uniform)
        {
            m_matrices[i] = matrix;
            return;
        }
    }
    m_matrixUniforms.push_back(uniform);
    m_matrices.push_back(matrix);
}

void Material::Bind()
{
    m_shaderProgram->Bind();
    for(size_t i = 0, textureUniformsLen = m_textureUniforms.size(); i < textureUniformsLen; i++)
    {
        glActiveTexture(GL_TEXTURE0 + i);
        glBindTexture(GL_TEXTURE_2D, m_textures[i]->GetGLTexture());
        glUniform1i(m_textureUniforms[i], i);
    }

    for(size_t i = 0, matrixUniformsLen = m_matrixUniforms.size(); i < matrixUniformsLen; i++)
    {
        glUniformMatrix4fv(m_matrixUniforms[i], 1, GL_FALSE, &(m_matrices[i][0][0]));
    }

}

void Material::UnBind()
{
    for(size_t i = 0, textureUniformsLen = m_textureUniforms.size(); i  < textureUniformsLen; i++)
    {
        glActiveTexture(GL_TEXTURE0 + i);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    m_shaderProgram->Unbind();
}