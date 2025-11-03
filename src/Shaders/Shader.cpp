
#include "Shader.h"

Shader::Shader(std::string filepath, GLenum shaderType)
{
    InitFromFile(filepath, shaderType);
}

Shader::~Shader()
{
    if (m_shader != 0)
    {
        glDeleteShader(m_shader);
    }
}


GLuint Shader::GetGLShader() 
{
    return m_shader;
}

GLenum Shader::GetGLShaderType()
{
    return m_type;
}
    
bool Shader::InitFromFile(std::string filePath, GLenum shaderType)
{
    std::ifstream file(filePath);

    if(!file.good())
    {
        std::cout << "can't read: " << filePath << std::endl;
        return false;
    }

    //gets end of file
    file.seekg(0, std::ios::end);

    //makes a string the size of the file
    std::string shaderCode;
    shaderCode.resize((size_t)file.tellg());

    //goes to beginning
    file.seekg(0, std::ios::beg);

    //read it
    file.read(&shaderCode[0], shaderCode.size());
    
    //close it
    file.close();

    return InitFromString(shaderCode, shaderType);
}
bool Shader::InitFromString(std::string shaderCode, GLenum shaderType)
{
    m_type = shaderType;
    m_shader = glCreateShader(shaderType);

    const char* shaderCodePointer = shaderCode.data();
    int shaderCodeLength = shaderCode.size();

    glShaderSource(m_shader, 1, &shaderCodePointer, &shaderCodeLength);
    glCompileShader(m_shader);

    GLint isCompiled;

    //check if fragment shader can cimpile
    //print error and delete if fails
    glGetShaderiv(m_shader, GL_COMPILE_STATUS, &isCompiled);

    if(!isCompiled)
    {
        char log[1024];
        glGetShaderInfoLog(m_shader, 1024, NULL, log);
        std::cout << "Shader compile found with error: " << std::endl << log << std::endl;

        glDeleteShader(m_shader);
        m_shader = 0;
        return false;
    }
    return true;
}
void Shader::IncRefCount()
{
    m_refCount++;
}

void Shader::DecRefCount()
{
    m_refCount--;
}