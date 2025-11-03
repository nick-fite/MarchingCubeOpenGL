#include "FPSController.h"
FPSController::FPSController()
{
    m_transform = Transform3D();
}

FPSController::~FPSController()
{

}

Transform3D FPSController::GetTransform()
{
    return m_transform;
}

void FPSController::Update(GLFWwindow* window, glm::vec2 viewportDimmensions, glm::vec2 mousePosition, float deltaTime)
{
    glm::vec2 mouseMovement = mousePosition - (viewportDimmensions / 2.0f);

    float yaw = m_transform.getRotation().y;
    yaw += (int)mouseMovement.x * .001f;

    float pitch = m_transform.getRotation().x;
    pitch -= (int)mouseMovement.y * .001f;

    float halfpi = 3.1416f/ 2.f;
    if(pitch < -halfpi) pitch = -halfpi;
    else if (pitch > halfpi) pitch = halfpi;

    m_transform.SetRotation(glm::vec3(pitch,yaw,0));

    glfwSetCursorPos(window, mousePosition.x - mouseMovement.x, mousePosition.y - mouseMovement.y);

    if(glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) 
    {
        m_transform.Translate(m_transform.GetForward() * 50.f * deltaTime);
    }
    if(glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) 
    {
        m_transform.Translate(m_transform.GetRight() * -50.f * deltaTime);
    }
    if(glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) 
    {
        m_transform.Translate(m_transform.GetForward() * -50.f * deltaTime);
    }
    if(glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) 
    {
        m_transform.Translate(m_transform.GetRight() * 50.f * deltaTime);
    }
}
