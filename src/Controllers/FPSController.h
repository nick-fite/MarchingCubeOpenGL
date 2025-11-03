#pragma once
#include "../Transform/Transform3D.h"
#include <GLFW/glfw3.h>

class FPSController {
private:
    Transform3D m_transform;

public:
    FPSController();
    ~FPSController();
    Transform3D GetTransform();
    void Update(GLFWwindow* window, glm::vec2 viewportDimmensions, glm::vec2 mousePosition, float deltaTime);
};