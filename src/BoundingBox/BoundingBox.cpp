#include "BoundingBox.h"
#include <gl/glew.h>

BoundingBox::BoundingBox(const glm::vec3& newMin, const glm::vec3& newMax): min(newMin), max(newMax)
{
   this->center = 
    glm::vec3(
        min.x + (max.x - min.x),
        min.y + (max.y - min.y),
        min.z + (max.z - min.z)
    ); 
}


BoundingBox BoundingBox::ApplyScale(Transform3D t)
{
    glm::vec3 newSize = t.GetMatrix() * glm::vec4(this->getSize(), 1.0);
    glm::vec3 newMin = this->GetCenter() - newSize / 2.0f;
    glm::vec3 newMax = this->GetCenter() + newSize / 2.0f;
    return BoundingBox(newMin, newMax);

}

glm::vec3 BoundingBox::GetCenter()
{
    return this->center;
}

glm::vec3 BoundingBox::getSize()
{
    return this->max - this->min;

}

bool BoundingBox::isInside(glm::vec3 point)
{
    return (
        min.x <= point.x && point.x <= max.x &&
        min.y <= point.y && point.y <= max.y &&
        min.z <= point.z && point.z <= max.z
    );
}

void BoundingBox::Draw() const {
    glColor3f(1.0f, 0.0f, 0.0f); // RGB: Red
    // Define the 8 corners of the bounding box
    glm::vec3 corners[8] = {
        min, 
        glm::vec3(max.x, min.y, min.z), 
        glm::vec3(max.x, max.y, min.z), 
        glm::vec3(min.x, max.y, min.z), 
        glm::vec3(min.x, min.y, max.z), 
        glm::vec3(max.x, min.y, max.z), 
        max,
        glm::vec3(min.x, max.y, max.z) 
    };

    // Define the edges of the bounding box (pairs of indices into the corners array)
    unsigned int edges[24] = {
        0, 1, 1, 2, 2, 3, 3, 0, // Front face
        4, 5, 5, 6, 6, 7, 7, 4, // Back face
        0, 4, 1, 5, 2, 6, 3, 7  // Connecting edges
    };

    glBegin(GL_LINES);
    for (int i = 0; i < 24; i += 2) {
        glVertex3fv(&corners[edges[i]].x);
        glVertex3fv(&corners[edges[i + 1]].x);
    }
    glEnd();
}
