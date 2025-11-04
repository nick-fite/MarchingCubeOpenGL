#include <glm/glm.hpp>
#include "../Transform/Transform3D.h"


class BoundingBox {
    glm::vec3 min, max, center;

public:
    BoundingBox(const glm::vec3& min, const glm::vec3& max);
    BoundingBox ApplyScale(Transform3D t);
    glm::vec3 GetCenter();
    glm::vec3 getSize();
    bool isInside(glm::vec3 point);
    void Draw() const;

};