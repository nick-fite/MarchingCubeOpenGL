#pragma once
#include <glm/gtc/matrix_transform.hpp>

//yes, header only files are not usually the best practice. However, given the nature of this tool,
//it makes more sense to me that it's just a header. No need to overcomplicate it.
class Transform3D {
private: 
    float scale;
    glm::vec3 rotation;
    glm::vec3 position;

    bool matrixDirty;
    bool inverseDirty;

    glm::mat4 rotMatrix;
    glm::mat4 matrix;
    glm::mat4 inverseMatrix;

public:
    Transform3D() {
        scale = 1;
        rotation = glm::vec3();
        position = glm::vec3();
        matrix = inverseMatrix = glm::mat4();
    };
    
    float getScale() const {
        return scale;
    };

    glm::vec3 getRotation() const {
        return rotation;
    };

    glm::vec3 getPosition() const {
        return position;
    };

    void SetScale(float scale) {
        this->scale = scale;
        matrixDirty = true;
    };

    void SetPosition(glm::vec3 pos) {
        position = pos;
        matrixDirty = inverseDirty =  true;
    };
    void SetRotation(glm::vec3 rot) {
        rotation = rot;
        matrixDirty = inverseDirty = true;

    };

    void RotateX(float rot){
        rotation.x += rot;
        matrixDirty = inverseDirty =  true;
    };
    void RotateY(float rot){
        rotation.y += rot;
        matrixDirty = inverseDirty =  true;
    };
    void RotateZ(float rot){
        rotation.z += rot;
        matrixDirty = inverseDirty =  true;
    };

    void Translate(glm::vec3 vec) {
        position += vec;
        matrixDirty = inverseDirty =  true;
    };

    glm::mat4 GetMatrix() {
        //recalculate matrix if it's changed
        if(matrixDirty)
        {
            glm::mat4 s = glm::mat4(
                scale, 0,0,0,
                0, scale, 0,0,
                0, 0, scale, 0,
                0, 0, 0, scale
            );

            glm::mat4 rx = glm::mat4(
                1, 0, 0, 0,
                0, cos(rotation.x), sin(rotation.x), 0,
                0, -sin(rotation.x), cos(rotation.x), 0,
                0,0,0,1
            );

            glm::mat4 ry = glm::mat4(
                cos(rotation.y), 0, sin(rotation.y), 0,
                0, 1, 0, 0,
                -sin(rotation.y), 0, cos(rotation.y), 0,
                0,0,0,1
            );

            glm::mat4 rz = glm::mat4(
                cos(rotation.z), sin(rotation.z), 0, 0,
                -sin(rotation.z), cos(rotation.z), 0,0,
                0,0,1,0,
                0,0,0,1
            );

            glm::mat4 t = glm::mat4(
                1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                position.x, position.y, position.z, 1
            );

            rotMatrix = ry * rx * rz;
            matrix = t * rotMatrix * s;

            matrixDirty = false;

        }
        return matrix;
    };


    glm::mat4 GetInverseMatrix() {
        if(inverseDirty) {
            glm::mat4 s = glm::mat4(
                1.f/scale, 0,0,0,
                0, 1.f/scale, 0,0,
                0, 0, 1.f/scale, 0,
                0, 0, 0, 1.f/scale
            );

            glm::mat4 rx = glm::mat4(
                1, 0, 0, 0,
                0, cos(rotation.x), -sin(rotation.x), 0,
                0, sin(rotation.x), cos(rotation.x), 0,
                0,0,0,1
            );

            glm::mat4 ry = glm::mat4(
                cos(rotation.y), 0, -sin(rotation.y), 0,
                0, 1, 0, 0,
                sin(rotation.y), 0, cos(rotation.y), 0,
                0,0,0,1
            );

            glm::mat4 rz = glm::mat4(
                cos(rotation.z), -sin(rotation.z), 0, 0,
                sin(rotation.z), cos(rotation.z), 0,0,
                0,0,1,0,
                0,0,0,1
            );

            glm::mat4 t = glm::mat4(
                1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                -position.x, -position.y, -position.z, 1
            );

            glm::mat4 composedRot = rz * rx * ry;
            inverseMatrix = s * composedRot * t;
            inverseDirty = false;
        }
        return inverseMatrix;
    };

    glm::vec3 GetUp() {
        if(matrixDirty) {
            GetMatrix();
        }
        glm::vec4 up = rotMatrix * glm::vec4(0,1,0,1);

        return glm::vec3(up);
    };
    
    glm::vec3 GetForward() {
        if(matrixDirty) {
            GetMatrix();
        }
        glm::vec4 forward = rotMatrix * glm::vec4(0,0,-1,-1);

        return glm::vec3(forward);
    };
    
    glm::vec3 GetRight() {
        if(matrixDirty) {
            GetMatrix();
        }
        glm::vec4 right = rotMatrix * glm::vec4(1,0,0,1);

        return glm::vec3(right);
    };
};