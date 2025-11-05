#include "Mesh.h"
#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>
#include <assimp/scene.h>
#include <assimp/DefaultLogger.hpp>
#include <assimp/LogStream.hpp>
#include <omp.h>

Mesh::Mesh(std::vector<MeshData> vertices, std::vector<unsigned short> indices)
{
    MakeMesh(vertices, indices);
}


Mesh::Mesh(std::string filePath, int meshNum)
{
    std::ifstream file(filePath);
    
    if(!file.good())
    {
        std::cout << "what the fuc\n"<< "bad file: " << filePath << std::endl;
        return;
    }

    Assimp::Importer importer;
    const aiScene* scene = NULL;
    unsigned int n = 0, t;

    for(int k = 0; k < 1; k++)
    {
        scene = importer.ReadFile(filePath, aiProcess_Triangulate | aiProcess_GenSmoothNormals |
                 aiProcess_FlipUVs | aiProcess_CalcTangentSpace |
                 aiProcess_GenBoundingBoxes);
    }

    //this what indicates how many meshes are being imported. If I take in multiple in the same fbx will this process it?
    struct aiMesh* mesh = scene->mMeshes[meshNum];

    // Convert aiVector3D to glm::vec3
    glm::vec3 minBounds(mesh->mAABB.mMin.x, mesh->mAABB.mMin.y, mesh->mAABB.mMin.z);
    glm::vec3 maxBounds(mesh->mAABB.mMax.x, mesh->mAABB.mMax.y, mesh->mAABB.mMax.z);
    innerBB = BoundingBox(minBounds, maxBounds);
    //outerBB = BoundingBox(minBounds, maxBounds);
    Transform3D t3d;
    t3d.SetScale(1.1f);
    outerBB = innerBB.ApplyScale(t3d);

    std::vector<float> positionsX;
    std::vector<float> positionsY;
    std::vector<float> positionsZ;

    for(t = 0; t < mesh->mNumVertices; ++t)
    {
        MeshData v;
        memcpy(&v.position, &mesh->mVertices[t], sizeof(glm::vec3));
        memcpy(&v.normal, &mesh->mNormals[t], sizeof(glm::vec3));
        memcpy(&v.texCoords, &mesh->mTextureCoords[0][t], sizeof(glm::vec2));

        m_vertices.push_back(v);

        positionsX.push_back(v.position.x);
        positionsY.push_back(v.position.y);
        positionsZ.push_back(v.position.z);
    }

    for(t = 0; t < mesh->mNumFaces; ++t)
    {
        const struct aiFace* face = &mesh->mFaces[t];

        for(int i = 0; i < 3; i++)
        {
            m_indices.push_back(face->mIndices[i]);
        }

        if(face->mNumIndices == 4)
        {
            m_indices.push_back(face->mIndices[0]);
            m_indices.push_back(face->mIndices[2]);
            m_indices.push_back(face->mIndices[3]);
        }
    }

    vertSizeX = innerBB.max.x - innerBB.min.x;
    vertSizeY = innerBB.max.y - innerBB.min.y;
    vertSizeZ = innerBB.max.z - innerBB.min.z;

    // Initialize minX, minY, minZ from the bounding box
    minX = innerBB.min.x;
    minY = innerBB.min.y;
    minZ = innerBB.min.z;

    // Store voxel sizes in class-level variables
    voxelSizeX = std::ceil(vertSizeX / voxelSize) * 4;
    voxelSizeY = std::ceil(vertSizeY / voxelSize) * 2;
    voxelSizeZ = std::ceil(vertSizeZ / voxelSize) * 4;

    // Debugging: Print initialized values
    std::cout << "Initialized minX: " << minX << ", minY: " << minY << ", minZ: " << minZ << std::endl;
    std::cout << "Stored voxelSizeX: " << voxelSizeX << ", voxelSizeY: " << voxelSizeY << ", voxelSizeZ: " << voxelSizeZ << std::endl;

    SizeX = voxelSizeX;
    SizeY = voxelSizeY;
    SizeZ = voxelSizeZ;

    glm::vec3 startPos = glm::vec3(voxelSizeX / 4 * voxelSize, voxelSizeY / 2 * voxelSize, 0);
    float surfaceProximityThreshold = voxelSize * 0.1f;

    std::cout << "Voxel grid size: " << SizeX << ", " << SizeY << ", " << SizeZ << std::endl;

    for(int x = 0; x <= SizeX; x++)
    {
        for(int y = 0; y <= SizeY; y++)
        {
            for(int z = 0; z <= SizeZ; z++)
            {
                glm::vec3 pos = glm::vec3(x, y, z);
                float posX = minX + ((float)x / voxelSizeX) * vertSizeX;
                float posY = minY + ((float)y / voxelSizeY) * vertSizeY;
                float posZ = minZ + ((float)z / voxelSizeZ) * vertSizeZ;

                float dist = closestTriangleDistance(glm::vec3(posX, posY, posZ));
                
                std::cout << "Generating voxel: " << x << ", " << y << ", " << z << " dist: " << dist << std::endl;
                bool isInside = isInsideMesh(pos);
                
                Voxels.push_back( isInside ? -dist : dist);

                //initSdf_TestData.push_back(MeshData(glm::vec3(posX, posY, posZ), glm::vec2(0,0), glm::vec3(0,0,0)));
                //Voxels.push_back(posX);
                //Voxels.push_back(posY);
                //Voxels.push_back(posZ);
            }
        }
    }
    
    GenerateMesh();

    std::vector<unsigned short> marchingCubesIndices;
    for(int i = 0; i < m_marchingCubesVertices.size(); i++) {
        marchingCubesIndices.push_back(i);
    }

    m_marchingCubesIndices = marchingCubesIndices;



    printf("verticies: %d, indicies: %d\n", (int)m_vertices.size(), (int)m_indices.size());
    MakeMesh(m_marchingCubesVertices, m_marchingCubesIndices);

}
void Mesh::MakeMesh(std::vector<MeshData> vertices, std::vector<unsigned short> indices)
{
    //vertices = vertices;
    //indices = indices;

    std::cout << "Making mesh with " << vertices.size() << " vertices and " << indices.size() << " indices." << std::endl;
    //vertex buffer
    glGenBuffers(1, &m_vertextBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, m_vertextBuffer);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(MeshData), &vertices[0], GL_STATIC_DRAW);

    std::cout << "Created vertex buffer with " << vertices.size() << " vertices." << std::endl;

    //indicies buffer
    glGenBuffers(1, &m_indexBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_indexBuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned short), &indices[0], GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

}

Mesh::~Mesh() 
{
    glDeleteBuffers(1, &m_vertextBuffer);
    glDeleteBuffers(1, &m_indexBuffer);
}

#define SetupAttribute(index, size, type, structure, element)\
    glVertexAttribPointer(index, size, type, 0, sizeof(structure), (void*)offsetof(structure, element));

void Mesh::DrawMesh()
{
    //bind buffers
    glBindBuffer(GL_ARRAY_BUFFER, m_vertextBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_indexBuffer);

    //setup attribute
    SetupAttribute(0,3,GL_FLOAT, MeshData, position);
    SetupAttribute(1,2,GL_FLOAT, MeshData, texCoords);
    SetupAttribute(2,3,GL_FLOAT, MeshData, normal);

    //turn on attributei
    for(int i = 0; i < 3; i++)
    {
        glEnableVertexAttribArray(i);
    }

    //draw it
    glDrawElements(GL_TRIANGLES, m_marchingCubesIndices.size(), GL_UNSIGNED_SHORT, (void*)0);

    //unbind buffer
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    //disable attribs
    for(int i = 0; i < 3; i++)
    {
        glDisableVertexAttribArray(i);
    }

    innerBB.Draw();
    outerBB.Draw();

}
float Mesh::closestTriangleDistance(glm::vec3 point)
{

    float max_float = std::numeric_limits<float>::max();

    for(int i = 0; i < m_indices.size(); i += 3)
    {
        glm::vec3 v0 = m_vertices[m_indices[i]].position;
        glm::vec3 v1 = m_vertices[m_indices[i + 1]].position;
        glm::vec3 v2 = m_vertices[m_indices[i + 2]].position;

        float dist = glm::distance(point, ClosestPointOnTriangleToPoint(v0, v1, v2, point));
        //std::cout << "voxel at " << point.x << ", " << point.y << ", " << point.z << " dist: " << dist << std::endl;
        max_float = std::min(max_float, dist);
    }

    return max_float;
}

glm::vec3 Mesh::ClosestPointOnTriangleToPoint( glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 point )
{
    glm::vec3 AB = b - a;
    glm::vec3 AC = c - a;
    glm::vec3 Normal = glm::normalize(glm::cross(AB, AC));
    glm::vec3 AP = point - a;
    float DistanceToPlane = glm::dot(AP, Normal);
    glm::vec3 ProjectedPoint = point - DistanceToPlane * Normal;

glm::vec3 Closest = ProjectedPoint; // Start with the projected point
    bool Inside = true; // Assume inside initially

    // Check if the point is "outside" any of the half-planes defined by the edges
    auto check_edge = [&](const glm::vec3& p1, const glm::vec3& p2) {
        glm::vec3 Edge = p2 - p1;
        glm::vec3 PointToEdgeStart = point - p1;
        // Project PointToEdgeStart onto Edge, clamped to the segment length
        float t = glm::clamp(glm::dot(PointToEdgeStart, Edge) / glm::dot(Edge, Edge), 0.0f, 1.0f);
        glm::vec3 ClosestOnEdge = p1 + t * Edge;
        return ClosestOnEdge;
    };
    
    // Check the three edges and compare distances to find the actual closest point if the projected point is outside
    glm::vec3 ClosestA = check_edge(a, b);
    glm::vec3 ClosestB = check_edge(b, c);
    glm::vec3 ClosestC = check_edge(c, a);

    float DistSqA = glm::distance(point, ClosestA);
    float DistSqB = glm::distance(point, ClosestB);
    float DistSqC = glm::distance(point, ClosestC);

    glm::vec3 v0 = b - a, v1 = c - a, v2 = point - a;
    float d00 = glm::dot(v0, v0);
    float d01 = glm::dot(v0, v1);
    float d11 = glm::dot(v1, v1);
    float d20 = glm::dot(v2, v0);
    float d21 = glm::dot(v2, v1);
    float denom = d00 * d11 - d01 * d01;

    // Projected point barycentric coords
    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;

    if (v >= 0 && w >= 0 && u >= 0) {
        // Point is inside the triangle's infinite prism, projected point is closest
        return ProjectedPoint;
    }

    glm::vec3 ClosestOnAB = check_edge(a, b);
    glm::vec3 ClosestOnBC = check_edge(b, c);
    glm::vec3 ClosestOnCA = check_edge(c, a);

    float DistSqAB = glm::distance(point, ClosestOnAB);
    float DistSqBC = glm::distance(point, ClosestOnBC);
    float DistSqCA = glm::distance(point, ClosestOnCA);

    if (DistSqAB <= DistSqBC && DistSqAB <= DistSqCA) {
        return ClosestOnAB;
    } else if (DistSqBC <= DistSqAB && DistSqBC <= DistSqCA) {
        return ClosestOnBC;
    } else {
        return ClosestOnCA;
    }
}

int Mesh::GetVoxelIndex(int X, int Y, int Z) const
{
    return Z * (SizeX + 1) * (SizeY + 1) + Y * (SizeX + 1) + X;
}

void Mesh::GenerateMesh()
{
    //for now, just testing it out
	if (true)
	{
		TriangleOrder[0] = 0;
		TriangleOrder[1] = 1;
		TriangleOrder[2] = 2;
	}
	else
	{
		TriangleOrder[0] = 2;
		TriangleOrder[1] = 1;
		TriangleOrder[2] = 0;
	}

	float Cube[8];
	for (int X = 0; X < SizeX; ++X)
	{
		for (int Y = 0; Y < SizeY; ++Y)
		{
			for (int Z = 0; Z < SizeZ; ++Z)
			{
                std::cout << "Marching cube at: " << X << ", " << Y << ", " << Z << std::endl;
				for (int i = 0; i < 8; ++i)
				{
					Cube[i] = Voxels[GetVoxelIndex(X + VertexOffset[i][0], Y + VertexOffset[i][1], Z + VertexOffset[i][2])];
				}
				March(X,Y,Z,Cube);
			}
		}
	}
}

void Mesh::March(int X, int Y, int Z, const float Cube[8])
{
	//UStaticMesh* StaticMesh = StaticMeshComponent ? StaticMeshComponent->GetStaticMesh() : nullptr;
	//if(!StaticMesh)
	//{
	//	return;
	//}
	//FBox boundingBox = StaticMesh->GetBoundingBox();
	//FVector meshMin = boundingBox.Min;
	//FVector meshMax = boundingBox.Max;
	const float meshWidth = maxX - minX;
	const float meshHeight = maxY - minY;

	int VertexMask = 0;
	glm::vec3 EdgeVertex[12];
	for (int i = 0; i < 8; ++i)
	{
		if (Cube[i] <= voxelSize * 0.5f)
				VertexMask |= (1 << i);
	}
	const int EdgeMask = CubeEdgeFlags[VertexMask];

	if (EdgeMask == 0 ) return;
	

	for (int i = 0; i < 12; ++i)
	{
		if ((EdgeMask & (1 << i)) != 0)
		{
			const float Offset = 0.5f;

			EdgeVertex[i].x = X + (VertexOffset[EdgeConnection[i][0]][0] + Offset * EdgeDirection[i][0]);
			EdgeVertex[i].y = Y + (VertexOffset[EdgeConnection[i][0]][1] + Offset * EdgeDirection[i][1]);
			EdgeVertex[i].z = Z + (VertexOffset[EdgeConnection[i][0]][2] + Offset * EdgeDirection[i][2]);
		}
	}

	for (int i = 0; i < 5; ++i)
	{
		if (TriangleConnectionTable[VertexMask][3*i] < 0) break;

        glm::vec3 V1 = glm::vec3(
            minX + (EdgeVertex[TriangleConnectionTable[VertexMask][3*i]].x / SizeX) * vertSizeX,
            minY + (EdgeVertex[TriangleConnectionTable[VertexMask][3*i]].y / SizeY) * vertSizeY, 
            minZ + (EdgeVertex[TriangleConnectionTable[VertexMask][3*i]].z / SizeZ) * vertSizeZ
        );

        glm::vec3 V2 = glm::vec3(
            minX + (EdgeVertex[TriangleConnectionTable[VertexMask][3*i + 1]].x / SizeX) * vertSizeX,
            minY + (EdgeVertex[TriangleConnectionTable[VertexMask][3*i + 1]].y / SizeY) * vertSizeY, 
            minZ + (EdgeVertex[TriangleConnectionTable[VertexMask][3*i + 1]].z / SizeZ) * vertSizeZ
        );

        glm::vec3 V3 = glm::vec3(
            minX + (EdgeVertex[TriangleConnectionTable[VertexMask][3*i + 2]].x / SizeX) * vertSizeX,
            minY + (EdgeVertex[TriangleConnectionTable[VertexMask][3*i + 2]].y / SizeY) * vertSizeY, 
            minZ + (EdgeVertex[TriangleConnectionTable[VertexMask][3*i + 2]].z / SizeZ) * vertSizeZ
        );

		glm::vec3 Normal = glm::cross(V2 - V1, V3 - V1);

		//FColor Color = FColor::MakeRandomColor();

		Normal = glm::normalize(Normal);

		//Vertices.Append({V1, V2, V3});
        Vertices.push_back(V1);
        Vertices.push_back(V2);
        Vertices.push_back(V3);

		glm::vec2 UV1 = glm::vec2(V1.x / meshWidth, (V1.y + V1.z) / (meshHeight + meshHeight));
		glm::vec2 UV2 = glm::vec2(V2.x / meshWidth, (V2.y + V2.z) / (meshHeight + meshHeight));
		glm::vec2 UV3 = glm::vec2(V3.x / meshWidth, (V3.y + V3.z) / (meshHeight + meshHeight));

        //UVs.Append({ UV1, UV2, UV3 });
        UVs.push_back(UV1);
        UVs.push_back(UV2);
        UVs.push_back(UV3);

		Triangles.push_back(vertexCount + TriangleOrder[0]);
		Triangles.push_back(vertexCount + TriangleOrder[1]);
		Triangles.push_back(vertexCount + TriangleOrder[2]);

		Normals.push_back(Normal);
		Normals.push_back(Normal);
		Normals.push_back(Normal);

		//Colors.Append({Color, Color, Color});
		vertexCount += 3;

        m_marchingCubesVertices.push_back(MeshData(V1, UV1, Normal));
        m_marchingCubesVertices.push_back(MeshData(V2, UV2, Normal));
        m_marchingCubesVertices.push_back(MeshData(V3, UV3, Normal));
	}
}


bool Mesh::isInsideMesh(const glm::vec3 p)
{
    std::vector<glm::vec3> directions = {
        glm::vec3(1.f, 0.f, 0.f),
        glm::vec3(0.f, 1.f, 0.f),
        glm::vec3(0.f, 0.f, 1.f),
        glm::vec3(-1.f, 0.f, 0.f),
        glm::vec3(0.f, -1.f, 0.f),
        glm::vec3(0.f, 0.f, -1.f),
        normalize(glm::vec3(1.f, 1.f, 1.f)),
        normalize(glm::vec3(-1.f, 1.f, 1.f)),
        normalize(glm::vec3(1.f, -1.f, 1.f)),
        normalize(glm::vec3(1.f, 1.f, -1.f))
    };

    int insideCount = 0;

    int totalRays = directions.size();
    const float epsilon = 0.0001f;

    for(const glm::vec3 dir : directions)
    {
        glm::vec3 RayStart = p + dir * epsilon;
        glm::vec3 RayEnd = p + dir * 10000.f;
        
        int hits = 0;
        bool hasHit = false;

        for(int i = 0; i < m_indices.size(); i += 3)
        {
            glm::vec3 A = m_vertices[m_indices[i]].position;
            glm::vec3 B = m_vertices[m_indices[i + 1]].position;
            glm::vec3 C = m_vertices[m_indices[i + 2]].position;

            if (glm::length(glm::cross(B - A, C - A)) < 0.0001f) {
                continue; // Skip degenerate triangles
            }

            glm::vec3 hitPoint;
            glm::vec3 normal;

            if(SegmentTriangleIntersection(RayStart, RayEnd, A,B,C, hitPoint, normal))
            {
                std::cout << "Ray hit triangle at: " << hitPoint.x << ", " << hitPoint.y << ", " << hitPoint.z << std::endl;
                hasHit = true;
                hits++;
            }
        }

        if(!hasHit)
        {
            glm::vec3 oppositeRayEnd = p - dir * 10000.f;
            for(int i = 0; i < m_indices.size(); i += 3)
            {
                glm::vec3 A = m_vertices[m_indices[i]].position;
                glm::vec3 B = m_vertices[m_indices[i + 1]].position;
                glm::vec3 C = m_vertices[m_indices[i + 2]].position;

                glm::vec3 hitPoint;
                glm::vec3 normal;

                if(SegmentTriangleIntersection(RayStart, oppositeRayEnd, A,B,C, hitPoint, normal))
                {
                    hits++;
                }
            }
        }

        if(hits % 2 == 1)
        {
            insideCount++;
        }
    }
    return insideCount > (totalRays / 2);
}


float Mesh::magnitude(glm::vec3 p)
{
    return sqrt(p.x*p.x+ p.y*p.y + p.z*p.z);
}

glm::vec3 Mesh::normalize(glm::vec3 p)
{
    float mag = magnitude(p);

    if(mag > 0)
    {
        return glm::vec3(p.x / mag, p.y / mag, p.z / mag);
    }
    else
    {
        return glm::vec3(0.0f, 0.0f, 0.0f);
    }

}


bool Mesh::SegmentTriangleIntersection(glm::vec3& start, glm::vec3& end, 
    glm::vec3& A, glm::vec3& B, glm::vec3& C, glm::vec3& outIntersectPoint, glm::vec3& outTriangleNormal)
{
    glm::vec3 edge1(B - A); 
    edge1 = normalize(edge1);
    glm::vec3 edge2(C - A); 
    edge2 = normalize(edge2);
    glm::vec3 triNormal = glm::cross(edge2, edge1);

    triNormal = normalize(triNormal);

    if(triNormal == glm::vec3(0,0,0))
    {
        return false;
    }

    bool bCollide = SegmentPlaneIntersection(start, end, triNormal, A, outIntersectPoint);
    if(!bCollide)
    {
        return false;
    }

    glm::vec3 baryCentricTriangle = computeBarycentricCoordinates(outIntersectPoint, A, B, C);
    if(baryCentricTriangle.x >= 0.0f &&baryCentricTriangle.y >= 0.0f &&baryCentricTriangle.z >= 0.0f )
    {
        outTriangleNormal = triNormal;
        return true;
    }
    return false;
    
}

bool Mesh::SegmentPlaneIntersection(
    const glm::vec3& StartPoint,
    const glm::vec3& EndPoint,
    const glm::vec3& PlaneNormal,
    const glm::vec3& PlanePoint,
    glm::vec3& OutIntersectionPoint)
{
    glm::vec3 SegmentDirection = EndPoint - StartPoint;
    glm::vec3 W = StartPoint - PlanePoint;

    float Denominator = glm::dot(PlaneNormal, SegmentDirection);
    float Numerator = -glm::dot(PlaneNormal, W);

    // Check if segment is parallel to the plane (denominator is close to zero)
    if (std::abs(Denominator) < 1e-6) { // Use a small epsilon for float comparison
        // If numerator is also zero, the segment lies within the plane.
        // For simplicity, this implementation returns false in this case
        // as it's often treated as a non-unique intersection in game math.
        return false;
    }

    // Compute the intersection parameter (t) for the infinite line
    float t = Numerator / Denominator;

    // Check if the intersection point lies within the segment (0 <= t <= 1)
    if (t >= 0.0f && t <= 1.0f) {
        OutIntersectionPoint = StartPoint + t * SegmentDirection;
        return true;
    }

    // Intersection occurs outside the segment bounds (on the infinite line)
    return false;
}
glm::vec3 Mesh::computeBarycentricCoordinates(const glm::vec3& p, const glm::vec3& a, const glm::vec3& b, const glm::vec3& c)
{
    glm::vec3 v0 = b - a;
    glm::vec3 v1 = c - a;
    glm::vec3 v2 = p - a;

    float d00 = glm::dot(v0, v0);
    float d01 = glm::dot(v0, v1);
    float d11 = glm::dot(v1, v1);
    float d20 = glm::dot(v2, v0);
    float d21 = glm::dot(v2, v1);

    float denom = d00 * d11 - d01 * d01;

    // Handle degenerate triangles (denom close to zero)
    if (glm::abs(denom) < 0.0001f) { 
        return glm::vec3(-1.0f); // Indicate an error or invalid state
    }

    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;

    return glm::vec3(u, v, w);


}
