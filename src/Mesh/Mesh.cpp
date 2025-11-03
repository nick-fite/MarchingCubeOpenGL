#include "Mesh.h"
#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>
#include <assimp/scene.h>
#include <assimp/DefaultLogger.hpp>
#include <assimp/LogStream.hpp>

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
        scene = importer.ReadFile(filePath, aiProcessPreset_TargetRealtime_Quality | aiProcess_PreTransformVertices);
    }

    //this what indicates how many meshes are being imported. If I take in multiple in the same fbx will this process it?
    const struct aiMesh* mesh = scene->mMeshes[meshNum];

    minY = std::numeric_limits<float>::max();
    maxY = std::numeric_limits<float>::lowest();
    minX = std::numeric_limits<float>::max();
    maxX = std::numeric_limits<float>::lowest();
    minZ = std::numeric_limits<float>::max();
    maxZ = std::numeric_limits<float>::lowest();

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

        minY = std::min(minY, v.position.y);
        maxY = std::max(maxY, v.position.y);
        minX = std::min(minX, v.position.x);
        maxX = std::max(maxX, v.position.x);
        minZ = std::min(minZ, v.position.z);
        maxZ = std::max(maxZ, v.position.z);

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

    vertSizeX = maxX - minX;
    vertSizeY = maxY - minY;
    vertSizeZ = maxZ - minZ;

    //for(t = 0; t < mesh->mNumVertices; ++t)
    //{
    //    m_vertices[t].position.y = (m_vertices[t].position.y - minY) / height;
    //}

    float voxelSizeX = std::ceil(vertSizeX / voxelSize) * 4;
    float voxelSizeY = std::ceil(vertSizeY / voxelSize) * 2;
    float voxelSizeZ = std::ceil(vertSizeZ / voxelSize) * 4;

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
                std::cout << "Generating voxel: " << x << ", " << y << ", " << z << std::endl;
                glm::vec3 pos = glm::vec3(x, y, z);
                float posX = minX + ((float)x / voxelSizeX) * vertSizeX;
                float posY = minY + ((float)y / voxelSizeY) * vertSizeY;
                float posZ = minZ + ((float)z / voxelSizeZ) * vertSizeZ;

                float dist = closestTriangleDistance(glm::vec3(posX, posY, posZ));

                Voxels.push_back(dist);

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
    MakeMesh(m_marchingCubesVertices, marchingCubesIndices);

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

    //turn on attribute
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
