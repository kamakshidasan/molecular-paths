#ifndef PROTEINRENDERER_H
#define PROTEINRENDERER_H

#include<vector>
#include<vertex.h>
#include<glslshader.h>

class ProteinRenderer
{
public:
    std::vector<Vertex>& vertList;

    ProteinRenderer(std::vector<Vertex>& vertexList);
    bool init();
    void render();
};

#endif // PROTEINRENDERER_H
