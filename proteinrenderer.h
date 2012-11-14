#ifndef PROTEINRENDERER_H
#define PROTEINRENDERER_H

#include<vector>
#include<vertex.h>
#include<glslshader.h>

class ProteinRenderer
{
private :
    GLuint sphereslistID;
    bool initialized;
public:
    std::vector<Vertex>& vertList;

    ProteinRenderer(std::vector<Vertex>& vertexList);
    bool init();
    void render();
    void makeDisplayList();
};

#endif // PROTEINRENDERER_H
