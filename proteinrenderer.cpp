#include "proteinrenderer.h"
#include "assert.h"
#include <iostream>

ProteinRenderer::ProteinRenderer(std::vector<Vertex>& vertexList)
    :vertList(vertexList)
{
}

using namespace shaders;

bool ProteinRenderer::init(){
    if(initSphereShader())return true;
    return false;
}

void ProteinRenderer::render(){
    glEnable(GL_COLOR_MATERIAL);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    GLSLShader sphereShader = (*getSphereShader());
    sphereShader.Use();
//    glColor4d(0.1, 0.6, 0.9, 1);
//    glColor4d(0.5, 0.1, 0.2, 1);
    glColor4d(0.35, 0, 0, 1);
    glBegin(GL_POINTS);
    glUniform1f(sphereShader("ug_add_radius"), (float) 0.0);
    glUniform1f(sphereShader("ug_alpha_value"), (float) 0.0);
    for(unsigned int i=0;i<vertList.size();i++){
        Vertex* v = &vertList[i];
        Vector3 pos = v->getCoordVector();
        glVertexAttrib1f(sphereShader["radius"], (float)v->Radius);
        glVertex3d(pos.X, pos.Y, pos.Z);
    }
    glEnd();
    sphereShader.UnUse();
    glDisable(GL_COLOR_MATERIAL);
}