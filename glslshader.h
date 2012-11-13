/***************************************************************************
 *   Copyright (C) 2010 by talha bin masood                                *                                                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef GLSLSHADER_H
#define GLSLSHADER_H
//A simple class for handling GLSL shader compilation
//Auhtor: Movania Muhammad Mobeen
#include <GL/glew.h>
#include <map>
#include <string>

using namespace std;

class GLSLShader
{
public:
        GLSLShader(void);
        ~GLSLShader(void);
        void LoadFromString(GLenum whichShader, const string& source);
        void LoadFromFile(GLenum whichShader, const string& filename);
        void CreateAndLinkProgram(GLuint geomIn, GLuint geomOut);
        void Use();
        void UnUse();
        void AddAttribute(const string& attribute);
        void AddUniform(const string& uniform);
        GLuint GetProgram() const;
        //An indexer that returns the location of the attribute/uniform
        GLuint operator[](const string& attribute);
        GLuint operator()(const string& uniform);
        //Program deletion
        void DeleteProgram();

        GLuint  _program;

private:
        enum ShaderType {VERTEX_SHADER, FRAGMENT_SHADER, GEOMETRY_SHADER};
        GLuint _shaders[3];//0->vertexshader, 1->fragmentshader, 2->geometryshader
        map<string,GLuint> _attributeList;
        map<string,GLuint> _uniformLocationList;
};
namespace shaders{
    bool initSphereShader();
    GLSLShader* getSphereShader();
    bool initCylinderShader();
    GLSLShader* getCylinderShader();
}
#endif // GLSLSHADER_H
