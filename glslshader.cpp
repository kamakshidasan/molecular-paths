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
//A simple class for handling GLSL shader compilation
//Author: Movania Muhammad Mobeen
//Last Modified: February 2, 2011

#include "glslshader.h"
#include <iostream>

GLSLShader::GLSLShader(void)
{
        _shaders[VERTEX_SHADER]=0;
        _shaders[FRAGMENT_SHADER]=0;
        _shaders[GEOMETRY_SHADER]=0;
        _attributeList.clear();
        _uniformLocationList.clear();
}

GLSLShader::~GLSLShader(void)
{
        _attributeList.clear();
        _uniformLocationList.clear();
}

void GLSLShader::LoadFromString(GLenum type, const string& source) {
        GLuint shader = glCreateShader (type);

        const char * ptmp = source.c_str();
        glShaderSource (shader, 1, &ptmp, NULL);

        //check whether the shader loads fine
        GLint status;
        glCompileShader (shader);
        glGetShaderiv (shader, GL_COMPILE_STATUS, &status);
        if (status == GL_FALSE) {
                GLint infoLogLength;
                glGetShaderiv (shader, GL_INFO_LOG_LENGTH, &infoLogLength);
                GLchar *infoLog= new GLchar[infoLogLength];
                glGetShaderInfoLog (shader, infoLogLength, NULL, infoLog);
                cerr<<"Compile log: "<<infoLog<<endl;
                delete [] infoLog;
        }
        switch(type){
        case GL_VERTEX_SHADER:
            _shaders[VERTEX_SHADER]=shader;
            break;
        case GL_FRAGMENT_SHADER:
            _shaders[FRAGMENT_SHADER]=shader;
            break;
        case GL_GEOMETRY_SHADER:
            _shaders[GEOMETRY_SHADER]=shader;
            break;
        }
}


void GLSLShader::CreateAndLinkProgram(GLuint geomIn, GLuint geomOut) {
        _program = glCreateProgram ();
        if (_shaders[VERTEX_SHADER] != 0) {
                glAttachShader (_program, _shaders[VERTEX_SHADER]);
        }
        if (_shaders[FRAGMENT_SHADER] != 0) {
                glAttachShader (_program, _shaders[FRAGMENT_SHADER]);
        }
        if (_shaders[GEOMETRY_SHADER] != 0) {
                glAttachShader (_program, _shaders[GEOMETRY_SHADER]);
        }

        glProgramParameteriEXT ( _program, GL_GEOMETRY_INPUT_TYPE_EXT, geomIn );
        glProgramParameteriEXT ( _program, GL_GEOMETRY_OUTPUT_TYPE_EXT, geomOut );

        int temp;
        glGetIntegerv ( GL_MAX_GEOMETRY_OUTPUT_VERTICES_EXT, &temp );
        glProgramParameteriEXT ( _program, GL_GEOMETRY_VERTICES_OUT_EXT, temp );

        //link and check whether the program links fine
        GLint status;
        glLinkProgram (_program);
        glGetProgramiv (_program, GL_LINK_STATUS, &status);
        if (status == GL_FALSE) {
                GLint infoLogLength;

                glGetProgramiv (_program, GL_INFO_LOG_LENGTH, &infoLogLength);
                GLchar *infoLog= new GLchar[infoLogLength];
                glGetProgramInfoLog (_program, infoLogLength, NULL, infoLog);
                cerr<<"Link log: "<<infoLog<<endl;
                delete [] infoLog;
        }

        glDeleteShader(_shaders[VERTEX_SHADER]);
        glDeleteShader(_shaders[FRAGMENT_SHADER]);
        glDeleteShader(_shaders[GEOMETRY_SHADER]);
}

void GLSLShader::Use() {
        glUseProgram(_program);
}

void GLSLShader::UnUse() {
        glUseProgram(0);
}

void GLSLShader::AddAttribute(const string& attribute) {
        _attributeList[attribute]= glGetAttribLocation(_program, attribute.c_str());
}

//An indexer that returns the location of the attribute
GLuint GLSLShader::operator [](const string& attribute) {
        return _attributeList[attribute];
}

void GLSLShader::AddUniform(const string& uniform) {
        _uniformLocationList[uniform] = glGetUniformLocation(_program, uniform.c_str());
}

GLuint GLSLShader::operator()(const string& uniform){
        return _uniformLocationList[uniform];
}
GLuint GLSLShader::GetProgram() const {
        return _program;
}
#include <fstream>
void GLSLShader::LoadFromFile(GLenum whichShader, const string& filename){
        ifstream fp;
        fp.open(filename.c_str(), ios_base::in);
        if(fp) {
                /*string line, buffer;
                while(getline(fp, line)) {
                        buffer.append(line);
                        buffer.append("\r\n");
                }               */
                string buffer(std::istreambuf_iterator<char>(fp), (std::istreambuf_iterator<char>()));
                //copy to source
                LoadFromString(whichShader, buffer);
        } else {
                cerr<<"Error loading shader: "<<filename<<endl;
        }
}

void GLSLShader::DeleteProgram(){
    glDeleteProgram(_program);
    _program=-1;
}

namespace shaders{
    static GLSLShader sphereShader;
    static GLSLShader cylinShader;

    static bool sphereShaderInitialized = false;
    static bool cylinShaderInitialized = false;

    bool initSphereShader(){
        if(sphereShaderInitialized)return true;
        GLenum err = glewInit();
        if (err != GLEW_OK){
            std::cerr << glewGetString(err) << std::endl;
            return false;
        }
        if (!GLEW_VERSION_2_1){
            std::cerr << "No support for OpenGL 2.1" << std::endl;
            return false;
        }
        sphereShader.LoadFromFile(GL_VERTEX_SHADER, "./shaders/sphere.vert");
        sphereShader.LoadFromFile(GL_GEOMETRY_SHADER, "./shaders/sphere.geom");
        sphereShader.LoadFromFile(GL_FRAGMENT_SHADER, "./shaders/sphere.frag");
        sphereShader.CreateAndLinkProgram(GL_POINTS, GL_QUADS);
        sphereShader.Use();
        sphereShader.AddAttribute("radius");
        sphereShader.AddUniform("ug_add_radius");
        sphereShader.AddUniform("ug_alpha_value");
        sphereShader.UnUse();
        sphereShaderInitialized = true;
        return true;
    }

    bool initCylinderShader(){
        if(cylinShaderInitialized)return true;
        GLenum err = glewInit();
        if (err != GLEW_OK){
            std::cerr << glewGetString(err) << std::endl;
            return false;
        }
        if (!GLEW_VERSION_2_1){
            std::cerr << "No support for OpenGL 2.1" << std::endl;
            return false;
        }
        cylinShader.LoadFromFile(GL_VERTEX_SHADER, "./shaders/cylinder_vert.glsl");
        cylinShader.LoadFromFile(GL_GEOMETRY_SHADER, "./shaders/cylinder_geom.glsl");
        cylinShader.LoadFromFile(GL_FRAGMENT_SHADER, "./shaders/cylinder_frag.glsl");
        cylinShader.CreateAndLinkProgram(GL_LINES, GL_QUADS);
        cylinShader.Use();
        cylinShader.AddAttribute("radius");
        cylinShader.UnUse();
        cylinShaderInitialized = true;
        return true;
    }

    GLSLShader* getSphereShader(){
        return &sphereShader;
    }

    GLSLShader* getCylinderShader(){
        return &cylinShader;
    }
}
