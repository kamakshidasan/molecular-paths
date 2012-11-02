/***************************************************************************
 *   Copyright (C) 2010 by raghavendra,,,                                  *
 *   raghavendra@incognito                                                 *
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
#include "pocketviewer.h"
#include "ui_mainwindow.h"
#include <iostream>

PocketViewer::PocketViewer(Processor* pr,QWidget* parent)
    : QGLViewer(parent)
{
    m_processor = pr;
    rank = 0;
    persistence = 0;
    currentPocNum = 0;
    alphaShape = false;
    wireFrame = false;
    allPockets = false;
    mouths = false;
    onlyPockets = false;
    onlyVoids = false;
    alphaSkinSurface = false;
    pocketSkinSurface = false;
    alphaSkinWireFrame = false;
    pocketSkinWireFrame = false;
    smoothShading = true;
    allOrIndivPocs = true;
    volumeEnabled = false;
    pocketWireframe = false;

    powerDiag = true;
    complementSpacePD = true;
    onlyInsideVerts = true;
    pruneIsolatedVerts = false;
    intersectEdges = false;

    cHull = false;
    cHullWF = false;
    cHullNorm = false;
}

void PocketViewer::init()
{

    glClearColor(0.0f,0.0f,0.0f,0.0f);
    glShadeModel(GL_SMOOTH);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

    glEnable(GL_LIGHTING);

    glLightModelfv(GL_LIGHT_MODEL_AMBIENT,LightMaterial::BrightLightModelAmb);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    // Lights
    /*glLightfv(GL_LIGHT0, GL_POSITION, LightMaterial::BrightTopLightPos);
    glLightfv(GL_LIGHT1, GL_POSITION, LightMaterial::BrightBottomLightPos);

    glLightfv(GL_LIGHT0, GL_DIFFUSE, LightMaterial::BrightTopLightDif);
    glLightfv(GL_LIGHT0, GL_AMBIENT, LightMaterial::BrightTopLightAmb);
    glLightfv(GL_LIGHT0, GL_SPECULAR, LightMaterial::BrightTopLightSpec);
    glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, LightMaterial::BrightTopLightDir);

    glLightfv(GL_LIGHT1, GL_DIFFUSE, LightMaterial::BrightBottomLightDif);
    glLightfv(GL_LIGHT1, GL_AMBIENT, LightMaterial::BrightBottomLightAmb);
    glLightfv(GL_LIGHT1, GL_SPECULAR, LightMaterial::BrightBottomLightSpec);
    glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, LightMaterial::BrightBottomLightDir);

    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);*/

    glEnable(GL_DEPTH_TEST);
    glDisable(GL_COLOR_MATERIAL);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(1.0);

    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

//    powerDiagram->makeDisplayList(m_processor->alcx->delcx, m_processor->vertexList, complementSpacePD);

    help();
}

QString PocketViewer::helpString() const
{
    QString text("<h2>P o c k e t V i e w e r</h2>");
      text += "A set of interfaces in the form of <b>checkboxes/radio buttons</b> are provided to show/hide <b>AlphaShapes</b> , <b>Pockets</b> and <b>SkinSurfaces</b>";
      text += "in both <b>solid</b> and <b>wireframe</b> modes.<br><br>";
      text += "<b>AlphaSlider</b> can be used to vary alpha rank and the rank can also be given as input as a text in";
      text += "the textbox adjacent to AlphaSlider.<br><br>";
      text += "Similarly persistence can also be specified using <b>PersistenceSlider</b> and the textbox adjacent.<br><br>";
      text += "To display/hide individual <b>pockets</b> and <b>mouths</b> one can enable/disable checkboxes provided for each";
      text += "pocket/mouth in the <b>Pocket Information/Mouth Information table</b> displayed in the right hand side of the";
      text += "viewer.<br><br>";
      text += "Use the mouse to move the camera around the object. ";
      text += "You can respectively revolve around, zoom and translate with the three mouse buttons. ";
      text += "Left and middle buttons pressed together rotate around the camera view direction axis<br><br>";
      text += "Pressing <b>Alt</b> and one of the function keys (<b>F1</b>..<b>F12</b>) defines a camera keyFrame. ";
      text += "Simply press the function key again to restore it. Several keyFrames define a ";
      text += "camera path. Paths are saved when you quit the application and restored at next start.<br><br>";
      text += "Press <b>F</b> to display the frame rate, <b>A</b> for the world axis, ";
      text += "<b>Alt+Return</b> for full screen mode and <b>Control+S</b> to save a snapshot. ";
      text += "See the <b>Keyboard</b> tab in this window for a complete shortcut list.<br><br>";
      text += "Double clicks automates single click actions: A left button double click aligns the closer axis with the camera (if close enough). ";
      text += "A middle button double click fits the zoom of the camera and the right button re-centers the scene.<br><br>";
      text += "A left button double click while holding right button pressed defines the camera <i>Revolve Around Point</i>. ";
      text += "See the <b>Mouse</b> tab and the documentation web pages for details.<br><br>";
      text += "Press <b>Escape</b> to exit the viewer.";
      return text;
}

bool remakePowerDiagramDL = true;

void PocketViewer::draw()
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glPushMatrix();

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, LightMaterial::MatAmb2);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, LightMaterial::MatDiff2);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, LightMaterial::MatSpec2);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, LightMaterial::MatShin2[1]);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, LightMaterial::MatEmission);

    //lights
    glLightfv(GL_LIGHT0, GL_POSITION, LightMaterial::BrightTopLightPos);
    glLightfv(GL_LIGHT1, GL_POSITION, LightMaterial::BrightBottomLightPos);

    glLightfv(GL_LIGHT0, GL_DIFFUSE, LightMaterial::BrightTopLightDif);
    glLightfv(GL_LIGHT0, GL_AMBIENT, LightMaterial::BrightTopLightAmb);
    glLightfv(GL_LIGHT0, GL_SPECULAR, LightMaterial::BrightTopLightSpec);
    glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, LightMaterial::BrightTopLightDir);

    glLightfv(GL_LIGHT1, GL_DIFFUSE, LightMaterial::BrightBottomLightDif);
    glLightfv(GL_LIGHT1, GL_AMBIENT, LightMaterial::BrightBottomLightAmb);
    glLightfv(GL_LIGHT1, GL_SPECULAR, LightMaterial::BrightBottomLightSpec);
    glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, LightMaterial::BrightBottomLightDir);

    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);

    //drawLight (GL_LIGHT0);
    //drawLight (GL_LIGHT1);

    if(m_processor->IsRenderable())
    {
        if(remakePowerDiagramDL){
            m_processor->powerDiagram->makeDisplayList(complementSpacePD, onlyInsideVerts, pruneIsolatedVerts, intersectEdges);
            remakePowerDiagramDL = false;
        }
        m_processor->Render(persistence,alphaShape,allPockets,onlyPockets,onlyVoids,pocketSkinSurface,mouths,allOrIndivPocs,currentPocNum,rank,wireFrame,alphaSkinSurface,
                            alphaSkinWireFrame,smoothShading,pocketSkinWireFrame,pocketWireframe, powerDiag, cHull, cHullWF, cHullNorm);
    }
    glPopMatrix();
}

void PocketViewer::setRank(int rank)
{
    this->rank = rank;
    remakePowerDiagramDL = true;
    m_processor->powerDiagram->constructGraph(true);
    updateGL();
}

void PocketViewer::setPersistence(int persistence)
{
    this->persistence = persistence;
    updateGL();
}

void PocketViewer::setAlphaShape()
{
    alphaShape = !alphaShape;
    updateGL();
}

void PocketViewer::setWireFrame()
{
    wireFrame = !wireFrame;
    updateGL();
}

void PocketViewer::setAllPockets()
{
    allPockets = !allPockets;
    updateGL();
}

void PocketViewer::setMouths()
{
    mouths = !mouths;
    updateGL();
}

void PocketViewer::setOnlyPockets()
{
    onlyPockets = !onlyPockets;
    updateGL();
}

void PocketViewer::setOnlyVoids()
{
    onlyVoids = !onlyVoids;
    updateGL();
}

void PocketViewer::setPowerDiag()
{
    powerDiag = !powerDiag;
    updateGL();
}

void PocketViewer::setComplementSpacePD()
{
    complementSpacePD = !complementSpacePD;
    remakePowerDiagramDL = true;
    updateGL();
}

void PocketViewer::setOnlyInsideVerts()
{
    onlyInsideVerts = !onlyInsideVerts;
    remakePowerDiagramDL = true;
    updateGL();
}

void PocketViewer::setPruneIsolatedVerts()
{
    pruneIsolatedVerts = !pruneIsolatedVerts;
    remakePowerDiagramDL = true;
    updateGL();
}

void PocketViewer::setIntersectEdges()
{
    intersectEdges = !intersectEdges;
    remakePowerDiagramDL = true;
    updateGL();
}

void PocketViewer::setCHull()
{
    cHull = !cHull;
    updateGL();
}

void PocketViewer::setCHullWF()
{
    cHullWF = !cHullWF;
    updateGL();
}

void PocketViewer::setCHullNorm()
{
    cHullNorm = !cHullNorm;
    updateGL();
}


void PocketViewer::setPocketSkinSurface()
{
    pocketSkinSurface = !pocketSkinSurface;
    updateGL();
}

void PocketViewer::setPocketSkinWireFrame()
{
    pocketSkinWireFrame = !pocketSkinWireFrame;
    updateGL();
}

void PocketViewer::setAlphaSkinSurface()
{
    alphaSkinSurface = !alphaSkinSurface;
    updateGL();
}

void PocketViewer::setAlphaSkinWireFrame()
{
    alphaSkinWireFrame = !alphaSkinWireFrame;
    updateGL();
}

void PocketViewer::setSmoothShading(bool flag)
{
    smoothShading = flag;
    updateGL();
}

void PocketViewer::setAllOrIndivPocs(bool flag)
{
    allOrIndivPocs = flag;
    updateGL();
}

void PocketViewer::setCurrentPocNum(int num)
{
    currentPocNum = num;
    updateGL();
}

void PocketViewer::setPocketWireFrame()
{
    pocketWireframe = !pocketWireframe;
    updateGL();
}
