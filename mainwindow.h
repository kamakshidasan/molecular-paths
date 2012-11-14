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
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <pocketviewer.h>
#include <processor.h>
#include <iostream>
#include <sys/time.h>

#include <QMainWindow>
#include <QFileDialog>
#include <QMessageBox>

namespace Ui
{
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    public:
        MainWindow(Processor *pr,QWidget *parent = 0);
        ~MainWindow();

    public slots:
        void about();
        void open(bool constantRadius, bool incrementRadius);
        void openAsSpheres();
        void openAsPoints();
        void openWithIncrementedRadius();
        void saveGraph();
        void savePath();
        void onAlphaValueZero();
        void onChangeFiltration();
        void onUndoChange();
        void onRankChange();
        void onEpsilonTextChanged();
        void onRankTextChanged();
        void onPersistenceChange();
        void onPersistenceTextChanged();
        void onCheckBoxAlphaShapeToggled();
        void onCheckBoxWireFrameToggled();
        void onCheckBoxAllPocketsToggled();
        void onCheckBoxMouthsToggled();
        void onCheckBoxOnlyPocketsToggled();
        void onCheckBoxOnlyVoidsToggled();

        void onCheckPowerDiagToggled();
        void onCheckComplementPDToggled();
        void onCheckInsideVerts();
        void onCheckPruneIsolatedVerts();
        void onCheckIntersectEdges();

        void onCheckCHull();
        void onCheckCHullWF();
        void onCheckCHullNorm();

        void onSTPathClick();
        void onEscapePathClickRepeated();
        void onEscapePathClickAll();
        void onEscapePathClick(bool repeated, int maxIter);
        void onShortestEscapePathClick();
        void onCheckShowPathToggled();
        void onCheckShowPathSkinToggled();
        void onCheckShowPathSkinWFToggled();
        void onCheckShowPathSpheresToggled();

        void onCheckSpaceFillToggled();
        void onCheckSpaceFillPDToggled();

        void onCheckBoxPocketSkinSurfaceToggled();
        void onCheckBoxPocketSkinWireFrameToggled();
        void onCheckBoxAlphaSkinSurfaceToggled();
        void onCheckBoxAlphaSkinWireFrameToggled();
        void onRadiobuttonFlatClicked();
        void onRadioButtonSmoothClicked();
        void onRadioButtonAllPocketsClicked();
        void onRadioButtonIndivPocketsClicked();
        void onSpinBoxPocketsValueChanged();
        void onCheckBoxVolumeToggled();
        void onCheckBoxPocketWireFrameToggled();

    private:
        Ui::MainWindow *ui;
        PocketViewer * m_viewer1;
        PocketViewer * m_viewer2;
        Processor *m_processor;
};

#endif // MAINWINDOW_H
