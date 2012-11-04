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
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "qcustomplot.h"

MainWindow::MainWindow(Processor *pr,QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    m_viewer1 = new PocketViewer(pr,ui->FirstViewerPlaceHolder);

    m_viewer1->resize(ui->FirstViewerPlaceHolder->size());

//    m_viewer2 = new PocketViewer(pr,ui->SecondViewerPlaceHolder);

//    m_viewer2->resize(ui->SecondViewerPlaceHolder->size());

    m_processor = pr;

    QButtonGroup *bg1 = new QButtonGroup(this);
    bg1->addButton(ui->radioButtonAllPockets);
    bg1->addButton(ui->radioButtonIndivPockets);

    QButtonGroup *bg2 = new QButtonGroup(this);
    bg2->addButton(ui->radioButtonFlat);
    bg2->addButton(ui->radioButtonSmooth);

    m_processor->table1 = ui->tableWidgetPocket;
    m_processor->table2 = ui->tableWidgetMouth;
    m_processor->totVol = ui->lineEditMolTotalVolume;
    m_processor->totSurf = ui->lineEditMolTotSurfArea;
    m_processor->pocVol = ui->lineEditPocketTotalVolume;
    m_processor->pocSurf = ui->lineEditPocTotSurfArea;
    m_processor->mouSurf = ui->lineEditMouthTotalSurfArea;
    m_processor->pspin = ui->spinBoxPockets;

    connect(ui->actionExit,SIGNAL(triggered()),this,SLOT(close()));
    connect(ui->actionAbout,SIGNAL(triggered()),this,SLOT(about()));
    connect(ui->actionOpen,SIGNAL(triggered()),this,SLOT(open()));
    connect(ui->actionAlphaZero,SIGNAL(triggered()),this,SLOT(onAlphaValueZero()));
    connect(ui->actionChangeFiltration,SIGNAL(triggered()),this,SLOT(onChangeFiltration()));
    connect(ui->actionUndoChange,SIGNAL(triggered()),this,SLOT(onUndoChange()));
    connect(ui->actionSave_Graph,SIGNAL(triggered()),this,SLOT(saveGraph()));

    connect(ui->lineEditEpsilonInput,SIGNAL(returnPressed()),this,SLOT(onEpsilonTextChanged()));

    connect(ui->alphaSlider,SIGNAL(sliderReleased()),this,SLOT(onRankChange()));
    connect(ui->lineEditAlphaRank,SIGNAL(returnPressed()),this,SLOT(onRankTextChanged()));
    connect(ui->persistenceSlider,SIGNAL(sliderReleased()),this,SLOT(onPersistenceChange()));
    connect(ui->lineEditPersistenceRank,SIGNAL(returnPressed()),this,SLOT(onPersistenceTextChanged()));
    connect(ui->checkBoxAlphaShape,SIGNAL(toggled(bool)),this,SLOT(onCheckBoxAlphaShapeToggled()));
    connect(ui->checkBoxWireFrame,SIGNAL(toggled(bool)),this,SLOT(onCheckBoxWireFrameToggled()));
    connect(ui->checkBoxAllPockets,SIGNAL(toggled(bool)),this,SLOT(onCheckBoxAllPocketsToggled()));
    connect(ui->checkBoxAllMouths,SIGNAL(toggled(bool)),this,SLOT(onCheckBoxMouthsToggled()));
    connect(ui->checkBoxOnlyPockets,SIGNAL(toggled(bool)),this,SLOT(onCheckBoxOnlyPocketsToggled()));
    connect(ui->checkBoxOnlyVoids,SIGNAL(toggled(bool)),this,SLOT(onCheckBoxOnlyVoidsToggled()));

    connect(ui->checkPowerDiag,SIGNAL(toggled(bool)),this,SLOT(onCheckPowerDiagToggled()));
    connect(ui->checkComplementSpacePD,SIGNAL(toggled(bool)),this,SLOT(onCheckComplementPDToggled()));
    connect(ui->checkInsideVerts,SIGNAL(toggled(bool)),this,SLOT(onCheckInsideVerts()));
    connect(ui->checkPruneIsolatedVerts,SIGNAL(toggled(bool)),this,SLOT(onCheckPruneIsolatedVerts()));
    connect(ui->checkEdgeIntersect,SIGNAL(toggled(bool)),this,SLOT(onCheckIntersectEdges()));


    connect(ui->checkCHull,SIGNAL(toggled(bool)),this,SLOT(onCheckCHull()));
    connect(ui->checkCHullWireFrame,SIGNAL(toggled(bool)),this,SLOT(onCheckCHullWF()));
    connect(ui->checkCHullNorm,SIGNAL(toggled(bool)),this,SLOT(onCheckCHullNorm()));

    connect(ui->buttonShortestSTPath,SIGNAL(clicked()),this,SLOT(onSTPathClick()));
    connect(ui->buttonShortestEscapePath,SIGNAL(clicked()),this,SLOT(onShortestEscapePathClick()));
    connect(ui->buttonShortestEscapePathAll,SIGNAL(clicked()),this,SLOT(onEscapePathClick()));

    connect(ui->checkBoxAlphaSkinSolid,SIGNAL(toggled(bool)),this,SLOT(onCheckBoxAlphaSkinSurfaceToggled()));
    connect(ui->checkBoxAlphaSkinWireFrame,SIGNAL(toggled(bool)),this,SLOT(onCheckBoxAlphaSkinWireFrameToggled()));
    connect(ui->checkBoxPocketSkinSolid,SIGNAL(toggled(bool)),this,SLOT(onCheckBoxPocketSkinSurfaceToggled()));
    connect(ui->checkBoxPocketSkinWireFrame,SIGNAL(toggled(bool)),this,SLOT(onCheckBoxPocketSkinWireFrameToggled()));

    connect(ui->radioButtonFlat,SIGNAL(clicked()),this,SLOT(onRadiobuttonFlatClicked()));
    connect(ui->radioButtonSmooth,SIGNAL(clicked()),this,SLOT(onRadioButtonSmoothClicked()));
    connect(ui->radioButtonAllPockets,SIGNAL(clicked()),this,SLOT(onRadioButtonAllPocketsClicked()));
    connect(ui->radioButtonIndivPockets,SIGNAL(clicked()),this,SLOT(onRadioButtonIndivPocketsClicked()));
    connect(ui->spinBoxPockets,SIGNAL(valueChanged(int)),this,SLOT(onSpinBoxPocketsValueChanged()));
    connect(ui->checkBoxVolumeEnabled,SIGNAL(toggled(bool)),this,SLOT(onCheckBoxVolumeToggled()));
    connect(ui->checkBoxPocketWireFrame,SIGNAL(toggled(bool)),this,SLOT(onCheckBoxPocketWireFrameToggled()));

    ui->labelMolTotalSurfaceArea->hide();
    ui->labelMolTotalVolume->hide();
    ui->labelMouthTotalSurfaceArea->hide();
    ui->labelPocketTotalSurfaceArea->hide();
    ui->labelPocketTotalVolume->hide();

    ui->lineEditMolTotalVolume->hide();
    ui->lineEditMolTotSurfArea->hide();
    ui->lineEditMouthTotalSurfArea->hide();
    ui->lineEditPocketTotalVolume->hide();
    ui->lineEditPocTotSurfArea->hide();

    ui->tableWidgetMouth->hide();
    ui->tableWidgetPocket->hide();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::about()
{
    QMessageBox::about(this,tr("About PocketViewer"),tr("<b>PocketViewer 0.1</b> displays voids and pockets " "in terms of both tetrahedra and " "skin surface." "For help see the help window"));
}

void MainWindow::open()
{
    double center[3],depth;
    QString FileName = QFileDialog::getOpenFileName(this,tr("Open File"),"/home/tbmasood/Desktop/raghavendra/ParsePDB");
    if(!FileName.isEmpty())
    {
        m_processor->read(FileName.toAscii().constData(),center,&depth);

        qglviewer::Vec v1(center[0],center[1],center[2]+1.5*depth);
        qglviewer::Vec v2(center[0],center[1],center[2]);
        m_viewer1->camera()->setSceneRadius(1.5*depth);
//        m_viewer2->camera()->setSceneRadius(1.5*depth);
        m_viewer1->camera()->setSceneCenter(v2);
//        m_viewer2->camera()->setSceneCenter(v2);

        m_viewer1->camera()->showEntireScene();
//        m_viewer2->camera()->showEntireScene();
        ui->alphaSlider->setMinimum(0);
        ui->alphaSlider->setMaximum(m_processor->getMaxRank());
        ui->lineEditAlphaRank->setValidator(new QIntValidator(0,m_processor->getMaxRank(),this));
        ui->persistenceSlider->setMinimum(0);
        ui->persistenceSlider->setMaximum(m_processor->getMaxPersistence()+1);
        ui->lineEditPersistenceRank->setValidator(new QIntValidator(0,m_processor->getMaxPersistence()+1,this));
    }
}

void MainWindow::saveGraph()
{
    QString FileName = QFileDialog::getSaveFileName(this,tr("Save Graph"),".");
    if(!FileName.isEmpty())
    {
        m_processor->powerDiagram->writeGraph(true, FileName.toAscii().constData());
    }
}

void MainWindow::onAlphaValueZero()
{
    int rank = m_processor->getRankForAlpha(0.0);
    float alpvalue;

    ui->alphaSlider->setValue(rank);

    QString r1,r2;

    QString &ran1 = r1.setNum(rank);
    ui->lineEditAlphaRank->setText(ran1);

    alpvalue = m_processor->getAlphaValue(rank);
    QString& ran2 = r2.setNum(alpvalue,'f',6);
    ui->lineEditAlphaValue->setText(ran2);

    m_processor->CalculateRelevant(rank);
    if(ui->checkBoxVolumeEnabled->isChecked())
    {
        m_processor->CalculateVolumes(rank);
    }

    m_viewer1->setRank(rank);
//    m_viewer2->setRank(rank);
}

void MainWindow::onRankChange()
{
    int rank = ui->alphaSlider->value();
    float alpvalue;

    QString r1,r2;

    QString& ran1 = r1.setNum(rank);
    ui->lineEditAlphaRank->setText(ran1);

    alpvalue = m_processor->getAlphaValue(rank);
    QString& ran2 = r2.setNum(alpvalue,'f',6);
    ui->lineEditAlphaValue->setText(ran2);

    //m_processor->CalculateEverythingFor(rank);
    m_processor->CalculateRelevant(rank);
    if(ui->checkBoxVolumeEnabled->isChecked())
    {
        m_processor->CalculateVolumes(rank);
    }

    m_viewer1->setRank(rank);
//    m_viewer2->setRank(rank);
}

void MainWindow::onRankTextChanged()
{
    QString r;
    int rank = atoi(ui->lineEditAlphaRank->text().toAscii().data());

    ui->alphaSlider->setValue(rank);
    float alpvalue = m_processor->getAlphaValue(rank);
    QString& ran = r.setNum(alpvalue,'f',5);
    ui->lineEditAlphaValue->setText(ran);

    //m_processor->CalculateEverythingFor(rank);
    m_processor->CalculateRelevant(rank);
    if(ui->checkBoxVolumeEnabled->isChecked())
    {
        m_processor->CalculateVolumes(rank);
    }

    m_viewer1->setRank(rank);
//    m_viewer2->setRank(rank);
}

void MainWindow::onPersistenceChange()
{
    int persistence;
    QString r1;

    persistence = ui->persistenceSlider->value();
    QString& ran1 = r1.setNum(persistence);
    ui->lineEditPersistenceRank->setText(ran1);

    m_viewer1->setPersistence(persistence);
//    m_viewer2->setPersistence(persistence);
}

void MainWindow::onPersistenceTextChanged()
{
    QString r;
    int persistence = atoi(ui->lineEditPersistenceRank->text().toAscii().data());

    ui->persistenceSlider->setValue(persistence);

    m_viewer1->setPersistence(persistence);
//    m_viewer2->setPersistence(persistence);
}

void MainWindow::onCheckBoxAlphaShapeToggled()
{
    m_viewer1->setAlphaShape();
//    m_viewer2->setAlphaShape();
}

void MainWindow::onCheckBoxWireFrameToggled()
{
    m_viewer1->setWireFrame();
//    m_viewer2->setWireFrame();
}

void MainWindow::onCheckBoxAllPocketsToggled()
{
    m_viewer1->setAllPockets();
//    m_viewer2->setAllPockets();
}

void MainWindow::onCheckBoxMouthsToggled()
{
    m_viewer1->setMouths();
//    m_viewer2->setMouths();
}

void MainWindow::onCheckBoxOnlyPocketsToggled()
{
    m_viewer1->setOnlyPockets();
//    m_viewer2->setOnlyPockets();
}

void MainWindow::onCheckBoxOnlyVoidsToggled()
{
    m_viewer1->setOnlyVoids();
//    m_viewer2->setOnlyVoids();
}

void MainWindow::onCheckPowerDiagToggled()
{
    m_viewer1->setPowerDiag();
//    m_viewer2->setOnlyVoids();
}

void MainWindow::onCheckComplementPDToggled()
{
    m_viewer1->setComplementSpacePD();
//    m_viewer2->setOnlyVoids();
}

void MainWindow::onCheckInsideVerts()
{
    m_viewer1->setOnlyInsideVerts();
//    m_viewer2->setOnlyVoids();
}

void MainWindow::onCheckPruneIsolatedVerts()
{
    m_viewer1->setPruneIsolatedVerts();
//    m_viewer2->setOnlyVoids();
}

void MainWindow::onCheckIntersectEdges()
{
    m_viewer1->setIntersectEdges();
//    m_viewer2->setOnlyVoids();
}


void MainWindow::onCheckCHull()
{
    m_viewer1->setCHull();
//    m_viewer2->setOnlyVoids();
}

void MainWindow::onCheckCHullWF()
{
    m_viewer1->setCHullWF();
//    m_viewer2->setOnlyVoids();
}

void MainWindow::onCheckCHullNorm()
{
    m_viewer1->setCHullNorm();
//    m_viewer2->setOnlyVoids();
}

void MainWindow::onSTPathClick(){
    int startIndex = ui->startCombo->currentText().toInt();
    int targetIndex = ui->targetCombo->currentText().toInt();
    QVector<double> X(100), Y(100);
    double length, minY, maxY;
    bool found = m_processor->powerDiagram->findShortestPath(startIndex, targetIndex, &X, &Y, &length, &minY, &maxY);
    if(!found){
        return;
    }

    int numGraphs = ui->graphWidget->graphCount();
    for(int i=0;i<numGraphs;i++){
        ui->graphWidget->removeGraph(0);
    }
    ui->graphWidget->addGraph();
    ui->graphWidget->graph(0)->setData(X, Y);
    // give the axes some labels:
    ui->graphWidget->xAxis->setLabel("Distance from Source");
    ui->graphWidget->yAxis->setLabel("Power Distance");
    // set axes ranges, so we see all data:
    ui->graphWidget->xAxis->setRange(0, length);
    ui->graphWidget->yAxis->setRange(minY, maxY);
    ui->graphWidget->replot();

    m_viewer1->updateGL();
}

void MainWindow::onShortestEscapePathClick(){
    int startIndex = ui->startCombo->currentText().toInt();
    QVector<double> X(100), Y(100);
    double length, minY, maxY;
    bool found = m_processor->powerDiagram->findShortestEscapePath(startIndex, &X, &Y, &length, &minY, &maxY);
    if(!found){
        return;
    }

    int numGraphs = ui->graphWidget->graphCount();
    for(int i=0;i<numGraphs;i++){
        ui->graphWidget->removeGraph(0);
    }
    ui->graphWidget->addGraph();
    ui->graphWidget->graph(0)->setData(X, Y);
    // give the axes some labels:
    ui->graphWidget->xAxis->setLabel("Distance from Source");
    ui->graphWidget->yAxis->setLabel("Power Distance");
    // set axes ranges, so we see all data:
    ui->graphWidget->xAxis->setRange(0, length);
    ui->graphWidget->yAxis->setRange(minY, maxY);
    ui->graphWidget->replot();

    m_viewer1->updateGL();
}

void MainWindow::onEscapePathClick(){
    int startIndex = ui->startCombo->currentText().toInt();
    std::vector<QVector<double> > Xs, Ys;
    std::vector<double> lengths, minYs, maxYs;
    int shortest = m_processor->powerDiagram->findShortestEscapePaths(startIndex, 100, &Xs, &Ys, &lengths, &minYs, &maxYs);
    if(shortest==-1){
        return;
    }
    int numGraphs = ui->graphWidget->graphCount();
    for(int i=0;i<numGraphs;i++){
        ui->graphWidget->removeGraph(0);
    }
    if(!Xs.empty()){
        double maxLength = lengths[0];
        double maxY = maxYs[0];
        double minY = minYs[0];
        for(int i=0;i<Xs.size();i++){
            ui->graphWidget->addGraph();
            ui->graphWidget->graph(i)->setData(Xs[i], Ys[i]);
            if(maxLength<lengths[i]){
                maxLength = lengths[i];
            }
            if(maxY<maxYs[i]){
                maxY = maxYs[i];
            }
            if(minY>minYs[i]){
                minY = minYs[i];
            }
            if(i==shortest){
                ui->graphWidget->graph(i)->setPen(QPen(Qt::red));
            }
        }
        // give the axes some labels:
        ui->graphWidget->xAxis->setLabel("Distance from Source");
        ui->graphWidget->yAxis->setLabel("Power Distance");
        // set axes ranges, so we see all data:
        ui->graphWidget->xAxis->setRange(0, maxLength);
        ui->graphWidget->yAxis->setRange(minY, maxY);
    }
    ui->graphWidget->replot();

    m_viewer1->updateGL();
}

void MainWindow::onCheckBoxPocketSkinSurfaceToggled()
{
    m_viewer1->setPocketSkinSurface();
//    m_viewer2->setPocketSkinSurface();
}

void MainWindow::onCheckBoxPocketSkinWireFrameToggled()
{
    m_viewer1->setPocketSkinWireFrame();
//    m_viewer2->setPocketSkinWireFrame();
}

void MainWindow::onCheckBoxAlphaSkinSurfaceToggled()
{
    m_viewer1->setAlphaSkinSurface();
//    m_viewer2->setAlphaSkinSurface();
}

void MainWindow::onCheckBoxAlphaSkinWireFrameToggled()
{
   m_viewer1->setAlphaSkinWireFrame();
//   m_viewer2->setAlphaSkinWireFrame();
}

void MainWindow::onRadiobuttonFlatClicked()
{
    m_viewer1->setSmoothShading(false);
//    m_viewer2->setSmoothShading(false);
}

void MainWindow::onRadioButtonSmoothClicked()
{
    m_viewer1->setSmoothShading(true);
//    m_viewer2->setSmoothShading(true);
}

void MainWindow::onRadioButtonAllPocketsClicked()
{
    if(ui->radioButtonAllPockets->isChecked())
    {
        if(ui->spinBoxPockets->isEnabled())
        {
            ui->spinBoxPockets->setEnabled(false);
        }
        m_viewer1->setAllOrIndivPocs(true);
//        m_viewer2->setAllOrIndivPocs(true);
    }
}

void MainWindow::onRadioButtonIndivPocketsClicked()
{
    if(ui->radioButtonIndivPockets->isChecked())
    {
        if(!ui->spinBoxPockets->isEnabled())
        {
            ui->spinBoxPockets->setEnabled(true);
        }
        m_viewer1->setAllOrIndivPocs(false);
//        m_viewer2->setAllOrIndivPocs(false);
    }
}

void MainWindow::onSpinBoxPocketsValueChanged()
{
    m_viewer1->setCurrentPocNum(ui->spinBoxPockets->value());
//    m_viewer2->setCurrentPocNum(ui->spinBoxPockets->value());
}

void MainWindow::onCheckBoxVolumeToggled()
{
    if(ui->checkBoxVolumeEnabled->isChecked())
    {
        int rank = ui->alphaSlider->value();

        ui->labelMolTotalSurfaceArea->show();
        ui->labelMolTotalVolume->show();
        ui->labelMouthTotalSurfaceArea->show();
        ui->labelPocketTotalSurfaceArea->show();
        ui->labelPocketTotalVolume->show();

        ui->lineEditMolTotalVolume->show();
        ui->lineEditMolTotSurfArea->show();
        ui->lineEditMouthTotalSurfArea->show();
        ui->lineEditPocketTotalVolume->show();
        ui->lineEditPocTotSurfArea->show();

        ui->tableWidgetMouth->show();
        ui->tableWidgetPocket->show();

        m_processor->CalculateVolumes(rank);
    }
    else
    {
        ui->labelMolTotalSurfaceArea->hide();
        ui->labelMolTotalVolume->hide();
        ui->labelMouthTotalSurfaceArea->hide();
        ui->labelPocketTotalSurfaceArea->hide();
        ui->labelPocketTotalVolume->hide();

        ui->lineEditMolTotalVolume->hide();
        ui->lineEditMolTotSurfArea->hide();
        ui->lineEditMouthTotalSurfArea->hide();
        ui->lineEditPocketTotalVolume->hide();
        ui->lineEditPocTotSurfArea->hide();

        ui->tableWidgetMouth->hide();
        ui->tableWidgetPocket->hide();
    }
}

void MainWindow::onCheckBoxPocketWireFrameToggled()
{
    m_viewer1->setPocketWireFrame();
//    m_viewer2->setPocketWireFrame();
}

void MainWindow::onChangeFiltration()
{
    int rank = ui->alphaSlider->value();

    m_processor->ModifyFiltration();

    m_processor->CalculateRelevant(rank);
    if(ui->checkBoxVolumeEnabled->isChecked())
    {
        m_processor->CalculateVolumes(rank);
    }

    m_viewer1->setRank(rank);
//    m_viewer2->setRank(rank);
}

void MainWindow::onUndoChange()
{
    int rank = ui->alphaSlider->value();

    m_processor->UndoModify();

    m_processor->CalculateRelevant(rank);
    if(ui->checkBoxVolumeEnabled->isChecked())
    {
        m_processor->CalculateVolumes(rank);
    }

    m_viewer1->setRank(rank);
//    m_viewer2->setRank(rank);
}

void MainWindow::onEpsilonTextChanged()
{
    double epsilon = atof(ui->lineEditEpsilonInput->text().toAscii().data());
    if(epsilon !=0)
    {
        if(epsilon<0) epsilon *= -1;
        double lowAlpha = -epsilon;
        //double highAlpha = epsilon;
        double highAlpha = 0.0;

        int lowRank = m_processor->getRankForAlpha(lowAlpha);
        int highRank = m_processor->getRankForAlpha(highAlpha);

        m_processor->ProcessEpsilonInterval(epsilon,lowRank,highRank);
    }
}
