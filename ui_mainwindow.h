/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created: Fri Nov 16 03:55:33 2012
**      by: Qt User Interface Compiler version 4.8.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QSlider>
#include <QtGui/QSpinBox>
#include <QtGui/QTabWidget>
#include <QtGui/QTableWidget>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionOpen;
    QAction *actionExit;
    QAction *actionAbout;
    QAction *actionAlphaZero;
    QAction *actionChangeFiltration;
    QAction *actionUndoChange;
    QAction *actionSave_Graph;
    QAction *actionOpen_With_incremented_Radius;
    QAction *actionOpen_As_Points;
    QAction *actionSave_Current_Path;
    QWidget *centralWidget;
    QWidget *FirstViewerPlaceHolder;
    QTabWidget *tabWidget_2;
    QWidget *tab_5;
    QGroupBox *groupBoxUserInput;
    QLineEdit *lineEditEpsilonInput;
    QGroupBox *groupBoxAlphaShape;
    QGridLayout *gridLayout_2;
    QCheckBox *checkBoxAlphaShape;
    QLabel *labelAlphaValue;
    QLineEdit *lineEditAlphaValue;
    QCheckBox *checkBoxWireFrame;
    QGroupBox *groupBoxPockets;
    QGridLayout *gridLayout;
    QCheckBox *checkBoxAllPockets;
    QCheckBox *checkBoxOnlyPockets;
    QCheckBox *checkBoxAllMouths;
    QCheckBox *checkBoxPocketWireFrame;
    QCheckBox *checkBoxOnlyVoids;
    QCheckBox *checkBoxVolumeEnabled;
    QLabel *labelAlphaRank;
    QSlider *alphaSlider;
    QLabel *labelPersistence;
    QSlider *persistenceSlider;
    QGroupBox *groupBoxSkinSurface;
    QVBoxLayout *verticalLayout_4;
    QCheckBox *checkBoxAlphaSkinSolid;
    QCheckBox *checkBoxPocketSkinSolid;
    QCheckBox *checkBoxAlphaSkinWireFrame;
    QCheckBox *checkBoxPocketSkinWireFrame;
    QRadioButton *radioButtonSmooth;
    QRadioButton *radioButtonFlat;
    QLineEdit *lineEditAlphaRank;
    QLineEdit *lineEditPersistenceRank;
    QGroupBox *groupBoxIndivPockets;
    QHBoxLayout *horizontalLayout_4;
    QRadioButton *radioButtonAllPockets;
    QRadioButton *radioButtonIndivPockets;
    QSpinBox *spinBoxPockets;
    QWidget *tab_6;
    QTabWidget *tabWidget;
    QWidget *tab;
    QGroupBox *groupBox;
    QVBoxLayout *verticalLayout_2;
    QCheckBox *checkPowerDiag;
    QCheckBox *checkComplementSpacePD;
    QCheckBox *checkInsideVerts;
    QCheckBox *checkPruneIsolatedVerts;
    QCheckBox *checkEdgeIntersect;
    QGroupBox *cHullBox;
    QVBoxLayout *verticalLayout;
    QCheckBox *checkCHull;
    QCheckBox *checkCHullWireFrame;
    QCheckBox *checkCHullNorm;
    QGroupBox *groupBoxProtein;
    QVBoxLayout *verticalLayout_5;
    QCheckBox *checkSpaceFill;
    QCheckBox *checkSpaceFillPD;
    QWidget *tab_3;
    QGroupBox *groupBox_2;
    QVBoxLayout *verticalLayout_3;
    QHBoxLayout *horizontalLayout;
    QLabel *startLabel;
    QLabel *startIndexLabel;
    QLabel *label_2;
    QLabel *targetIndexLabel;
    QHBoxLayout *horizontalLayout_6;
    QPushButton *buttonShortestEscapePath;
    QPushButton *buttonShortestEscapePathAll;
    QHBoxLayout *horizontalLayout_7;
    QLabel *label;
    QLineEdit *maxIterDijkstra;
    QPushButton *buttonShortestEscapePathRepeated;
    QHBoxLayout *horizontalLayout_8;
    QPushButton *buttonShortestSTPath;
    QHBoxLayout *horizontalLayout_9;
    QCheckBox *checkShowPath;
    QCheckBox *checkShowPathSpheres;
    QHBoxLayout *horizontalLayout_2;
    QCheckBox *checkShowPathSkin;
    QCheckBox *checkShowPathSkinWF;
    QCustomPlot *graphWidget;
    QWidget *tab_7;
    QHBoxLayout *horizontalLayout_3;
    QGroupBox *groupBoxMoleclueProperties;
    QLabel *labelMolTotalVolume;
    QLabel *labelMolTotalSurfaceArea;
    QLineEdit *lineEditMolTotalVolume;
    QLineEdit *lineEditMolTotSurfArea;
    QGroupBox *groupBoxPocketProperties;
    QLabel *labelPocketTotalVolume;
    QLabel *labelPocketTotalSurfaceArea;
    QLineEdit *lineEditPocketTotalVolume;
    QLineEdit *lineEditPocTotSurfArea;
    QTableWidget *tableWidgetPocket;
    QGroupBox *groupBoxMouthProperties;
    QLabel *labelMouthTotalSurfaceArea;
    QLineEdit *lineEditMouthTotalSurfArea;
    QTableWidget *tableWidgetMouth;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QMenu *menuHelp;
    QMenu *menu_Tools;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(1270, 939);
        MainWindow->setMinimumSize(QSize(1270, 721));
        MainWindow->setMaximumSize(QSize(1270, 1000));
        actionOpen = new QAction(MainWindow);
        actionOpen->setObjectName(QString::fromUtf8("actionOpen"));
        actionExit = new QAction(MainWindow);
        actionExit->setObjectName(QString::fromUtf8("actionExit"));
        actionAbout = new QAction(MainWindow);
        actionAbout->setObjectName(QString::fromUtf8("actionAbout"));
        actionAlphaZero = new QAction(MainWindow);
        actionAlphaZero->setObjectName(QString::fromUtf8("actionAlphaZero"));
        actionAlphaZero->setCheckable(false);
        actionChangeFiltration = new QAction(MainWindow);
        actionChangeFiltration->setObjectName(QString::fromUtf8("actionChangeFiltration"));
        actionUndoChange = new QAction(MainWindow);
        actionUndoChange->setObjectName(QString::fromUtf8("actionUndoChange"));
        actionSave_Graph = new QAction(MainWindow);
        actionSave_Graph->setObjectName(QString::fromUtf8("actionSave_Graph"));
        actionOpen_With_incremented_Radius = new QAction(MainWindow);
        actionOpen_With_incremented_Radius->setObjectName(QString::fromUtf8("actionOpen_With_incremented_Radius"));
        actionOpen_As_Points = new QAction(MainWindow);
        actionOpen_As_Points->setObjectName(QString::fromUtf8("actionOpen_As_Points"));
        actionSave_Current_Path = new QAction(MainWindow);
        actionSave_Current_Path->setObjectName(QString::fromUtf8("actionSave_Current_Path"));
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        FirstViewerPlaceHolder = new QWidget(centralWidget);
        FirstViewerPlaceHolder->setObjectName(QString::fromUtf8("FirstViewerPlaceHolder"));
        FirstViewerPlaceHolder->setGeometry(QRect(10, 300, 1251, 611));
        FirstViewerPlaceHolder->setMinimumSize(QSize(440, 440));
        FirstViewerPlaceHolder->setMaximumSize(QSize(1300, 800));
        tabWidget_2 = new QTabWidget(centralWidget);
        tabWidget_2->setObjectName(QString::fromUtf8("tabWidget_2"));
        tabWidget_2->setGeometry(QRect(10, 0, 1251, 291));
        tab_5 = new QWidget();
        tab_5->setObjectName(QString::fromUtf8("tab_5"));
        groupBoxUserInput = new QGroupBox(tab_5);
        groupBoxUserInput->setObjectName(QString::fromUtf8("groupBoxUserInput"));
        groupBoxUserInput->setGeometry(QRect(690, 70, 131, 71));
        lineEditEpsilonInput = new QLineEdit(groupBoxUserInput);
        lineEditEpsilonInput->setObjectName(QString::fromUtf8("lineEditEpsilonInput"));
        lineEditEpsilonInput->setGeometry(QRect(0, 30, 131, 25));
        groupBoxAlphaShape = new QGroupBox(tab_5);
        groupBoxAlphaShape->setObjectName(QString::fromUtf8("groupBoxAlphaShape"));
        groupBoxAlphaShape->setGeometry(QRect(10, 70, 321, 85));
        gridLayout_2 = new QGridLayout(groupBoxAlphaShape);
        gridLayout_2->setSpacing(6);
        gridLayout_2->setContentsMargins(11, 11, 11, 11);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        checkBoxAlphaShape = new QCheckBox(groupBoxAlphaShape);
        checkBoxAlphaShape->setObjectName(QString::fromUtf8("checkBoxAlphaShape"));

        gridLayout_2->addWidget(checkBoxAlphaShape, 0, 0, 1, 1);

        labelAlphaValue = new QLabel(groupBoxAlphaShape);
        labelAlphaValue->setObjectName(QString::fromUtf8("labelAlphaValue"));

        gridLayout_2->addWidget(labelAlphaValue, 0, 1, 1, 1);

        lineEditAlphaValue = new QLineEdit(groupBoxAlphaShape);
        lineEditAlphaValue->setObjectName(QString::fromUtf8("lineEditAlphaValue"));

        gridLayout_2->addWidget(lineEditAlphaValue, 1, 1, 1, 1);

        checkBoxWireFrame = new QCheckBox(groupBoxAlphaShape);
        checkBoxWireFrame->setObjectName(QString::fromUtf8("checkBoxWireFrame"));

        gridLayout_2->addWidget(checkBoxWireFrame, 1, 0, 1, 1);

        groupBoxPockets = new QGroupBox(tab_5);
        groupBoxPockets->setObjectName(QString::fromUtf8("groupBoxPockets"));
        groupBoxPockets->setGeometry(QRect(360, 70, 309, 105));
        gridLayout = new QGridLayout(groupBoxPockets);
        gridLayout->setSpacing(6);
        gridLayout->setContentsMargins(11, 11, 11, 11);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        checkBoxAllPockets = new QCheckBox(groupBoxPockets);
        checkBoxAllPockets->setObjectName(QString::fromUtf8("checkBoxAllPockets"));

        gridLayout->addWidget(checkBoxAllPockets, 0, 0, 1, 1);

        checkBoxOnlyPockets = new QCheckBox(groupBoxPockets);
        checkBoxOnlyPockets->setObjectName(QString::fromUtf8("checkBoxOnlyPockets"));

        gridLayout->addWidget(checkBoxOnlyPockets, 0, 1, 1, 1);

        checkBoxAllMouths = new QCheckBox(groupBoxPockets);
        checkBoxAllMouths->setObjectName(QString::fromUtf8("checkBoxAllMouths"));

        gridLayout->addWidget(checkBoxAllMouths, 1, 0, 1, 1);

        checkBoxPocketWireFrame = new QCheckBox(groupBoxPockets);
        checkBoxPocketWireFrame->setObjectName(QString::fromUtf8("checkBoxPocketWireFrame"));

        gridLayout->addWidget(checkBoxPocketWireFrame, 2, 0, 1, 1);

        checkBoxOnlyVoids = new QCheckBox(groupBoxPockets);
        checkBoxOnlyVoids->setObjectName(QString::fromUtf8("checkBoxOnlyVoids"));

        gridLayout->addWidget(checkBoxOnlyVoids, 1, 1, 1, 1);

        checkBoxVolumeEnabled = new QCheckBox(groupBoxPockets);
        checkBoxVolumeEnabled->setObjectName(QString::fromUtf8("checkBoxVolumeEnabled"));

        gridLayout->addWidget(checkBoxVolumeEnabled, 2, 1, 1, 1);

        labelAlphaRank = new QLabel(tab_5);
        labelAlphaRank->setObjectName(QString::fromUtf8("labelAlphaRank"));
        labelAlphaRank->setGeometry(QRect(10, 10, 75, 20));
        alphaSlider = new QSlider(tab_5);
        alphaSlider->setObjectName(QString::fromUtf8("alphaSlider"));
        alphaSlider->setGeometry(QRect(90, 10, 991, 16));
        alphaSlider->setOrientation(Qt::Horizontal);
        alphaSlider->setTickPosition(QSlider::TicksBelow);
        labelPersistence = new QLabel(tab_5);
        labelPersistence->setObjectName(QString::fromUtf8("labelPersistence"));
        labelPersistence->setGeometry(QRect(10, 40, 75, 20));
        persistenceSlider = new QSlider(tab_5);
        persistenceSlider->setObjectName(QString::fromUtf8("persistenceSlider"));
        persistenceSlider->setGeometry(QRect(90, 40, 991, 16));
        persistenceSlider->setOrientation(Qt::Horizontal);
        persistenceSlider->setTickPosition(QSlider::TicksBelow);
        groupBoxSkinSurface = new QGroupBox(tab_5);
        groupBoxSkinSurface->setObjectName(QString::fromUtf8("groupBoxSkinSurface"));
        groupBoxSkinSurface->setGeometry(QRect(860, 70, 188, 180));
        verticalLayout_4 = new QVBoxLayout(groupBoxSkinSurface);
        verticalLayout_4->setSpacing(6);
        verticalLayout_4->setContentsMargins(11, 11, 11, 11);
        verticalLayout_4->setObjectName(QString::fromUtf8("verticalLayout_4"));
        checkBoxAlphaSkinSolid = new QCheckBox(groupBoxSkinSurface);
        checkBoxAlphaSkinSolid->setObjectName(QString::fromUtf8("checkBoxAlphaSkinSolid"));

        verticalLayout_4->addWidget(checkBoxAlphaSkinSolid);

        checkBoxPocketSkinSolid = new QCheckBox(groupBoxSkinSurface);
        checkBoxPocketSkinSolid->setObjectName(QString::fromUtf8("checkBoxPocketSkinSolid"));

        verticalLayout_4->addWidget(checkBoxPocketSkinSolid);

        checkBoxAlphaSkinWireFrame = new QCheckBox(groupBoxSkinSurface);
        checkBoxAlphaSkinWireFrame->setObjectName(QString::fromUtf8("checkBoxAlphaSkinWireFrame"));

        verticalLayout_4->addWidget(checkBoxAlphaSkinWireFrame);

        checkBoxPocketSkinWireFrame = new QCheckBox(groupBoxSkinSurface);
        checkBoxPocketSkinWireFrame->setObjectName(QString::fromUtf8("checkBoxPocketSkinWireFrame"));

        verticalLayout_4->addWidget(checkBoxPocketSkinWireFrame);

        radioButtonSmooth = new QRadioButton(groupBoxSkinSurface);
        radioButtonSmooth->setObjectName(QString::fromUtf8("radioButtonSmooth"));

        verticalLayout_4->addWidget(radioButtonSmooth);

        radioButtonFlat = new QRadioButton(groupBoxSkinSurface);
        radioButtonFlat->setObjectName(QString::fromUtf8("radioButtonFlat"));

        verticalLayout_4->addWidget(radioButtonFlat);

        lineEditAlphaRank = new QLineEdit(tab_5);
        lineEditAlphaRank->setObjectName(QString::fromUtf8("lineEditAlphaRank"));
        lineEditAlphaRank->setGeometry(QRect(1090, 10, 151, 25));
        lineEditPersistenceRank = new QLineEdit(tab_5);
        lineEditPersistenceRank->setObjectName(QString::fromUtf8("lineEditPersistenceRank"));
        lineEditPersistenceRank->setGeometry(QRect(1090, 40, 151, 25));
        groupBoxIndivPockets = new QGroupBox(tab_5);
        groupBoxIndivPockets->setObjectName(QString::fromUtf8("groupBoxIndivPockets"));
        groupBoxIndivPockets->setGeometry(QRect(10, 180, 329, 60));
        horizontalLayout_4 = new QHBoxLayout(groupBoxIndivPockets);
        horizontalLayout_4->setSpacing(6);
        horizontalLayout_4->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        radioButtonAllPockets = new QRadioButton(groupBoxIndivPockets);
        radioButtonAllPockets->setObjectName(QString::fromUtf8("radioButtonAllPockets"));

        horizontalLayout_4->addWidget(radioButtonAllPockets);

        radioButtonIndivPockets = new QRadioButton(groupBoxIndivPockets);
        radioButtonIndivPockets->setObjectName(QString::fromUtf8("radioButtonIndivPockets"));

        horizontalLayout_4->addWidget(radioButtonIndivPockets);

        spinBoxPockets = new QSpinBox(groupBoxIndivPockets);
        spinBoxPockets->setObjectName(QString::fromUtf8("spinBoxPockets"));

        horizontalLayout_4->addWidget(spinBoxPockets);

        tabWidget_2->addTab(tab_5, QString());
        tab_6 = new QWidget();
        tab_6->setObjectName(QString::fromUtf8("tab_6"));
        tabWidget = new QTabWidget(tab_6);
        tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
        tabWidget->setGeometry(QRect(10, 0, 1231, 251));
        tabWidget->setTabPosition(QTabWidget::South);
        tabWidget->setDocumentMode(false);
        tab = new QWidget();
        tab->setObjectName(QString::fromUtf8("tab"));
        groupBox = new QGroupBox(tab);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        groupBox->setGeometry(QRect(10, 10, 199, 155));
        verticalLayout_2 = new QVBoxLayout(groupBox);
        verticalLayout_2->setSpacing(6);
        verticalLayout_2->setContentsMargins(11, 11, 11, 11);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        checkPowerDiag = new QCheckBox(groupBox);
        checkPowerDiag->setObjectName(QString::fromUtf8("checkPowerDiag"));
        checkPowerDiag->setChecked(true);

        verticalLayout_2->addWidget(checkPowerDiag);

        checkComplementSpacePD = new QCheckBox(groupBox);
        checkComplementSpacePD->setObjectName(QString::fromUtf8("checkComplementSpacePD"));
        checkComplementSpacePD->setChecked(true);

        verticalLayout_2->addWidget(checkComplementSpacePD);

        checkInsideVerts = new QCheckBox(groupBox);
        checkInsideVerts->setObjectName(QString::fromUtf8("checkInsideVerts"));
        checkInsideVerts->setChecked(true);

        verticalLayout_2->addWidget(checkInsideVerts);

        checkPruneIsolatedVerts = new QCheckBox(groupBox);
        checkPruneIsolatedVerts->setObjectName(QString::fromUtf8("checkPruneIsolatedVerts"));

        verticalLayout_2->addWidget(checkPruneIsolatedVerts);

        checkEdgeIntersect = new QCheckBox(groupBox);
        checkEdgeIntersect->setObjectName(QString::fromUtf8("checkEdgeIntersect"));

        verticalLayout_2->addWidget(checkEdgeIntersect);

        checkPowerDiag->raise();
        checkInsideVerts->raise();
        checkPruneIsolatedVerts->raise();
        checkEdgeIntersect->raise();
        checkComplementSpacePD->raise();
        cHullBox = new QGroupBox(tab);
        cHullBox->setObjectName(QString::fromUtf8("cHullBox"));
        cHullBox->setGeometry(QRect(230, 10, 192, 105));
        verticalLayout = new QVBoxLayout(cHullBox);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        checkCHull = new QCheckBox(cHullBox);
        checkCHull->setObjectName(QString::fromUtf8("checkCHull"));

        verticalLayout->addWidget(checkCHull);

        checkCHullWireFrame = new QCheckBox(cHullBox);
        checkCHullWireFrame->setObjectName(QString::fromUtf8("checkCHullWireFrame"));

        verticalLayout->addWidget(checkCHullWireFrame);

        checkCHullNorm = new QCheckBox(cHullBox);
        checkCHullNorm->setObjectName(QString::fromUtf8("checkCHullNorm"));

        verticalLayout->addWidget(checkCHullNorm);

        groupBoxProtein = new QGroupBox(tab);
        groupBoxProtein->setObjectName(QString::fromUtf8("groupBoxProtein"));
        groupBoxProtein->setGeometry(QRect(430, 10, 221, 81));
        verticalLayout_5 = new QVBoxLayout(groupBoxProtein);
        verticalLayout_5->setSpacing(6);
        verticalLayout_5->setContentsMargins(11, 11, 11, 11);
        verticalLayout_5->setObjectName(QString::fromUtf8("verticalLayout_5"));
        checkSpaceFill = new QCheckBox(groupBoxProtein);
        checkSpaceFill->setObjectName(QString::fromUtf8("checkSpaceFill"));
        checkSpaceFill->setChecked(true);

        verticalLayout_5->addWidget(checkSpaceFill);

        checkSpaceFillPD = new QCheckBox(groupBoxProtein);
        checkSpaceFillPD->setObjectName(QString::fromUtf8("checkSpaceFillPD"));

        verticalLayout_5->addWidget(checkSpaceFillPD);

        tabWidget->addTab(tab, QString());
        tab_3 = new QWidget();
        tab_3->setObjectName(QString::fromUtf8("tab_3"));
        groupBox_2 = new QGroupBox(tab_3);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        groupBox_2->setGeometry(QRect(2, 9, 491, 211));
        verticalLayout_3 = new QVBoxLayout(groupBox_2);
        verticalLayout_3->setSpacing(6);
        verticalLayout_3->setContentsMargins(11, 11, 11, 11);
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setSpacing(6);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        startLabel = new QLabel(groupBox_2);
        startLabel->setObjectName(QString::fromUtf8("startLabel"));
        startLabel->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout->addWidget(startLabel);

        startIndexLabel = new QLabel(groupBox_2);
        startIndexLabel->setObjectName(QString::fromUtf8("startIndexLabel"));

        horizontalLayout->addWidget(startIndexLabel);

        label_2 = new QLabel(groupBox_2);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout->addWidget(label_2);

        targetIndexLabel = new QLabel(groupBox_2);
        targetIndexLabel->setObjectName(QString::fromUtf8("targetIndexLabel"));

        horizontalLayout->addWidget(targetIndexLabel);


        verticalLayout_3->addLayout(horizontalLayout);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setSpacing(6);
        horizontalLayout_6->setObjectName(QString::fromUtf8("horizontalLayout_6"));
        buttonShortestEscapePath = new QPushButton(groupBox_2);
        buttonShortestEscapePath->setObjectName(QString::fromUtf8("buttonShortestEscapePath"));

        horizontalLayout_6->addWidget(buttonShortestEscapePath);

        buttonShortestEscapePathAll = new QPushButton(groupBox_2);
        buttonShortestEscapePathAll->setObjectName(QString::fromUtf8("buttonShortestEscapePathAll"));

        horizontalLayout_6->addWidget(buttonShortestEscapePathAll);


        verticalLayout_3->addLayout(horizontalLayout_6);

        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setSpacing(6);
        horizontalLayout_7->setObjectName(QString::fromUtf8("horizontalLayout_7"));
        label = new QLabel(groupBox_2);
        label->setObjectName(QString::fromUtf8("label"));

        horizontalLayout_7->addWidget(label);

        maxIterDijkstra = new QLineEdit(groupBox_2);
        maxIterDijkstra->setObjectName(QString::fromUtf8("maxIterDijkstra"));

        horizontalLayout_7->addWidget(maxIterDijkstra);

        buttonShortestEscapePathRepeated = new QPushButton(groupBox_2);
        buttonShortestEscapePathRepeated->setObjectName(QString::fromUtf8("buttonShortestEscapePathRepeated"));

        horizontalLayout_7->addWidget(buttonShortestEscapePathRepeated);


        verticalLayout_3->addLayout(horizontalLayout_7);

        horizontalLayout_8 = new QHBoxLayout();
        horizontalLayout_8->setSpacing(6);
        horizontalLayout_8->setObjectName(QString::fromUtf8("horizontalLayout_8"));
        buttonShortestSTPath = new QPushButton(groupBox_2);
        buttonShortestSTPath->setObjectName(QString::fromUtf8("buttonShortestSTPath"));

        horizontalLayout_8->addWidget(buttonShortestSTPath);


        verticalLayout_3->addLayout(horizontalLayout_8);

        horizontalLayout_9 = new QHBoxLayout();
        horizontalLayout_9->setSpacing(6);
        horizontalLayout_9->setObjectName(QString::fromUtf8("horizontalLayout_9"));
        checkShowPath = new QCheckBox(groupBox_2);
        checkShowPath->setObjectName(QString::fromUtf8("checkShowPath"));
        checkShowPath->setChecked(true);

        horizontalLayout_9->addWidget(checkShowPath);

        checkShowPathSpheres = new QCheckBox(groupBox_2);
        checkShowPathSpheres->setObjectName(QString::fromUtf8("checkShowPathSpheres"));

        horizontalLayout_9->addWidget(checkShowPathSpheres);


        verticalLayout_3->addLayout(horizontalLayout_9);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        checkShowPathSkin = new QCheckBox(groupBox_2);
        checkShowPathSkin->setObjectName(QString::fromUtf8("checkShowPathSkin"));
        checkShowPathSkin->setChecked(true);

        horizontalLayout_2->addWidget(checkShowPathSkin);

        checkShowPathSkinWF = new QCheckBox(groupBox_2);
        checkShowPathSkinWF->setObjectName(QString::fromUtf8("checkShowPathSkinWF"));
        checkShowPathSkinWF->setChecked(true);

        horizontalLayout_2->addWidget(checkShowPathSkinWF);


        verticalLayout_3->addLayout(horizontalLayout_2);

        graphWidget = new QCustomPlot(tab_3);
        graphWidget->setObjectName(QString::fromUtf8("graphWidget"));
        graphWidget->setGeometry(QRect(510, 0, 711, 221));
        tabWidget->addTab(tab_3, QString());
        tabWidget_2->addTab(tab_6, QString());
        tab_7 = new QWidget();
        tab_7->setObjectName(QString::fromUtf8("tab_7"));
        horizontalLayout_3 = new QHBoxLayout(tab_7);
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        groupBoxMoleclueProperties = new QGroupBox(tab_7);
        groupBoxMoleclueProperties->setObjectName(QString::fromUtf8("groupBoxMoleclueProperties"));
        labelMolTotalVolume = new QLabel(groupBoxMoleclueProperties);
        labelMolTotalVolume->setObjectName(QString::fromUtf8("labelMolTotalVolume"));
        labelMolTotalVolume->setEnabled(true);
        labelMolTotalVolume->setGeometry(QRect(10, 30, 90, 20));
        labelMolTotalSurfaceArea = new QLabel(groupBoxMoleclueProperties);
        labelMolTotalSurfaceArea->setObjectName(QString::fromUtf8("labelMolTotalSurfaceArea"));
        labelMolTotalSurfaceArea->setGeometry(QRect(10, 70, 120, 20));
        lineEditMolTotalVolume = new QLineEdit(groupBoxMoleclueProperties);
        lineEditMolTotalVolume->setObjectName(QString::fromUtf8("lineEditMolTotalVolume"));
        lineEditMolTotalVolume->setGeometry(QRect(140, 30, 180, 25));
        lineEditMolTotSurfArea = new QLineEdit(groupBoxMoleclueProperties);
        lineEditMolTotSurfArea->setObjectName(QString::fromUtf8("lineEditMolTotSurfArea"));
        lineEditMolTotSurfArea->setGeometry(QRect(140, 70, 180, 25));

        horizontalLayout_3->addWidget(groupBoxMoleclueProperties);

        groupBoxPocketProperties = new QGroupBox(tab_7);
        groupBoxPocketProperties->setObjectName(QString::fromUtf8("groupBoxPocketProperties"));
        labelPocketTotalVolume = new QLabel(groupBoxPocketProperties);
        labelPocketTotalVolume->setObjectName(QString::fromUtf8("labelPocketTotalVolume"));
        labelPocketTotalVolume->setGeometry(QRect(10, 30, 90, 20));
        labelPocketTotalSurfaceArea = new QLabel(groupBoxPocketProperties);
        labelPocketTotalSurfaceArea->setObjectName(QString::fromUtf8("labelPocketTotalSurfaceArea"));
        labelPocketTotalSurfaceArea->setGeometry(QRect(10, 60, 120, 20));
        lineEditPocketTotalVolume = new QLineEdit(groupBoxPocketProperties);
        lineEditPocketTotalVolume->setObjectName(QString::fromUtf8("lineEditPocketTotalVolume"));
        lineEditPocketTotalVolume->setGeometry(QRect(140, 20, 180, 25));
        lineEditPocTotSurfArea = new QLineEdit(groupBoxPocketProperties);
        lineEditPocTotSurfArea->setObjectName(QString::fromUtf8("lineEditPocTotSurfArea"));
        lineEditPocTotSurfArea->setGeometry(QRect(140, 60, 180, 25));
        tableWidgetPocket = new QTableWidget(groupBoxPocketProperties);
        if (tableWidgetPocket->columnCount() < 5)
            tableWidgetPocket->setColumnCount(5);
        QTableWidgetItem *__qtablewidgetitem = new QTableWidgetItem();
        tableWidgetPocket->setHorizontalHeaderItem(0, __qtablewidgetitem);
        QTableWidgetItem *__qtablewidgetitem1 = new QTableWidgetItem();
        tableWidgetPocket->setHorizontalHeaderItem(1, __qtablewidgetitem1);
        QTableWidgetItem *__qtablewidgetitem2 = new QTableWidgetItem();
        tableWidgetPocket->setHorizontalHeaderItem(2, __qtablewidgetitem2);
        QTableWidgetItem *__qtablewidgetitem3 = new QTableWidgetItem();
        tableWidgetPocket->setHorizontalHeaderItem(3, __qtablewidgetitem3);
        QTableWidgetItem *__qtablewidgetitem4 = new QTableWidgetItem();
        tableWidgetPocket->setHorizontalHeaderItem(4, __qtablewidgetitem4);
        tableWidgetPocket->setObjectName(QString::fromUtf8("tableWidgetPocket"));
        tableWidgetPocket->setGeometry(QRect(10, 90, 391, 151));
        tableWidgetPocket->setColumnCount(5);

        horizontalLayout_3->addWidget(groupBoxPocketProperties);

        groupBoxMouthProperties = new QGroupBox(tab_7);
        groupBoxMouthProperties->setObjectName(QString::fromUtf8("groupBoxMouthProperties"));
        labelMouthTotalSurfaceArea = new QLabel(groupBoxMouthProperties);
        labelMouthTotalSurfaceArea->setObjectName(QString::fromUtf8("labelMouthTotalSurfaceArea"));
        labelMouthTotalSurfaceArea->setGeometry(QRect(10, 30, 120, 20));
        lineEditMouthTotalSurfArea = new QLineEdit(groupBoxMouthProperties);
        lineEditMouthTotalSurfArea->setObjectName(QString::fromUtf8("lineEditMouthTotalSurfArea"));
        lineEditMouthTotalSurfArea->setGeometry(QRect(140, 30, 180, 25));
        tableWidgetMouth = new QTableWidget(groupBoxMouthProperties);
        if (tableWidgetMouth->columnCount() < 4)
            tableWidgetMouth->setColumnCount(4);
        QTableWidgetItem *__qtablewidgetitem5 = new QTableWidgetItem();
        tableWidgetMouth->setHorizontalHeaderItem(0, __qtablewidgetitem5);
        QTableWidgetItem *__qtablewidgetitem6 = new QTableWidgetItem();
        tableWidgetMouth->setHorizontalHeaderItem(1, __qtablewidgetitem6);
        QTableWidgetItem *__qtablewidgetitem7 = new QTableWidgetItem();
        tableWidgetMouth->setHorizontalHeaderItem(2, __qtablewidgetitem7);
        QTableWidgetItem *__qtablewidgetitem8 = new QTableWidgetItem();
        tableWidgetMouth->setHorizontalHeaderItem(3, __qtablewidgetitem8);
        tableWidgetMouth->setObjectName(QString::fromUtf8("tableWidgetMouth"));
        tableWidgetMouth->setGeometry(QRect(10, 89, 391, 151));
        tableWidgetMouth->setColumnCount(4);

        horizontalLayout_3->addWidget(groupBoxMouthProperties);

        tabWidget_2->addTab(tab_7, QString());
        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1270, 22));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuHelp = new QMenu(menuBar);
        menuHelp->setObjectName(QString::fromUtf8("menuHelp"));
        menu_Tools = new QMenu(menuBar);
        menu_Tools->setObjectName(QString::fromUtf8("menu_Tools"));
        MainWindow->setMenuBar(menuBar);
#ifndef QT_NO_SHORTCUT
        labelAlphaValue->setBuddy(lineEditAlphaValue);
        labelAlphaRank->setBuddy(alphaSlider);
        labelPersistence->setBuddy(persistenceSlider);
        labelMolTotalVolume->setBuddy(lineEditMolTotalVolume);
        labelMolTotalSurfaceArea->setBuddy(lineEditMolTotSurfArea);
        labelPocketTotalVolume->setBuddy(lineEditPocketTotalVolume);
        labelPocketTotalSurfaceArea->setBuddy(lineEditPocTotSurfArea);
        labelMouthTotalSurfaceArea->setBuddy(lineEditMouthTotalSurfArea);
#endif // QT_NO_SHORTCUT

        menuBar->addAction(menuFile->menuAction());
        menuBar->addAction(menu_Tools->menuAction());
        menuBar->addAction(menuHelp->menuAction());
        menuFile->addAction(actionOpen);
        menuFile->addAction(actionOpen_With_incremented_Radius);
        menuFile->addAction(actionOpen_As_Points);
        menuFile->addSeparator();
        menuFile->addAction(actionSave_Graph);
        menuFile->addAction(actionSave_Current_Path);
        menuFile->addSeparator();
        menuFile->addAction(actionExit);
        menuHelp->addAction(actionAbout);
        menu_Tools->addAction(actionAlphaZero);
        menu_Tools->addAction(actionChangeFiltration);
        menu_Tools->addAction(actionUndoChange);

        retranslateUi(MainWindow);

        tabWidget_2->setCurrentIndex(1);
        tabWidget->setCurrentIndex(1);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "Pocket Viewer", 0, QApplication::UnicodeUTF8));
        actionOpen->setText(QApplication::translate("MainWindow", "&Open", 0, QApplication::UnicodeUTF8));
        actionExit->setText(QApplication::translate("MainWindow", "&Exit", 0, QApplication::UnicodeUTF8));
        actionAbout->setText(QApplication::translate("MainWindow", "&About", 0, QApplication::UnicodeUTF8));
        actionAlphaZero->setText(QApplication::translate("MainWindow", "&AlphaValue = 0.0", 0, QApplication::UnicodeUTF8));
        actionChangeFiltration->setText(QApplication::translate("MainWindow", "&Change Filtration", 0, QApplication::UnicodeUTF8));
        actionUndoChange->setText(QApplication::translate("MainWindow", "&Undo Change", 0, QApplication::UnicodeUTF8));
        actionSave_Graph->setText(QApplication::translate("MainWindow", "Save Graph", 0, QApplication::UnicodeUTF8));
        actionOpen_With_incremented_Radius->setText(QApplication::translate("MainWindow", "Open (With incremented Radius)", 0, QApplication::UnicodeUTF8));
        actionOpen_As_Points->setText(QApplication::translate("MainWindow", "Open As Points", 0, QApplication::UnicodeUTF8));
        actionSave_Current_Path->setText(QApplication::translate("MainWindow", "Save Current Path(s)", 0, QApplication::UnicodeUTF8));
        groupBoxUserInput->setTitle(QApplication::translate("MainWindow", "Epsilon Input", 0, QApplication::UnicodeUTF8));
        groupBoxAlphaShape->setTitle(QApplication::translate("MainWindow", "AlphaShape", 0, QApplication::UnicodeUTF8));
        checkBoxAlphaShape->setText(QApplication::translate("MainWindow", "AlphaShape", 0, QApplication::UnicodeUTF8));
        labelAlphaValue->setText(QApplication::translate("MainWindow", "Alpha Value :", 0, QApplication::UnicodeUTF8));
        checkBoxWireFrame->setText(QApplication::translate("MainWindow", "WireFrame", 0, QApplication::UnicodeUTF8));
        groupBoxPockets->setTitle(QApplication::translate("MainWindow", "Pockets", 0, QApplication::UnicodeUTF8));
        checkBoxAllPockets->setText(QApplication::translate("MainWindow", "All Pockets", 0, QApplication::UnicodeUTF8));
        checkBoxOnlyPockets->setText(QApplication::translate("MainWindow", "Only Pockets", 0, QApplication::UnicodeUTF8));
        checkBoxAllMouths->setText(QApplication::translate("MainWindow", "Mouths", 0, QApplication::UnicodeUTF8));
        checkBoxPocketWireFrame->setText(QApplication::translate("MainWindow", "PocketWireFrame", 0, QApplication::UnicodeUTF8));
        checkBoxOnlyVoids->setText(QApplication::translate("MainWindow", "Only Voids", 0, QApplication::UnicodeUTF8));
        checkBoxVolumeEnabled->setText(QApplication::translate("MainWindow", "Display Properties", 0, QApplication::UnicodeUTF8));
        labelAlphaRank->setText(QApplication::translate("MainWindow", "Alpha Rank", 0, QApplication::UnicodeUTF8));
        labelPersistence->setText(QApplication::translate("MainWindow", "Persistence", 0, QApplication::UnicodeUTF8));
        groupBoxSkinSurface->setTitle(QApplication::translate("MainWindow", "SkinSurface", 0, QApplication::UnicodeUTF8));
        checkBoxAlphaSkinSolid->setText(QApplication::translate("MainWindow", "AlphaSkinSurface", 0, QApplication::UnicodeUTF8));
        checkBoxPocketSkinSolid->setText(QApplication::translate("MainWindow", "PocketSkinSurface", 0, QApplication::UnicodeUTF8));
        checkBoxAlphaSkinWireFrame->setText(QApplication::translate("MainWindow", "AlphaSkinWireFrame", 0, QApplication::UnicodeUTF8));
        checkBoxPocketSkinWireFrame->setText(QApplication::translate("MainWindow", "PocketSkinWireFrame", 0, QApplication::UnicodeUTF8));
        radioButtonSmooth->setText(QApplication::translate("MainWindow", "Smooth Shading", 0, QApplication::UnicodeUTF8));
        radioButtonFlat->setText(QApplication::translate("MainWindow", "Flat Shading", 0, QApplication::UnicodeUTF8));
        groupBoxIndivPockets->setTitle(QApplication::translate("MainWindow", "Individual Pockets", 0, QApplication::UnicodeUTF8));
        radioButtonAllPockets->setText(QApplication::translate("MainWindow", "All Pockets", 0, QApplication::UnicodeUTF8));
        radioButtonIndivPockets->setText(QApplication::translate("MainWindow", "Individual Pockets", 0, QApplication::UnicodeUTF8));
        tabWidget_2->setTabText(tabWidget_2->indexOf(tab_5), QApplication::translate("MainWindow", "Alpha Shape and Pockets Controls", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QApplication::translate("MainWindow", "PowerDiagram", 0, QApplication::UnicodeUTF8));
        checkPowerDiag->setText(QApplication::translate("MainWindow", "Show Power Diagram", 0, QApplication::UnicodeUTF8));
        checkComplementSpacePD->setText(QApplication::translate("MainWindow", "Complement Space PD", 0, QApplication::UnicodeUTF8));
        checkInsideVerts->setText(QApplication::translate("MainWindow", "Only Inside Vertices", 0, QApplication::UnicodeUTF8));
        checkPruneIsolatedVerts->setText(QApplication::translate("MainWindow", "Prune Isolated Vertices", 0, QApplication::UnicodeUTF8));
        checkEdgeIntersect->setText(QApplication::translate("MainWindow", "Intersect Edges", 0, QApplication::UnicodeUTF8));
        cHullBox->setTitle(QApplication::translate("MainWindow", "Convex Hull", 0, QApplication::UnicodeUTF8));
        checkCHull->setText(QApplication::translate("MainWindow", "Show Convex Hull", 0, QApplication::UnicodeUTF8));
        checkCHullWireFrame->setText(QApplication::translate("MainWindow", "Convex Hull Wireframe", 0, QApplication::UnicodeUTF8));
        checkCHullNorm->setText(QApplication::translate("MainWindow", "Show Hull Normals", 0, QApplication::UnicodeUTF8));
        groupBoxProtein->setTitle(QApplication::translate("MainWindow", "Space Fill Visualization", 0, QApplication::UnicodeUTF8));
        checkSpaceFill->setText(QApplication::translate("MainWindow", "Show Space Fill Molecule", 0, QApplication::UnicodeUTF8));
        checkSpaceFillPD->setText(QApplication::translate("MainWindow", "Show Space Fill PD", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab), QApplication::translate("MainWindow", "Power Diagram", 0, QApplication::UnicodeUTF8));
        groupBox_2->setTitle(QApplication::translate("MainWindow", "Shortest Escape Route", 0, QApplication::UnicodeUTF8));
        startLabel->setText(QApplication::translate("MainWindow", "Start Node:", 0, QApplication::UnicodeUTF8));
        startIndexLabel->setText(QApplication::translate("MainWindow", "-1", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("MainWindow", "Target Node:", 0, QApplication::UnicodeUTF8));
        targetIndexLabel->setText(QApplication::translate("MainWindow", "-1", 0, QApplication::UnicodeUTF8));
        buttonShortestEscapePath->setText(QApplication::translate("MainWindow", "Find Shortest Escape Path", 0, QApplication::UnicodeUTF8));
        buttonShortestEscapePathAll->setText(QApplication::translate("MainWindow", "Find Shortest Escape Path (All)", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("MainWindow", "Max Iterations:", 0, QApplication::UnicodeUTF8));
        maxIterDijkstra->setText(QApplication::translate("MainWindow", "10", 0, QApplication::UnicodeUTF8));
        buttonShortestEscapePathRepeated->setText(QApplication::translate("MainWindow", "Escape paths (Repeated Dijkstra)", 0, QApplication::UnicodeUTF8));
        buttonShortestSTPath->setText(QApplication::translate("MainWindow", "Find Shortest Source-Target path", 0, QApplication::UnicodeUTF8));
        checkShowPath->setText(QApplication::translate("MainWindow", "Show Path(s)", 0, QApplication::UnicodeUTF8));
        checkShowPathSpheres->setText(QApplication::translate("MainWindow", "Show Path Spheres", 0, QApplication::UnicodeUTF8));
        checkShowPathSkin->setText(QApplication::translate("MainWindow", "Show Path Skin", 0, QApplication::UnicodeUTF8));
        checkShowPathSkinWF->setText(QApplication::translate("MainWindow", "Show Path Skin Wireframe", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab_3), QApplication::translate("MainWindow", "Shortest Path", 0, QApplication::UnicodeUTF8));
        tabWidget_2->setTabText(tabWidget_2->indexOf(tab_6), QApplication::translate("MainWindow", "Power Diagram Controls", 0, QApplication::UnicodeUTF8));
        groupBoxMoleclueProperties->setTitle(QApplication::translate("MainWindow", "Molecular Properties", 0, QApplication::UnicodeUTF8));
        labelMolTotalVolume->setText(QApplication::translate("MainWindow", "Total Volume", 0, QApplication::UnicodeUTF8));
        labelMolTotalSurfaceArea->setText(QApplication::translate("MainWindow", "Total Surface Area", 0, QApplication::UnicodeUTF8));
        groupBoxPocketProperties->setTitle(QApplication::translate("MainWindow", "Pocket Properties", 0, QApplication::UnicodeUTF8));
        labelPocketTotalVolume->setText(QApplication::translate("MainWindow", "Total Volume", 0, QApplication::UnicodeUTF8));
        labelPocketTotalSurfaceArea->setText(QApplication::translate("MainWindow", "Total Surface Area", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem = tableWidgetPocket->horizontalHeaderItem(0);
        ___qtablewidgetitem->setText(QApplication::translate("MainWindow", "Identifier", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem1 = tableWidgetPocket->horizontalHeaderItem(1);
        ___qtablewidgetitem1->setText(QApplication::translate("MainWindow", "Select", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem2 = tableWidgetPocket->horizontalHeaderItem(2);
        ___qtablewidgetitem2->setText(QApplication::translate("MainWindow", "No.of Mouths", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem3 = tableWidgetPocket->horizontalHeaderItem(3);
        ___qtablewidgetitem3->setText(QApplication::translate("MainWindow", "Volume", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem4 = tableWidgetPocket->horizontalHeaderItem(4);
        ___qtablewidgetitem4->setText(QApplication::translate("MainWindow", "Surface Area", 0, QApplication::UnicodeUTF8));
        groupBoxMouthProperties->setTitle(QApplication::translate("MainWindow", "Mouth Properties", 0, QApplication::UnicodeUTF8));
        labelMouthTotalSurfaceArea->setText(QApplication::translate("MainWindow", "Total Surface Area", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem5 = tableWidgetMouth->horizontalHeaderItem(0);
        ___qtablewidgetitem5->setText(QApplication::translate("MainWindow", "Identifier", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem6 = tableWidgetMouth->horizontalHeaderItem(1);
        ___qtablewidgetitem6->setText(QApplication::translate("MainWindow", "Select", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem7 = tableWidgetMouth->horizontalHeaderItem(2);
        ___qtablewidgetitem7->setText(QApplication::translate("MainWindow", "Parent Pocket", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem8 = tableWidgetMouth->horizontalHeaderItem(3);
        ___qtablewidgetitem8->setText(QApplication::translate("MainWindow", "Surface Area", 0, QApplication::UnicodeUTF8));
        tabWidget_2->setTabText(tabWidget_2->indexOf(tab_7), QApplication::translate("MainWindow", "Pocket and Mouth Properties", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("MainWindow", "&File", 0, QApplication::UnicodeUTF8));
        menuHelp->setTitle(QApplication::translate("MainWindow", "&Help", 0, QApplication::UnicodeUTF8));
        menu_Tools->setTitle(QApplication::translate("MainWindow", "&Tools", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
