/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created: Sat Nov 3 04:07:38 2012
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
#include <QtGui/QComboBox>
#include <QtGui/QFrame>
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
    QWidget *centralWidget;
    QSlider *alphaSlider;
    QLabel *labelAlphaRank;
    QLabel *labelPersistence;
    QSlider *persistenceSlider;
    QFrame *lineVOne;
    QFrame *lineHOne;
    QFrame *lineHTwo;
    QLineEdit *lineEditAlphaRank;
    QLineEdit *lineEditPersistenceRank;
    QGroupBox *groupBoxAlphaShape;
    QCheckBox *checkBoxAlphaShape;
    QCheckBox *checkBoxWireFrame;
    QLabel *labelAlphaValue;
    QLineEdit *lineEditAlphaValue;
    QGroupBox *groupBoxPockets;
    QCheckBox *checkBoxAllPockets;
    QCheckBox *checkBoxOnlyPockets;
    QCheckBox *checkBoxOnlyVoids;
    QCheckBox *checkBoxAllMouths;
    QCheckBox *checkBoxPocketWireFrame;
    QCheckBox *checkBoxVolumeEnabled;
    QGroupBox *groupBoxSkinSurface;
    QCheckBox *checkBoxPocketSkinSolid;
    QRadioButton *radioButtonFlat;
    QRadioButton *radioButtonSmooth;
    QCheckBox *checkBoxAlphaSkinSolid;
    QCheckBox *checkBoxAlphaSkinWireFrame;
    QCheckBox *checkBoxPocketSkinWireFrame;
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
    QWidget *FirstViewerPlaceHolder;
    QWidget *IndivPocMenuPlaceHolder;
    QGroupBox *groupBoxIndivPockets;
    QRadioButton *radioButtonAllPockets;
    QRadioButton *radioButtonIndivPockets;
    QSpinBox *spinBoxPockets;
    QGroupBox *groupBoxUserInput;
    QLineEdit *lineEditEpsilonInput;
    QTabWidget *tabWidget;
    QWidget *tab;
    QGroupBox *groupBox;
    QVBoxLayout *verticalLayout_2;
    QCheckBox *checkPowerDiag;
    QCheckBox *checkComplementSpacePD;
    QCheckBox *checkInsideVerts;
    QCheckBox *checkPruneIsolatedVerts;
    QCheckBox *checkEdgeIntersect;
    QWidget *tab_2;
    QGroupBox *cHullBox;
    QVBoxLayout *verticalLayout;
    QCheckBox *checkCHull;
    QCheckBox *checkCHullWireFrame;
    QCheckBox *checkCHullNorm;
    QWidget *tab_3;
    QGroupBox *groupBox_2;
    QVBoxLayout *verticalLayout_3;
    QHBoxLayout *horizontalLayout;
    QLabel *startLabel;
    QComboBox *startCombo;
    QPushButton *buttonShortestEscapePath;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_2;
    QComboBox *targetCombo;
    QPushButton *buttonShortestSTPath;
    QWidget *tab_4;
    QCustomPlot *graphWidget;
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
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        alphaSlider = new QSlider(centralWidget);
        alphaSlider->setObjectName(QString::fromUtf8("alphaSlider"));
        alphaSlider->setGeometry(QRect(90, 10, 660, 16));
        alphaSlider->setOrientation(Qt::Horizontal);
        alphaSlider->setTickPosition(QSlider::TicksBelow);
        labelAlphaRank = new QLabel(centralWidget);
        labelAlphaRank->setObjectName(QString::fromUtf8("labelAlphaRank"));
        labelAlphaRank->setGeometry(QRect(10, 10, 75, 20));
        labelPersistence = new QLabel(centralWidget);
        labelPersistence->setObjectName(QString::fromUtf8("labelPersistence"));
        labelPersistence->setGeometry(QRect(10, 40, 75, 20));
        persistenceSlider = new QSlider(centralWidget);
        persistenceSlider->setObjectName(QString::fromUtf8("persistenceSlider"));
        persistenceSlider->setGeometry(QRect(90, 40, 660, 16));
        persistenceSlider->setOrientation(Qt::Horizontal);
        persistenceSlider->setTickPosition(QSlider::TicksBelow);
        lineVOne = new QFrame(centralWidget);
        lineVOne->setObjectName(QString::fromUtf8("lineVOne"));
        lineVOne->setGeometry(QRect(920, 0, 3, 911));
        lineVOne->setMinimumSize(QSize(3, 700));
        lineVOne->setMaximumSize(QSize(3, 1000));
        lineVOne->setFrameShape(QFrame::VLine);
        lineVOne->setFrameShadow(QFrame::Sunken);
        lineHOne = new QFrame(centralWidget);
        lineHOne->setObjectName(QString::fromUtf8("lineHOne"));
        lineHOne->setGeometry(QRect(920, 130, 350, 3));
        lineHOne->setMinimumSize(QSize(350, 3));
        lineHOne->setMaximumSize(QSize(350, 3));
        lineHOne->setFrameShape(QFrame::HLine);
        lineHOne->setFrameShadow(QFrame::Sunken);
        lineHTwo = new QFrame(centralWidget);
        lineHTwo->setObjectName(QString::fromUtf8("lineHTwo"));
        lineHTwo->setGeometry(QRect(920, 420, 350, 3));
        lineHTwo->setMinimumSize(QSize(350, 3));
        lineHTwo->setMaximumSize(QSize(350, 3));
        lineHTwo->setFrameShape(QFrame::HLine);
        lineHTwo->setFrameShadow(QFrame::Sunken);
        lineEditAlphaRank = new QLineEdit(centralWidget);
        lineEditAlphaRank->setObjectName(QString::fromUtf8("lineEditAlphaRank"));
        lineEditAlphaRank->setGeometry(QRect(760, 10, 150, 25));
        lineEditPersistenceRank = new QLineEdit(centralWidget);
        lineEditPersistenceRank->setObjectName(QString::fromUtf8("lineEditPersistenceRank"));
        lineEditPersistenceRank->setGeometry(QRect(760, 50, 150, 25));
        groupBoxAlphaShape = new QGroupBox(centralWidget);
        groupBoxAlphaShape->setObjectName(QString::fromUtf8("groupBoxAlphaShape"));
        groupBoxAlphaShape->setGeometry(QRect(10, 80, 221, 70));
        checkBoxAlphaShape = new QCheckBox(groupBoxAlphaShape);
        checkBoxAlphaShape->setObjectName(QString::fromUtf8("checkBoxAlphaShape"));
        checkBoxAlphaShape->setGeometry(QRect(0, 20, 111, 24));
        checkBoxWireFrame = new QCheckBox(groupBoxAlphaShape);
        checkBoxWireFrame->setObjectName(QString::fromUtf8("checkBoxWireFrame"));
        checkBoxWireFrame->setGeometry(QRect(0, 40, 95, 24));
        labelAlphaValue = new QLabel(groupBoxAlphaShape);
        labelAlphaValue->setObjectName(QString::fromUtf8("labelAlphaValue"));
        labelAlphaValue->setGeometry(QRect(110, 20, 90, 20));
        lineEditAlphaValue = new QLineEdit(groupBoxAlphaShape);
        lineEditAlphaValue->setObjectName(QString::fromUtf8("lineEditAlphaValue"));
        lineEditAlphaValue->setGeometry(QRect(110, 40, 101, 25));
        groupBoxPockets = new QGroupBox(centralWidget);
        groupBoxPockets->setObjectName(QString::fromUtf8("groupBoxPockets"));
        groupBoxPockets->setGeometry(QRect(240, 70, 291, 81));
        checkBoxAllPockets = new QCheckBox(groupBoxPockets);
        checkBoxAllPockets->setObjectName(QString::fromUtf8("checkBoxAllPockets"));
        checkBoxAllPockets->setGeometry(QRect(0, 20, 101, 24));
        checkBoxOnlyPockets = new QCheckBox(groupBoxPockets);
        checkBoxOnlyPockets->setObjectName(QString::fromUtf8("checkBoxOnlyPockets"));
        checkBoxOnlyPockets->setGeometry(QRect(150, 20, 111, 24));
        checkBoxOnlyVoids = new QCheckBox(groupBoxPockets);
        checkBoxOnlyVoids->setObjectName(QString::fromUtf8("checkBoxOnlyVoids"));
        checkBoxOnlyVoids->setGeometry(QRect(150, 40, 105, 24));
        checkBoxAllMouths = new QCheckBox(groupBoxPockets);
        checkBoxAllMouths->setObjectName(QString::fromUtf8("checkBoxAllMouths"));
        checkBoxAllMouths->setGeometry(QRect(0, 40, 90, 24));
        checkBoxPocketWireFrame = new QCheckBox(groupBoxPockets);
        checkBoxPocketWireFrame->setObjectName(QString::fromUtf8("checkBoxPocketWireFrame"));
        checkBoxPocketWireFrame->setGeometry(QRect(0, 60, 141, 24));
        checkBoxVolumeEnabled = new QCheckBox(groupBoxPockets);
        checkBoxVolumeEnabled->setObjectName(QString::fromUtf8("checkBoxVolumeEnabled"));
        checkBoxVolumeEnabled->setGeometry(QRect(150, 60, 141, 24));
        groupBoxSkinSurface = new QGroupBox(centralWidget);
        groupBoxSkinSurface->setObjectName(QString::fromUtf8("groupBoxSkinSurface"));
        groupBoxSkinSurface->setGeometry(QRect(730, 90, 181, 151));
        checkBoxPocketSkinSolid = new QCheckBox(groupBoxSkinSurface);
        checkBoxPocketSkinSolid->setObjectName(QString::fromUtf8("checkBoxPocketSkinSolid"));
        checkBoxPocketSkinSolid->setGeometry(QRect(10, 40, 171, 24));
        radioButtonFlat = new QRadioButton(groupBoxSkinSurface);
        radioButtonFlat->setObjectName(QString::fromUtf8("radioButtonFlat"));
        radioButtonFlat->setGeometry(QRect(10, 120, 171, 24));
        radioButtonSmooth = new QRadioButton(groupBoxSkinSurface);
        radioButtonSmooth->setObjectName(QString::fromUtf8("radioButtonSmooth"));
        radioButtonSmooth->setGeometry(QRect(10, 100, 171, 24));
        checkBoxAlphaSkinSolid = new QCheckBox(groupBoxSkinSurface);
        checkBoxAlphaSkinSolid->setObjectName(QString::fromUtf8("checkBoxAlphaSkinSolid"));
        checkBoxAlphaSkinSolid->setGeometry(QRect(10, 20, 171, 24));
        checkBoxAlphaSkinWireFrame = new QCheckBox(groupBoxSkinSurface);
        checkBoxAlphaSkinWireFrame->setObjectName(QString::fromUtf8("checkBoxAlphaSkinWireFrame"));
        checkBoxAlphaSkinWireFrame->setGeometry(QRect(10, 60, 171, 24));
        checkBoxPocketSkinWireFrame = new QCheckBox(groupBoxSkinSurface);
        checkBoxPocketSkinWireFrame->setObjectName(QString::fromUtf8("checkBoxPocketSkinWireFrame"));
        checkBoxPocketSkinWireFrame->setGeometry(QRect(10, 80, 171, 24));
        groupBoxMoleclueProperties = new QGroupBox(centralWidget);
        groupBoxMoleclueProperties->setObjectName(QString::fromUtf8("groupBoxMoleclueProperties"));
        groupBoxMoleclueProperties->setGeometry(QRect(930, 10, 330, 125));
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
        groupBoxPocketProperties = new QGroupBox(centralWidget);
        groupBoxPocketProperties->setObjectName(QString::fromUtf8("groupBoxPocketProperties"));
        groupBoxPocketProperties->setGeometry(QRect(930, 140, 340, 271));
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
        tableWidgetPocket->setGeometry(QRect(10, 100, 320, 161));
        tableWidgetPocket->setColumnCount(5);
        groupBoxMouthProperties = new QGroupBox(centralWidget);
        groupBoxMouthProperties->setObjectName(QString::fromUtf8("groupBoxMouthProperties"));
        groupBoxMouthProperties->setGeometry(QRect(930, 420, 340, 240));
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
        tableWidgetMouth->setGeometry(QRect(10, 70, 320, 160));
        tableWidgetMouth->setColumnCount(4);
        labelMouthTotalSurfaceArea->raise();
        lineEditMouthTotalSurfArea->raise();
        tableWidgetMouth->raise();
        tableWidgetPocket->raise();
        FirstViewerPlaceHolder = new QWidget(centralWidget);
        FirstViewerPlaceHolder->setObjectName(QString::fromUtf8("FirstViewerPlaceHolder"));
        FirstViewerPlaceHolder->setGeometry(QRect(10, 240, 900, 671));
        FirstViewerPlaceHolder->setMinimumSize(QSize(440, 440));
        FirstViewerPlaceHolder->setMaximumSize(QSize(1000, 800));
        IndivPocMenuPlaceHolder = new QWidget(centralWidget);
        IndivPocMenuPlaceHolder->setObjectName(QString::fromUtf8("IndivPocMenuPlaceHolder"));
        IndivPocMenuPlaceHolder->setGeometry(QRect(10, 160, 331, 61));
        groupBoxIndivPockets = new QGroupBox(IndivPocMenuPlaceHolder);
        groupBoxIndivPockets->setObjectName(QString::fromUtf8("groupBoxIndivPockets"));
        groupBoxIndivPockets->setGeometry(QRect(0, 10, 321, 51));
        radioButtonAllPockets = new QRadioButton(groupBoxIndivPockets);
        radioButtonAllPockets->setObjectName(QString::fromUtf8("radioButtonAllPockets"));
        radioButtonAllPockets->setGeometry(QRect(10, 20, 103, 24));
        radioButtonIndivPockets = new QRadioButton(groupBoxIndivPockets);
        radioButtonIndivPockets->setObjectName(QString::fromUtf8("radioButtonIndivPockets"));
        radioButtonIndivPockets->setGeometry(QRect(110, 20, 150, 24));
        spinBoxPockets = new QSpinBox(groupBoxIndivPockets);
        spinBoxPockets->setObjectName(QString::fromUtf8("spinBoxPockets"));
        spinBoxPockets->setGeometry(QRect(260, 20, 54, 25));
        groupBoxUserInput = new QGroupBox(centralWidget);
        groupBoxUserInput->setObjectName(QString::fromUtf8("groupBoxUserInput"));
        groupBoxUserInput->setGeometry(QRect(360, 160, 131, 71));
        lineEditEpsilonInput = new QLineEdit(groupBoxUserInput);
        lineEditEpsilonInput->setObjectName(QString::fromUtf8("lineEditEpsilonInput"));
        lineEditEpsilonInput->setGeometry(QRect(0, 30, 131, 25));
        tabWidget = new QTabWidget(centralWidget);
        tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
        tabWidget->setGeometry(QRect(930, 670, 331, 241));
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
        tabWidget->addTab(tab, QString());
        tab_2 = new QWidget();
        tab_2->setObjectName(QString::fromUtf8("tab_2"));
        cHullBox = new QGroupBox(tab_2);
        cHullBox->setObjectName(QString::fromUtf8("cHullBox"));
        cHullBox->setGeometry(QRect(10, 10, 192, 105));
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

        tabWidget->addTab(tab_2, QString());
        tab_3 = new QWidget();
        tab_3->setObjectName(QString::fromUtf8("tab_3"));
        groupBox_2 = new QGroupBox(tab_3);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        groupBox_2->setGeometry(QRect(10, 9, 250, 160));
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

        startCombo = new QComboBox(groupBox_2);
        startCombo->setObjectName(QString::fromUtf8("startCombo"));
        startCombo->setEditable(true);

        horizontalLayout->addWidget(startCombo);


        verticalLayout_3->addLayout(horizontalLayout);

        buttonShortestEscapePath = new QPushButton(groupBox_2);
        buttonShortestEscapePath->setObjectName(QString::fromUtf8("buttonShortestEscapePath"));

        verticalLayout_3->addWidget(buttonShortestEscapePath);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        label_2 = new QLabel(groupBox_2);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_2->addWidget(label_2);

        targetCombo = new QComboBox(groupBox_2);
        targetCombo->setObjectName(QString::fromUtf8("targetCombo"));
        targetCombo->setEditable(true);

        horizontalLayout_2->addWidget(targetCombo);


        verticalLayout_3->addLayout(horizontalLayout_2);

        buttonShortestSTPath = new QPushButton(groupBox_2);
        buttonShortestSTPath->setObjectName(QString::fromUtf8("buttonShortestSTPath"));

        verticalLayout_3->addWidget(buttonShortestSTPath);

        tabWidget->addTab(tab_3, QString());
        tab_4 = new QWidget();
        tab_4->setObjectName(QString::fromUtf8("tab_4"));
        graphWidget = new QCustomPlot(tab_4);
        graphWidget->setObjectName(QString::fromUtf8("graphWidget"));
        graphWidget->setGeometry(QRect(0, 0, 331, 211));
        tabWidget->addTab(tab_4, QString());
        MainWindow->setCentralWidget(centralWidget);
        alphaSlider->raise();
        labelAlphaRank->raise();
        labelPersistence->raise();
        persistenceSlider->raise();
        lineVOne->raise();
        lineHOne->raise();
        lineHTwo->raise();
        lineEditAlphaRank->raise();
        lineEditPersistenceRank->raise();
        groupBoxAlphaShape->raise();
        groupBoxPockets->raise();
        groupBoxSkinSurface->raise();
        groupBoxMoleclueProperties->raise();
        groupBoxPocketProperties->raise();
        groupBoxMouthProperties->raise();
        IndivPocMenuPlaceHolder->raise();
        groupBoxUserInput->raise();
        tabWidget->raise();
        FirstViewerPlaceHolder->raise();
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
        labelAlphaRank->setBuddy(alphaSlider);
        labelPersistence->setBuddy(persistenceSlider);
        labelAlphaValue->setBuddy(lineEditAlphaValue);
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
        menuFile->addAction(actionExit);
        menuFile->addAction(actionSave_Graph);
        menuHelp->addAction(actionAbout);
        menu_Tools->addAction(actionAlphaZero);
        menu_Tools->addAction(actionChangeFiltration);
        menu_Tools->addAction(actionUndoChange);

        retranslateUi(MainWindow);

        tabWidget->setCurrentIndex(3);


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
        labelAlphaRank->setText(QApplication::translate("MainWindow", "Alpha Rank", 0, QApplication::UnicodeUTF8));
        labelPersistence->setText(QApplication::translate("MainWindow", "Persistence", 0, QApplication::UnicodeUTF8));
        groupBoxAlphaShape->setTitle(QApplication::translate("MainWindow", "AlphaShape", 0, QApplication::UnicodeUTF8));
        checkBoxAlphaShape->setText(QApplication::translate("MainWindow", "AlphaShape", 0, QApplication::UnicodeUTF8));
        checkBoxWireFrame->setText(QApplication::translate("MainWindow", "WireFrame", 0, QApplication::UnicodeUTF8));
        labelAlphaValue->setText(QApplication::translate("MainWindow", "Alpha Value :", 0, QApplication::UnicodeUTF8));
        groupBoxPockets->setTitle(QApplication::translate("MainWindow", "Pockets", 0, QApplication::UnicodeUTF8));
        checkBoxAllPockets->setText(QApplication::translate("MainWindow", "All Pockets", 0, QApplication::UnicodeUTF8));
        checkBoxOnlyPockets->setText(QApplication::translate("MainWindow", "Only Pockets", 0, QApplication::UnicodeUTF8));
        checkBoxOnlyVoids->setText(QApplication::translate("MainWindow", "Only Voids", 0, QApplication::UnicodeUTF8));
        checkBoxAllMouths->setText(QApplication::translate("MainWindow", "Mouths", 0, QApplication::UnicodeUTF8));
        checkBoxPocketWireFrame->setText(QApplication::translate("MainWindow", "PocketWireFrame", 0, QApplication::UnicodeUTF8));
        checkBoxVolumeEnabled->setText(QApplication::translate("MainWindow", "Display Properties", 0, QApplication::UnicodeUTF8));
        groupBoxSkinSurface->setTitle(QApplication::translate("MainWindow", "SkinSurface", 0, QApplication::UnicodeUTF8));
        checkBoxPocketSkinSolid->setText(QApplication::translate("MainWindow", "PocketSkinSurface", 0, QApplication::UnicodeUTF8));
        radioButtonFlat->setText(QApplication::translate("MainWindow", "Flat Shading", 0, QApplication::UnicodeUTF8));
        radioButtonSmooth->setText(QApplication::translate("MainWindow", "Smooth Shading", 0, QApplication::UnicodeUTF8));
        checkBoxAlphaSkinSolid->setText(QApplication::translate("MainWindow", "AlphaSkinSurface", 0, QApplication::UnicodeUTF8));
        checkBoxAlphaSkinWireFrame->setText(QApplication::translate("MainWindow", "AlphaSkinWireFrame", 0, QApplication::UnicodeUTF8));
        checkBoxPocketSkinWireFrame->setText(QApplication::translate("MainWindow", "PocketSkinWireFrame", 0, QApplication::UnicodeUTF8));
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
        groupBoxIndivPockets->setTitle(QApplication::translate("MainWindow", "Individual Pockets", 0, QApplication::UnicodeUTF8));
        radioButtonAllPockets->setText(QApplication::translate("MainWindow", "All Pockets", 0, QApplication::UnicodeUTF8));
        radioButtonIndivPockets->setText(QApplication::translate("MainWindow", "Individual Pockets", 0, QApplication::UnicodeUTF8));
        groupBoxUserInput->setTitle(QApplication::translate("MainWindow", "Epsilon Input", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QApplication::translate("MainWindow", "PowerDiagram", 0, QApplication::UnicodeUTF8));
        checkPowerDiag->setText(QApplication::translate("MainWindow", "Show Power Diagram", 0, QApplication::UnicodeUTF8));
        checkComplementSpacePD->setText(QApplication::translate("MainWindow", "Complement Space PD", 0, QApplication::UnicodeUTF8));
        checkInsideVerts->setText(QApplication::translate("MainWindow", "Only Inside Vertices", 0, QApplication::UnicodeUTF8));
        checkPruneIsolatedVerts->setText(QApplication::translate("MainWindow", "Prune Isolated Vertices", 0, QApplication::UnicodeUTF8));
        checkEdgeIntersect->setText(QApplication::translate("MainWindow", "Intersect Edges", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab), QApplication::translate("MainWindow", "Power Diagram", 0, QApplication::UnicodeUTF8));
        cHullBox->setTitle(QApplication::translate("MainWindow", "Convex Hull", 0, QApplication::UnicodeUTF8));
        checkCHull->setText(QApplication::translate("MainWindow", "Show Convex Hull", 0, QApplication::UnicodeUTF8));
        checkCHullWireFrame->setText(QApplication::translate("MainWindow", "Convex Hull Wireframe", 0, QApplication::UnicodeUTF8));
        checkCHullNorm->setText(QApplication::translate("MainWindow", "Show Hull Normals", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab_2), QApplication::translate("MainWindow", "Convex Hull", 0, QApplication::UnicodeUTF8));
        groupBox_2->setTitle(QApplication::translate("MainWindow", "Shortest Escape Route", 0, QApplication::UnicodeUTF8));
        startLabel->setText(QApplication::translate("MainWindow", "Start Node:", 0, QApplication::UnicodeUTF8));
        buttonShortestEscapePath->setText(QApplication::translate("MainWindow", "Find Shortest Escape Path", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("MainWindow", "Target Node:", 0, QApplication::UnicodeUTF8));
        buttonShortestSTPath->setText(QApplication::translate("MainWindow", "Find Shortest Source-Target path", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab_3), QApplication::translate("MainWindow", "Shortest Path", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab_4), QApplication::translate("MainWindow", "Path Profile", 0, QApplication::UnicodeUTF8));
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
