# -------------------------------------------------
# Project created by QtCreator 2010-04-27T02:30:29
# -------------------------------------------------
QT += opengl \
    xml
INCLUDEPATH *= /usr/include/QGLViewer
LIBS *= -L/usr/lib \
    -lqglviewer-qt4 \
    -lgmp
TARGET = PocketViewer
TEMPLATE = app
SOURCES += main.cpp \
    pocketviewer.cpp \
    mainwindow.cpp \
    processor.cpp \
    vertex.cpp \
    vector3.cpp \
    edge.cpp \
    triangle.cpp \
    tetrahedron.cpp \
    sos.cpp \
    alfmasternode.cpp \
    linktrig.cpp \
    tlistnode.cpp \
    node.cpp \
    deluanaycomplex.cpp \
    alphacomplex.cpp \
    mheapnode.cpp \
    size.cpp \
    disjointset.cpp \
    depth.cpp \
    deluanayflow.cpp \
    lightmaterial.cpp \
    flow.cpp \
    filereader.cpp \
    pocket.cpp \
    mouth.cpp \
    spacefillmeasure.cpp \
    volume.cpp \
    trigvol.cpp \
    skinsurface.cpp \
    rankmap.cpp \
    simplexmasterlistmap.cpp \
    metric.cpp \
    powerdiagram.cpp \
    graph.cpp \
    qcustomplot.cpp
HEADERS += pocketviewer.h \
    mainwindow.h \
    processor.h \
    vertex.h \
    vector3.h \
    edge.h \
    triangle.h \
    tetrahedron.h \
    sos.h \
    alfmasternode.h \
    linktrig.h \
    tlistnode.h \
    node.h \
    deluanaycomplex.h \
    alphacomplex.h \
    mheapnode.h \
    size.h \
    disjointset.h \
    depth.h \
    deluanayflow.h \
    lightmaterial.h \
    flow.h \
    filereader.h \
    pocket.h \
    mouth.h \
    spacefillmeasure.h \
    volume.h \
    trigvol.h \
    skinsurface.h \
    rankmap.h \
    simplexmasterlistmap.h \
    metric.h \
    powerdiagram.h \
    graph.h \
    qcustomplot.h
FORMS += mainwindow.ui
