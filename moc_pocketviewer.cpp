/****************************************************************************
** Meta object code from reading C++ file 'pocketviewer.h'
**
** Created: Tue Nov 6 17:05:05 2012
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "pocketviewer.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'pocketviewer.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_PocketViewer[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

       0        // eod
};

static const char qt_meta_stringdata_PocketViewer[] = {
    "PocketViewer\0"
};

void PocketViewer::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    Q_UNUSED(_o);
    Q_UNUSED(_id);
    Q_UNUSED(_c);
    Q_UNUSED(_a);
}

const QMetaObjectExtraData PocketViewer::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject PocketViewer::staticMetaObject = {
    { &QGLViewer::staticMetaObject, qt_meta_stringdata_PocketViewer,
      qt_meta_data_PocketViewer, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &PocketViewer::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *PocketViewer::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *PocketViewer::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_PocketViewer))
        return static_cast<void*>(const_cast< PocketViewer*>(this));
    if (!strcmp(_clname, "Processor"))
        return static_cast< Processor*>(const_cast< PocketViewer*>(this));
    return QGLViewer::qt_metacast(_clname);
}

int PocketViewer::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGLViewer::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    return _id;
}
QT_END_MOC_NAMESPACE
