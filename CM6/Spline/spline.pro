QT       += core gui printsupport

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    gauss.cpp \
    main.cpp \
    mainwindow.cpp \
    qcustomplot.cpp \
    sndwin.cpp \
    spline.cpp \
    trdwin.cpp

HEADERS += \
    gauss.h \
    head.h \
    mainwindow.h \
    qcustomplot.h \
    sndwin.h \
    spline.h \
    trdwin.h

FORMS += \
    mainwindow.ui \
    sndwin.ui \
    trdwin.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

DISTFILES += \
    nodes.txt \
    points.txt
