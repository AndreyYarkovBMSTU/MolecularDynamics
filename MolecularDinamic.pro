#-------------------------------------------------
#
# Project created by QtCreator 2018-09-27T14:28:45
#
#-------------------------------------------------

QT       += core gui widgets

TARGET = MolecularDinamic
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++14
SOURCES += \
        main.cpp \
        mainwindow.cpp \
    thermostat/thermostat.cpp \
    particle/particle.cpp \
    particle/particlesystem.cpp \
    particle/state.cpp \
    numeralequations/numeralequations.cpp \
    moleculardinamic/moleculardinamic.cpp \
    moleculardinamic/potential.cpp \
    material/material.cpp \
    environment/environment.cpp \
    electrostatics/externalfields.cpp \
    electrostatics/interaction.cpp \
    electrostatics/methods.cpp \
    namespase/mathematics.cpp \
    namespase/phys.cpp

HEADERS += \
        mainwindow.h \
    electrostatics/externalfields.h \
    electrostatics/interaction.h \
    electrostatics/methods.h \
    environment/environment.h \
    material/material.h \
    moleculardinamic/moleculardinamic.h \
    moleculardinamic/potential.h \
    namespase/mathematics.h \
    namespase/phys.h \
    numeralequations/numeralequations.h \
    particle/particle.h \
    particle/particlesystem.h \
    particle/state.h \
    thermostat/thermostat.h \
    header.h \
    properties.h \
    object.h

FORMS += \
        mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
