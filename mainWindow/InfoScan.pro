QT +=       core gui \
            printsupport \
            xml\
            widgets\
            concurrent widgets

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

QT += webenginewidgets

requires(qtConfig(listwidget))

qtHaveModule(printsupport): QT += printsupport


# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

include(FlatTabWidget/FlatTabWidget.pri);


SOURCES += \
    configuration.cpp \
    contact.cpp \
    dataanalysis.cpp \
    dataexample.cpp \
    fastp2.cpp \
    fastpexample.cpp \
    fastpresult.cpp \
    functionanalysis.cpp \
    functionresult.cpp \
    hisat.cpp \
    hisatexample.cpp \
    hisatresult.cpp \
    lncexample.cpp \
    main.cpp \
    mainpanel.cpp \
    mainwindow.cpp \
    sclncrnafun.cpp \
    spatialexample.cpp \
    spatialscan.cpp \
    tutroial.cpp

HEADERS += \
    configuration.h \
    contact.h \
    dataanalysis.h \
    dataexample.h \
    fastp2.h \
    fastpexample.h \
    fastpresult.h \
    functionanalysis.h \
    functionresult.h \
    hisat.h \
    hisatexample.h \
    hisatresult.h \
    lncexample.h \
    mainpanel.h \
    mainwindow.h \
    sclncrnafun.h \
    spatialexample.h \
    spatialscan.h \
    tutroial.h

FORMS += \
    configuration.ui \
    contact.ui \
    dataanalysis.ui \
    dataexample.ui \
    fastp2.ui \
    fastpexample.ui \
    fastpresult.ui \
    fileview.ui \
    functionanalysis.ui \
    functionresult.ui \
    hisat.ui \
    hisatexample.ui \
    hisatresult.ui \
    lncexample.ui \
    mainpanel.ui \
    mainwindow.ui \
    sclncrnafun.ui \
    spatialexample.ui \
    spatialscan.ui \
    tutroial.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

DISTFILES += \
    ../build-InfoScan-Desktop-Debug/qss/AMOLED.qss \
    ../build-InfoScan-Desktop-Debug/qss/Adaptic.qss \
    ../build-InfoScan-Desktop-Debug/qss/Aqua.qss \
    ../build-InfoScan-Desktop-Debug/qss/Chatbee.qss \
    ../build-InfoScan-Desktop-Debug/qss/Combinear.qss \
    ../build-InfoScan-Desktop-Debug/qss/DeepBox.qss \
    ../build-InfoScan-Desktop-Debug/qss/Diffnes.qss \
    ../build-InfoScan-Desktop-Debug/qss/Diplaytap.qss \
    ../build-InfoScan-Desktop-Debug/qss/Dtor.qss \
    ../build-InfoScan-Desktop-Debug/qss/EasyCode.qss \
    ../build-InfoScan-Desktop-Debug/qss/Eclippy.qss \
    ../build-InfoScan-Desktop-Debug/qss/Fibers.qss \
    ../build-InfoScan-Desktop-Debug/qss/Fibrary.qss \
    ../build-InfoScan-Desktop-Debug/qss/Filmovio.qss \
    ../build-InfoScan-Desktop-Debug/qss/Genetive.qss \
    ../build-InfoScan-Desktop-Debug/qss/Geoo.qss \
    ../build-InfoScan-Desktop-Debug/qss/Gravira.qss \
    ../build-InfoScan-Desktop-Debug/qss/HackBook.qss \
    ../build-InfoScan-Desktop-Debug/qss/Hookmark.qss \
    ../build-InfoScan-Desktop-Debug/qss/Incrypt.qss \
    ../build-InfoScan-Desktop-Debug/qss/Integrid.qss \
    ../build-InfoScan-Desktop-Debug/qss/Irrorater.qss \
    ../build-InfoScan-Desktop-Debug/qss/MailSy.qss \
    ../build-InfoScan-Desktop-Debug/qss/ManjaroMix.qss \
    ../build-InfoScan-Desktop-Debug/qss/MaterialDark.qss \
    ../build-InfoScan-Desktop-Debug/qss/Medize.qss \
    ../build-InfoScan-Desktop-Debug/qss/MyWatch.qss \
    ../build-InfoScan-Desktop-Debug/qss/NeonButtons.qss \
    ../build-InfoScan-Desktop-Debug/qss/Obit.qss \
    ../build-InfoScan-Desktop-Debug/qss/Perstfic.qss \
    ../build-InfoScan-Desktop-Debug/qss/Photoxo.qss \
    ../build-InfoScan-Desktop-Debug/qss/PicPax.qss \
    ../build-InfoScan-Desktop-Debug/qss/Remover.qss \
    ../build-InfoScan-Desktop-Debug/qss/Scalcula.qss \
    ../build-InfoScan-Desktop-Debug/qss/SpyBot.qss \
    ../build-InfoScan-Desktop-Debug/qss/SyNet.qss \
    ../build-InfoScan-Desktop-Debug/qss/TCobra.qss \
    ../build-InfoScan-Desktop-Debug/qss/Takezo.qss \
    ../build-InfoScan-Desktop-Debug/qss/Ubuntu.qss \
    ../build-InfoScan-Desktop-Debug/qss/VisualScript.qss \
    ../build-InfoScan-Desktop-Debug/qss/Webmas.qss \
    ../build-InfoScan-Desktop-Debug/qss/WordOffice.qss \
    ../build-InfoScan-Desktop-Debug/qss/Wstartpage.qss \
    FlatTabWidget/3rdparty/LICENSE \
    FlatTabWidget/3rdparty/README.md \
    FlatTabWidget/FlatTabWidget.pri \
    FlatTabWidget/LICENSE \
    FlatTabWidget/README.md
RESOURCES += \
    images.qrc \
    qdarkstyle/dark/darkstyle.qrc \
    qdarkstyle/light/lightstyle.qrc

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../miniconda3/envs/QT/lib/release/ -licuuc
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../miniconda3/envs/QT/lib/debug/ -licuuc
else:unix: LIBS += -L$$PWD/../../../miniconda3/envs/QT/lib/ -licuuc

INCLUDEPATH += $$PWD/../../../miniconda3/envs/QT/include
DEPENDPATH += $$PWD/../../../miniconda3/envs/QT/include
