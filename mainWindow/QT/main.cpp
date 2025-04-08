#include "mainpanel.h"
#include "ui_mainpanel.h"
#include "fastp2.h"
#include "ui_fastp2.h"
#include "sclncrnafun.h"
#include "ui_sclncrnafun.h"
#include "dataanalysis.h"
#include "ui_dataanalysis.h"
#include "hisat.h"
#include "ui_hisat.h"
#include "tutroial.h"
#include "ui_tutroial.h"
#include "contact.h"
#include "ui_contact.h"
#include "functionanalysis.h"
#include "ui_functionanalysis.h"
#include "spatialscan.h"
#include "ui_spatialscan.h"

#include <flattabwidget.h>
#include "mainwindow.h"
#include <QApplication>
#include <QCommandLineParser>
#include <QCommandLineOption>
#include <QSplashScreen>
#include <QDesktopWidget>
#include <QFileDialog>
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    QCoreApplication::setOrganizationName("Sun Yat-Sen University");
    QCoreApplication::setApplicationName("InfoScan");
    QCoreApplication::setApplicationVersion(QT_VERSION_STR);
    QCommandLineParser parser;
    parser.setApplicationDescription(QCoreApplication::applicationName());
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addPositionalArgument("file", "The file to open.");
    parser.process(a);
//    Q_INIT_RESOURCE(InfoScan);

    QSplashScreen *splash = new QSplashScreen;
    splash->setPixmap(QPixmap("images/logo2.jpg"));
    splash->show();

    //让对话框延迟一段时间显示
    for(int i=0;i<5000;i++)
    {
        splash->repaint();
    }


    FlatTabWidget* w = new FlatTabWidget();
    w->setContentsMargins(9,9,9,9);
    w->addPage("Introduction",new class Mainpanel);
    w->addPage("InfoUpload",new class Fastp2);
    w->addPage("InfoAssembly",new class Hisat);
    w->addPage("NovelScan",new class SClncRNAFun);
    w->addPage("CellInfo",new class DataAnalysis);
    w->addPage("FuncScan",new class FunctionAnalysis);
    w->addPage("SpatialScan",new class SpatialScan);
//    w->addPage("Tutorial",new class Tutroial);
//    w->addPage("Contact",new class Contact);

    QGridLayout  *layout = new QGridLayout();
    layout->addWidget(w);
    w->setLayout(layout);

    MainWindow mainWin;
    mainWin.setWindowState(Qt::WindowMaximized);
    mainWin.resize(1000,800);
    mainWin.setCentralWidget(w);
    mainWin.show();
    splash->finish(&mainWin);
    delete splash;

    return a.exec();
}
