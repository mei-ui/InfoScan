#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>
#include <QFile>
#include <QToolBar>
#include <QDesktopServices>
#include <QUrl>
#include "configuration.h"
#include <QWebEngineView>
#include <QWebEngineSettings>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    createActions();
    createStatusBar();
    setWindowTitle(tr("InfoScan"));
    setWindowIcon(QIcon("./images/logo2.jpg"));
    setUnifiedTitleAndToolBarOnMac(true);
 //   QFile styleFile( "./qss/Ubuntu.qss" );
  //  styleFile.open( QFile::ReadOnly );
 //   QString style(styleFile.readAll());
 //   qApp->setStyleSheet(style);

    QFile f(":qdarkstyle/light/lightstyle.qss");

    if (!f.exists())   {
        printf("Unable to set stylesheet, file not found\n");
    }
    else   {
        f.open(QFile::ReadOnly | QFile::Text);
        QTextStream ts(&f);
        qApp->setStyleSheet(ts.readAll());
    }
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::appsetStyleSheet1()
{
    QFile styleFile( "qss/Fibrary.qss" );
    styleFile.open( QFile::ReadOnly );
    QString style(styleFile.readAll());
    qApp->setStyleSheet(style);

}
void MainWindow::appsetStyleSheet2()
{
    QFile styleFile( "qss/Scalcula.qss" );
    styleFile.open( QFile::ReadOnly );
    QString style(styleFile.readAll());
    qApp->setStyleSheet(style);

}
void MainWindow::appsetStyleSheet3()
{
    QFile styleFile( "qss/Medize.qss" );
    styleFile.open( QFile::ReadOnly );
    QString style(styleFile.readAll());
    qApp->setStyleSheet(style);

}

void MainWindow::about()
{
   QMessageBox::about(this, tr("About InfoScan"),
            tr("The <b>InfoScan</b> program performs genome-wise identification of lncRNA."
               "Using <b>InfoScan</b>'s program, you can"
               "choose suitable parameters to run and "
               "visualize, and explore single cell data."));
}
void MainWindow::sraexplorer()
{
    QString link = "https://sra-explorer.info/";
    QWebEngineView *LiveView = new QWebEngineView();
    LiveView->settings()->setAttribute(QWebEngineSettings::PluginsEnabled, true);
    LiveView->setAttribute(Qt::WA_DeleteOnClose);
    LiveView->load(QUrl(link));
    LiveView->resize(1024, 768);
    LiveView->show();
 //   QDesktopServices::openUrl(QUrl(link));
}

void MainWindow::configuration()
{
    class Configuration configuration;
    configuration.setWindowModality(Qt::WindowModal);
    configuration.setWindowIcon(QIcon(":/images/logo.png"));
    configuration.setWindowTitle("Softwares checking and configuration");
    configuration.exec();
}
void MainWindow::createStatusBar()
{
    statusBar()->showMessage(tr("Ready"));
}

void MainWindow::createActions()
{

//Tools:
//    addToolBarBreak()
    QToolBar *ToolsBar = addToolBar(tr("Tools"));
    ToolsBar->setFixedHeight(100);
    ToolsBar->setIconSize(QSize(30, 30));
    ToolsBar->setToolButtonStyle( Qt::ToolButtonTextUnderIcon );




//View:

//    viewMenu = menuBar()->addMenu(tr("&View"));

//    menuBar()->addSeparator();


//Global:
    addToolBarBreak();
    QMenu *globalMenu = menuBar()->addMenu(tr("&Global"));


//    globalMenu->addMenu(tr("themes"));
    QMenu* themes = globalMenu->addMenu(tr("&Themes"));

    QAction *configurationAct = globalMenu->addAction(tr("&Configuration"), this, &MainWindow::configuration);
    configurationAct->setStatusTip(tr("Show the application's About box"));

    const QIcon themesIcon = QIcon::fromTheme("document-open", QIcon(":/images/themelogo.png"));
    QAction *themesAct1 = new QAction(themesIcon, tr("&Combinear"), this);
    themesAct1->setStatusTip(tr("Change Combinear theme"));
    connect(themesAct1, &QAction::triggered, this, &MainWindow::appsetStyleSheet1);
    themes->addAction(themesAct1);

    QAction *themesAct2 = new QAction(themesIcon, tr("&WordOffice"), this);
    themesAct2->setStatusTip(tr("Change WordOffice theme"));
    connect(themesAct2, &QAction::triggered, this, &MainWindow::appsetStyleSheet2);
    themes->addAction(themesAct2);

    QAction *themesAct3 = new QAction(themesIcon, tr("&NeonButtons"), this);
    themesAct3->setStatusTip(tr("Change NeonButtons theme"));
    connect(themesAct3, &QAction::triggered, this, &MainWindow::appsetStyleSheet3);
    themes->addAction(themesAct3);



//Help:

    QMenu *helpMenu = menuBar()->addMenu(tr("&Help"));

    QAction *aboutAct = helpMenu->addAction(tr("&About"), this, &MainWindow::about);
    aboutAct->setStatusTip(tr("Show the application's About box"));

    QAction *aboutQtAct = helpMenu->addAction(tr("About &Qt"), qApp, &QApplication::aboutQt);
    aboutQtAct->setStatusTip(tr("Show the PseudoToolBox library's About box"));

    QAction *sraexplorertAct = helpMenu->addAction(tr("Link to &sraexplorer"), this, &MainWindow::sraexplorer);
    sraexplorertAct->setStatusTip(tr("Website link of SRA-Explorer"));
}
