#include "mainpanel.h"
#include "ui_mainpanel.h"
#include <QDebug>
#include <QProcess>
#include <QDesktopServices>
#include <QWebEngineSettings>
#include <QUrl>

Mainpanel::Mainpanel(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Mainpanel)
{
    ui->setupUi(this);
}

Mainpanel::~Mainpanel()
{
    delete ui;
}

void Mainpanel::on_check_clicked()
{
    configuration = new Configuration(this);
    configuration->show();
}

void Mainpanel::on_check_2_clicked()
{
    QWebEngineView *LiveView = new QWebEngineView();
    LiveView->settings()->setAttribute(QWebEngineSettings::PluginsEnabled, true);
    LiveView->setAttribute(Qt::WA_DeleteOnClose);
 //  QString start_path = QCoreApplication::applicationDirPath()+"/html/index.html";
  //  qDebug() << start_path << endl;
    LiveView->load(QUrl("https://infoscan-docs.readthedocs.io/en/latest/index.html"));
    LiveView->resize(1024, 768);
    LiveView->show();
}
