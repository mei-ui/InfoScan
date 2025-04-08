#include "configuration.h"
#include "ui_configuration.h"
#include <QProcess>
#include <QDebug>
#include <QFile>
#include <QFileInfo>
#include <QProgressDialog>
#include <QtConcurrent>
#include <QMessageBox>
using namespace QtConcurrent;

Configuration::Configuration(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Configuration)
{
    ui->setupUi(this);
}

Configuration::~Configuration()
{
    delete ui;
}
void Configuration::on_pushButton_clicked()
{

    QStringList config_command;
    config_command<<"0_software_test.sh"
                 <<"-y"
                <<"snakemake/InfoScan.yaml";
    QString config_command_join = config_command.join(" ");
    config_command_join = QCoreApplication::applicationDirPath()+"/snakemake/script/"+config_command_join;
    qDebug() << "config_command_join:" << config_command_join;

    QEventLoop loop;
    QProgressDialog dialog;
    dialog.setWindowModality(Qt::WindowModal);
    dialog.setWindowTitle("Configuration: checking and installation");
    dialog.setWindowIcon(QIcon(":/images/logo2.jpg"));
    dialog.setMinimumWidth(800);
    dialog.setRange(0,0);
    dialog.setLabelText(QString("Configuration: checking required softwares and tools for InfoScan... (%1 threads)...").arg(QThread::idealThreadCount()));
    dialog.show();


    QProcess p;
    p.start("bash",QStringList()<<"-c"<<config_command_join);
    connect(&p,SIGNAL(finished(int,QProcess::ExitStatus)),&loop,SLOT(quit()));
    loop.exec();

    QString getStr=QString(p.readAllStandardOutput());
    ui->textEdit->setText(getStr);
    QMessageBox::information(this,"Configuration log: \n","All checking and installation are done!");

    QString getStr_Error=QString(p.readAllStandardError());
    if(getStr_Error!=""){
    QMessageBox::information(this,"Configuration Error: \n",getStr_Error);}
}
