#include "fastp2.h"
#include "ui_fastp2.h"
#include "fastpexample.h"
#include "ui_fastpexample.h"
#include <QDebug>
#include <QProcess>
#include <QDesktopServices>
#include <QWebEngineSettings>
#include <QUrl>
#include "fastp2.h"

Fastpexample::Fastpexample(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Fastpexample)
{
    ui->setupUi(this);
    connect(ui->lineEdit,SIGNAL(returnPressed()),this,SLOT(showCurrentDirFiles()));
    connect(ui->listWidget,SIGNAL(itemDoubleClicked(QListWidgetItem*)),this,SLOT(showNextDirFiles(QListWidgetItem*)));

}

Fastpexample::~Fastpexample()
{
    delete ui;
}
void Fastpexample::showNextDirFiles(QListWidgetItem *item)
{
    ///获取鼠标双击的文件名字
    QString strName = item->text();
    QDir dir;
    //设置路径为当前目录路径
    dir.setPath(ui->lineEdit->text());
    if(!dir.cd(strName))         //如果进入失败，则为普通文件，创建新的进程打开对应文件
    {
        QStringList arguments;
        QString start_path = QDir::currentPath()+"/"+ui->lineEdit->text()+"/" + strName;
        arguments << start_path;
        qDebug()<<"arguments"<<arguments;
//        QProcess* process = new QProcess;

        // qDebug() << start_path << Qt::endl;
//        QDesktopServices :: openUrl(QUrl(start_path));

        QWebEngineView *LiveView = new QWebEngineView();
        LiveView->settings()->setAttribute(QWebEngineSettings::PluginsEnabled, true);
        LiveView->setAttribute(Qt::WA_DeleteOnClose);
        LiveView->load(QUrl("file://"+start_path));
        LiveView->resize(1024, 768);
        LiveView->show();

//       process->start("vim",arguments);  //开启新进程打开文件

        return;}

    //重新设置路径
    dir.cd(strName);
    //更新当前显示路径,并显示当前目录下所有文件
    ui->lineEdit->setText(dir.absolutePath());
    showCurrentDirFiles();

}

void Fastpexample::showCurrentDirFiles()
{
    //获取当前输入的目录
    QDir currentDir(ui->lineEdit->text());
    qDebug()<<"currentdir"<<currentDir;
    QStringList fileList;
    fileList<<"*.html";
    QFileInfoList infoList = currentDir.entryInfoList(fileList,QDir::AllEntries,QDir::DirsFirst);
    //在QListWidget里显示文件列表
    showFileInfoList(infoList);
}
void Fastpexample::showFileInfoList(QFileInfoList pInfoList)
{
    ui->listWidget->clear();
    for(int i=0;i<pInfoList.size();i++)
    {
        QFileInfo tmpInfo = pInfoList.at(i);
        QString pFileName = tmpInfo.fileName();
        QListWidgetItem *tmpItem = new QListWidgetItem(pFileName);
        if(tmpInfo.isDir())
            tmpItem->setIcon(*getItemPropertyIcon(1));
        else
            tmpItem->setIcon(*getItemPropertyIcon(2));
        ui->listWidget->addItem(tmpItem);
    }
}
QIcon * Fastpexample::getItemPropertyIcon(int iType)
{
    QDir dir;
    QString path = dir.currentPath();
    qDebug() <<"path"<< path ;
    switch(iType)
    {
    case 1:
        return new QIcon(path+"/images/open.png");
        break;
    case 2:
        return new QIcon(path+"/images/new.png");
        break;
    }
    return NULL;
}


