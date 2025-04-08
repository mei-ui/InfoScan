#include "spatialscan.h"
#include "ui_spatialscan.h"
#include <QDebug>
#include <QFileDialog>
#include <QMessageBox>
#include <QDesktopWidget>
#include <QProcess>
#include <QProgressDialog>
#include <QtConcurrent>

using namespace QtConcurrent;

SpatialScan::SpatialScan(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SpatialScan)
{
    ui->setupUi(this);
    connect(ui->scRNA_data, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->stRNA_data, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));

    connect(ui->celltype, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->gene, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
}

SpatialScan::~SpatialScan()
{
    delete ui;
}
void SpatialScan::checkLineEdits()
{
    bool ok_1 =!ui->scRNA_data->text().isEmpty()
            && !ui->stRNA_data->text().isEmpty()
            && !ui->celltype->text().isEmpty()
            && !ui->gene->text().isEmpty();
    ui->start->setEnabled(ok_1);

}
QString convertAbsolutePathToRelative(const QString& absolute_path)
{
    QString current_path = QDir::currentPath(); // 当前路径，如 C:\windows的形式
    QString tmp_str = absolute_path;  // 绝对路径，如 C:\windows\system32\cmd.exe的形式
    // if 语句首先判断绝对路径是否是以当前路径开头，以及绝对路径是否大于相对路径，如果满足条件，调用remove函数从绝对路径中剔除当前路径，剩下的即为相对路径，如 `system32\cmd.exe`
    if (absolute_path.startsWith(current_path) && absolute_path.size() > current_path.size())
    {
        tmp_str.remove(0, current_path.size() + 1);
    }
    return tmp_str;  // 返回相对路径
}
void SpatialScan::on_browse_scRNA_clicked()
{


    QFileDialog dialog(this);
    dialog.setAcceptMode(QFileDialog::AcceptSave);
    QString FilePath = "11_data_analysis/celltype";
    FilePath = QCoreApplication::applicationDirPath()+"/snakemake/result/"+FilePath;
    QString meta = dialog.getOpenFileName(this,tr("Locate rds"),FilePath,tr("Find scRNA celltype result File (*.rds)"));

    if (!meta.isEmpty()) {
        ui->scRNA_data->setText(meta);
    }
}
void SpatialScan::on_browse_stRNA_clicked()
{
    QFileDialog dialog(this);
    dialog.setAcceptMode(QFileDialog::AcceptSave);
    QString meta = dialog.getOpenFileName(this,tr("Locate rds"),"rds",tr("Find GeneSetEnrichment result File (*.rds)"));

    if (!meta.isEmpty()) {
        ui->stRNA_data->setText(meta);
    }
}

void SpatialScan::on_start_clicked()
{
    QString scRNA_data=ui->scRNA_data->text();
    QString stRNA_data=ui->stRNA_data->text();
    QString celltype=ui->celltype->text();
    QString gene=ui->gene->text();
    QStringList SpatialScan_config;
    SpatialScan_config<<"SC: \""+scRNA_data+"\""
               <<"ST: \""+stRNA_data+"\""
                 <<"Celltype: \""+celltype+"\""
                   <<"Gene: \""+gene+"\"";
    QString SpatialScan_config_join = SpatialScan_config.join(" \n");
    QString qPath("snakemake/config/SpatialScan_config.yml");
    QFile qFile(qPath);
      if (qFile.open(QIODevice::WriteOnly)) {
        QTextStream out(&qFile); out << SpatialScan_config_join;
        qFile.close();
      }

      QStringList SpatialScan_command;
      SpatialScan_command<<"SpatialScan.sh"
                <<"-I"
                 <<"snakemake/config/SpatialScan_config.yml";


      QString SpatialScan_command_join = SpatialScan_command.join(" \n");
      qDebug()<<"SpatialScan command:"<<SpatialScan_command_join;

//      qDebug() << QCoreApplication::applicationDirPath()+"/snakemake/script/";
      SpatialScan_command_join = QCoreApplication::applicationDirPath()+"/snakemake/script/"+SpatialScan_command_join;
      qDebug() << "SpatialScan_command_join:" << SpatialScan_command_join;

      QEventLoop loop;
      QProgressDialog dialog;
      dialog.setWindowModality(Qt::WindowModal);
      dialog.setWindowTitle("Spatial analysis");
      dialog.setWindowIcon(QIcon(":/images/qt-logo.png"));
      dialog.setMinimumWidth(800);
      dialog.setRange(0,0);
      dialog.setLabelText(QString("Begin Spatial analysis (%1 threads)...").arg(QThread::idealThreadCount()));
      dialog.show();

      QProcess p;
      p.start("bash",QStringList()<<"-c"<<SpatialScan_command_join);
      connect(&p,SIGNAL(finished(int,QProcess::ExitStatus)),&loop,SLOT(quit()));
      loop.exec();
/*      QGraphicsScene* originalScene =new QGraphicsScene(this);
      ui->graphicsView->setScene(originalScene);
      QStringList plot_result;
      plot_result<<lncRNA_id+".png";
      qDebug() << "plot_result:" << plot_result;
      QString  plot_result_join = plot_result.join(" ");
      QPixmap *pix=new QPixmap(plot_result_join);
      originalScene->addPixmap(*pix);
 //     ui->graphicsView->resize(pix->width()+500,pix->height()+500);
      ui->graphicsView->show();*/
      //ui->lncRNAFinder_usage->setText(getStr);
      QMessageBox::information(this,"SpatialScan log: \n","A SpatialScan job is done!");

      QString getStr_Error=QString(p.readAllStandardError());
      if(getStr_Error!=""){
      QMessageBox::information(this,"SpatialScan Error: \n",getStr_Error);
      }

}

void SpatialScan::on_example_clicked()
{
    spatialexample = new SpatialExample(this);
    spatialexample->show();
}

