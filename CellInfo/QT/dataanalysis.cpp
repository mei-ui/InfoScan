#include "dataanalysis.h"
#include "ui_dataanalysis.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QTextStream>
#include <QDebug>
#include <QProcess>
#include <QDesktopWidget>
#include <QProgressDialog>
#include <QtConcurrent>

DataAnalysis::DataAnalysis(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DataAnalysis)
{
    ui->setupUi(this);
    connect(ui->line_phsatCon, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->line_specie, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->SingleR, SIGNAL(toggled(bool)), this, SLOT(checkLineEdits()));
    connect(ui->Celltypist, SIGNAL(toggled(bool)), this, SLOT(checkLineEdits()));
    connect(ui->SingleR_comboBox,SIGNAL(currentTextChanged(QString)),this,SLOT(display_other_items(QString)));
    connect(ui->Celltypist_comboBox,SIGNAL(currentTextChanged(QString)),this,SLOT(display_other_items(QString)));
    connect(ui->line_metadata, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->line_group, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->comboBox_genome,SIGNAL(currentTextChanged(QString)),this,SLOT(display_other_items(QString)));
    connect(ui->spinBox_mt, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->spinBox_ercc, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->spinBox_max, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->spinBox_min, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->spinBox_pc, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->doubleSpinBox_resolution, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->doubleSpinBox_conservation, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->spinBox_fpkm, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->line_pathway, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->checkBox_group,SIGNAL(toggled(bool)), this, SLOT(checkLineEdits()));
    
}

DataAnalysis::~DataAnalysis()
{
    delete ui;
}

void DataAnalysis::checkLineEdits()
{
    bool ok_align = !ui->line_phsatCon->text().isEmpty()
    && !ui->line_specie->text().isEmpty()
    && !ui->SingleR->isChecked() || !ui->Celltypist->isChecked();
    ui->start->setEnabled(ok_align);
/* 是否提供metadata信息*/
    if(ui->checkBox_group->isChecked()){
        ui->line_metadata->setEnabled(true);
        ui->line_group->setEnabled(true);
        ui->browse_metadata->setEnabled(true);
    }else if(!ui->checkBox_group->isChecked()){
        ui->line_metadata->setEnabled(false);
        ui->line_group->setEnabled(false);
        ui->browse_metadata->setEnabled(false);};
/* 选择细胞类型鉴定的方法 */
    if(ui->SingleR->isChecked()){
        ui->SingleR_comboBox->setEnabled(true);
        ui->Celltypist_comboBox->setEnabled(false);}
    if(!ui->SingleR->isChecked()){
        ui->SingleR_comboBox->setEnabled(false);}
    if(ui->Celltypist->isChecked()){
        ui->Celltypist_comboBox->setEnabled(true);
        ui->SingleR_comboBox->setEnabled(false);}
    if(!ui->Celltypist->isChecked()){
        ui->Celltypist_comboBox->setEnabled(false);}

}
void DataAnalysis::display_other_items(QString)
{

    ui->comboBox_genome->setItemText(0,"Human");
    ui->comboBox_genome->setItemText(1,"Mouse");
    if(ui->comboBox_genome->currentText().contains("Human"))
    {
        ui->line_phsatCon->setText("genome/hg38/hg38.60way.phastCons.bw");
        ui->line_specie->setText("homo sapiens");
        ui->line_pathway->setText("snakemake/script/pathway/hg38_gene2File.txt");
    }

    if(ui->comboBox_genome->currentText().contains("Mouse"))
    {
        ui->line_phsatCon->setText("genome/mm10/mm10.60way.phastCons.bw");
        ui->line_specie->setText("mus musculus");
        ui->line_pathway->setText("snakemake/script/pathway/mm10_gene2File.txt");

    }
}

void DataAnalysis::on_browse_metadata_clicked()
{
    QFileDialog dialog(this);
    dialog.setAcceptMode(QFileDialog::AcceptSave);
    QString meta = dialog.getOpenFileName(this,tr("Locate metadata"),"text|txt",tr("Find metadata File (*.txt)"));

    if (!meta.isEmpty()) {
        ui->line_metadata->setText(meta);
    }
}

void DataAnalysis::on_config_clicked()
{
   dataexample = new DataExample(this);
   dataexample->show();
}

void DataAnalysis::loadFile(const QString &fileName)
{
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot read file %1:\n%2.")
                             .arg(QDir::toNativeSeparators(fileName), file.errorString()));
        return;
    }

    QTextStream in(&file);
    QString config_text=in.readAll();
    config_text.remove(QRegExp("[\n\t\r\":]"));
    QStringList config_list = config_text.split(' ');
    qDebug()<<"config_list: "<<config_list;


}
void DataAnalysis::on_start_clicked()
{

    QString GENOME= ui->comboBox_genome ->currentText();
    QString Phsatcon=ui->line_phsatCon->text();
    QString specie= ui->line_specie->text();
    QString check1=ui->checkBox_group->text();
    if(ui->checkBox_group->isChecked()){check1="true";}else{check1="false";}
    QString metadata=ui->line_metadata->text();
    QString pathway=ui->line_pathway->text();
    QString mt=ui->spinBox_mt->text();
    QString ercc=ui->spinBox_ercc->text();
    QString maxGene=ui->spinBox_max->text();
    QString minGene=ui->spinBox_min->text();
    QString resolution=ui->doubleSpinBox_resolution->text();
    QString pca=ui->spinBox_pc->text();
    QString group=ui->line_group->text();
    QString conservation=ui->doubleSpinBox_conservation->text();
    QString fpkm=ui->spinBox_fpkm->text();
    QString singleR= ui->SingleR_comboBox ->currentText();
    QString celltypist= ui->Celltypist_comboBox ->currentText();
    if(singleR=="HumanPrimaryCellAtlasData - (The Human Primary Cell Atlas (Mabbott et al. 2013) (Human) )"){singleR="snakemake/script/ref_Human_all.RData";}


    QStringList DataAnalysis_config;
    DataAnalysis_config<<"Genome: \""+GENOME+"\""
               <<"minGene: \""+minGene+"\""
              <<"maxGene: \""+maxGene+"\""
                 <<"pctMT: \""+mt+"\""
                <<"PCAnum: \""+pca+"\""
                 <<"metadata: \""+metadata+"\""
                  <<"phastcons: \""+Phsatcon+"\""
                    <<"fpkm_filter: \""+fpkm+"\""
                      <<"conservation_filter: \""+conservation+"\""
                    <<"G0file: \""+pathway+"\""
                      <<"specie: \""+specie+"\""
                         <<"resolution: \""+resolution+"\""
                           <<"Celltypist: \""+celltypist+"\""
                           <<"pctERCC: \""+ercc+"\""
                          <<"SingleR: \""+singleR+"\""
                          <<"group_chose: \""+group+"\"";

    QString DataAnalysis_config_join = DataAnalysis_config.join(" \n");
    qDebug()<<"DataAnalysis config:"<<DataAnalysis_config_join;

    QString qPath("snakemake/config/DataAnalysis_config.yml");
    QFile qFile(qPath);
      if (qFile.open(QIODevice::WriteOnly)) {
        QTextStream out(&qFile); out << DataAnalysis_config_join;
        qFile.close();
      }

      if(ui->checkBox_group->isChecked()){
      QStringList DataAnalysis_command;
      DataAnalysis_command<<"scDataAnalysis.sh"
                <<"-c"
                 <<"snakemake/config/DataAnalysis_config.yml"
                <<"-i"
                <<"snakemake/logs/CellInfo.log";

      QString DataAnalysis_command_join = DataAnalysis_command.join(" ");
      qDebug()<<"lncRNAFinder command:"<<DataAnalysis_command_join;

//      qDebug() << QCoreApplication::applicationDirPath()+"/snakemake/script/";
      DataAnalysis_command_join = QCoreApplication::applicationDirPath()+"/snakemake/script/"+DataAnalysis_command_join;
      qDebug() << "lncRNAFinder_command_join:" << DataAnalysis_command_join;

      QEventLoop loop;
      QProgressDialog dialog;
      dialog.setWindowModality(Qt::WindowModal);
      dialog.setWindowTitle("CellInfo: for scRNA data analysis");
      dialog.setWindowIcon(QIcon(":/images/qt-logo.png"));
      dialog.setMinimumWidth(800);
      dialog.setRange(0,0);
      dialog.setLabelText(QString("Begin CellInfo analysis (%1 threads)...").arg(QThread::idealThreadCount()));
      dialog.show();

      QProcess p;
      p.start("bash",QStringList()<<"-c"<<DataAnalysis_command_join);
      connect(&p,SIGNAL(finished(int,QProcess::ExitStatus)),&loop,SLOT(quit()));
      loop.exec();
      //QString getStr=QString(p.readAllStandardOutput());
      //ui->lncRNAFinder_usage->setText(getStr);
      QMessageBox::information(this,"CellInfo log: \n","A CellInfo job is done!");

      QString getStr_Error=QString(p.readAllStandardError());
      if(getStr_Error!=""){
      QMessageBox::information(this,"CellInfo Error: \n",getStr_Error);
      }}

      if(!ui->checkBox_group->isChecked()){
      QStringList DataAnalysis_command;
      DataAnalysis_command<<"scDataAnalysis2.sh"
                <<"-c"
                 <<"snakemake/config/DataAnalysis_config.yml"
                <<"-i"
                <<"snakemake/logs/CellInfo.log";

      QString DataAnalysis_command_join = DataAnalysis_command.join(" ");
      qDebug()<<"lncRNAFinder command:"<<DataAnalysis_command_join;

//      qDebug() << QCoreApplication::applicationDirPath()+"/snakemake/script/";
      DataAnalysis_command_join = QCoreApplication::applicationDirPath()+"/snakemake/script/"+DataAnalysis_command_join;
      qDebug() << "lncRNAFinder_command_join:" << DataAnalysis_command_join;

      QEventLoop loop;
      QProgressDialog dialog;
      dialog.setWindowModality(Qt::WindowModal);
      dialog.setWindowTitle("CellInfo: for scRNA data analysis");
      dialog.setWindowIcon(QIcon("images/qt-logo.png"));
      dialog.setMinimumWidth(800);
      dialog.setRange(0,0);
      dialog.setLabelText(QString("Begin CellInfo analysis (%1 threads)...").arg(QThread::idealThreadCount()));
      dialog.show();

      QProcess p;
      p.start("bash",QStringList()<<"-c"<<DataAnalysis_command_join);
      connect(&p,SIGNAL(finished(int,QProcess::ExitStatus)),&loop,SLOT(quit()));
      loop.exec();
      //QString getStr=QString(p.readAllStandardOutput());
      //ui->lncRNAFinder_usage->setText(getStr);
      QMessageBox::information(this,"CellInfo log: \n","A CellInfo job is done!");

      QString getStr_Error=QString(p.readAllStandardError());
      if(getStr_Error!=""){
      QMessageBox::information(this,"CellInfo Error: \n",getStr_Error);
      }}


}

