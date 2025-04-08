#include "fastp2.h"
#include "ui_fastp2.h"
#include <QDebug>
#include <QFileDialog>
#include <QMessageBox>
#include <QScreen>
#include <QProcess>
#include <QProgressDialog>
#include <QtConcurrent>

using namespace QtConcurrent;

Fastp2::Fastp2(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Fastp2)
{
    ui->setupUi(this);
    connect(ui->R1, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->R2, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));

    connect(ui->q_value, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->thread_value, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));

    connect(ui->data_dir, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->checkBox_adapter,SIGNAL(toggled(bool)), this, SLOT(checkLineEdits()));
    connect(ui->checkBox_example,SIGNAL(toggled(bool)), this, SLOT(checkLineEdits()));
    connect(ui->line_example,SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->line_bam,SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->line_gtf,SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->comboBox,SIGNAL(currentTextChanged(QString)),this,SLOT(display_other_items(QString)));
    connect(ui->comboBox_2,SIGNAL(currentTextChanged(QString)),this,SLOT(display_other_items(QString)));
    connect(ui->line_output,SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->line_output2,SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));

/*    QVBoxLayout* mainLayout = new QVBoxLayout;
    mainLayout->addWidget(ui->scrollArea);
    setLayout(mainLayout);

*/

}

Fastp2::~Fastp2()
{
    delete ui;
}
void Fastp2::checkLineEdits()
{
    bool ok_1 =!ui->q_value->text().isEmpty()
            && !ui->thread_value->text().isEmpty()
            && !ui->line_output->text().isEmpty()
            && !ui->data_dir->text().isEmpty();
    ui->start->setEnabled(ok_1);
    bool ok2 =!ui->line_bam->text().isEmpty()
            && !ui->line_output2->text().isEmpty();

    ui->begin->setEnabled(ok2);

   if(ui->checkBox_adapter->isChecked()){
        ui->R1->setEnabled(true);
        ui->R2->setEnabled(true);

    }else if (!ui->checkBox_adapter->isChecked()) {
        ui->R1->setEnabled(false);
        ui->R2->setEnabled(false);
    }

    if(ui->checkBox_example->isChecked()){
        ui->line_example->setEnabled(true);
        QString example_dir=ui->line_example->text();
        ui->line_bam->setText(example_dir);
    }else if (!ui->checkBox_example->isChecked()) {
        ui->line_example->setEnabled(false);

    }

}
void Fastp2::display_other_items(QString)
{

    ui->comboBox->setItemText(0,"Human");
    ui->comboBox->setItemText(1,"Mouse");
    if(ui->comboBox->currentText().contains("Human"))
    {
        ui->line_gtf->setText("hg38/gencode.annotation.gtf");
    }

    if(ui->comboBox->currentText().contains("Mouse"))
    {
        ui->line_gtf->setText("mm10/gencode.vM25.annotation.gtf");
    }
    ui->comboBox_2->setItemText(0,"paired-end");
    ui->comboBox_2->setItemText(1,"single-end");
}
void Fastp2::on_browse_data_clicked()
{
    QFileDialog dialog(this);
    QString data_path = dialog.getExistingDirectory(this,
                                tr("Find directory"), QDir::currentPath());
        if (!data_path.isEmpty()) {
            ui->data_dir->setText(data_path);
        }
}




void Fastp2::loadFile(const QString &fileName)
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
    config_text.remove(QRegularExpression("[\n\t\r\":]"));
    QStringList config_list = config_text.split(' ');
    qDebug()<<"config_list: "<<config_list;

    if(config_list[1].contains("true")){ui->checkBox_adapter->setChecked(true); ui->R1->setText(config_list[3]);ui->R2->setText(config_list[5]);}else {ui->checkBox_adapter->setChecked(false);}

    ui->data_dir->setText(config_list[7]);
    ui->q_value->setValue(config_list[10].toInt());
    ui->thread_value->setValue(config_list[12].toInt());


}

void Fastp2::on_start_clicked()
{
    QString Reads1=ui->R1->text();
    QString Reads2=ui->R2->text();
    QString Reads=ui->R1->text();

    QString qvalue=ui->q_value->text();
    QString threadval=ui->thread_value->text();
    QString data_path=ui->data_dir->text();
    QString library_preparation=ui->comboBox_2->currentText();
    QString out=ui->line_output->text();

    qDebug()<<"reads1 adapter <reads1.fastq>:"<<Reads1;
    qDebug()<<"reads2 adapter <reads2.fastq>:"<<Reads2;

    qDebug()<<"quality_phred:"<<qvalue;
    qDebug()<<"thread:"<<threadval;

    qDebug()<<"data dir<path>:"<<data_path;

//  checkBox:
    QString adapter_input = ui->checkBox_adapter->text();
    if(ui->checkBox_adapter->isChecked()){adapter_input = "true";}else {adapter_input = "false";}

    // 获取应用Bundle的资源目录路径
       QDir resourcesDir(QApplication::applicationDirPath() + "/../Resources");
       QString scriptPath = resourcesDir.absoluteFilePath("snakemake/");

       // 检查脚本文件是否存在
       if (!QFileInfo::exists(scriptPath)) {
           QMessageBox::critical(this, "Error",
               QString("Script not found at:\n%1\n\nApplication dir: %2")
                   .arg(scriptPath)
                   .arg(QApplication::applicationDirPath()));
           return;
       }

    QStringList fastp_config;
    fastp_config<<"adapter_input: \""+adapter_input+"\""
            <<"adapter_read1: \""+Reads1+"\""
            <<"adapter_read2: \""+Reads2+"\""
              <<"adapter_read: \""+Reads+"\""
              <<"data_dir: \""+data_path+"\""
              <<"quality_phred : \""+qvalue+"\""
                <<"library_preparation : \""+library_preparation+"\""
               <<"thread: \""+threadval+"\""
                <<"outputdir: \""+out+"\""
                  <<"appDir: \""+scriptPath+"\"";


    QString fastp_config_join = fastp_config.join(" \n");
    qDebug()<<"fastp config:"<<fastp_config_join;

    QString qPath(out+"/snakemake/config/upload_fastq_config.yml");
    QString dirName = out+"/snakemake/config/";
    QDir dir(dirName);
    if(!dir.exists())
    {
        dir.mkdir(dirName);
        qDebug()<<"Document created successfully";
    }
    QFile qFile(qPath);
      if (qFile.open(QIODevice::WriteOnly)) {
        QTextStream out(&qFile); out << fastp_config_join;
        qFile.close();
      }


      QStringList fastp_command;
      fastp_command<<"cp.sh"
                 <<data_path
               <<scriptPath
                 <<out
              <<out+"/snakemake/logs/upload_fastq.log";

      QString fastp_command_join = fastp_command.join(" ");
      qDebug()<<"fastp command:"<<fastp_command_join;

//      qDebug() << QCoreApplication::applicationDirPath()+"/snakemake/script/";
      fastp_command_join = scriptPath+"script/"+fastp_command_join;
      qDebug() << "fastp_command_join:" << fastp_command_join;

      QEventLoop loop;
      QProgressDialog dialog;
      dialog.setWindowModality(Qt::WindowModal);
      dialog.setWindowTitle("fastp: for removing adapter and low quality reads");
      dialog.setWindowIcon(QIcon(":/images/qt-logo.png"));
      dialog.setMinimumWidth(800);
      dialog.setRange(0,0);
      dialog.setLabelText(QString("Begin Quality Control ...").arg(QThread::idealThreadCount()));
      dialog.show();

      QProcess p;
      p.start("bash",QStringList()<<"-c"<<fastp_command_join);
      connect(&p,SIGNAL(finished(int,QProcess::ExitStatus)),&loop,SLOT(quit()));
      loop.exec();
      //QString getStr=QString(p.readAllStandardOutput());
      //ui->fastp_usage->setText(getStr);
      QMessageBox::information(this,"fastp log: \n","A fastp job is done!");

      QString getStr_Error=QString(p.readAllStandardError());
      if(getStr_Error!=""){
      QMessageBox::information(this,"fastp Error: \n",getStr_Error);
      }
}

void Fastp2::on_next_clicked()
{
    hisat = new Hisat(this);
    hisat->show();
}

void Fastp2::on_config_clicked()
{
    fastpexample = new Fastpexample(this);
    fastpexample -> show();
}

void Fastp2::on_next_lnc_clicked()
{
    sclncrnafun = new SClncRNAFun(this);
    sclncrnafun -> show();
}

void Fastp2::on_browse_bam_clicked()
{
    QFileDialog dialog(this);
    QString data_path = dialog.getExistingDirectory(this,
                                tr("Find directory"), QDir::currentPath());
        if (!data_path.isEmpty()) {
            ui->line_bam->setText(data_path);
        }
}

void Fastp2::on_begin_clicked()
{
    QString outdir = ui->line_bam->text();
    QString bam=ui->line_bam->text();
    QString gtf=ui->line_gtf->text();
    QString out=ui->line_output2->text();
    QDir resourcesDir(QApplication::applicationDirPath() + "/../Resources");
    QString scriptPath = resourcesDir.absoluteFilePath("snakemake/");
    QDir resourcesDir2(QApplication::applicationDirPath() + "/../../../");
    QString scriptPath2 = resourcesDir2.absoluteFilePath("genome/");
    QStringList fastp_config;
    fastp_config<<"bam_path: \""+bam+"\""
                <<"GTFfile: \""+scriptPath2+gtf+"\""
               <<"outputdir: \""+out+"\""
                 <<"appDir: \""+scriptPath+"\"";


    QString fastp_config_join = fastp_config.join(" \n");
    qDebug()<<"fastp config:"<<fastp_config_join;

    QString qPath(out+"/snakemake/config/upload_bam_config.yml");
    QString dirName = out+"/snakemake/config/";
    QDir dir(dirName);
    if(!dir.exists())
    {
        dir.mkdir(dirName);
        qDebug()<<"Document created successfully";
    }
    QFile qFile(qPath);
      if (qFile.open(QIODevice::WriteOnly)) {
        QTextStream out(&qFile); out << fastp_config_join;
        qFile.close();
      }

      QString logPath(out+"/snakemake/logs/upload_bam.log");
      QString logdirName = out+"/snakemake/logs/";
      QDir logdir(logdirName);
      if(!logdir.exists())
      {
          logdir.mkdir(logdirName);
          qDebug()<<"Document created successfully";
      }
      QFile logFile(logPath);
        if (logFile.open(QIODevice::WriteOnly)) {
          QTextStream out2(&logFile); out2 << "#upload bam file log";
          logFile.close();
        }
        // 获取应用Bundle的资源目录路径
      QStringList fastp_command;
      fastp_command<<"cp_bam.sh"
                 <<bam
                <<scriptPath
               <<out
              <<out+"/snakemake/logs/upload_bam.log";

      QString fastp_command_join = fastp_command.join(" ");
      qDebug()<<"fastp command:"<<fastp_command_join;

//      qDebug() << QCoreApplication::applicationDirPath()+"/snakemake/script/";

      fastp_command_join = scriptPath+"script/"+fastp_command_join;
      qDebug() << "fastp_command_join:" << fastp_command_join;

      QEventLoop loop;
      QProgressDialog dialog;
      dialog.setWindowModality(Qt::WindowModal);
      dialog.setWindowTitle("InfoUpload: Used to upload data,assemble and merge transcripts ");
      dialog.setWindowIcon(QIcon(":/images/qt-logo.png"));
      dialog.setMinimumWidth(800);
      dialog.setRange(0,0);
      dialog.setLabelText(QString("Begin InfoUpload Job ...").arg(QThread::idealThreadCount()));
      dialog.show();

      QProcess p;
      p.start("bash",QStringList()<<"-c"<<fastp_command_join);
      connect(&p,SIGNAL(finished(int,QProcess::ExitStatus)),&loop,SLOT(quit()));
      loop.exec();
      //QString getStr=QString(p.readAllStandardOutput());
      //ui->fastp_usage->setText(getStr);
      QMessageBox::information(this,"InfoUploa log: \n","A upload data job is done!");

      QString getStr_Error=QString(p.readAllStandardError());
      if(getStr_Error!=""){
      QMessageBox::information(this,"InfoUpload Error: \n",getStr_Error);
      }
}

void Fastp2::on_example_result_clicked()
{
    fastpresult = new FastpResult(this);
    fastpresult->show();
}





void Fastp2::on_browse_output2_clicked()
{
    QFileDialog dialog(this);
    QString data_path = dialog.getExistingDirectory(this,
                                tr("Find directory"), QDir::currentPath());
        if (!data_path.isEmpty()) {
            ui->line_output2->setText(data_path);
        }
}

void Fastp2::on_browse_output_clicked()
{
    QFileDialog dialog(this);
    QString data_path = dialog.getExistingDirectory(this,
                                tr("Find directory"), QDir::currentPath());
        if (!data_path.isEmpty()) {
            ui->line_output->setText(data_path);
        }
}
