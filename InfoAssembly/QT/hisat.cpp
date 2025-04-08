#include "hisat.h"
#include "ui_hisat.h"
#include <QDialog>
#include <QDebug>
#include <QFileDialog>
#include <QMessageBox>
#include <QProcess>
#include <QScreen>
#include <QProgressDialog>
#include <QtConcurrent>
#include <QRegularExpression>

Hisat::Hisat(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Hisat)
{
    ui->setupUi(this);
    //  build:
 /*       connect(ui->line_fasta, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
        connect(ui->line_gtf, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
        connect(ui->line_snp, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
        connect(ui->line_index, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
        connect(ui->build_thread, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));*/
    //  align:
        connect(ui->align_thread, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
        connect(ui->comboBox_genome,SIGNAL(currentTextChanged(QString)),this,SLOT(display_other_items(QString)));
        connect(ui->checkBox_params2,SIGNAL(toggled(bool)), this, SLOT(checkLineEdits()));
        connect(ui->checkBox_params3,SIGNAL(toggled(bool)), this, SLOT(checkLineEdits()));
        connect(ui->checkBox_params4,SIGNAL(toggled(bool)), this, SLOT(checkLineEdits()));
        connect(ui->checkBox,SIGNAL(toggled(bool)), this, SLOT(checkLineEdits()));
        connect(ui->checkBox_2,SIGNAL(toggled(bool)), this, SLOT(checkLineEdits()));
        connect(ui->checkBox_3,SIGNAL(toggled(bool)), this, SLOT(checkLineEdits()));
        connect(ui->checkBox_4,SIGNAL(toggled(bool)), this, SLOT(checkLineEdits()));
        connect(ui->checkBox_5,SIGNAL(toggled(bool)), this, SLOT(checkLineEdits()));
        connect(ui->line_fast, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
        connect(ui->line_sensitive, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
        connect(ui->line_verysensitive, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
        connect(ui->line_pseudo, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
        connect(ui->line_qc, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
        connect(ui->hisat_gtf, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
        connect(ui->hisat_fasta, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
        connect(ui->hisat_index, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
        connect(ui->line_out, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
  /*      QVBoxLayout* mainLayout = new QVBoxLayout;
        mainLayout->addWidget(ui->scrollArea);
        setLayout(mainLayout);
*/
}

Hisat::~Hisat()
{
    delete ui;
}

void Hisat::checkLineEdits()
{
  /*  bool ok_build = !ui->line_fasta->text().isEmpty()
    && !ui->line_gtf->text().isEmpty()
    && !ui->line_snp->text().isEmpty()
    && !ui->build_thread->text().isEmpty()
    && !ui->line_index->text().isEmpty();
    ui->build->setEnabled(ok_build);*/

    bool ok_align = !ui->align_thread->text().isEmpty()
            && !ui->line_out->text().isEmpty();
    ui->start->setEnabled(ok_align);
    if(ui->checkBox_4->isChecked()){
        ui->min_intron->setEnabled(true);
    }else if(!ui->checkBox_4->isChecked()){
        ui->min_intron->setEnabled(false);}

    if(ui->checkBox_5->isChecked()){
        ui->max_intron->setEnabled(true);
    }else if(!ui->checkBox_5->isChecked()){
        ui->max_intron->setEnabled(false);}

    if(ui->checkBox_params2->isChecked()){
        ui->checkBox->setEnabled(true);
        if(ui->checkBox->isChecked()){
         ui->line_fast->setEnabled(true);
        }else if(!ui->checkBox->isChecked()){
         ui->line_fast->setEnabled(false);
        }
        ui->checkBox_2->setEnabled(true);
        if(ui->checkBox_2->isChecked()){
         ui->line_sensitive->setEnabled(true);
        }else if(!ui->checkBox_2->isChecked()){
         ui->line_sensitive->setEnabled(false);
        }
        ui->checkBox_3->setEnabled(true);
        if(ui->checkBox_3->isChecked()){
         ui->line_verysensitive->setEnabled(true);
        }else if(!ui->checkBox_3->isChecked()){
         ui->line_verysensitive->setEnabled(false);
        }

    }else if(!ui->checkBox_params2->isChecked()){
        ui->checkBox->setEnabled(false);
        ui->checkBox_2->setEnabled(false);
        ui->checkBox_3->setEnabled(false);
        ui->line_fast->setEnabled(false);
        ui->line_sensitive->setEnabled(false);
        ui->line_verysensitive->setEnabled(false);}

    if(ui->checkBox_params4->isChecked()){
        ui->line_pseudo->setEnabled(true);
    }else if(!ui->checkBox_params4->isChecked()){
        ui->line_pseudo->setEnabled(false);}

    if(ui->checkBox_params3->isChecked()){
        ui->line_qc->setEnabled(true);
    }else if(!ui->checkBox_params4->isChecked()){
        ui->line_qc->setEnabled(false);}

}

void Hisat::on_next_clicked()
{
    sclncrnafun = new SClncRNAFun(this);
    sclncrnafun->show();
}


void Hisat::display_other_items(QString)
{

    ui->comboBox_genome->setItemText(0,"Human");
    ui->comboBox_genome->setItemText(1,"Mouse");
    if(ui->comboBox_genome->currentText().contains("Human"))
    {
        ui->hisat_gtf->setText("hg38/gencode.annotation.gtf");
        ui->hisat_fasta->setText("hg38/hg38.fa");
        ui->hisat_index->setText("hg38/index/hg38_hisat2");
    }

    if(ui->comboBox_genome->currentText().contains("Mouse"))
    {
        ui->hisat_gtf->setText("mm10/gencode.vM25.annotation.gtf");
        ui->hisat_fasta->setText("mm10/mm10.fa");
        ui->hisat_index->setText("mm10/index/mm10_hisat2");
    }


}


/*void Hisat::on_build_clicked()
{
//   Basic-Generating genome indexes:
   QString BuildThread=ui->line_fasta->text();
    QString GenomeFastaFiles=ui->line_fasta->text();
    QString GTFfile=ui->line_gtf->text();
    QString SNPfile=ui->line_snp->text();
    QString GenomeIndex=ui->line_index->text();

    qDebug()<<"NumberOfThreads:"<<BuildThread;
    qDebug()<<"/path/to/fasta...:"<<GenomeFastaFiles;
    qDebug()<<"/path/to/annotations.gtf:"<<GTFfile;
    qDebug()<<"/path/to/snp:"<<SNPfile;
    qDebug()<<"/path/to/index:"<<GenomeIndex;

    QStringList build_command;
    build_command<<"snakemake/script/hisat2_build.sh"
               <<BuildThread
             <<GenomeFastaFiles
            <<GTFfile
           <<SNPfile
          <<GenomeIndex;

    QString build_command_join = build_command.join(" ");
    qDebug()<<"huisat2 build command:"<<build_command_join;

      build_command_join = QCoreApplication::applicationDirPath()+"/snakemake/script/"+build_command_join;
      qDebug() << "build_command_join:" << build_command_join;

      QEventLoop loop;

      QProcess *hisat2_build_pro = new QProcess;
      build_command_join = "\""+build_command_join+"\"";
      qDebug() << "build_command_join:" << build_command_join;
      QString hisat2_build_wrap = "nohup bash "+build_command_join;
      qDebug() << hisat2_build_wrap;
      hisat2_build_pro->startDetached(hisat2_build_wrap);
      hisat2_build_pro->waitForStarted();
      hisat2_build_pro->waitForFinished();

      loop.exec();

}*/
void Hisat::on_config_clicked()
{
    QString out=ui->line_out->text();
    hisatexample = new Hisatexample(this);
    // 添加参数传递 ↓
    hisatexample->setCurrentPath(out); // 传递out参数
    hisatexample -> show();
}

void Hisat::loadFile(const QString &fileName)
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


    ui->align_thread->setValue(config_list[12].toInt());
    ui->comboBox_genome->setCurrentText(config_list[14]);
    ui->hisat_gtf->setText(config_list[16]);




}

void Hisat::on_start_clicked()
{
    //build
/*    QString BuildThread=ui->line_fasta->text();
    QString GenomeFastaFiles=ui->line_fasta->text();
    QString GTFfile=ui->line_gtf->text();
    QString SNPfile=ui->line_snp->text();
    QString GenomeIndex=ui->line_index->text();*/
//  align
    QString check1=ui->checkBox_params2->text();
    if(ui->checkBox_params2->isChecked()){check1="true";}else{check1="false";}

    QString hisat2Thread=ui->align_thread->text();
    QString genome=ui->comboBox_genome->currentText();
    QString gtf_file=ui->hisat_gtf->displayText();
    QString genome_fa=ui->hisat_fasta->displayText();
    QString index=ui->hisat_index->displayText();
    QString params1=ui->line_fast->text();
    if(ui->checkBox->isChecked()){params1="--fast";}else{params1="";}
    QString params2=ui->line_sensitive->text();
    if(ui->checkBox_2->isChecked()){params2="--sensitive";}else{params2="";}
    QString params3=ui->line_verysensitive->displayText();
    if(ui->checkBox_3->isChecked()){params3="--very-sensitive";}else{params3="";}
    QString params4=ui->min_intron->text();
    if(ui->checkBox_4->isChecked()){params4=ui->min_intron->text();}else{params4="20";}
    QString params5=ui->max_intron->text();
    if(ui->checkBox_5->isChecked()){params5=ui->max_intron->text();}else{params5="500000";}
    QString params6=ui->line_pseudo->text();
    if(ui->checkBox_params4->isChecked()){params6=ui->line_pseudo->text();}else{params6="";}
    QString params7=ui->line_qc->text();
    if(ui->checkBox_params3->isChecked()){params7=ui->line_qc->text();}else{params7="";}
    QString out=ui->line_out->text();



    qDebug()<<"hisat run thread:"<<hisat2Thread;
    // 获取应用Bundle的资源目录路径
       QDir resourcesDir(QApplication::applicationDirPath() + "/../Resources");
       QString scriptPath = resourcesDir.absoluteFilePath("snakemake/");
       QDir resourcesDir2(QApplication::applicationDirPath() + "/../../../");
       QString scriptPath2 = resourcesDir2.absoluteFilePath("genome/");
        qDebug()<<"genome dir:"<<scriptPath2;
    QStringList hisat2_config;
    hisat2_config<<"Genome: \""+genome+"\""
            <<"GenomeFastaFiles: \""+scriptPath2+genome_fa+"\""
            <<"GTFfile: \""+scriptPath2+gtf_file+"\""
              <<"GenomeIndex: \""+scriptPath2+index+"\""
               <<"thread: \""+hisat2Thread+"\""
              <<"params_fast: \""+params1+"\""
                <<"params_sensitive: \""+params2+"\""
               <<"min_intron: \""+params4+"\""
                <<"max_intron: \""+params5+"\""
                  <<"params_verysensitive: \""+params3+"\""
                   <<"avoid pseudogenes: \""+params6+"\""
                     <<"outputDir: \""+out+"\""
                     <<"appDir: \""+scriptPath+"\""
                      <<"params_qc: \""+params7+"\"";


    QString hisat2_config_join = hisat2_config.join(" \n");
    qDebug()<<"hisat2 config:"<<hisat2_config_join;

    QString qPath(out+"/snakemake/config/hisat2_config.yml");
    QString dirName = out+"/snakemake/config/";
    QDir dir(dirName);
    if(!dir.exists())
    {
        dir.mkdir(dirName);
        qDebug()<<"Document created successfully";
    }
    QFile qFile(qPath);
      if (qFile.open(QIODevice::WriteOnly)) {
        QTextStream out(&qFile); out << hisat2_config_join;
        qFile.close();
      }
      QString logPath(out+"/snakemake/logs/InfoAssembly.log");
      QString logdirName = out+"/snakemake/logs/";
      QDir logdir(logdirName);
      if(!logdir.exists())
      {
          logdir.mkdir(logdirName);
          qDebug()<<"Document created successfully";
      }
      QFile logFile(logPath);
        if (logFile.open(QIODevice::WriteOnly)) {
          QTextStream out2(&logFile); out2 << "#InfoAssembly log";
          logFile.close();
        }


      QStringList hisat2_command;
      hisat2_command<<"hisat2.sh"
                <<"-c"
                 <<out+"/snakemake/config/hisat2_config.yml"
                <<"-o"
                 <<out
                <<"-t"
                 <<scriptPath
                <<"-i"
                <<out+"/snakemake/logs/InfoAssembly.log";

      QString hisat2_command_join = hisat2_command.join(" ");
      qDebug()<<"hisat2 command:"<<hisat2_command_join;

//      qDebug() << QCoreApplication::applicationDirPath()+"/snakemake/script/";

      hisat2_command_join = scriptPath+"script/"+hisat2_command_join;
      qDebug() << "hisat2_command_join:" << hisat2_command_join;

      QEventLoop loop;
      QProgressDialog dialog;
      dialog.setWindowModality(Qt::WindowModal);
      dialog.setWindowTitle("InfoAssembly: for transcript alignment,assembly,merge");
      dialog.setWindowIcon(QIcon("InfoAssembly:/images/qt-logo.png"));
      dialog.setMinimumWidth(800);
      dialog.setRange(0,0);
      dialog.setLabelText(QString("Begin InfoAssembly Job ...").arg(QThread::idealThreadCount()));
      dialog.show();

      QProcess p;
      p.start("bash",QStringList()<<"-c"<<hisat2_command_join);
      connect(&p,SIGNAL(finished(int,QProcess::ExitStatus)),&loop,SLOT(quit()));
      loop.exec();
      //QString getStr=QString(p.readAllStandardOutput());
      //ui->hisat2_usage->setText(getStr);
      QMessageBox::information(this,"InfoAssembly log: \n","A InfoAssembly job is done!");

      QString getStr_Error=QString(p.readAllStandardError());
      if(getStr_Error!=""){
      QMessageBox::information(this,"InfoAssembly Error: \n",getStr_Error);
      }
}

void Hisat::on_example_result_clicked()
{
    hisatresult = new HisatResult(this);
    hisatresult->show();
}

void Hisat::on_browse_data_clicked()
{
    QString out=ui->line_out->text();
    QFileDialog dialog(this);
    QString data_path = dialog.getExistingDirectory(this,
                                tr("Find directory"), out);
        if (!data_path.isEmpty()) {
            ui->line_out->setText(data_path);
        }
}
