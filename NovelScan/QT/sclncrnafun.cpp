#include "sclncrnafun.h"
#include "ui_sclncrnafun.h"
#include "Htmlreportviewer.h"
#include <QFileDialog>
#include <QDialog>
#include <QMessageBox>
#include <QTextStream>
#include <QDebug>
#include <QProcess>
#include <QScreen>
#include <QProgressDialog>
#include <QWebEngineSettings>
#include <QtConcurrent>
#include <QFileInfo>

SClncRNAFun::SClncRNAFun(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SClncRNAFun)
{
    ui->setupUi(this);
    connect(ui->spinBox_fpkm, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->spinBox_length, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->spinBox_exon, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->line_fa,SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->line_gtf,SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->line_ERCC,SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->line_protein,SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->line_lncRNA,SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->line_model,SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->checkBox_exon,SIGNAL(toggled(bool)), this, SLOT(checkLineEdits()));
    connect(ui->line_pseudo,SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->line_other,SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->line_phastcons,SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->comboBox_genome,SIGNAL(currentTextChanged(QString)),this,SLOT(display_other_items(QString)));
    connect(ui->comboBox_cutoff,SIGNAL(currentTextChanged(QString)),this,SLOT(display_other_items(QString)));
    connect(ui->line_out,SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));

 /*   QVBoxLayout* mainLayout = new QVBoxLayout;
    mainLayout->addWidget(ui->scrollArea);
    setLayout(mainLayout);
*/
}

SClncRNAFun::~SClncRNAFun()
{
    delete ui;
}
void SClncRNAFun::checkLineEdits()
{
    bool start_ok1 = !ui->spinBox_length->text().isEmpty()
            && !ui->line_out->text().isEmpty();
    ui->start->setEnabled(start_ok1);

    if(ui->checkBox_exon->isChecked()){
        ui->spinBox_exon->setEnabled(true);
    }else if(!ui->checkBox_exon->isChecked()){
        ui->spinBox_exon->setEnabled(false);
    }
}

void SClncRNAFun::display_other_items(QString)
{

    ui->comboBox_genome->setItemText(0,"Human");
    ui->comboBox_genome->setItemText(1,"Mouse");
    if(ui->comboBox_genome->currentText().contains("Human"))
    {
        ui->comboBox_cutoff->setCurrentText("0.364");
        ui->line_gtf->setText("hg38/gencode.annotation.gtf");
        ui->line_fa->setText("hg38/hg38.fa");
        ui->line_ERCC->setText("ERCC92.gtf");
        ui->line_protein->setText("hg38/GENCODE_protein_coding_gene.gtf");
        ui->line_lncRNA->setText("hg38/GENCODE_long_noncoding_RNAs.gtf");
        ui->line_model->setText("hg38/human_model.hdf5");
        ui->line_pseudo->setText("hg38/GENCODE_pseudogene.gtf");
        ui->line_other->setText("hg38/GENCODE_other_rnas.gtf");
        ui->line_phastcons->setText("hg38/hg38.phastCons.bw");
    }

    if(ui->comboBox_genome->currentText().contains("Mouse"))
    {
        ui->comboBox_cutoff->setCurrentText("0.44");
        ui->line_gtf->setText("mm10/gencode.vM25.annotation.gtf");
        ui->line_fa->setText("mm10/mm10.fa");
        ui->line_ERCC->setText("ERCC92.gtf");
        ui->line_protein->setText("mm10/GENCODE_protein_coding_gene.gtf");
        ui->line_lncRNA->setText("mm10/GENCODE_long_noncoding_RNAs.gtf");
        ui->line_model->setText("mm10/mouse_model.hdf5");
        ui->line_pseudo->setText("mm10/GENCODE_pseudogene.gtf");
        ui->line_other->setText("mm10/GENCODE_other_rnas.gtf");
        ui->line_phastcons->setText("mm10/mm10.60way.phastCons.bw");
    }


}



void SClncRNAFun::on_config_clicked()
{
   lncexample = new LncExample(this);
   lncexample->show();
}

void SClncRNAFun::loadFile(const QString &fileName)
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
}
void SClncRNAFun::on_next_clicked()
{
    dataanalysis = new DataAnalysis(this);
    dataanalysis->show();
}
void SClncRNAFun::on_start_clicked()
{

    QString GENOME= ui->comboBox_genome ->currentText();
    QString GENOME2=GENOME.toLower();
    QString CPATcutoff= ui->comboBox_cutoff->currentText();
    QString fpkm=ui->spinBox_fpkm->text();
    QString length=ui->spinBox_length->text();
    QString exon=ui->spinBox_exon->text();
    QString gtf=ui->line_gtf->text();
    QString fa=ui->line_fa->text();
    QString ercc=ui->line_ERCC->text();
    QString protein=ui->line_protein->text();
    QString lncRNA=ui->line_lncRNA->text();
    QString model=ui->line_model->text();
    QString GENCODE_pseudogene=ui->line_pseudo->text();
    QString GENCODE_other_rna=ui->line_other->text();
    QString phastcons=ui->line_phastcons->text();
    QString out=ui->line_out->text();

    if(ui->checkBox_exon->isChecked()){exon=ui->spinBox_exon->text();}else{exon="0";}

    QString CPATlogitModel=GENOME+"_logitModel.RData";
    QString CPATHexamer=GENOME+"_Hexamer.tsv";


    qDebug()<<"genome:"<<GENOME2;
    qDebug()<<"CPATcutoff:"<<CPATcutoff;
    // 获取应用Bundle的资源目录路径
       QDir resourcesDir(QApplication::applicationDirPath() + "/../Resources/");
       QString scriptPath = resourcesDir.absoluteFilePath("snakemake/");
       QDir resourcesDir2(QApplication::applicationDirPath() + "/../../../");
       QString scriptPath2 = resourcesDir2.absoluteFilePath("genome/");

    QStringList lncRNAFinder_config;
    lncRNAFinder_config<<"CPATlogitModel: \""+scriptPath+"CPAT/dat/"+CPATlogitModel+"\""
               <<"CPATHexamer: \""+scriptPath+"CPAT/dat/"+CPATHexamer+"\""
              <<"genome: \""+GENOME2+"\""
                 <<"CPATcutoff: \""+CPATcutoff+"\""
                <<"fpkm_filter: \""+fpkm+"\""
                 <<"length_filter: \""+length+"\""
                  <<"exon_filter: \""+exon+"\""
                    <<"gtf: \""+scriptPath2+gtf+"\""
                      <<"fa: \""+scriptPath2+fa+"\""
                    <<"ercc: \""+scriptPath2+ercc+"\""
                      <<"protein: \""+scriptPath2+protein+"\""
                        <<"outdir: \""+out+"\""
                          <<"appDir: \""+scriptPath+"\""
                        <<"lncRNA: \""+scriptPath2+lncRNA+"\""
                            <<"GENCODE_pseudogene: \""+scriptPath2+GENCODE_pseudogene+"\""
                              <<"GENCODE_other_rna: \""+scriptPath2+GENCODE_other_rna+"\""
                                <<"phastcons: \""+scriptPath2+phastcons+"\""
                          <<"model: \""+scriptPath2+model+"\"";

    QString lncRNAFinder_config_join = lncRNAFinder_config.join(" \n");
    qDebug()<<"lncRNAFinder config:"<<lncRNAFinder_config_join;

    QString qPath(out+"/snakemake/config/NovelScan_config.yml");
    QString dirName = out+"/snakemake/config/";
    QDir dir(dirName);
    if(!dir.exists())
    {
        dir.mkdir(dirName);
        qDebug()<<"Document created successfully";
    }
    QFile qFile(qPath);
      if (qFile.open(QIODevice::WriteOnly)) {
        QTextStream out(&qFile); out << lncRNAFinder_config_join;
        qFile.close();
      }
      QString logPath(out+"/snakemake/logs/NovelScan.log");
      QString logdirName = out+"/snakemake/logs/";
      QDir logdir(logdirName);
      if(!logdir.exists())
      {
          logdir.mkdir(logdirName);
          qDebug()<<"Document created successfully";
      }
      QFile logFile(logPath);
        if (logFile.open(QIODevice::WriteOnly)) {
          QTextStream out2(&logFile); out2 << "#NovelScan log";
          logFile.close();
        }
          QStringList lncRNAFinder_command;
          lncRNAFinder_command<<"lncRNAFinder.sh"
                    <<"-c"
                     <<out+"/snakemake/config/NovelScan_config.yml"
                    <<"-t"
                     <<scriptPath
                    <<"-o"
                     <<out
                    <<"-i"
                    <<out+"/snakemake/logs/NovelScan.log";

          QString lncRNAFinder_command_join = lncRNAFinder_command.join(" ");
          qDebug()<<"lncRNAFinder command:"<<lncRNAFinder_command_join;

    //      qDebug() << QCoreApplication::applicationDirPath()+"/snakemake/script/";
          lncRNAFinder_command_join =scriptPath+"script/"+lncRNAFinder_command_join;
          qDebug() << "lncRNAFinder_command_join:" << lncRNAFinder_command_join;

          QEventLoop loop;
          QProgressDialog dialog;
          dialog.setWindowModality(Qt::WindowModal);
          dialog.setWindowTitle("NovelScan: for novel lncRNA identity");
          dialog.setWindowIcon(QIcon(":/images/qt-logo.png"));
          dialog.setMinimumWidth(800);
          dialog.setRange(0,0);
          dialog.setLabelText(QString("Begin novel lncRNA Find ...").arg(QThread::idealThreadCount()));
          dialog.show();

          QProcess p;
          p.start("bash",QStringList()<<"-c"<<lncRNAFinder_command_join);
          connect(&p,SIGNAL(finished(int,QProcess::ExitStatus)),&loop,SLOT(quit()));
          loop.exec();
          //QString getStr=QString(p.readAllStandardOutput());
          //ui->lncRNAFinder_usage->setText(getStr);
          QMessageBox::information(this,"NovelScan log: \n","A NovelScan job is done!");

          QString getStr_Error=QString(p.readAllStandardError());
          if(getStr_Error!=""){
          QMessageBox::information(this,"NovelScan Error: \n",getStr_Error);}


}

/* void SClncRNAFun::on_config_2_clicked()
{
    QString out = ui->line_out->text();
    QString start_path = out + "/snakemake/result/NovelScan_report.html";

    HtmlReportViewer *viewer = new HtmlReportViewer(start_path, this);
    viewer->show();
} */
void SClncRNAFun::on_config_2_clicked()
{
    QWebEngineView *LiveView = new QWebEngineView();
    QString out = ui->line_out->text();

    // 确保生成的 HTML 文件路径正确
    QString start_path = out + "/snakemake/result/NovelScan_report.html";
    qDebug() << "NovelScan Report" << start_path;

    // 启用插件以支持 PDF 显示
    //LiveView->settings()->setAttribute(QWebEngineSettings::PluginsEnabled, true);
    LiveView->setAttribute(Qt::WA_DeleteOnClose);

    // 确保启用滚动条，允许用户滚动查看内容
    LiveView->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    // 加载 HTML 文件
    LiveView->load(QUrl("file://" + start_path));
    LiveView->resize(1024, 768);
    LiveView->show();
}


void SClncRNAFun::on_browse_data_clicked()
{
    QFileDialog dialog(this);
    QString data_path = dialog.getExistingDirectory(this,
                                tr("Find directory"), QDir::currentPath());
        if (!data_path.isEmpty()) {
            ui->line_out->setText(data_path);
        }
}
