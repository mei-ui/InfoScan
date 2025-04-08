#include "functionanalysis.h"
#include "ui_functionanalysis.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QTextStream>
#include <QDebug>
#include <QProcess>
#include <QDesktopWidget>
#include <QProgressDialog>
#include <QtConcurrent>
#include <QGraphicsScene>

FunctionAnalysis::FunctionAnalysis(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::FunctionAnalysis)
{
    ui->setupUi(this);
    connect(ui->line_id, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->line_pathway, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->line_geneset, SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
    connect(ui->comboBox,SIGNAL(currentTextChanged(QString)),this,SLOT(display_other_items(QString)));
    connect(ui->line_lnc,SIGNAL(textChanged(QString)), this, SLOT(checkLineEdits()));
}


FunctionAnalysis::~FunctionAnalysis()
{
    delete ui;
}
void FunctionAnalysis::checkLineEdits()
{
    bool ok_align = !ui->line_id->text().isEmpty()
    && !ui->line_geneset->text().isEmpty();
    ui->begin->setEnabled(ok_align);

    bool ok_align2 = !ui->line_lnc->text().isEmpty();
    ui->start->setEnabled(ok_align2);

}
void FunctionAnalysis::on_browse_lncRNA_clicked()
{
    QFileDialog dialog(this);
    dialog.setAcceptMode(QFileDialog::AcceptSave);
    QString FilePath = "11_data_analysis/coexpression";
    FilePath = QCoreApplication::applicationDirPath()+"/snakemake/result/"+FilePath;
    QString meta = dialog.getOpenFileName(this,tr("Locate txt"),FilePath,tr("Find GeneSetEnrichment result File (*.txt)"));

    if (!meta.isEmpty()) {
        ui->line_lnc->setText(meta);
    }
}
void FunctionAnalysis::on_start_clicked()
{
    QString lncRNA_id=ui->line_lnc->text();
    QString outputfile=lncRNA_id+".pdf";
    QStringList FunctionAnalysis_config;
    FunctionAnalysis_config<<"lncRNA_id: \""+lncRNA_id+"\""
               <<"outputfile: \""+outputfile+"\"";
    QString FunctionAnalysis_config_join = FunctionAnalysis_config.join(" \n");
    QString qPath("snakemake/config/FunctionPlot_config.yml");
    QFile qFile(qPath);
      if (qFile.open(QIODevice::WriteOnly)) {
        QTextStream out(&qFile); out << FunctionAnalysis_config_join;
        qFile.close();
      }
      QStringList FunAnalysis_command;
      FunAnalysis_command<<"function_plot.sh"
                <<"-I"
                 <<lncRNA_id
                <<"-O"
               <<outputfile;

      QString FunAnalysis_command_join = FunAnalysis_command.join(" \n");
      qDebug()<<"lncRNAFinder command:"<<FunAnalysis_command_join;

//      qDebug() << QCoreApplication::applicationDirPath()+"/snakemake/script/";
      FunAnalysis_command_join = QCoreApplication::applicationDirPath()+"/snakemake/script/"+FunAnalysis_command_join;
      qDebug() << "FunAnalysis_command_join:" << FunAnalysis_command_join;

      QEventLoop loop;
      QProgressDialog dialog;
      dialog.setWindowModality(Qt::WindowModal);
      dialog.setWindowTitle("genesetEnrichment Analysis");
      dialog.setWindowIcon(QIcon(":/images/qt-logo.png"));
      dialog.setMinimumWidth(800);
      dialog.setRange(0,0);
      dialog.setLabelText(QString("Begin Function analysis (%1 threads)...").arg(QThread::idealThreadCount()));
      dialog.show();

      QProcess p;
      p.start("bash",QStringList()<<"-c"<<FunAnalysis_command_join);
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
      QMessageBox::information(this,"FuncScan log: \n","A FunctionAnalysis job is done!");

      QString getStr_Error=QString(p.readAllStandardError());
      if(getStr_Error!=""){
      QMessageBox::information(this,"FuncScan Error: \n",getStr_Error);
      }

}

void FunctionAnalysis::display_other_items(QString)
{

    ui->comboBox->setItemText(0,"Human");
    ui->comboBox->setItemText(1,"Mouse");
    if(ui->comboBox->currentText().contains("Human"))
    {
        ui->line_pathway->setText("snakemake/script/pathway/hg38_gene2File.txt");
    }

    if(ui->comboBox->currentText().contains("Mouse"))
    {

        ui->line_pathway->setText("snakemake/script/pathway/mm10_gene2File.txt");

    }
}

void FunctionAnalysis::on_browse_geneset_clicked()
{
    QFileDialog dialog(this);
    dialog.setAcceptMode(QFileDialog::AcceptSave);
    QString FilePath = "pathway";
    FilePath = QCoreApplication::applicationDirPath()+"/snakemake/script/"+FilePath;
    QString meta = dialog.getOpenFileName(this,tr("Locate GeneSet"),"gmt",tr("Find GeneSet File (*.gmt)"));

    if (!meta.isEmpty()) {
        ui->line_geneset->setText(meta);
    }
}
void FunctionAnalysis::on_begin_clicked()
{
    QString Genome=ui->comboBox->currentText();
    QString pathway=ui->line_pathway->text();
    QString gene_id=ui->line_id->text();
    QString gene_set=ui->line_geneset->text();

    QStringList FunctionAnalysis_config;
    FunctionAnalysis_config<<"Genome: \""+Genome+"\""
               <<"gene_id: \""+gene_id+"\""
              <<"gene_set: \""+gene_set+"\""
                    <<"G0file: \""+pathway+"\"";

    QString DataAnalysis_config_join = FunctionAnalysis_config.join(" \n");
    qDebug()<<"FunctionAnalysisconfig config:"<<DataAnalysis_config_join;

    QString qPath("snakemake/config/FunctionAnalysis_config.yml");
    QFile qFile(qPath);
      if (qFile.open(QIODevice::WriteOnly)) {
        QTextStream out(&qFile); out << DataAnalysis_config_join;
        qFile.close();
      }
      QStringList FunAnalysis_command;
      FunAnalysis_command<<"FunctionAnalysis.sh"
                <<"-c"
                 <<"snakemake/config/FunctionAnalysis_config.yml"
                <<"-i"
                  <<"snakemake/logs/FuncScan.log";
      QString FunAnalysis_command_join = FunAnalysis_command.join(" ");
      qDebug()<<"lncRNAFinder command:"<<FunAnalysis_command_join;

//      qDebug() << QCoreApplication::applicationDirPath()+"/snakemake/script/";
      FunAnalysis_command_join = QCoreApplication::applicationDirPath()+"/snakemake/script/"+FunAnalysis_command_join;
      qDebug() << "FunAnalysis_command_join:" << FunAnalysis_command_join;

      QEventLoop loop;
      QProgressDialog dialog;
      dialog.setWindowModality(Qt::WindowModal);
      dialog.setWindowTitle("genesetEnrichment Analysis");
      dialog.setWindowIcon(QIcon(":/images/qt-logo.png"));
      dialog.setMinimumWidth(800);
      dialog.setRange(0,0);
      dialog.setLabelText(QString("Begin Function analysis (%1 threads)...").arg(QThread::idealThreadCount()));
      dialog.show();

      QProcess p;
      p.start("bash",QStringList()<<"-c"<<FunAnalysis_command_join);
      connect(&p,SIGNAL(finished(int,QProcess::ExitStatus)),&loop,SLOT(quit()));
      loop.exec();
      //QString getStr=QString(p.readAllStandardOutput());
      //ui->lncRNAFinder_usage->setText(getStr);
      QMessageBox::information(this,"FuncScan log: \n","A FunctionAnalysis job is done!");

      QString getStr_Error=QString(p.readAllStandardError());
      if(getStr_Error!=""){
      QMessageBox::information(this,"FuncScan Error: \n",getStr_Error);
      }
}

void FunctionAnalysis::on_example_result_clicked()
{
    functionresult = new FunctionResult(this);
    functionresult -> show();
}
