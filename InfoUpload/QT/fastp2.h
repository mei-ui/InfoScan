#ifndef FASTP2_H
#define FASTP2_H

#include <QDialog>
#include "hisat.h"
#include "fastpexample.h"
#include "sclncrnafun.h"
#include "fastpresult.h"
#include <QWidget>

extern QString outdir;

namespace Ui {
class Fastp2;
}

class Fastp2 : public QDialog
{
    Q_OBJECT

public:
    explicit Fastp2(QWidget *parent = nullptr);
    ~Fastp2();

public:
    void loadFile(const QString &fileName);

signals:
    void display(int number);

private slots:
    void checkLineEdits();
    void on_start_clicked();
    void on_browse_data_clicked();
    void on_browse_bam_clicked();
    void on_begin_clicked();
    void on_next_lnc_clicked();
    void on_next_clicked();
    void on_config_clicked();
    void on_example_result_clicked();
    void display_other_items(QString);


    void on_browse_output2_clicked();

    void on_browse_output_clicked();

private:
    Ui::Fastp2 *ui;
    Hisat *hisat;
    Fastpexample *fastpexample;
    SClncRNAFun *sclncrnafun;
    FastpResult *fastpresult;
    QString fileName;
    QFile *file;
};

#endif // FASTP2_H
