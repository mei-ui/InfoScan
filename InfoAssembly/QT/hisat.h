#ifndef HISAT_H
#define HISAT_H

#include <QDialog>
#include "sclncrnafun.h"
#include "hisatexample.h"
#include "hisatresult.h"
QT_BEGIN_NAMESPACE
class QFile;
QT_END_NAMESPACE

namespace Ui {
class Hisat;
}

class Hisat : public QDialog
{
    Q_OBJECT

public:
    explicit Hisat(QWidget *parent = nullptr);
    ~Hisat();

public:
    void loadFile(const QString &fileName);

signals:
    void display(int number);

private slots:
    void checkLineEdits();
    void on_config_clicked();
    void on_start_clicked();
 //   void on_build_clicked();
    void on_browse_data_clicked();
    void display_other_items(QString);
    void on_next_clicked();
    void on_example_result_clicked();

private:
    Ui::Hisat *ui;
    SClncRNAFun *sclncrnafun;
    Hisatexample *hisatexample;
    HisatResult *hisatresult;
    QString fileName;
    QFile *file;
};

#endif // HISAT_H

