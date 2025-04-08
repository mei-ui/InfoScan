#ifndef SCLNCRNAFUN_H
#define SCLNCRNAFUN_H

#include <QDialog>
#include <QProcess>
#include <QScrollBar>
#include "dataanalysis.h"
#include "lncexample.h"

QT_BEGIN_NAMESPACE
class QFile;
QT_END_NAMESPACE

namespace Ui {
class SClncRNAFun;
}

class SClncRNAFun : public QDialog
{
    Q_OBJECT

public:
    explicit SClncRNAFun(QWidget *parent = nullptr);
    ~SClncRNAFun();

public:
    void loadFile(const QString &fileName);
signals:
    void display(int number);
private slots:
    void checkLineEdits();
    void on_config_clicked();
    void on_start_clicked();
    void display_other_items(QString);
    void on_next_clicked();

    void on_config_2_clicked();
    void on_browse_data_clicked();

private:
    Ui::SClncRNAFun *ui;
    LncExample *lncexample;
    DataAnalysis *dataanalysis;
    QString fileName;
    QFile *file;
};

#endif // SCLNCRNAFUN_H
