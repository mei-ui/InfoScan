#ifndef DATAANALYSIS_H
#define DATAANALYSIS_H

#include <QDialog>
#include <QProcess>
#include <QScrollBar>
#include "dataexample.h"

QT_BEGIN_NAMESPACE
class QFile;
QT_END_NAMESPACE

namespace Ui {
class DataAnalysis;
}

class DataAnalysis : public QDialog
{
    Q_OBJECT

public:
    explicit DataAnalysis(QWidget *parent = nullptr);
    ~DataAnalysis();
public:
    void loadFile(const QString &fileName);

private slots:
    void checkLineEdits();

    void on_config_clicked();

    void on_start_clicked();

    void display_other_items(QString);

    void on_browse_metadata_clicked();

private:
    Ui::DataAnalysis *ui;
    DataExample *dataexample;
    QString fileName;
    QFile *file;
};
#endif // DATAANALYSIS_H
