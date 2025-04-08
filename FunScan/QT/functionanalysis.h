#ifndef FUNCTIONANALYSIS_H
#define FUNCTIONANALYSIS_H

#include <QDialog>
#include <QProcess>
#include <QScrollBar>
#include "functionresult.h"
QT_BEGIN_NAMESPACE
class QFile;
QT_END_NAMESPACE
namespace Ui {
class FunctionAnalysis;
}

class FunctionAnalysis : public QDialog
{
    Q_OBJECT

public:
    explicit FunctionAnalysis(QWidget *parent = nullptr);
    ~FunctionAnalysis();
private slots:
    void checkLineEdits();

    void on_begin_clicked();
    void on_start_clicked();

    void display_other_items(QString);
    void on_browse_lncRNA_clicked();
    void on_browse_geneset_clicked();
    void on_example_result_clicked();

private:
    Ui::FunctionAnalysis *ui;
    FunctionResult *functionresult;
    QString fileName;
    QFile *file;
};

#endif // FUNCTIONANALYSIS_H
