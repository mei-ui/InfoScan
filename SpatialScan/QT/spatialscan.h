#ifndef SPATIALSCAN_H
#define SPATIALSCAN_H

#include <QDialog>
#include "spatialexample.h"

namespace Ui {
class SpatialScan;
}

class SpatialScan : public QDialog
{
    Q_OBJECT

public:
    explicit SpatialScan(QWidget *parent = nullptr);
    ~SpatialScan();

public:
    void loadFile(const QString &fileName);

signals:
    void display(int number);

private slots:
    void checkLineEdits();
    void on_start_clicked();
    void on_browse_scRNA_clicked();
    void on_browse_stRNA_clicked();

    void on_example_clicked();

private:
    Ui::SpatialScan *ui;
    SpatialExample *spatialexample;
    QString fileName;
};

#endif // SPATIALSCAN_H
