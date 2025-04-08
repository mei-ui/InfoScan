#ifndef FASTPEXAMPLE_H
#define FASTPEXAMPLE_H

#include <QDialog>
#include <QIcon>
#include <QDir>
#include "ui_fastpexample.h"
#include <QListWidgetItem>
#include <QWebEngineView>
#include <QWebChannel>
#include <QWebEngineProfile>

namespace Ui {
class Fastpexample;
}

class Fastpexample : public QDialog
{
    Q_OBJECT

public:
    explicit Fastpexample(QWidget *parent = nullptr);
    ~Fastpexample();    

public slots:
    void showCurrentDirFiles();
    void showNextDirFiles(QListWidgetItem *item);
    void showFileInfoList(QFileInfoList pInfoList);
    QIcon *getItemPropertyIcon(int iType);

private:
    Ui::Fastpexample *ui;
};

#endif // FASTPEXAMPLE_H
