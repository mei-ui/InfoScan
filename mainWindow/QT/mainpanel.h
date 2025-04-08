#ifndef MAINPANEL_H
#define MAINPANEL_H

#include <QDialog>
#include <QIcon>
#include <QDir>
#include "configuration.h"
#include <QListWidgetItem>
#include <QWebEngineView>
#include <QWebChannel>
#include <QWebEngineProfile>
namespace Ui {
class Mainpanel;
}

class Mainpanel : public QDialog
{
    Q_OBJECT

public:
    explicit Mainpanel(QWidget *parent = nullptr);
    ~Mainpanel();

private slots:
    void on_check_clicked();

    void on_check_2_clicked();

private:
    Ui::Mainpanel *ui;
    Configuration *configuration;
};

#endif // MAINPANEL_H
