#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
private slots:
    void about();
    void sraexplorer();
    void configuration();
    void appsetStyleSheet1();
    void appsetStyleSheet2();
    void appsetStyleSheet3();
private:
    Ui::MainWindow *ui;
private:
    void createActions();
    void createStatusBar();
};
#endif // MAINWINDOW_H
