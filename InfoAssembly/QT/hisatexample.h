#ifndef HISATEXAMPLE_H
#define HISATEXAMPLE_H
#include <QListWidgetItem>
#include <QWebEngineView>
#include <QWebChannel>
#include <QWebEngineProfile>
#include <QDialog>
#include <QIcon>
#include <QDir>

namespace Ui {
class Hisatexample;
}

class Hisatexample : public QDialog
{
    Q_OBJECT

public:
    explicit Hisatexample(QWidget *parent = nullptr);
    ~Hisatexample();
    void setCurrentPath(const QString& path); // 添加公共方法
public slots:
    void showCurrentDirFiles();
    void showNextDirFiles(QListWidgetItem *item);
    void showFileInfoList(QFileInfoList pInfoList);
    QIcon *getItemPropertyIcon(int iType);
private:
    Ui::Hisatexample *ui;
};

#endif // HISATEXAMPLE_H
