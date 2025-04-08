#ifndef DATAEXAMPLE_H
#define DATAEXAMPLE_H

#include <QDialog>

namespace Ui {
class DataExample;
}

class DataExample : public QDialog
{
    Q_OBJECT

public:
    explicit DataExample(QWidget *parent = nullptr);
    ~DataExample();

private:
    Ui::DataExample *ui;
};

#endif // DATAEXAMPLE_H
