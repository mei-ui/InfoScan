#ifndef TUTROIAL_H
#define TUTROIAL_H

#include <QDialog>

namespace Ui {
class Tutroial;
}

class Tutroial : public QDialog
{
    Q_OBJECT

public:
    explicit Tutroial(QWidget *parent = nullptr);
    ~Tutroial();

private:
    Ui::Tutroial *ui;
};

#endif // TUTROIAL_H
