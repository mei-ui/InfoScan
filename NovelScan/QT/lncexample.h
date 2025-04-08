#ifndef LNCEXAMPLE_H
#define LNCEXAMPLE_H

#include <QDialog>

namespace Ui {
class LncExample;
}

class LncExample : public QDialog
{
    Q_OBJECT

public:
    explicit LncExample(QWidget *parent = nullptr);
    ~LncExample();

private:
    Ui::LncExample *ui;
};

#endif // LNCEXAMPLE_H
