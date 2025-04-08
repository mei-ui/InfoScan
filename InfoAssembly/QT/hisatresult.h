#ifndef HISATRESULT_H
#define HISATRESULT_H

#include <QDialog>

namespace Ui {
class HisatResult;
}

class HisatResult : public QDialog
{
    Q_OBJECT

public:
    explicit HisatResult(QWidget *parent = nullptr);
    ~HisatResult();

private:
    Ui::HisatResult *ui;
};

#endif // HISATRESULT_H
