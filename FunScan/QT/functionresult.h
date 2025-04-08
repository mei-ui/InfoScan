#ifndef FUNCTIONRESULT_H
#define FUNCTIONRESULT_H

#include <QDialog>

namespace Ui {
class FunctionResult;
}

class FunctionResult : public QDialog
{
    Q_OBJECT

public:
    explicit FunctionResult(QWidget *parent = nullptr);
    ~FunctionResult();

private:
    Ui::FunctionResult *ui;
};

#endif // FUNCTIONRESULT_H
