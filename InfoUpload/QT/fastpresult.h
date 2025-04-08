#ifndef FASTPRESULT_H
#define FASTPRESULT_H

#include <QDialog>

namespace Ui {
class FastpResult;
}

class FastpResult : public QDialog
{
    Q_OBJECT

public:
    explicit FastpResult(QWidget *parent = nullptr);
    ~FastpResult();

private:
    Ui::FastpResult *ui;
};

#endif // FASTPRESULT_H
