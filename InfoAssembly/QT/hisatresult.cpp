#include "hisatresult.h"
#include "ui_hisatresult.h"

HisatResult::HisatResult(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::HisatResult)
{
    ui->setupUi(this);
}

HisatResult::~HisatResult()
{
    delete ui;
}
