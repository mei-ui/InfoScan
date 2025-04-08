#include "fastpresult.h"
#include "ui_fastpresult.h"

FastpResult::FastpResult(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::FastpResult)
{
    ui->setupUi(this);
}

FastpResult::~FastpResult()
{
    delete ui;
}
