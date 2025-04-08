#include "functionresult.h"
#include "ui_functionresult.h"

FunctionResult::FunctionResult(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::FunctionResult)
{
    ui->setupUi(this);
}

FunctionResult::~FunctionResult()
{
    delete ui;
}
