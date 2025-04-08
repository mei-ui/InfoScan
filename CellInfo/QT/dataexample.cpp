#include "dataexample.h"
#include "ui_dataexample.h"

DataExample::DataExample(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DataExample)
{
    ui->setupUi(this);
}

DataExample::~DataExample()
{
    delete ui;
}
