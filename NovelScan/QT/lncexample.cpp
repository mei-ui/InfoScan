#include "lncexample.h"
#include "ui_lncexample.h"

LncExample::LncExample(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::LncExample)
{
    ui->setupUi(this);
}

LncExample::~LncExample()
{
    delete ui;
}
