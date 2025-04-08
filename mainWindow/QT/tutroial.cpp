#include "tutroial.h"
#include "ui_tutroial.h"

Tutroial::Tutroial(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Tutroial)
{
    ui->setupUi(this);
}

Tutroial::~Tutroial()
{
    delete ui;
}
