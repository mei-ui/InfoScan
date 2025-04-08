#include "spatialexample.h"
#include "ui_spatialexample.h"

SpatialExample::SpatialExample(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SpatialExample)
{
    ui->setupUi(this);
}

SpatialExample::~SpatialExample()
{
    delete ui;
}
