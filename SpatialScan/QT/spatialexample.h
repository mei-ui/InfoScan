#ifndef SPATIALEXAMPLE_H
#define SPATIALEXAMPLE_H

#include <QDialog>

namespace Ui {
class SpatialExample;
}

class SpatialExample : public QDialog
{
    Q_OBJECT

public:
    explicit SpatialExample(QWidget *parent = nullptr);
    ~SpatialExample();

private:
    Ui::SpatialExample *ui;
};

#endif // SPATIALEXAMPLE_H
