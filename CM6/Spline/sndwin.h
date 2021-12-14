#ifndef SNDWIN_H
#define SNDWIN_H

#include <QWidget>
#include "spline.h"
#include "qcustomplot.h"

namespace Ui {
class SndWin;
}

class SndWin : public QWidget
{
    Q_OBJECT

public:
    explicit SndWin(QWidget *parent = nullptr, spline **_s = nullptr);
    ~SndWin();

private slots:
    void on_lineEdit_textEdited(const QString &arg1);

    void on_pushButton_clicked();

private:
    Ui::SndWin *ui;
    spline *s;
};

#endif // SNDWIN_H
