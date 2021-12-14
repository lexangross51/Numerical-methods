#ifndef TRDWIN_H
#define TRDWIN_H

#include <QWidget>
#include "spline.h"

namespace Ui {
class TrdWin;
}

class TrdWin : public QWidget
{
    Q_OBJECT

public:
    explicit TrdWin(QWidget *parent = nullptr, spline **_s = nullptr);
    ~TrdWin();

private slots:
    void on_lineEdit_textEdited(const QString &arg1);

    void on_pushButton_clicked();

private:
    Ui::TrdWin *ui;
    spline *s;
};

#endif // TRDWIN_H
