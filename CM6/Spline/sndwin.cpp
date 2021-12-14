#include "sndwin.h"
#include "ui_sndwin.h"

SndWin::SndWin(QWidget *parent, spline **_s) :
    QWidget(parent),
    ui(new Ui::SndWin)
{
    ui->setupUi(this);
    this->setWindowTitle("Задать точки");
    ui->tableWidget->setColumnCount(2);
    ui->tableWidget->setHorizontalHeaderLabels(QStringList() << "X" << "Y");
    ui->tableWidget->setColumnWidth(0, 4);
    ui->tableWidget->setColumnWidth(1, 4);

    s = *_s;
}

SndWin::~SndWin()
{
    delete ui;
}

void SndWin::on_lineEdit_textEdited(const QString &arg1)
{
    s->n = arg1.toInt();
    s->points.clear();
    s->points.resize(s->n);
    ui->tableWidget->setRowCount(s->n);
}

void SndWin::on_pushButton_clicked()
{
    for (size_t i = 0; i < s->n; i++)
    {
        s->points[i].x = ui->tableWidget->item(i, 0)->data(Qt::DisplayRole).toDouble();
        s->points[i].y = ui->tableWidget->item(i, 1)->data(Qt::DisplayRole).toDouble();
    }
    this->close();
}

