#include "trdwin.h"
#include "ui_trdwin.h"

TrdWin::TrdWin(QWidget *parent, spline **_s) :
    QWidget(parent),
    ui(new Ui::TrdWin)
{
    ui->setupUi(this);
    s = *_s;
    this->setWindowTitle("Задать узлы");
    ui->tableWidget->setColumnCount(1);
    ui->tableWidget->setColumnWidth(0, 4);
    ui->tableWidget->setColumnWidth(1, 4);
}

TrdWin::~TrdWin()
{
    delete ui;
}

void TrdWin::on_lineEdit_textEdited(const QString &arg1)
{
    s->m = arg1.toDouble();

    s->nodes.clear();
    s->nodes.resize(s->m);
    ui->tableWidget->setRowCount(s->m);
}


void TrdWin::on_pushButton_clicked()
{
    for (size_t i = 0; i < s->m; i++)
        s->nodes[i] = ui->tableWidget->item(i, 0)->data(Qt::DisplayRole).toDouble();

    this->close();
}

