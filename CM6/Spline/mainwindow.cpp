#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    this->setWindowTitle("Построение сплайнов");
    s = new spline();
    sw = new SndWin(0, &s);
    tw = new TrdWin(0, &s);
}

MainWindow::~MainWindow()
{
    delete ui;
}

// считать узлы
void MainWindow::on_pushButton_clicked()
{
    s->read_nodes("nodes.txt");
}

// задать узлы
//void MainWindow::on_pushButton_2_clicked()
//{
//   tw->show();
//}

// считать точки
void MainWindow::on_pushButton_3_clicked()
{
    s->read_points("points.txt");
}

// задать точки
//void MainWindow::on_pushButton_4_clicked()
//{
//    sw->show();
//}

// построить сплайн
void MainWindow::on_pushButton_5_clicked()
{
    s->slae->size = 2 * s->m;
    s->slae->A.clear();
    s->slae->b.clear();
    s->slae->q.clear();

    s->slae->A.resize(s->slae->size);
    for (size_t i = 0; i < s->slae->size; i++)
        s->slae->A[i].resize(s->slae->size);
    s->slae->b.resize(s->slae->size);
    s->slae->q.resize(s->slae->size);

    s->splineVals.clear();
    s->build();

    QVector<double> x(s->splineVals.size());
    QVector<double> y(s->splineVals.size());

    for (auto i = 0; i < x.size(); i++)
    {
        x[i] = s->splineVals[i].x;
        y[i] = s->splineVals[i].y;
    }

    ui->widget->addGraph();
    //ui->widget->graph(2)->setScatterStyle(QCPScatterStyle::ss);
    ui->widget->graph(2)->setPen(QColor(Qt::green));

    //ui->widget->clearGraphs();
    ui->widget->graph(2)->setData(x, y);
    ui->widget->replot();
    ui->widget->update();
}

// альфа
void MainWindow::on_lineEdit_textEdited(const QString &arg1)
{
    s->alpha_reg = arg1.toDouble();
}

// бета
void MainWindow::on_lineEdit_2_textEdited(const QString &arg1)
{
    s->beta_reg = arg1.toDouble();
}

// Отобразить узлы
void MainWindow::on_pushButton_6_clicked()
{
    QVector<double> x(s->m), y(s->m);
    for (auto i = 0; i < x.size(); i++)
    {
        x[i] = s->nodes[i];
        y[i] = 0;
    }

    ui->widget->addGraph();
    ui->widget->graph(0)->setScatterStyle(QCPScatterStyle::ssCrossCircle);
    ui->widget->graph(0)->setPen(QColor(Qt::red));
    ui->widget->graph(0)->setLineStyle(QCPGraph::lsNone);

    ui->widget->xAxis->setRange(x[0]-1, x[x.size() - 1] + 1);
    ui->widget->yAxis->setRange(-0.5, 5);

    ui->widget->graph(0)->setData(x, y);
    ui->widget->replot();
    ui->widget->update();
}

// Отобразить точки
void MainWindow::on_pushButton_7_clicked()
{
    QVector<double> x(s->n), y(s->n);
    for (auto i = 0; i < x.size(); i++)
    {
        x[i] = s->points[i].x;
        y[i] = s->points[i].y;
    }

    ui->widget->addGraph();
    ui->widget->graph(1)->setScatterStyle(QCPScatterStyle::ssCrossCircle);
    ui->widget->graph(1)->setPen(QColor(Qt::blue));
    ui->widget->graph(1)->setLineStyle(QCPGraph::lsNone);

    double minY, maxY;
    minY = y[0];
    maxY = y[0];

    for (auto i : y)
    {
        if (i < minY)
            minY = i;
        else if (i > maxY)
            maxY = i;
    }
    if (minY < 0)
        ui->widget->yAxis->setRange(minY - 0.5, maxY + 1);
    else
        ui->widget->yAxis->setRange(minY - 0.5, maxY + 1);

    ui->widget->graph(1)->setData(x, y);
    ui->widget->replot();
    ui->widget->update();
}

