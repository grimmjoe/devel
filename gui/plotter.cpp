#include "plotter.h"
#include "qcustomplot.h"
#include <QtGui>

Plotter::Plotter(QWidget* p, tMatrix theM, double beg, double end, const QString& title)
	: QWidget(p)
	, m_matrix(theM.getNumRows(), theM.getNumCols(), tFunctionPtr(nullptr))
	, m_begin(beg)
	, m_end(end)
	, m_delta(0.05)
{
	m_matrix = theM;
	this->setWindowTitle(title);
	plotters.resize(theM.getNumRows(), std::vector<QCustomPlot*>(theM.getNumCols(), 0));
	createWidgets();
}

void Plotter::createWidgets()
{
	m_table = new QTableWidget(this);
	updateTableWithFunction(m_matrix);
	m_current_plotter = 0;
	m_plotter_layout = new QHBoxLayout();
	QWidget* w = new QWidget();
	w->setLayout(m_plotter_layout);
	QSplitter* sp = new QSplitter(Qt::Horizontal, this);
	sp->addWidget(m_table);
	sp->addWidget(w);
	QHBoxLayout* main = new QHBoxLayout();
	main->addWidget(sp);
	this->setLayout(main);

	m_table->setSelectionBehavior(QAbstractItemView::SelectItems);
	m_table->setSelectionMode(QAbstractItemView::SingleSelection);
	connect(m_table, SIGNAL(itemSelectionChanged()), this, SLOT(on_select()));
	m_table->setCurrentCell(0, 0, QItemSelectionModel::Select);
}

void Plotter::on_select()
{
	if (m_current_plotter)
		m_current_plotter->hide();
	QList<QTableWidgetItem*> selected = m_table->selectedItems();
	assert (selected.size() == 1);
	QTableWidgetItem* s = selected.back();
	if (!plotters[s->row()][s->column()])
	{
		plotters[s->row()][s->column()] = new QCustomPlot(this);
		m_current_plotter = plotters[s->row()][s->column()];
		setupPlotter(m_current_plotter, m_matrix[s->row()+1][s->column()+1]);
		//m_current_plotter->setMinimumHeight(100);
		//m_current_plotter->setMinimumWidth(100);
		m_plotter_layout->addWidget(m_current_plotter);
	}
	else
	{
		m_current_plotter = plotters[s->row()][s->column()];
		m_current_plotter->show();
	}
}

void Plotter::setupPlotter(QCustomPlot* plotter, tFunctionPtr func)
{
	plotter->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes |
					QCP::iSelectLegend | QCP::iSelectPlottables);
	plotter->xAxis->setRange(m_begin, m_end);
	//plotter->yAxis->setRange(-5, 5);
	plotter->axisRect()->setupFullAxesBox();

	plotter->plotLayout()->insertRow(0);
	plotter->plotLayout()->addElement(0, 0, new QCPPlotTitle(plotter, this->windowTitle()));

	plotter->xAxis->setLabel("Parameter");
	plotter->yAxis->setLabel("Element");
	plotter->legend->setVisible(true);
	QFont legendFont = font();
	legendFont.setPointSize(10);
	plotter->legend->setFont(legendFont);
	plotter->legend->setSelectedFont(legendFont);
	plotter->legend->setSelectableParts(QCPLegend::spItems); // legend box shall not be selectable, only legend items

	QVector<double> x;
	QVector<double> y;
	x.reserve((m_end-m_begin)/m_delta);
	y.reserve((m_end-m_begin)/m_delta);
	for (double beg = m_begin; beg <= m_end; beg += m_delta)
	{
		x.push_back(beg);
		y.push_back((*func)(beg));
	}
	plotter->addGraph();
	plotter->graph()->setName(this->windowTitle());
	plotter->graph()->setData(x, y);
	plotter->graph()->setLineStyle((QCPGraph::LineStyle)1);
	//plotter->graph()->setLineStyle((QCPGraph::LineStyle)(rand()%5+1));
	//if (rand()%100 > 75)
	//		plotter->graph()->setScatterStyle(QCPScatterStyle((QCPScatterStyle::ScatterShape)(rand()%9+1)));
	QPen graphPen;
	graphPen.setColor(QColor(rand()%245+10, rand()%245+10, rand()%245+10));
	graphPen.setWidthF(5);
	plotter->graph()->setPen(graphPen);
	plotter->replot();
}

void Plotter::updateTableWithFunction(const tMatrix& mat)
{
	int m = mat.getNumRows();
	int n = mat.getNumCols();
	m_table->setRowCount(m);
	m_table->setColumnCount(n);
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			m_table->setItem(i, j,
				new QTableWidgetItem(QString::fromStdString(mat[i+1][j+1]->toString())));
	//m_table->resizeRowsToContents();
	//m_table->resizeColumnsToContents();
}
