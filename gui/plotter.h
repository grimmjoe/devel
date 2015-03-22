#ifndef PLOTTER_H__
#define PLOTTER_H__

#include "core/inc/apparatus.h"
#include "core/inc/matrix.h"
#include "core/inc/parser.h"
#include <string>
#include <QtGui>

class QTableWidget;
class QCustomPlot;

typedef core::parser<double> tParser;
typedef core::apparatus<double> tApp;
typedef typename tApp::tMatrixDiscrete tDiscrete;
typedef typename tApp::tMatrixDiscretes tDiscretes;
typedef typename tApp::tFunction tFunction;
typedef typename tApp::tFunctionPtr tFunctionPtr;
typedef typename tApp::tFuncMatrix tMatrix;
typedef typename tApp::tDiffInfo tDiffInfo;

class Plotter : public QWidget
{
	Q_OBJECT
public:
	Plotter(QWidget* p, tMatrix theM, double beg, double end, const QString& title);
protected:
	tMatrix m_matrix;
	double m_begin;
	double m_end;

	std::vector<std::vector<QCustomPlot*> > plotters;
	double m_delta;
protected:
	void createWidgets();
	void updateTableWithFunction(const tMatrix& mat);
	void setupPlotter(QCustomPlot* plotter, tFunctionPtr func);
protected:
	QTableWidget* m_table;
	QCustomPlot* m_current_plotter;
	QBoxLayout* m_plotter_layout;
protected slots:
	void on_select();
};

#endif // PLOTTER_H__
