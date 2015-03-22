#include "main_window.h"
#include "wizards.h"
#include "plotter.h"

#include <QtGui>

#include <iostream>
#include <cassert>

static const char* sToolName = "GIMT";
static const char* sTitle = "Generalized Inverse Matrices Tool";

MainWindow::MainWindow(QWidget* p)
	: QMainWindow(p)
	, m_di(0, 0, 0)
	, m_app(m_di)
	, m_func_matrix(1, 1, nullptr)
	, m_binv_restore(1, 1, nullptr)
	, m_qinv_restore(1, 1, nullptr)
	, m_bqinv_restore(1, 1, nullptr)
	, m_drazin_recursive_restore(1, 1, nullptr)
	, m_drazin_skeleton_restore(1, 1, nullptr)
	, m_drazin_canonical_restore(1, 1, nullptr)
{
	// GIM = Generalized Inverse Matrix
	this->setObjectName(sToolName);
	this->setWindowTitle(sTitle);
	createWidgets();
	createActions();
	createMenus();
	createToolBar();

	disableAllActions();
	enableBasicActions();
	QRect screenGeometry = QApplication::desktop()->screenGeometry();
	int x = (screenGeometry.width()-this->width()) / 2;
	int y = (screenGeometry.height()-this->height()) / 2;
	this->move(x, y);

	m_rank = 0;
	m_epsilon = 0.000001;
	if (m_epsilon != 0)
	{
		m_app.setComparator(std::shared_ptr<core::comparator<double> >(new core::comparator<double>(m_epsilon)));
	}

	is_binv = false;
	is_qinv = false;
	is_bqinv = false;
	is_drazin_recursive = false;
	is_drazin_skeleton = false;
	is_drazin_canonical = false;
}

void MainWindow::createWidgets()
{
	m_matrix_widget = new QTableWidget();
	m_central_widget = new QTabWidget();
	m_central_widget->setTabsClosable(true);
	this->setCentralWidget(m_central_widget);
	m_left_dock = new QDockWidget(tr("Matrix View"), this);
	m_left_dock->setAllowedAreas(Qt::LeftDockWidgetArea);
	m_left_dock->setWidget(m_matrix_widget);
	this->addDockWidget(Qt::LeftDockWidgetArea, m_left_dock);
}

void MainWindow::createActions()
{
	// File menu
	createFileMenuActions();

	// View menu
	createViewMenuActions();

	// Tools
	createToolsMenuActions();

	// Help
	createHelpMenuActions();

}

void MainWindow::createFileMenuActions()
{
	m_new_action = new QAction(QIcon(":/images/new.png"), tr("&New"), this);
	m_new_action->setShortcuts(QKeySequence::New);
	m_new_action->setStatusTip("New session");
	connect(m_new_action, SIGNAL(triggered()), this, SLOT(on_new()));

	m_open_action = new QAction(QIcon(":/images/open.png"), tr("&Open"), this);
	m_open_action->setShortcuts(QKeySequence::Open);
	m_open_action->setStatusTip("Open session");
	connect(m_open_action, SIGNAL(triggered()), this, SLOT(on_open()));

	m_save_action = new QAction(QIcon(":/images/save.png"), tr("&Save"), this);
	m_save_action->setShortcuts(QKeySequence::Save);
	m_save_action->setStatusTip("Save session");
	connect(m_save_action, SIGNAL(triggered()), this, SLOT(on_save()));

	m_save_as_action = new QAction(QIcon(":/images/save.png"), tr("Save &As"), this);
	m_save_as_action->setShortcuts(QKeySequence::SaveAs);
	m_save_as_action->setStatusTip("Save session as...");
	connect(m_save_as_action, SIGNAL(triggered()), this, SLOT(on_saveAs()));
	
	m_close_action = new QAction(QIcon(":/images/close.png"), tr("&Close"), this);
	m_close_action->setShortcuts(QKeySequence::Close);
	m_close_action->setStatusTip("Close session");
	connect(m_close_action, SIGNAL(triggered()), this, SLOT(on_close()));

	m_quit_action = new QAction(QIcon(":/images/exit.png"), tr("&Exit"), this);
	m_quit_action->setShortcuts(QKeySequence::Quit);
	m_quit_action->setStatusTip("Exit the tool");
	connect(m_quit_action, SIGNAL(triggered()), this, SLOT(on_quit()));
}

void MainWindow::on_new()
{
	NewWizard wizard(this);
	if (!wizard.exec())
		return;
	
	int numRows = wizard.field("new.numRows").toInt();
	int numCols = wizard.field("new.numCols").toInt();
	int numDiscs = wizard.field("new.numDiscs").toInt();
	double center = wizard.field("new.center").toDouble();
	double scale = wizard.field("new.scale").toDouble();
	m_di.tv = center;
	m_di.H = scale;
	m_di.K = numDiscs;
	m_app.setDiffInfo(m_di);
	tParser p;
	m_func_matrix.resize(numRows, numCols, tFunctionPtr(nullptr));
	std::vector<std::vector<QString> > theMat = wizard.getMatrix();
	for (int i = 0; i < numRows; ++i)
		for (int j = 0; j < numCols; ++j)
			m_func_matrix[i+1][j+1] = p.parse(theMat[i][j].toStdString());
	std::cout << "Finished parsing\n";
	for (int i = 1; i <= numRows; ++i)
	{
		for (int j = 1; j <= numCols; ++j)
			std::cout << m_func_matrix[i][j]->toString() << " ";
		std::cout << std::endl;
	}
	clearMatrixInfo();
	clearCentralWidget();
	updateMatrixInfo();
	std::cout << "Updated matrix info\n";
	updateToolsActions();
	std::cout << "Updated tools actions\n";
	updateMatrixWidget();
	std::cout << "Updated matrix widget\n";
	updateCentralWidget();
}

void MainWindow::clearCentralWidget()
{
	m_central_widget->clear();
	//int c = m_central_widget->count();
	//for (int i = 0; i <= c; ++i)
	//	m_central_widget->removeTab(i);
}

void MainWindow::updateCentralWidget()
{
	QTableWidget* table = new QTableWidget();
	int index = m_central_widget->addTab(table, "A(t) - Discretes");
	m_central_widget->setCurrentIndex(index);
	updateTableWithDiscretes(table, m_matrix_discretes, "A");
}

void MainWindow::updateTableWithDiscretes(QTableWidget* table, tDiscretes& discretes, const QString& prefix)
{
	assert (m_di.K >= 0);
	int m = discretes[0].getNumRows();
	int n = discretes[0].getNumCols();
	int K = m_di.K;
	if (m_di.K != 0)
		--K;
	std::cout << "K = " << K << std::endl;
	if (K >= 3)
	{
		table->setColumnCount(3*std::max(2, n)+2);
		table->setRowCount((K/3+1)*(1+m)+K/3);
	}
	else
	{
		table->setColumnCount((K+1)*std::max(2, n)+K);
		table->setRowCount((1+m)+1);
	}
	std::cout << "Num of rows = " << table->rowCount() << std::endl;
	std::cout << "Num of cols = " << table->columnCount() << std::endl;
	int index = 0;
	int col = 0;
	for (int k = 0; k <= K; ++k)
	{
		index = (k/3)*(1+m+1);
		col = (k % 3)*(1+n);
		std::cout << "k = " << k << std::endl;
		std::cout << "index = " << index << std::endl;
		std::cout << "col = " << col << std::endl;
		int origCol = col;
		table->setItem(index, col, new QTableWidgetItem(QString("%1[%2]").arg(prefix).arg(k)));
		++index;
		col = origCol;
		for (int i = 0; i < m; ++i)
		{
			for (int j = 0; j < n; ++j)
				table->setItem(index+i, col+j, new QTableWidgetItem(QString("%1").arg(discretes[k][i+1][j+1])));
		}
	}
	table->resizeRowsToContents();
	table->resizeColumnsToContents();
}

void MainWindow::clearMatrixInfo()
{
	m_matrix_discretes.clear();
	m_binv_discretes.clear();
	m_qinv_discretes.clear();
	m_bqinv_discretes.clear();
	m_drazin_recursive.clear();
	m_drazin_canonical.clear();
	m_drazin_skeleton.clear();
}

void MainWindow::updateMatrixInfo()
{
	m_app.applyDiffTrans(m_func_matrix, m_matrix_discretes);
	std::cout << "Got discretes\n";
	m_rank = m_app.getRank(m_matrix_discretes);
	std::cout << "r = " << m_rank << std::endl;
}

void MainWindow::updateToolsActions()
{
	m_b_inverse_action->setEnabled(m_app.isBInvertible(m_func_matrix, m_rank));
	m_q_inverse_action->setEnabled(m_app.isQInvertible(m_func_matrix, m_rank));
	m_bq_inverse_action->setEnabled(m_app.isBQInvertible(m_func_matrix, m_rank));

	m_drazin_recursive_inverse_action->setEnabled(m_app.isDrazinInvertible(m_func_matrix));
	m_drazin_skeleton_inverse_action->setEnabled(m_app.isDrazinInvertible(m_func_matrix));
	m_drazin_canonical_inverse_action->setEnabled(m_app.isDrazinInvertible(m_func_matrix));

	setChecks();
	m_restore_action->setEnabled(false);
	m_plot_action->setEnabled(false);
}

void MainWindow::updateMatrixWidget()
{
	int m = m_func_matrix.getNumRows();
	int n = m_func_matrix.getNumCols();
	m_matrix_widget->setRowCount(m+8);
	int index = 0;
	m_matrix_widget->setColumnCount(std::max(2, n));
	m_matrix_widget->setItem(index++, 0, new QTableWidgetItem("A(t)"));
	for (int i = 1; i <= m; ++i)
		for (int j = 0; j < n; ++j)
			m_matrix_widget->setItem(i, j, new QTableWidgetItem(QString::fromStdString(m_func_matrix[i][j+1]->toString())));
	
	index += m;
	++index;
	m_matrix_widget->setItem(index, 0, new QTableWidgetItem("Rank:"));
	m_matrix_widget->setItem(index, 1, new QTableWidgetItem(QString("%1").arg(m_rank)));
	++index;
	m_matrix_widget->setItem(index++, 0, new QTableWidgetItem("Diff. Trans."));
	m_matrix_widget->setItem(index, 0, new QTableWidgetItem("K:"));
	m_matrix_widget->setItem(index, 1, new QTableWidgetItem(QString("%1").arg(m_di.K)));
	++index;
	m_matrix_widget->setItem(index, 0, new QTableWidgetItem("tv:"));
	m_matrix_widget->setItem(index, 1, new QTableWidgetItem(QString("%1").arg(m_di.tv)));
	++index;
	m_matrix_widget->setItem(index, 0, new QTableWidgetItem("H:"));
	m_matrix_widget->setItem(index, 1, new QTableWidgetItem(QString("%1").arg(m_di.H)));
	m_matrix_widget->resizeRowsToContents();
	m_matrix_widget->resizeColumnsToContents();
}

void MainWindow::on_open()
{
}

void MainWindow::on_save()
{
}

void MainWindow::on_saveAs()
{
}

void MainWindow::on_close()
{
}

void MainWindow::on_quit()
{
}

void MainWindow::createViewMenuActions()
{
}

void MainWindow::createToolsMenuActions()
{
	// B, Q, BQ
	m_b_inverse_action = new QAction(QIcon(":/images/b_inv.png"), tr("(B)-Inverse"), this);
	m_b_inverse_action->setShortcut(tr("Ctrl+Shift+B"));
	m_b_inverse_action->setStatusTip("(B)-Inverse the matrix");
	connect(m_b_inverse_action, SIGNAL(triggered()), this, SLOT(on_b_inverse()));

	m_q_inverse_action = new QAction(QIcon(":/images/q_inv.png"), tr("(Q)-Inverse"), this);
	m_q_inverse_action->setShortcut(tr("Ctrl+Shift+Q"));
	m_q_inverse_action->setStatusTip("(Q)-Inverse the matrix");
	connect(m_q_inverse_action, SIGNAL(triggered()), this, SLOT(on_q_inverse()));

	m_bq_inverse_action = new QAction(QIcon(":/images/bq_inv.png"), tr("(BQ)-Inverse"), this);
	m_bq_inverse_action->setShortcut(tr("Ctrl+Alt+B"));
	m_bq_inverse_action->setStatusTip("(BQ)-Inverse the matrix");
	connect(m_bq_inverse_action, SIGNAL(triggered()), this, SLOT(on_bq_inverse()));

	m_check_gen_action = new QAction(QIcon(":/images/check1.png"), tr("Check B, Q and BQ"), this);
	m_check_gen_action->setShortcut(tr("Ctrl+C"));
	m_check_gen_action->setStatusTip("Check B, Q and BQ inverses");
	connect(m_check_gen_action, SIGNAL(triggered()), this, SLOT(on_check_gen()));

	// Drazin
	m_drazin_recursive_inverse_action = new QAction(QIcon(":/images/drazin1.png"), tr("Drazin Inverse Recursive"), this);
	m_drazin_recursive_inverse_action->setShortcut(tr("Ctrl+D"));
	m_drazin_recursive_inverse_action->setStatusTip("Drazin inverse of the matrix");
	connect(m_drazin_recursive_inverse_action, SIGNAL(triggered()), this, SLOT(on_drazin_recursive()));

	m_drazin_skeleton_inverse_action = new QAction(QIcon(":/images/drazin2.png"), tr("Drazin Inverse Skeleton"), this);
	m_drazin_skeleton_inverse_action->setShortcut(tr("Shift+D"));
	m_drazin_skeleton_inverse_action->setStatusTip("Drazin inverse of the matrix");
	connect(m_drazin_skeleton_inverse_action, SIGNAL(triggered()), this, SLOT(on_drazin_skeleton()));

	m_drazin_canonical_inverse_action = new QAction(QIcon(":/images/drazin3.png"), tr("Drazin Inverse Canonical"), this);
	m_drazin_canonical_inverse_action->setShortcut(tr("Shift+Ctrl+D"));
	m_drazin_canonical_inverse_action->setStatusTip("Drazin inverse of the matrix");
	connect(m_drazin_canonical_inverse_action, SIGNAL(triggered()), this, SLOT(on_drazin_canonical()));

	m_check_drazin_action = new QAction(QIcon(":/images/check2.png"), tr("Check Drazin Inverse"), this);
	m_check_drazin_action->setShortcut(tr("Ctrl+Shift+C"));
	m_check_drazin_action->setStatusTip("Check Drazin Inverse");
	connect(m_check_drazin_action, SIGNAL(triggered()), this, SLOT(on_check_drazin()));

	// Diff trans
	m_restore_action = new QAction(QIcon(":/images/restore.png"), tr("Restore the original"), this);
	m_restore_action->setShortcut(tr("Ctrl+R"));
	m_restore_action->setStatusTip("Restore the original");
	connect(m_restore_action, SIGNAL(triggered()), this, SLOT(on_restore()));

	m_compare_action = new QAction(QIcon(":/images/compare.png"), tr("Compare with the original"), this);
	m_compare_action->setShortcut(tr("Shift+C"));
	m_compare_action->setStatusTip("compare with the original");
	connect(m_compare_action, SIGNAL(triggered()), this, SLOT(on_compare()));

	// Other
	m_plot_action = new QAction(QIcon(":/images/plot.png"), tr("Plot"), this);
	m_plot_action->setShortcut(tr("Shift+P"));
	m_plot_action->setStatusTip("Plot the result");
	connect(m_plot_action, SIGNAL(triggered()), this, SLOT(on_plot()));

	m_options_action = new QAction(QIcon(":/images/options.png"), tr("&Options..."), this);
	m_options_action->setShortcut(tr("Shift+O"));
	m_options_action->setStatusTip("Display Options");
	connect(m_options_action, SIGNAL(triggered()), this, SLOT(on_options()));

}

void MainWindow::on_b_inverse()
{
	m_app.getBInverse(m_matrix_discretes, m_rank, m_binv_discretes);
	QTableWidget* table = new QTableWidget();
	int index = m_central_widget->addTab(table, "(B)-Inverse");
	m_central_widget->setCurrentIndex(index);
	updateTableWithDiscretes(table, m_binv_discretes, "Binv");
	m_restore_action->setEnabled(true);
}

void MainWindow::on_q_inverse()
{
	m_app.getQInverse(m_matrix_discretes, m_rank, m_qinv_discretes);
	QTableWidget* table = new QTableWidget();
	int index = m_central_widget->addTab(table, "(Q)-Inverse");
	m_central_widget->setCurrentIndex(index);
	updateTableWithDiscretes(table, m_qinv_discretes, "Qinv");
	m_restore_action->setEnabled(true);
}

void MainWindow::on_bq_inverse()
{
	m_app.getBQInverse(m_matrix_discretes, m_bqinv_discretes);
	QTableWidget* table = new QTableWidget();
	int index = m_central_widget->addTab(table, "(B, Q)-Inverse");
	m_central_widget->setCurrentIndex(index);
	updateTableWithDiscretes(table, m_bqinv_discretes, "BQinv");
	m_restore_action->setEnabled(true);
}

void MainWindow::on_check_gen()
{
	bool ok = true;
	if (m_binv_discretes.size() != 0)
	{
		if(!m_app.checkB_Q_BQ_Inverse(m_matrix_discretes, m_binv_discretes, m_di.K>0? m_di.K-1:0))
		{
			ok = false;
			QMessageBox::critical(this, QString("Check B-inverse"), tr("Main equality isn't held"), QMessageBox::Ok);

		}
	}

	if (m_qinv_discretes.size() != 0)
	{
		if(!m_app.checkB_Q_BQ_Inverse(m_matrix_discretes, m_qinv_discretes, m_di.K>0? m_di.K-1:0))
		{
			ok = false;
			QMessageBox::critical(this, QString("Check Q-inverse"), tr("Main equality isn't held"), QMessageBox::Ok);

		}
	}

	if (m_bqinv_discretes.size() != 0)
	{
		if(!m_app.checkB_Q_BQ_Inverse(m_matrix_discretes, m_bqinv_discretes, m_di.K>0? m_di.K-1:0))
		{
			ok = false;
			QMessageBox::critical(this, QString("Check (B,Q)-inverse"), tr("Main equality isn't held"), QMessageBox::Ok);

		}
	}

	if (ok)
	{
		QMessageBox::information(this, QString("Checks"), tr("Checks for all calculated generalized inverses passed successfully"), QMessageBox::Ok);
	}
}


void MainWindow::on_drazin_recursive()
{
	m_app.getDrazinInverseRecursive(m_matrix_discretes, m_drazin_recursive);
	QTableWidget* table = new QTableWidget();
	int index = m_central_widget->addTab(table, "Drazin Inverse Recursive");
	m_central_widget->setCurrentIndex(index);
	updateTableWithDiscretes(table, m_drazin_recursive, "AD");
	m_restore_action->setEnabled(true);
}

void MainWindow::on_drazin_skeleton()
{
	m_app.getDrazinInverseRecursive(m_matrix_discretes, m_drazin_skeleton);
	QTableWidget* table = new QTableWidget();
	int index = m_central_widget->addTab(table, "Drazin Inverse Skeleton");
	m_central_widget->setCurrentIndex(index);
	updateTableWithDiscretes(table, m_drazin_skeleton, "AD");
	m_restore_action->setEnabled(true);
}

void MainWindow::on_drazin_canonical()
{
	m_app.getDrazinInverseRecursive(m_matrix_discretes, m_drazin_canonical);
	QTableWidget* table = new QTableWidget();
	int index = m_central_widget->addTab(table, "Drazin Inverse Canonical");
	m_central_widget->setCurrentIndex(index);
	updateTableWithDiscretes(table, m_drazin_canonical, "AD");
	m_restore_action->setEnabled(true);
}

void MainWindow::on_check_drazin()
{
	bool ok = true;
	if (m_drazin_recursive.size() != 0)
	{
		if(!m_app.checkDrazinInverse(m_matrix_discretes, m_drazin_recursive, m_di.K>0? m_di.K-1:0))
		{
			ok = false;
			QMessageBox::critical(this, QString("Check Drazin-recursive"), tr("Main equalities aren't held"), QMessageBox::Ok);

		}
	}

	if (m_drazin_skeleton.size() != 0)
	{
		if(!m_app.checkDrazinInverse(m_matrix_discretes, m_drazin_skeleton, m_di.K>0? m_di.K-1:0))
		{
			ok = false;
			QMessageBox::critical(this, QString("Check Drazin-skeleton"), tr("Main equalities aren't held"), QMessageBox::Ok);

		}
	}

	if (m_drazin_canonical.size() != 0)
	{
		if(!m_app.checkDrazinInverse(m_matrix_discretes, m_drazin_canonical, m_di.K>0? m_di.K-1:0))
		{
			ok = false;
			QMessageBox::critical(this, QString("Check Drazin-canonical"), tr("Main equalities aren't held"), QMessageBox::Ok);

		}
	}

	if (ok)
	{
		QMessageBox::information(this, QString("Checks"), tr("Checks for all calculated Drazn inverses passed successfully"), QMessageBox::Ok);
	}
}

void MainWindow::on_restore()
{
	RestoreWizard wizard(this,
		(m_binv_discretes.size() != 0),
		(m_qinv_discretes.size() != 0),
		(m_bqinv_discretes.size() != 0),
		(m_drazin_recursive.size() != 0),
		(m_drazin_skeleton.size() != 0),
		(m_drazin_canonical.size() != 0)
		);
	if (!wizard.exec())
		return;
	is_binv = wizard.field("restore.b").toBool();
	is_qinv = wizard.field("restore.q").toBool();
	is_bqinv = wizard.field("restore.bq").toBool();
	is_drazin_recursive = wizard.field("restore.dr").toBool();
	is_drazin_skeleton = wizard.field("restore.ds").toBool();
	is_drazin_canonical = wizard.field("restore.dc").toBool();
	bool taylor_single = wizard.field("restore.taylor.single").toBool();
	bool pade = wizard.field("restore.pade").toBool();
	if (pade)
	{
		m_pade_m = wizard.field("restore.pade.m").toInt();
		m_pade_n = wizard.field("restore.pade.n").toInt();
	}

	if (is_binv)
		restore(m_binv_discretes, pade ? eRestorePade : (taylor_single ? eRestoreTaylorSingle : eRestoreTaylorMulti), "(B)-Inverse Restored", m_binv_restore);

	if (is_qinv)
		restore(m_qinv_discretes, pade ? eRestorePade : (taylor_single ? eRestoreTaylorSingle : eRestoreTaylorMulti), "(Q)-Inverse Restored", m_qinv_restore);

	if (is_bqinv)
		restore(m_bqinv_discretes, pade ? eRestorePade : (taylor_single ? eRestoreTaylorSingle : eRestoreTaylorMulti), "(B,Q)-Inverse Restored", m_bqinv_restore);

	if (is_drazin_recursive)
		restore(m_drazin_recursive, pade ? eRestorePade : (taylor_single ? eRestoreTaylorSingle : eRestoreTaylorMulti), "Drazin Recursive Inverse Restored", m_drazin_recursive_restore);

	if (is_drazin_skeleton)
		restore(m_drazin_skeleton, pade ? eRestorePade : (taylor_single ? eRestoreTaylorSingle : eRestoreTaylorMulti), "Drazin Skeleton Inverse Restored", m_drazin_skeleton_restore);

	if (is_drazin_canonical)
		restore(m_drazin_canonical, pade ? eRestorePade : (taylor_single ? eRestoreTaylorSingle : eRestoreTaylorMulti), "Drazin Canonical Inverse Restored", m_drazin_canonical_restore);
	
	m_plot_action->setEnabled(true);
}

void MainWindow::restore(const tDiscretes& discretes, eRestType rt, const QString& title, tMatrix& origMatrix)
{
	switch (rt)
	{
	case eRestorePade:
		return restorePade(discretes, title, origMatrix);
	case eRestoreTaylorSingle:
		return restoreTaylorSingle(discretes, title, origMatrix);
	case eRestoreTaylorMulti:
		return restoreTaylorMulti(discretes, title, origMatrix);
	default:
		assert (! "Unknown restoration type");
	}
}

void MainWindow::restorePade(const tDiscretes& discretes, const QString& title, tMatrix& origMatrix)
{
}

void MainWindow::restoreTaylorSingle(const tDiscretes& discretes, const QString& title, tMatrix& origMatrix)
{
	assert(m_di.K >= 0);
	origMatrix.resize(discretes[0].getNumCols(), discretes[0].getNumRows(), tFunctionPtr(nullptr));
	m_app.restoreTaylorSingle(discretes, origMatrix, m_di.K > 0 ? m_di.K-1 : 0);
	QTableWidget* table = new QTableWidget();
	int index = m_central_widget->addTab(table, title);
	m_central_widget->setCurrentIndex(index);
	updateTableWithFuncMatrix(table, origMatrix);
}

void MainWindow::restoreTaylorMulti(const tDiscretes& discretes, const QString& title, tMatrix& origMatrix)
{
}

void MainWindow::updateTableWithFuncMatrix(QTableWidget* table, const tMatrix& mat)
{
	int m = mat.getNumRows();
	int n = mat.getNumCols();
	table->setRowCount(m);
	table->setColumnCount(n);
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			table->setItem(i, j,
				new QTableWidgetItem(QString::fromStdString(mat[i+1][j+1]->toString())));
	table->resizeRowsToContents();
	table->resizeColumnsToContents();
}

void MainWindow::on_compare()
{
}

void MainWindow::on_plot()
{
	PlotWizard wizard(this, is_binv, is_qinv, is_bqinv,
				is_drazin_recursive, is_drazin_skeleton, is_drazin_canonical);
	if (!wizard.exec())
		return;
	double xbeg = wizard.field("plot.xbeg").toDouble();
	double xend = wizard.field("plot.xend").toDouble();
	if (wizard.field("plot.b").toBool())
	{
		plot(m_binv_restore, xbeg, xend, "(Q)-inverse Restored");
	}
	if (wizard.field("plot.q").toBool())
	{
		plot(m_qinv_restore, xbeg, xend, "(B)-inverse Restored");
	}
	if (wizard.field("plot.bq").toBool())
	{
		plot(m_bqinv_restore, xbeg, xend, "(B, Q)-Inverse Restored");
	}
	if (wizard.field("plot.dr").toBool())
	{
		plot(m_drazin_recursive_restore, xbeg, xend, "Drazin Recursive Inverse Restored");
	}
	if (wizard.field("plot.ds").toBool())
	{
		plot(m_drazin_skeleton_restore, xbeg, xend, "Drazin Skeleton Inverse Restored");
	}
	if (wizard.field("plot.dc").toBool())
	{
		plot(m_drazin_canonical_restore, xbeg, xend, "Drazin Canonical Inverse Restored");
	}
}

void MainWindow::plot(const tMatrix& theMatrix, double begin, double end, const QString& title)
{
	std::cout << "Creating the plotter\n";
	Plotter* plotter = new Plotter(0, theMatrix, begin, end, title);
	std::cout << "Created the plotter\n";
	plotter->show();
	std::cout << "Called show()\n";
}

void MainWindow::on_options()
{
}

void MainWindow::createHelpMenuActions()
{
	m_about_action = new QAction(tr("&About"), this);
	m_about_action->setStatusTip(tr("Show an information about this product"));
	connect(m_about_action, SIGNAL(triggered()), this, SLOT(on_about()));
	
	m_about_qt_action = new QAction(tr("About Qt"), this);
	m_about_qt_action->setStatusTip(tr("Show information about Qt"));
	connect(m_about_qt_action, SIGNAL(triggered()), qApp, SLOT(aboutQt()));
}

void MainWindow::on_about()
{
	QMessageBox::about(this, QString("About %1").arg(sToolName),
					QString("<h2>%1</h2>"
							"<p>Copyright &copy; A.H.A. Inc."
							"<p>This is an application that allows the users"
							" to maniuplate generalized (B, Q, and BQ) and Drazin inverse "
							" matrices.\n The whole concept is based on differential trans"
							"formations.\nIt was developed"
							" using C++ with it's toolkit Qt.\n\nDeveloped by:"
							"Hamlet Aslanyan      A.H.A. Inc.").arg(sTitle));

}

void MainWindow::disableAllActions()
{
	// File
	m_new_action->setEnabled(false);
	m_open_action->setEnabled(false);
	m_save_action->setEnabled(false);
	m_save_as_action->setEnabled(false);
	m_close_action->setEnabled(false);
	m_quit_action->setEnabled(false);

	// view

	// Tools
	m_b_inverse_action->setEnabled(false);
	m_q_inverse_action->setEnabled(false);
	m_bq_inverse_action->setEnabled(false);
	m_drazin_recursive_inverse_action->setEnabled(false);
	m_drazin_skeleton_inverse_action->setEnabled(false);
	m_drazin_canonical_inverse_action->setEnabled(false);
	m_restore_action->setEnabled(false);
	m_compare_action->setEnabled(false);
	m_plot_action->setEnabled(false);
	m_options_action->setEnabled(false);
	setChecks();

	m_about_action->setEnabled(false);
	m_about_qt_action->setEnabled(false);
}

void MainWindow::setChecks()
{
	m_check_gen_action->setEnabled(
		m_b_inverse_action->isEnabled() ||
		m_q_inverse_action->isEnabled() ||
		m_bq_inverse_action->isEnabled());

	m_check_drazin_action->setEnabled(
		m_drazin_recursive_inverse_action->isEnabled() ||
		m_drazin_skeleton_inverse_action->isEnabled() ||
		m_drazin_canonical_inverse_action->isEnabled());
}

void MainWindow::enableBasicActions()
{
	// File
	m_new_action->setEnabled(true);
	m_quit_action->setEnabled(true);

	// Help
	m_about_action->setEnabled(true);
	m_about_qt_action->setEnabled(true);
}

void MainWindow::createMenus()
{
	// File menu
	m_file_menu = this->menuBar()->addMenu(tr("&File"));
	m_file_menu->addAction(m_new_action);
	m_file_menu->addAction(m_open_action);
	m_file_menu->addAction(m_save_action);
	m_file_menu->addAction(m_save_as_action);
	m_file_menu->addAction(m_close_action);
	m_file_menu->addSeparator();
	m_file_menu->addAction(m_quit_action);

	// View Menu
	m_view_menu = this->menuBar()->addMenu(tr("&View"));

	// Tools Menu
	m_tools_menu = this->menuBar()->addMenu(tr("&Tools"));
	QMenu* sub_menu = m_tools_menu->addMenu("&Generalized Inverse");
	sub_menu->addAction(m_b_inverse_action);
	sub_menu->addAction(m_q_inverse_action);
	sub_menu->addAction(m_bq_inverse_action);
	sub_menu->addAction(m_check_gen_action);

	sub_menu = m_tools_menu->addMenu("&Drazin Inverse");
	sub_menu->addAction(m_drazin_recursive_inverse_action);
	sub_menu->addAction(m_drazin_skeleton_inverse_action);
	sub_menu->addAction(m_drazin_canonical_inverse_action);
	sub_menu->addAction(m_check_drazin_action);

	sub_menu = m_tools_menu->addMenu("D&ifferential Transformations");
	sub_menu->addAction(m_restore_action);
	sub_menu->addAction(m_compare_action);

	m_tools_menu->addSeparator();
	m_tools_menu->addAction(m_plot_action);
	m_tools_menu->addAction(m_options_action);

	// Help menu
	m_help_menu = this->menuBar()->addMenu(tr("&Help"));
	m_help_menu->addAction(m_about_action);
	m_help_menu->addAction(m_about_qt_action);
}

void MainWindow::createToolBar()
{
	// File toolbar
	m_file_tool = this->addToolBar(tr("File"));
	m_file_tool->addAction(m_new_action);
	m_file_tool->addAction(m_open_action);
	m_file_tool->addAction(m_save_action);
	m_file_tool->addAction(m_save_as_action);

	// Tools toolbar
	m_tools_tool = this->addToolBar(tr("Tools"));
	m_tools_tool->addAction(m_b_inverse_action);
	m_tools_tool->addAction(m_q_inverse_action);
	m_tools_tool->addAction(m_bq_inverse_action);
	m_tools_tool->addAction(m_check_gen_action);
	m_tools_tool->addSeparator();
	m_tools_tool->addAction(m_drazin_recursive_inverse_action);
	m_tools_tool->addAction(m_drazin_skeleton_inverse_action);
	m_tools_tool->addAction(m_drazin_canonical_inverse_action);
	m_tools_tool->addAction(m_check_drazin_action);
	m_tools_tool->addSeparator();
	m_tools_tool->addAction(m_restore_action);
	m_tools_tool->addAction(m_compare_action);
	m_tools_tool->addSeparator();
	m_tools_tool->addAction(m_plot_action);
}

