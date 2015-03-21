#include "main_window.h"
#include "wizards.h"

#include <QtGui>

#include <iostream>

static const char* sToolName = "GIMT";
static const char* sTitle = "Generalized Inverse Matrices Tool";

MainWindow::MainWindow(QWidget* p)
	: QMainWindow(p)
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
}

void MainWindow::createWidgets()
{
	m_matrix_widget = new QWidget();
	m_central_widget = new QWidget();
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
	std::cout << "numRows = " << numRows << std::endl;
	std::cout << "numCols = " << numCols << std::endl;
	std::cout << "Matrix =:\n";
	for (int i = 1; i <= numRows; ++i)
	{
		for (int j = 1; j <= numCols; ++j)
			std::cout << wizard.getMatrix()->item(i, j)->text().toStdString() << " ";

		std::cout << std::endl;
	}
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
}

void MainWindow::on_q_inverse()
{
}

void MainWindow::on_bq_inverse()
{
}

void MainWindow::on_check_gen()
{
}

void MainWindow::on_drazin_recursive()
{
}

void MainWindow::on_drazin_skeleton()
{
}

void MainWindow::on_drazin_canonical()
{
}

void MainWindow::on_check_drazin()
{
}

void MainWindow::on_restore()
{
}

void MainWindow::on_compare()
{
}

void MainWindow::on_plot()
{
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

