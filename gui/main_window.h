#ifndef MAIN_WINDOW_H__
#define MAIN_WINDOW_H__

#include <QMainWindow>
#include <QTableWidget>
#include <QTabWidget>
#include "core/inc/apparatus.h"
#include "core/inc/matrix.h"
#include "core/inc/parser.h"
#include <string>

typedef core::parser<double> tParser;
typedef core::apparatus<double> tApp;
typedef typename tApp::tMatrixDiscrete tDiscrete;
typedef typename tApp::tMatrixDiscretes tDiscretes;
typedef typename tApp::tFunction tFunction;
typedef typename tApp::tFunctionPtr tFunctionPtr;
typedef typename tApp::tFuncMatrix tMatrix;
typedef typename tApp::tDiffInfo tDiffInfo;

class MainWindow : public QMainWindow
{
	Q_OBJECT
public:
	MainWindow(QWidget* p = 0);

	enum eRestType
	{
		eRestoreTaylorSingle = 0,
		eRestoreTaylorMulti = 1,
		eRestorePade = 2
	};

protected:
	// Creation

	void createWidgets();

	void createActions();
	void createFileMenuActions();
	void createViewMenuActions();
	void createToolsMenuActions();
	void createHelpMenuActions();

	void createMenus();

	void createToolBar();

protected:
	// Manipulation
	void disableAllActions();
	void setChecks();
	void enableBasicActions();

	void clearMatrixInfo();
	void clearCentralWidget();
	void updateCentralWidget();
	void updateTableWithDiscretes(QTableWidget* table, tDiscretes& discretes, const QString& prefix);
	void updateTableWithFuncMatrix(QTableWidget* table, const tMatrix& mat);
	void updateMatrixInfo();
	void updateToolsActions();
	void updateMatrixWidget();
private:
	// Widgets
	QTableWidget* m_matrix_widget;
	QTabWidget* m_central_widget;
	QDockWidget* m_left_dock;

private:
	// Actions
	QAction* m_new_action;
	QAction* m_open_action;
	QAction* m_save_action;
	QAction* m_save_as_action;
	QAction* m_close_action;
	QAction* m_quit_action;

	// Tools
	QAction* m_b_inverse_action;
	QAction* m_q_inverse_action;
	QAction* m_bq_inverse_action;
	QAction* m_check_gen_action;

	QAction* m_drazin_recursive_inverse_action;
	QAction* m_drazin_skeleton_inverse_action;
	QAction* m_drazin_canonical_inverse_action;
	QAction* m_check_drazin_action;

	QAction* m_restore_action;
	QAction* m_compare_action;

	QAction* m_plot_action;
	QAction* m_options_action;

	// Help
	QAction* m_about_action;
	QAction* m_about_qt_action;

private:
	// Menu
	QMenu* m_file_menu;
	QMenu* m_view_menu;
	QMenu* m_tools_menu;
	QMenu* m_help_menu;

	// ToolBar
	QToolBar* m_file_tool;
	QToolBar* m_tools_tool;

protected:
	tDiffInfo m_di;
	tApp m_app;
	tMatrix m_func_matrix;
	tMatrix m_binv_restore;
	tMatrix m_qinv_restore;
	tMatrix m_bqinv_restore;
	tMatrix m_drazin_recursive_restore;
	tMatrix m_drazin_skeleton_restore;
	tMatrix m_drazin_canonical_restore;
	int m_rank;
	tDiscretes m_matrix_discretes;
	tDiscretes m_binv_discretes;
	tDiscretes m_qinv_discretes;
	tDiscretes m_bqinv_discretes;

	tDiscretes m_drazin_recursive;
	tDiscretes m_drazin_skeleton;
	tDiscretes m_drazin_canonical;

	double m_epsilon;
	int m_pade_m;
	int m_pade_n;

	bool is_binv;
	bool is_qinv;
	bool is_bqinv;
	bool is_drazin_recursive;
	bool is_drazin_skeleton;
	bool is_drazin_canonical;

protected:
	void restore(const tDiscretes&, eRestType, const QString&, tMatrix& origMatrix);
	void restorePade(const tDiscretes&, const QString&, tMatrix& origMatrix);
	void restoreTaylorSingle(const tDiscretes&, const QString&, tMatrix& origMatrix);
	void restoreTaylorMulti(const tDiscretes&, const QString&, tMatrix& origMatrix);

	void plot(const tMatrix& theMatrix, double begin, double end, const QString& title);

public slots:
	// slots
	void on_new();
	void on_open();
	void on_save();
	void on_saveAs();
	void on_close();
	void on_quit();

	// tools
	void on_b_inverse();
	void on_q_inverse();
	void on_bq_inverse();
	void on_check_gen();

	void on_drazin_recursive();
	void on_drazin_skeleton();
	void on_drazin_canonical();
	void on_check_drazin();

	void on_restore();
	void on_compare();

	void on_plot();
	void on_options();

	// help
	void on_about();
};

#endif // MAIN_WINDOW_H__
