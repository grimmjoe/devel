#ifndef MAIN_WINDOW_H__
#define MAIN_WINDOW_H__

#include <QMainWindow>

class MainWindow : public QMainWindow
{
	Q_OBJECT
public:
	MainWindow(QWidget* p = 0);

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
private:
	// Widgets
	QWidget* m_matrix_widget;
	QWidget* m_central_widget;
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
