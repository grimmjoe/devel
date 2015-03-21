#include <QApplication>
#include "main_window.h"

int main(int argc, char* argv[])
{
	QApplication app(argc, argv);
	MainWindow* mw = new MainWindow();
	mw->show();

	app.exec();
	return 0;
}
