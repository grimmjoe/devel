#include <QtGui>
#include <QApplication>
#include <QWidget>

int main(int argc, char* argv[])
{
	QApplication app(argc, argv);
	QWidget* w = new QWidget();
	w->show();

	app.exec();
	return 0;
}
