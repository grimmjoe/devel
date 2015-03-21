#include "wizards.h"
#include <QtGui>

#include <iostream>

NewWizard::NewWizard(QWidget* p)
	: QWizard(p)
{
	this->setPage(Page_New, new NewPage);
	this->setPage(Page_Matrix, new MatrixPage);
	this->setStartId(Page_New);

	#ifndef Q_WS_MAC
	setWizardStyle(ModernStyle);
	#endif  
	setOption(HaveHelpButton, true);
	setPixmap(QWizard::LogoPixmap, QPixmap(":/images/logo.png"));
	connect(this, SIGNAL(helpRequested()), this, SLOT(showHelp()));
	setWindowTitle(tr("New Session Wizard"));
}

QTableWidget* NewWizard::getMatrix() const
{
	return m_table;
}

void NewWizard::setMatrix(QTableWidget* t) const
{
	m_table = t;
}

void NewWizard::showHelp()
{
	static QString lastMessage;
	QString message;

	switch(currentId()) {
		case Page_New:
				message = tr("You must specify the dimensions of the matrix and "
							 "the configuration of differential transformations");
				break;
		case Page_Matrix:
				message = tr("You must specify each element of the matrix based on "
								" the agreed syntax");
				break;
		default:
				message = tr("This help wasn't a help at all.");
				break;
	}

	if (lastMessage == message) {
			message = tr("Sorry, I've already helped you as much as I could. "
							"Maybe you should ask a human.?");
	}

	QMessageBox::information(this, tr("New Session Wizard Help"), message);
	lastMessage = message;
}

NewPage::NewPage(QWidget* p)
	: QWizardPage(p)
{
	setTitle(tr("New Session"));
    setPixmap(QWizard::WatermarkPixmap, QPixmap(":/images/watermark.jpg"));

	qApp->setStyleSheet("QGroupBox { border: 1px solid gray; border-radius: 9px; margin-top: 0.5em; } QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 3px 0 3px; }");
	
	QGroupBox* matDim = new QGroupBox(tr("&Matrix Dimension"));
	QHBoxLayout* md = new QHBoxLayout();
	QLabel* numRows = new QLabel(tr("Number of &Rows:"));
	QSpinBox* r = new QSpinBox();
	r->setRange(1, 1000);
	numRows->setBuddy(r);
	QLabel* numCols = new QLabel(tr("Number of &Columns:"));
	QSpinBox* c = new QSpinBox();
	c->setRange(1, 1000);
	numCols->setBuddy(c);
	md->addWidget(numRows);
	md->addWidget(r);
	md->addWidget(numCols);
	md->addWidget(c);
	md->addStretch();
	matDim->setLayout(md);

	registerField("new.numRows*", r);
	registerField("new.numCols*", c);

	QGroupBox* diffTrans = new QGroupBox("&Differential Transformations");
	QGridLayout* gl = new QGridLayout();
	QLabel* numDiscs = new QLabel(tr("Number of &Discretes:"));
	QSpinBox* nd = new QSpinBox();
	nd->setRange(0, 1000);
	numDiscs->setBuddy(nd);
	gl->addWidget(numDiscs, 1, 1, Qt::AlignLeft);
	gl->addWidget(nd, 1, 2, Qt::AlignRight);

	QLabel* centre = new QLabel(tr("Center of &approximation:"));
	QLineEdit* ca = new QLineEdit();
	ca->setAlignment(Qt::AlignRight);
	ca->setValidator(new QDoubleValidator());
	centre->setBuddy(ca);
	gl->addWidget(centre, 2, 1, Qt::AlignLeft);
	gl->addWidget(ca, 2, 2, Qt::AlignRight);

	QLabel* H = new QLabel(tr("&Scale:"));
	QLineEdit* s = new QLineEdit();
	s->setAlignment(Qt::AlignRight);
	s->setValidator(new QDoubleValidator());
	H->setBuddy(s);
	gl->addWidget(H, 3, 1, Qt::AlignLeft);
	gl->addWidget(s, 3, 2, Qt::AlignRight);

	diffTrans->setLayout(gl);

	registerField("new.numDiscs*", nd);
	registerField("new.center*", ca);
	registerField("new.scale*", ca);

	QVBoxLayout* theLayout = new QVBoxLayout();
	theLayout->addWidget(matDim);
	theLayout->addWidget(diffTrans);
	this->setLayout(theLayout);
}

int NewPage::nextId() const
{
	return NewWizard::Page_Matrix;
}

MatrixPage::MatrixPage(QWidget* p)
	: QWizardPage(p)
{
	setTitle(tr("New Session"));
    setPixmap(QWizard::WatermarkPixmap, QPixmap(":/images/watermark.jpg"));
}

void MatrixPage::initializePage()
{

	int m = this->wizard()->field("new.numRows").toInt();
	int n = this->wizard()->field("new.numCols").toInt();
	std::cout << "m = " << m << ", n = " << n << std::endl;
	QTableWidget* tw = new QTableWidget();
	tw->setRowCount(m);
	tw->setColumnCount(m);
	QVBoxLayout* theLayout = new QVBoxLayout();
	theLayout->addWidget(tw);
	this->setLayout(theLayout);
	NewWizard* theW = dynamic_cast<NewWizard*>(this->wizard());
	if (theW != 0)
		theW->setMatrix(tw);

	return QWizardPage::initializePage();
}

int MatrixPage::nextId() const
{
	return -1;
}






