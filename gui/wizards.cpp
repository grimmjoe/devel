#include "wizards.h"
#include <QtGui>

#include <iostream>

NewWizard::NewWizard(QWidget* p)
	: QWizard(p)
{
	this->setPage(Page_New, new NewPage(this));
	m_matrix_page = new MatrixPage(this);
	this->setPage(Page_Matrix, m_matrix_page);
	this->setStartId(Page_New);

	#ifndef Q_WS_MAC
	setWizardStyle(ModernStyle);
	#endif  
	setOption(HaveHelpButton, true);
	setPixmap(QWizard::LogoPixmap, QPixmap(":/images/logo.png"));
	connect(this, SIGNAL(helpRequested()), this, SLOT(showHelp()));
	setWindowTitle(tr("New Session Wizard"));
	vec.clear();
}

std::vector<std::vector<QString> > NewWizard::getMatrix() const
{
	return vec;
}

void NewWizard::accept()
{
	std::cout << "Accept\n";
	vec.clear();
	if (m_matrix_page)
	{
		vec = m_matrix_page->getValues();
	}
	std::cout << "Accepted\n";
	return QWizard::accept();
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

	registerField("new.numRows", r);
	registerField("new.numCols", c);

	QGroupBox* diffTrans = new QGroupBox("&Differential Transformations");
	QGridLayout* gl = new QGridLayout();
	QLabel* numDiscs = new QLabel(tr("Number of &Discretes:"));
	QSpinBox* nd = new QSpinBox();
	nd->setRange(0, 1000);
	nd->setValue(1);
	numDiscs->setBuddy(nd);
	gl->addWidget(numDiscs, 1, 1, Qt::AlignLeft);
	gl->addWidget(nd, 1, 2, Qt::AlignRight);

	QLabel* centre = new QLabel(tr("Center of &approximation:"));
	QLineEdit* ca = new QLineEdit();
	ca->setAlignment(Qt::AlignRight);
	ca->setValidator(new QDoubleValidator());
	ca->setText("1");
	centre->setBuddy(ca);
	gl->addWidget(centre, 2, 1, Qt::AlignLeft);
	gl->addWidget(ca, 2, 2, Qt::AlignRight);

	QLabel* H = new QLabel(tr("&Scale:"));
	QLineEdit* s = new QLineEdit();
	s->setAlignment(Qt::AlignRight);
	s->setValidator(new QDoubleValidator());
	H->setBuddy(s);
	s->setText("1");
	gl->addWidget(H, 3, 1, Qt::AlignLeft);
	gl->addWidget(s, 3, 2, Qt::AlignRight);

	diffTrans->setLayout(gl);

	registerField("new.numDiscs", nd);
	registerField("new.center", ca);
	registerField("new.scale", s);

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
	m_table = 0;
	vec.clear();
}

void MatrixPage::initializePage()
{

	int m = this->wizard()->field("new.numRows").toInt();
	int n = this->wizard()->field("new.numCols").toInt();
	std::cout << "m = " << m << ", n = " << n << std::endl;
	if (m_table == 0)
	{
		m_table = new QTableWidget(this);
		connect(m_table, SIGNAL(cellChanged(int, int)), this, SLOT(on_cellChanged(int, int)));
		m_table->setRowCount(m);
		m_table->setColumnCount(n);
		vec.resize(m, std::vector<QString>(n, ""));
		//m_table->setRowHidden(m, true);
		//m_table->setColumnHidden(n, true);
		QVBoxLayout* theLayout = new QVBoxLayout();
		theLayout->addWidget(m_table);
		this->setLayout(theLayout);
	}
	else
	{
		m_table->setRowCount(m);
		m_table->setColumnCount(n);
		vec.resize(m, std::vector<QString>(n, ""));
		//m_table->setRowHidden(m, true);
		//m_table->setColumnHidden(n, true);
	}

	return QWizardPage::initializePage();
}

void MatrixPage::on_cellChanged(int r, int c)
{
	std::cout << m_table->item(r, c)->text().toStdString() << std::endl;
	//this->setField(QString("new.table[%1, %2]").arg(r).arg(c), m_table->item(r, c)->text());
	std::cout << "r = " << r << ", c = " << c << std::endl;
	vec[r][c] = m_table->item(r, c)->text();
}

QTableWidget* MatrixPage::getTable()
{
	return m_table;
}

std::vector<std::vector<QString> > MatrixPage::getValues()
{
	return vec;
}

int MatrixPage::nextId() const
{
	return -1;
}

RestoreWizard::RestoreWizard(QWidget* p, bool _b, bool _q, bool _bq, bool _dr, bool _ds, bool _dc)
	: QWizard(p)
{
	b = _b;
	q = _q;
	bq = _bq;
	dr = _dr;
	ds = _ds;
	dc = _dc;
	PageRestore* pr = new PageRestore(this, b, q, bq, dr, ds, dc);
	this->setPage(Page_Restore, pr);
	this->setPage(Page_Pade, new PagePade(this));
	this->setPage(Page_Conclusion, new PageConclusion(this));
	this->setStartId(Page_Restore);

	#ifndef Q_WS_MAC
	setWizardStyle(ModernStyle);
	#endif  
	setOption(HaveHelpButton, true);
	setPixmap(QWizard::LogoPixmap, QPixmap(":/images/logo.png"));
	connect(this, SIGNAL(helpRequested()), this, SLOT(showHelp()));
	setWindowTitle(tr("Restoration Wizard"));
}

void RestoreWizard::showHelp()
{
	static QString lastMessage;
	QString message;

	switch(currentId()) {
		case Page_Restore:
				message = tr("You must choose which discretes to restore and how");
				break;
		case Page_Pade:
				message = tr("You must specify values for Pade restoration");
				break;
		default:
				message = tr("This help wasn't a help at all.");
				break;
	}

	if (lastMessage == message) {
			message = tr("Sorry, I've already helped you as much as I could. "
							"Maybe you should ask a human.?");
	}

	QMessageBox::information(this, tr("Restoration Wizard Help"), message);
	lastMessage = message;
}

PageRestore::PageRestore(QWidget* p, bool _b, bool _q, bool _bq, bool _dr, bool _ds, bool _dc)
	: QWizardPage(p)
{
	b = _b;
	q = _q;
	bq = _bq;
	dr = _dr;
	ds = _ds;
	dc = _dc;
	setTitle(tr("Restoration Page"));
    setPixmap(QWizard::WatermarkPixmap, QPixmap(":/images/watermark.jpg"));

	qApp->setStyleSheet("QGroupBox { border: 1px solid gray; border-radius: 9px; margin-top: 0.5em; } QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 3px 0 3px; }");

	QCheckBox* bc = new QCheckBox("(B)-Inverse");
	QCheckBox* qc = new QCheckBox("(Q)-Inverse");
	QCheckBox* bqc = new QCheckBox("(B, Q)-Inverse");
	QCheckBox* drc = new QCheckBox("Drazin Recursive Inverse");
	QCheckBox* dsc = new QCheckBox("Drazin Skeleton Inverse");
	QCheckBox* dcc = new QCheckBox("Drazin Canonical Inverse");
	bc->hide();
	qc->hide();
	bqc->hide();
	drc->hide();
	dsc->hide();
	dcc->hide();
	if (b)
	{
		bc->show();
		bc->setCheckState(Qt::Checked);
	}

	if (q)
	{
		qc->show();
		qc->setCheckState(Qt::Checked);
	}

	if (bq)
	{
		bqc->show();
		bqc->setCheckState(Qt::Checked);
	}

	if (dr)
	{
		drc->show();
		drc->setCheckState(Qt::Checked);
	}

	if (ds)
	{
		dsc->show();
		dsc->setCheckState(Qt::Checked);
	}

	if (dc)
	{
		dcc->show();
		dcc->setCheckState(Qt::Checked);
	}

	
	QGroupBox* inv = new QGroupBox(tr("&Inverses to Restore"));
	QHBoxLayout* md = new QHBoxLayout();
	QVBoxLayout* g1 = new QVBoxLayout();
	g1->addWidget(bc);
	g1->addWidget(qc);
	g1->addWidget(bqc);
	QVBoxLayout* d1 = new QVBoxLayout();
	d1->addWidget(drc);
	d1->addWidget(dsc);
	d1->addWidget(dcc);

	md->addLayout(g1);
	md->addLayout(d1);
	inv->setLayout(md);

	registerField("restore.b", bc);
	registerField("restore.q", qc);
	registerField("restore.bq", bqc);
	registerField("restore.dr", drc);
	registerField("restore.ds", dsc);
	registerField("restore.dc", dcc);

	QGroupBox* restoreType = new QGroupBox("&Type of Restoration");
	taylor_single = new QRadioButton("Single-Point Taylor");
	pade = new QRadioButton("Pade");
	taylor_single->setChecked(true);
	QVBoxLayout* vb = new QVBoxLayout();
	vb->addWidget(taylor_single);
	vb->addWidget(pade);
	restoreType->setLayout(vb);

	registerField("restore.taylor.single", taylor_single);
	registerField("restore.pade", pade);

	QVBoxLayout* theLayout = new QVBoxLayout();
	theLayout->addWidget(inv);
	theLayout->addWidget(restoreType);
	this->setLayout(theLayout);
}

int PageRestore::nextId() const
{
	if (pade->isChecked())
		return RestoreWizard::Page_Pade;
	return RestoreWizard::Page_Conclusion;
}


PagePade::PagePade(QWidget* p)
	: QWizardPage(p)
{
	setTitle(tr("Pade Restoration Input"));
    setPixmap(QWizard::WatermarkPixmap, QPixmap(":/images/watermark.jpg"));

	qApp->setStyleSheet("QGroupBox { border: 1px solid gray; border-radius: 9px; margin-top: 0.5em; } QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 3px 0 3px; }");

	QGroupBox* matDim = new QGroupBox(tr("&Pade Restoration Parameters"));
	QHBoxLayout* md = new QHBoxLayout();
	QLabel* numRows = new QLabel(tr("Degree of the &numerator(m):"));
	QSpinBox* r = new QSpinBox();
	r->setRange(1, 1000);
	numRows->setBuddy(r);
	QLabel* numCols = new QLabel(tr("Degree of the &denumerator(n):"));
	QSpinBox* c = new QSpinBox();
	c->setRange(0, 1000);
	numCols->setBuddy(c);
	md->addWidget(numRows);
	md->addWidget(r);
	md->addWidget(numCols);
	md->addWidget(c);
	md->addStretch();
	matDim->setLayout(md);

	registerField("restore.pade.m", r);
	registerField("restore.pade.n", c);

	QVBoxLayout* theLayout = new QVBoxLayout();
	theLayout->addWidget(matDim);
	this->setLayout(theLayout);
}

int PagePade::nextId() const
{
	return RestoreWizard::Page_Conclusion;
}

PageConclusion::PageConclusion(QWidget* p)
	: QWizardPage(p)
{
	setTitle(tr("Restoration Conclusion"));
    setPixmap(QWizard::WatermarkPixmap, QPixmap(":/images/watermark.jpg"));

	qApp->setStyleSheet("QGroupBox { border: 1px solid gray; border-radius: 9px; margin-top: 0.5em; } QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 3px 0 3px; }");

	QVBoxLayout* theLayout = new QVBoxLayout();
	QLabel* theLabel = new QLabel("Press Finish to exit the wizard");
	theLayout->addWidget(theLabel);
	this->setLayout(theLayout);
}

int PageConclusion::nextId() const
{
	return -1;
}

PlotWizard::PlotWizard(QWidget* p, bool _b, bool _q, bool _bq, bool _dr, bool _ds, bool _dc)
	: QWizard(p)
{
	b = _b;
	q = _q;
	bq = _bq;
	dr = _dr;
	ds = _ds;
	dc = _dc;
	PageChoose* pr = new PageChoose(this, b, q, bq, dr, ds, dc);
	this->setPage(Page_Choose, pr);
	this->setPage(Page_Range, new PageRange(this));
	this->setStartId(Page_Choose);

	#ifndef Q_WS_MAC
	setWizardStyle(ModernStyle);
	#endif  
	setOption(HaveHelpButton, true);
	setPixmap(QWizard::LogoPixmap, QPixmap(":/images/logo.png"));
	connect(this, SIGNAL(helpRequested()), this, SLOT(showHelp()));
	setWindowTitle(tr("Plot Wizard"));
}

void PlotWizard::showHelp()
{
	static QString lastMessage;
	QString message;

	switch(currentId()) {
		case Page_Choose:
				message = tr("You must choose which discretes to plot and how");
				break;
		case Page_Range:
				message = tr("You must specify values for the plotting");
				break;
		default:
				message = tr("This help wasn't a help at all.");
				break;
	}

	if (lastMessage == message) {
			message = tr("Sorry, I've already helped you as much as I could. "
							"Maybe you should ask a human.?");
	}

	QMessageBox::information(this, tr("Plot Wizard Help"), message);
	lastMessage = message;
}

PageChoose::PageChoose(QWidget* p, bool _b, bool _q, bool _bq, bool _dr, bool _ds, bool _dc)
	: QWizardPage(p)
{
	b = _b;
	q = _q;
	bq = _bq;
	dr = _dr;
	ds = _ds;
	dc = _dc;
	setTitle(tr("Plot Choice Page"));
    setPixmap(QWizard::WatermarkPixmap, QPixmap(":/images/watermark.jpg"));

	qApp->setStyleSheet("QGroupBox { border: 1px solid gray; border-radius: 9px; margin-top: 0.5em; } QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 3px 0 3px; }");

	QCheckBox* bc = new QCheckBox("(B)-Inverse");
	QCheckBox* qc = new QCheckBox("(Q)-Inverse");
	QCheckBox* bqc = new QCheckBox("(B, Q)-Inverse");
	QCheckBox* drc = new QCheckBox("Drazin Recursive Inverse");
	QCheckBox* dsc = new QCheckBox("Drazin Skeleton Inverse");
	QCheckBox* dcc = new QCheckBox("Drazin Canonical Inverse");
	bc->hide();
	qc->hide();
	bqc->hide();
	drc->hide();
	dsc->hide();
	dcc->hide();
	if (b)
	{
		bc->show();
		bc->setCheckState(Qt::Checked);
	}

	if (q)
	{
		qc->show();
		qc->setCheckState(Qt::Checked);
	}

	if (bq)
	{
		bqc->show();
		bqc->setCheckState(Qt::Checked);
	}

	if (dr)
	{
		drc->show();
		drc->setCheckState(Qt::Checked);
	}

	if (ds)
	{
		dsc->show();
		dsc->setCheckState(Qt::Checked);
	}

	if (dc)
	{
		dcc->show();
		dcc->setCheckState(Qt::Checked);
	}

	
	QGroupBox* inv = new QGroupBox(tr("&Inverses to Plot"));
	QHBoxLayout* md = new QHBoxLayout();
	QVBoxLayout* g1 = new QVBoxLayout();
	g1->addWidget(bc);
	g1->addWidget(qc);
	g1->addWidget(bqc);
	QVBoxLayout* d1 = new QVBoxLayout();
	d1->addWidget(drc);
	d1->addWidget(dsc);
	d1->addWidget(dcc);

	md->addLayout(g1);
	md->addLayout(d1);
	inv->setLayout(md);

	registerField("plot.b", bc);
	registerField("plot.q", qc);
	registerField("plot.bq", bqc);
	registerField("plot.dr", drc);
	registerField("plot.ds", dsc);
	registerField("plot.dc", dcc);

	QVBoxLayout* theLayout = new QVBoxLayout();
	theLayout->addWidget(inv);
	this->setLayout(theLayout);
}

int PageChoose::nextId() const
{
	return PlotWizard::Page_Range;
}

PageRange::PageRange(QWidget* p)
	: QWizardPage(p)
{
	setTitle(tr("Plot Range Input"));
    setPixmap(QWizard::WatermarkPixmap, QPixmap(":/images/watermark.jpg"));

	qApp->setStyleSheet("QGroupBox { border: 1px solid gray; border-radius: 9px; margin-top: 0.5em; } QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 3px 0 3px; }");

	QGroupBox* matDim = new QGroupBox(tr("&Plot Parameters"));
	QHBoxLayout* md = new QHBoxLayout();
	QLabel* numRows = new QLabel(tr("x &begin:"));
	QLineEdit* r = new QLineEdit();
	r->setAlignment(Qt::AlignRight);
	r->setValidator(new QDoubleValidator());
	r->setText("-10");
	numRows->setBuddy(r);
	QLabel* numCols = new QLabel(tr("x &end:"));
	QLineEdit* c = new QLineEdit();
	c->setAlignment(Qt::AlignRight);
	c->setValidator(new QDoubleValidator());
	c->setText("10");
	numRows->setBuddy(r);
	numCols->setBuddy(c);
	md->addWidget(numRows);
	md->addWidget(r);
	md->addWidget(numCols);
	md->addWidget(c);
	md->addStretch();
	matDim->setLayout(md);

	registerField("plot.xbeg", r);
	registerField("plot.xend", c);

	QVBoxLayout* theLayout = new QVBoxLayout();
	theLayout->addWidget(matDim);
	this->setLayout(theLayout);
}

int PageRange::nextId() const
{
	return -1;
}



