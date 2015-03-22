#ifndef __WIZARDS_H__
#define __WIZARDS_H__

#include <QWizard>
#include <QWizardPage>
#include <vector>
#include <QString>

class QTableWidget;
class MatrixPage;
class QRadioButton;

class NewWizard : public QWizard
{   
	Q_OBJECT

	public:
		NewWizard(QWidget *parent = 0);
		enum { Page_New, Page_Matrix };

		std::vector<std::vector<QString> > getMatrix() const;
	private slots:
			void showHelp();
	private:
		mutable MatrixPage* m_matrix_page;
		std::vector<std::vector<QString> > vec;
	public slots:
		void accept();
};

class NewPage : public QWizardPage
{
	Q_OBJECT
	public:
		NewPage(QWidget* parent = 0);
		int nextId() const;
	private:
};

class MatrixPage : public QWizardPage
{
	Q_OBJECT
	public:
		MatrixPage(QWidget* parent = 0);
		int nextId() const;
		virtual void initializePage();
		QTableWidget* getTable();
		std::vector<std::vector<QString> > getValues();
	protected:
		QTableWidget* m_table;
		std::vector<std::vector<QString> > vec;
	public slots:
		void on_cellChanged(int, int);
};

class RestoreWizard : public QWizard
{   
	Q_OBJECT

	public:
		RestoreWizard(QWidget* parent, bool _b, bool _q, bool _bq, bool _dr, bool _ds, bool _dc);
		enum { Page_Restore, Page_Pade, Page_Conclusion };
	private slots:
			void showHelp();
	protected:
		bool b;
		bool q;
		bool bq;
		bool dr;
		bool ds;
		bool dc;
};

class PageRestore : public QWizardPage
{
	Q_OBJECT
	public:
		PageRestore(QWidget* parent, bool _b, bool _q, bool _bq, bool _dr, bool _ds, bool _dc);
		int nextId() const;
	protected:
		bool b;
		bool q;
		bool bq;
		bool dr;
		bool ds;
		bool dc;
	protected:
		QRadioButton* taylor_single;
		QRadioButton* pade;
};

class PagePade : public QWizardPage
{
	Q_OBJECT
	public:
		PagePade(QWidget* parent = 0);
		int nextId() const;
};

class PageConclusion : public QWizardPage
{
	Q_OBJECT
	public:
		PageConclusion(QWidget* parent = 0);
		int nextId() const;
};

class PlotWizard : public QWizard
{   
	Q_OBJECT

	public:
		PlotWizard(QWidget* parent, bool _b, bool _q, bool _bq, bool _dr, bool _ds, bool _dc);
		enum { Page_Choose, Page_Range };
	private slots:
			void showHelp();
	protected:
		bool b;
		bool q;
		bool bq;
		bool dr;
		bool ds;
		bool dc;
};

class PageChoose : public QWizardPage
{
	Q_OBJECT
	public:
		PageChoose(QWidget* parent, bool _b, bool _q, bool _bq, bool _dr, bool _ds, bool _dc);
		int nextId() const;
	protected:
		bool b;
		bool q;
		bool bq;
		bool dr;
		bool ds;
		bool dc;
};

class PageRange : public QWizardPage
{
	Q_OBJECT
	public:
		PageRange(QWidget* parent = 0);
		int nextId() const;
};

#endif // __WIZARDS_H__
