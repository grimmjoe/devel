#ifndef __WIZARDS_H__
#define __WIZARDS_H__

#include <QWizard>
#include <QWizardPage>

class QTableWidget;

class NewWizard : public QWizard
{   
	Q_OBJECT

	public:
		NewWizard(QWidget *parent = 0);
		enum { Page_New, Page_Matrix };

		void setMatrix(QTableWidget* t) const;
		QTableWidget* getMatrix() const;
	private slots:
			void showHelp();
	private:
		mutable QTableWidget* m_table;
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
};

#endif // __WIZARDS_H__
