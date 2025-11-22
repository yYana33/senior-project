#pragma once

#include <QMainWindow>
#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QTextEdit>
#include <QMenuBar>
#include <QToolBar>
#include <QStatusBar>
#include <QAction>
#include <QFileDialog>
#include <QMessageBox>
#include <QInputDialog>

class SequenceCanvas;
class ResultsPanel;
class Controller;

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    MainWindow(QWidget* parent = nullptr);
    ~MainWindow();

private slots:
    void loadFile();

    void searchSequence();
    void findORFs();
    void alignSequences();

    void showAbout();

private:
    void setupUI();
    void setupMenus();
    void setupToolbar();

    QTextEdit* sequenceDisplay;
    QTextEdit* resultsDisplay;

    QWidget* centralWidget;
    QVBoxLayout* mainLayout;
    QHBoxLayout* contentLayout;

    QAction* loadAction;
    QAction* searchAction;
    QAction* orfAction;
    QAction* alignAction;
    QAction* exitAction;
    QAction* aboutAction;

    QString currentSequence;
    QString currentHeader;
};
