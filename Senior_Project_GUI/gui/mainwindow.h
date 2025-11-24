#pragma once
#include <QMainWindow>
#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QMenuBar>
#include <QToolBar>
#include <QStatusBar>
#include <QAction>
#include <QFileDialog>
#include <QMessageBox>
#include <QTextEdit>
#include <QInputDialog>
#include "Controller.h"

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    MainWindow(QWidget* parent = nullptr);
    ~MainWindow();

private slots:
    void loadFile();
    void loadTestSequence();
    void loadSecondFile();
    void searchSequence();
    void findORFs();
    void alignSequences();
    void showAbout();

private:
    void setupUI();
    void setupMenus();
    void setupToolbar();
    void updateDisplay();

    QWidget* centralWidget;
    QVBoxLayout* mainLayout;
    QHBoxLayout* contentLayout;
    QTextEdit* sequenceDisplay;
    QTextEdit* resultsDisplay;

    QAction* loadAction;
    QAction* loadSecondAction;
    QAction* searchAction;
    QAction* orfAction;
    QAction* alignAction;
    QAction* exitAction;
    QAction* aboutAction;
    QAction* testAction;

    Controller controller;
};

