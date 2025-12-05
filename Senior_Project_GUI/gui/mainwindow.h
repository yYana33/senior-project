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
#include "SequenceVisualizer.h"
#include "GCContentCircle.h"

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
    void onViewModeChanged(SequenceVisualizer::VisualizationMode mode);
    void onZoomLevelChanged(double zoom);

private:
    void setupUI();
    void setupMenus();
    void setupToolbar();
    void setupViewToolbar();
    void updateDisplay();
    void clearSearchResults();
    void clearORFResults();
    void clearAlignmentResults();
    void clearAllResults();


    QAction* clearSearchAction;
    QAction* clearORFAction;
    QAction* clearAlignmentAction;
    QAction* clearAllAction;

    QWidget* centralWidget;
    QVBoxLayout* mainLayout;
    QHBoxLayout* contentLayout;
    QTextEdit* resultsDisplay;

    QAction* loadAction;
    QAction* loadSecondAction;
    QAction* searchAction;
    QAction* orfAction;
    QAction* alignAction;
    QAction* exitAction;
    QAction* aboutAction;
    QAction* testAction;

    Controller* controller;
    GCContentCircle* gcCircle;

    QToolBar* viewToolBar;
    QAction* detailedViewAction;
    QAction* featureViewAction;
    QAction* autoViewAction;
    QAction* zoomInAction;
    QAction* zoomOutAction;
    QAction* zoomResetAction;
    SequenceVisualizer* sequenceVisualizer;
    void updateVisualization();
};

